CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsdti(y,th,bi,ani,dir,det,n1,n2,n3,h,rho,lambda,
     1                  thnew,mask)
C
C   y        -  observed diffusion tensor data
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(A)^(1/3) 
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(6,n1,n2,n3),th(6,n1,n2,n3),thnew(6,n1,n2,n3),h,rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),dir(3,n1,n2,n3),
     2       det(n1,n2,n3)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k
      real*8 wij,adist,sw,swy(6),h3,thi(6),bii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,dtidist,sij,diri(3),anii,deti,z,sew
      external adist,dtidist
      logical aws,mask(n1,n2,n3) 
      aws=lambda.lt.1e20
      h3=h*h*h
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               sw=0.d0
               deti=det(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
               DO k=1,6
                  swy(k)=0.d0
                  thi(k)=th(k,i1,i2,i3)/deti
               END DO
C        if(thi(1).le.0.d0.or.thi(4).le.0.d0.or.thi(6).le.0.d0) THEN
C        call intpr("i1",2,i1,1)
C        call intpr("i2",2,i2,1)
C        call intpr("i3",2,i3,1)
C        call intpr("mask",4,mask(i1,i2,i3),1)
C         call dblepr("y",1,y(1,i1,i2,i3),6)
C         call dblepr("thi0",4,th(1,i1,i2,i3),6)
C         call dblepr("deti",4,deti,1)
C         call dblepr("bii",3,bii,1)
C         call dblepr("rho",3,rho,1)
C        call dblepr("y",1,y(1,i1,i2,i3),6)
C     END IF
               thi(1)=thi(1)+rho*bii
               thi(4)=thi(4)+rho*bii
               thi(6)=thi(6)+rho*bii
C               call dblepr("thi0",4,thi,6)
               call eigen3(thi,ew,ev,ierr)
               if(ierr.ne.0) THEN
               call intpr("ierr",4,ierr,1)
                  thi(1)=1
                  thi(2)=0
                  thi(3)=0
                  thi(4)=1
                  thi(5)=0
                  thi(6)=1
               ELSE
                  sew=ew(1)*ew(2)*ew(3)
                  sew=dexp(dlog(sew)/3.d0)
                  ew(1)=ew(1)/sew
                  ew(2)=ew(2)/sew
                  ew(3)=ew(3)/sew
                  thi(1)=ev(1,1)*ev(1,1)*ew(1)+ev(1,2)*ev(1,2)*ew(2)+
     1                   ev(1,3)*ev(1,3)*ew(3)
                  thi(2)=ev(1,1)*ev(2,1)*ew(1)+ev(1,2)*ev(2,2)*ew(2)+
     1                   ev(1,3)*ev(2,3)*ew(3)
                  thi(3)=ev(1,1)*ev(3,1)*ew(1)+ev(1,2)*ev(3,2)*ew(2)+
     1                   ev(1,3)*ev(3,3)*ew(3)
                  thi(4)=ev(2,1)*ev(2,1)*ew(1)+ev(2,2)*ev(2,2)*ew(2)+
     1                   ev(2,3)*ev(2,3)*ew(3)
                  thi(5)=ev(2,1)*ev(3,1)*ew(1)+ev(2,2)*ev(3,2)*ew(2)+
     1                   ev(2,3)*ev(3,3)*ew(3)
                  thi(6)=ev(3,1)*ev(3,1)*ew(1)+ev(3,2)*ev(3,2)*ew(2)+
     1                   ev(3,3)*ev(3,3)*ew(3)
               END IF
C               call intpr("i3",2,i3,1)
C               call dblepr("thi",3,thi,6)
               anii=ani(i1,i2,i3)
               DO k=1,3
                  diri(k)=dir(k,i1,i2,i3)
               END DO
               call rangex(thi,h,j1a,j1e)
C               call intpr("j1e",3,j1e,1)
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey(thi,j1,h,j2a,j2e)
C               call intpr("j2e",3,j2e,1)
                 DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez(thi,j1,j2,h,j3a,j3e)
C               call intpr("j3e",3,j3e,1)
                      DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE 
                        wij=adist(thi,j1,j2,j3)
                        if(wij.lt.0.d0) call dblepr("wij",3,wij,1)
C     triangular location kernel
                        if(wij.gt.h3) CYCLE
C                           call dblepr("h3",2,h3,1)
C                           call dblepr("outside",7,wij,1)
C                           CYCLE
                        wij = 1.d0 - wij/h3
                        IF(aws) THEN
                        sij=dtidist(diri,dir(1,jj1,jj2,jj3),anii)*bii
                           if(sij.lt.0.d0) THEN
C                              call dblepr("sij",3,sij,1)
                              sij=0.d0
                           END IF
                           if(sij.gt.1.d0) CYCLE
                           wij=wij*(1.d0-sij)
                        END IF
                        sw=sw+wij
                        DO k=1,6
                           swy(k)=swy(k)+wij*y(k,jj1,jj2,jj3)
                        END DO
                     END DO
                  END DO
               END DO
               if(sw.eq.0) THEN
                  call intpr("sw=0 in i1",9,i1,1)
                  call intpr("i2",2,i2,1)
                  call intpr("i3",2,i3,1)
               END IF
               bi(i1,i2,i3)=sw
               DO k=1,6
                  thnew(k,i1,i2,i3)=swy(k)/sw
               END DO
               call eigen3(thnew(1,i1,i2,i3),ew,ev,ierr)
               if(ierr.ne.0) THEN
                  ani(i1,i2,i3)=0.d0
                  DO k=1,3
                     dir(k,i1,i2,i3)=0.d0
                  END DO
                  det(i1,i2,i3)=1
               ELSE
                  mew=(ew(1)+ew(2)+ew(3))/3.d0
                  z1=ew(1)-mew
                  z2=ew(2)-mew
                  z3=ew(3)-mew
                  z=3.d0*(z1*z1+z2*z2+z3*z3)
                  z1=ew(1)
                  z2=ew(2)
                  z3=ew(3)
                  mew=2.d0*(z1*z1+z2*z2+z3*z3)
                  if(mew.le.1d-20) mew=1.d0
                  ani(i1,i2,i3)=dsqrt(z/mew)
                  DO k=1,3
                     dir(k,i1,i2,i3)=ev(k,1)
                  END DO
                  z=ew(1)*ew(2)*ew(3)
                  IF(z.le.1d-20) THEN
                     det(i1,i2,i3)=0.d0
                  ELSE
                     det(i1,i2,i3)=dexp(dlog(z)/3)
                  END IF
               ENDIF
               call rchkusr()
            END DO
         END DO
C         call intpr("i1",2,i1,1)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsdti2(y,th,bi,ani,andir,det,bcov,sigma2,n1,n2,n3,h,
     1                  zext,rho,lambda,thnew,mask)
C
C   y        -  observed diffusion tensor data
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(A)^(1/3) 
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(6,n1,n2,n3),th(6,n1,n2,n3),thnew(6,n1,n2,n3),h,rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),andir(3,n1,n2,n3),
     2       det(n1,n2,n3),bcov(6,6),sigma2(n1,n2,n3),zext
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k
      real*8 wij,adist,sw,swy(6),h3,thi(6),bii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,dtidist2,sij,anii,deti,z,sew
      external adist,dtidist2
      logical aws,mask(n1,n2,n3) 
      aws=lambda.lt.1e20
      h3=h*h*h
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               sw=0.d0
               deti=det(i1,i2,i3)
               bii=bi(i1,i2,i3)
               DO k=1,6
                  swy(k)=0.d0
                  thi(k)=th(k,i1,i2,i3)/deti
               END DO
               thi(1)=thi(1)+rho*bii
               thi(4)=thi(4)+rho*bii
               thi(6)=thi(6)+rho*bii
               call eigen3(thi,ew,ev,ierr)
               if(ierr.ne.0) THEN
                  call intpr("ierr",4,ierr,1)
                  thi(1)=1
                  thi(2)=0
                  thi(3)=0
                  thi(4)=1
                  thi(5)=0
                  thi(6)=1
               ELSE
                  sew=ew(1)*ew(2)*ew(3)
                  sew=dexp(dlog(sew)/3.d0)
                  ew(1)=ew(1)/sew
                  ew(2)=ew(2)/sew
                  ew(3)=ew(3)/sew
                  thi(1)=ev(1,1)*ev(1,1)*ew(1)+ev(1,2)*ev(1,2)*ew(2)+
     1                   ev(1,3)*ev(1,3)*ew(3)
                  thi(2)=ev(1,1)*ev(2,1)*ew(1)+ev(1,2)*ev(2,2)*ew(2)+
     1                   ev(1,3)*ev(2,3)*ew(3)
                  thi(3)=ev(1,1)*ev(3,1)*ew(1)+ev(1,2)*ev(3,2)*ew(2)+
     1                   ev(1,3)*ev(3,3)*ew(3)
                  thi(4)=ev(2,1)*ev(2,1)*ew(1)+ev(2,2)*ev(2,2)*ew(2)+
     1                   ev(2,3)*ev(2,3)*ew(3)
                  thi(5)=ev(2,1)*ev(3,1)*ew(1)+ev(2,2)*ev(3,2)*ew(2)+
     1                   ev(2,3)*ev(3,3)*ew(3)
                  thi(6)=ev(3,1)*ev(3,1)*ew(1)+ev(3,2)*ev(3,2)*ew(2)+
     1                   ev(3,3)*ev(3,3)*ew(3)
               END IF
               anii=ani(i1,i2,i3)
               call rangex(thi,h,j1a,j1e)
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey(thi,j1,h,j2a,j2e)
                 DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez(thi,j1,j2,h,j3a,j3e,zext)
                      DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE 
                        wij=adist(thi,j1,j2,j3,zext)
C     triangular location kernel
                        if(wij.gt.h3) CYCLE
                        wij = (1.d0 - wij/h3)
                        IF(aws) THEN
                        sij=dtidist2(th(1,i1,i2,i3),
     1                      th(1,jj1,jj2,jj3),bcov)*bii/lambda
                           if(sij.le.0.d0.and.j3.ne.0) THEN
                              call dblepr("sij",3,sij,1)
                              call dblepr("bii",3,bii,1)
                              call dblepr("lam",3,lambda,1)
                              call dblepr("sig",3,sigma2,1)
                              sij=0.d0
                           END IF
                           if(sij.gt.1.d0) CYCLE
                           wij=wij*(1.d0-sij)/sigma2(jj1,jj2,jj3)
                        END IF
                        sw=sw+wij
                        DO k=1,6
                           swy(k)=swy(k)+wij*y(k,jj1,jj2,jj3)
                        END DO
                     END DO
                  END DO
               END DO
               if(sw.eq.0) THEN
                  call intpr("sw=0 in i1",9,i1,1)
                  call intpr("i2",2,i2,1)
                  call intpr("i3",2,i3,1)
               END IF
               bi(i1,i2,i3)=sw
               DO k=1,6
                  thnew(k,i1,i2,i3)=swy(k)/sw
               END DO
               call eigen3(thnew(1,i1,i2,i3),ew,ev,ierr)
               if(ierr.ne.0) THEN
                  ani(i1,i2,i3)=0.d0
                  det(i1,i2,i3)=1
                  DO k=1,3
                     andir(k,i1,i2,i3)=0.d0
                  END DO
               ELSE
                  mew=(ew(1)+ew(2)+ew(3))/3.d0
                  z1=ew(1)-mew
                  z2=ew(2)-mew
                  z3=ew(3)-mew
                  z=3.d0*(z1*z1+z2*z2+z3*z3)
                  z1=ew(1)
                  z2=ew(2)
                  z3=ew(3)
                  mew=2.d0*(z1*z1+z2*z2+z3*z3)
                  if(mew.le.1d-20) mew=1.d0
                  ani(i1,i2,i3)=dsqrt(z/mew)
                  z=ew(1)*ew(2)*ew(3)
                  DO k=1,3
                     andir(k,i1,i2,i3)=ev(k,3)
                  END DO
                  IF(z.le.1d-20) THEN
                     det(i1,i2,i3)=0.d0
                  ELSE
                     det(i1,i2,i3)=dexp(dlog(z)/3)
                  END IF
               ENDIF
               call rchkusr()
            END DO
         END DO
C         call intpr("i1",2,i1,1)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize anisotropy index and direction of main anisotropy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine initdti(th,n1,n2,n3,ani,dir,det,mask)
C
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
C   ani      -  anisotropy index (computed internally)
C   dir      -  direction of main anisotropy (computed internally)
C   det      -  det(A)^(1/3) 
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 th(6,n1,n2,n3),ani(n1,n2,n3),dir(3,n1,n2,n3),det(n1,n2,n3)
      integer i1,i2,i3,ierr,k
      real*8 ew(3),ev(3,3),mew,z1,z2,z3,z
      logical mask(n1,n2,n3)
C  compute anisotropy index and direction of main anisotropy (nneded in statistical penalty)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               mask(i1,i2,i3)=.TRUE.
               call eigen3(th(1,i1,i2,i3),ew,ev,ierr)
               if(ierr.ne.0) THEN
                  ani(i1,i2,i3)=0
                  DO k=1,3
                     dir(k,i1,i2,i3)=0
                  END DO
                  mask(i1,i2,i3)=.FALSE.
               ELSE
                  mew=(ew(1)+ew(2)+ew(3))/3.d0
                  z1=ew(1)-mew
                  z2=ew(2)-mew
                  z3=ew(3)-mew
                  z=3.d0*(z1*z1+z2*z2+z3*z3)
                  z1=ew(1)
                  z2=ew(2)
                  z3=ew(3)
                  mew=2.d0*(z1*z1+z2*z2+z3*z3)
                  if(mew.le.1d-20) mew=1.d0
                  ani(i1,i2,i3)=dsqrt(z/mew)
                  DO k=1,3
                     dir(k,i1,i2,i3)=ev(k,3)
                  END DO
                  z=ew(1)*ew(2)*ew(3)
                  IF(z.le.1d-20) THEN
                     det(i1,i2,i3)=0.d0
                     mask(i1,i2,i3)=.FALSE.
                  ELSE
                     det(i1,i2,i3)=dexp(dlog(z)/3)
                  END IF
                  IF(ew(3).le.0.d0) mask(i1,i2,i3)=.FALSE.
                  IF(ew(2).le.0.d0) mask(i1,i2,i3)=.FALSE.
                  IF(ew(1).le.0.d0) mask(i1,i2,i3)=.FALSE.
               ENDIF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize anisotropy index and direction of main anisotropy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine projdt(th,n1,n2,n3,thnew,mask)
C
C   th       -  observed diffusion tensor data
C   thnew    -  projected tensor data
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 th(6,n1,n2,n3),thnew(6,n1,n2,n3)
      integer i1,i2,i3,ierr,k
      real*8 ew(3),ev(3,3),z1,z2,z3
      logical mask(n1,n2,n3)
C  compute anisotropy index and direction of main anisotropy (nneded in statistical penalty)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               mask(i1,i2,i3)=.TRUE.
               call eigen3(th(1,i1,i2,i3),ew,ev,ierr)
               IF(ierr.ne.0) THEN
                  DO k=1,6
                     thnew(k,i1,i2,i3)=0.d0
                  END DO
                  mask(i1,i2,i3)=.FALSE.
               ELSE IF(dmin1(ew(1),ew(2),ew(3)).lt.1d-10) THEN
                  z1=dmax1(ew(1),1d-10)
                  z2=dmax1(ew(2),1d-10)
                  z3=dmax1(ew(3),1d-10)
                  thnew(1,i1,i2,i3)=z1*ev(1,1)*ev(1,1)+
     1                     z2*ev(1,2)*ev(1,2)+z2*ev(1,3)*ev(1,3)
                  thnew(2,i1,i2,i3)=z1*ev(1,1)*ev(2,1)+
     1                     z2*ev(1,2)*ev(2,2)+z2*ev(1,3)*ev(2,3)
                  thnew(3,i1,i2,i3)=z1*ev(1,1)*ev(3,1)+
     1                     z2*ev(1,2)*ev(3,2)+z2*ev(1,3)*ev(3,3)
                  thnew(4,i1,i2,i3)=z1*ev(2,1)*ev(2,1)+
     1                     z2*ev(2,2)*ev(2,2)+z2*ev(2,3)*ev(2,3)
                  thnew(5,i1,i2,i3)=z1*ev(2,1)*ev(3,1)+
     1                     z2*ev(2,2)*ev(3,2)+z2*ev(2,3)*ev(3,3)
                  thnew(6,i1,i2,i3)=z1*ev(3,1)*ev(3,1)+
     1                     z2*ev(3,2)*ev(3,2)+z2*ev(3,3)*ev(3,3)
               ELSE
                  DO k=1,6
                     thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
                  END DO
               ENDIF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3(y,lambda,theta,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(6),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-50
      l=1
      DO i=1,3
         DO j=i,3
	    a(i,j)=y(l)
            l=l+1
         END DO
      END DO
      lwork=104
      liwork=50
      call dsyevr('V','A','U',n,a,n,vl,vu,1,n,eps,m,lambda,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
C      if(lambda(3).lt.dmax1(lambda(2),lambda(1)))
C     1       call dblepr("reverse ev",10,lambda,3)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute anisotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function adist(a,x,y,z,zext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer x,y,z
      real*8 a(6),zz,zext
      zz=z*zext
      adist=a(1)*x*x+a(4)*y*y+a(6)*zz*zz+
     1               2.d0*(a(2)*x*y+a(3)*x*zz+a(5)*y*zz)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute statistical penalty for dti
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function dtidist(diri,dirj,ani)
      implicit logical (a-z)
      real*8 diri(3),dirj(3),ani,z
      z=(diri(1)*dirj(1)+diri(2)*dirj(2)+diri(3)*dirj(3))
      dtidist=ani*ani*(1.d0-z*z)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute statistical penalty for dti
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function dtidist2(th1,th2,Bcov)
      implicit logical (a-z)
      real*8 th1(6),th2(6),Bcov(6,6),z,zd,zd2
      integer i,j
      z=0.d0
      DO i=1,6
         zd=th1(i)-th2(i)
         z=z+zd*zd*Bcov(i,i)
         if(i.lt.6) THEN
            DO j=i+1,6
               zd2=th1(j)-th2(j)
               z=z+2.d0*zd*zd2*Bcov(i,j)
            END DO
         END IF
      END DO
      dtidist2=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangex(a,h,ia,ie)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ia,ie
      real*8 a(6),h
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      t=p5*p4/p2/p3
      z=(p1-p44/p2-p55/p3+2.d0*p6*t)*s-p66*p55/p3-p66*p44/p2+
     1                                             2.d0*p66*p6*t
      if(dabs(z).le.1e-40) THEN
         call dblepr("denominator in rangex",21,z,1)
         z=1
      END IF
      z=s/z
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ia=-z*h
      ie=z*h
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangey(a,ix,h,ja,je)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ja,je,ix
      real*8 a(6),h
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44,x
      x=ix/h
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      if(s.le.1e-10) THEN
         call dblepr("denominator in rangey",21,s,1)
         call dblepr("a",1,a,6)
         call dblepr("z",1,z,1)
         z=1
      END IF
      t=-(p4*p3 - p6*p5)* x
      z=(p44*p3*p3-2.d0*p4*p3*p6*p5+p66*p55-p1*p3*s+p55*s)*x*x+p3*s
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ja=(t-z)*h/s
      je=(t+z)*h/s
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangez(a,ix,jy,h,ka,ke,zext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ka,ke,ix,jy
      real*8 a(6),h
      real*8 p1,p2,p3,p4,p5,p6,z,t,x,y,zext
      x=ix/h
      y=jy/h
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      z=(p6*y + p5*x)
      t= - z/p3
      z=(-z*t - p2*y*y - p1*x*x - 2.d0*p4*x*y + 1.d0)/p3
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ka=(t-z)*h/zext
      ke=(t+z)*h/zext
      RETURN
      END
