CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsrice(s0,si,th,bi,ani,andir,det,bcov,solvebtb,
     1 sigma2,sigma2h,n1,n2,n3,ngrad,h,zext,rho,lambda,thnew,s0new,
     2 sinew,sigma2n,swsi,eps,wwij,nij,s0ij,siij,sivar,ls0i)
C
C   y        -  observed diffusion tensor data
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(A)
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,nij
      real*8 s0(n1,n2,n3),th(6,n1,n2,n3),thnew(6,n1,n2,n3),h,rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),andir(3,n1,n2,n3),
     2       det(n1,n2,n3),bcov(6,6),sigma2,zext,wwij(nij),
     3       si(n1,n2,n3,ngrad),solvebtb(6,ngrad),swsi(ngrad),
     4       s0new(n1,n2,n3),sigma2h(n1,n2,n3),sigma2n(n1,n2,n3),
     5       s0ij(nij),siij(nij,ngrad),sinew(n1,n2,n3,ngrad),s0var,
     6       sivar(ngrad),ls0i(ngrad)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,l,iwij
      real*8 wij,adist,sw,sws0,h3,thi(6),bii,sqrbii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,dtidist2,sij,deti,z,sew,eps,eps3,ss2,sw0
      external adist,dtidist2
      logical aws
      aws=lambda.lt.1e20
      h3=h*h*h
      eps3=eps*eps*eps
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               sw=0.d0
               sw0=0.d0
               sws0=0.d0
               ss2=0.d0
               DO k=1,ngrad
                  swsi(k)=0.d0
               END DO
               deti=exp(log(det(i1,i2,i3))/3)
               bii=bi(i1,i2,i3)
               sqrbii=sigma2h(i1,i2,i3)/sqrt(bii)
               DO k=1,6
                  thi(k)=th(k,i1,i2,i3)/deti
               END DO
               thi(1)=thi(1)+rho*sqrbii
               thi(4)=thi(4)+rho*sqrbii
               thi(6)=thi(6)+rho*sqrbii
C  this is scale invariant sice sqrbii scales with sqrt(sigma2) (standard deviation)
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
                  sew=exp(log(sew)/3.d0)
                  ew(1)=ew(1)/sew
                  ew(2)=ew(2)/sew
                  ew(3)=ew(3)/sew
                  thi(1)=ev(1,1)*ev(1,1)/ew(1)+ev(1,2)*ev(1,2)/ew(2)+
     1                   ev(1,3)*ev(1,3)/ew(3)
                  thi(2)=ev(1,1)*ev(2,1)/ew(1)+ev(1,2)*ev(2,2)/ew(2)+
     1                   ev(1,3)*ev(2,3)/ew(3)
                  thi(3)=ev(1,1)*ev(3,1)/ew(1)+ev(1,2)*ev(3,2)/ew(2)+
     1                   ev(1,3)*ev(3,3)/ew(3)
                  thi(4)=ev(2,1)*ev(2,1)/ew(1)+ev(2,2)*ev(2,2)/ew(2)+
     1                   ev(2,3)*ev(2,3)/ew(3)
                  thi(5)=ev(2,1)*ev(3,1)/ew(1)+ev(2,2)*ev(3,2)/ew(2)+
     1                   ev(2,3)*ev(3,3)/ew(3)
                  thi(6)=ev(3,1)*ev(3,1)/ew(1)+ev(3,2)*ev(3,2)/ew(2)+
     1                   ev(3,3)*ev(3,3)/ew(3)
               END IF
               call rangex(thi,h,j1a,j1e)
               iwij=0
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
                        wij=adist(thi,j1,j2,j3,zext)
C     triangular location kernel
                        if(wij.gt.h3) CYCLE
                        wij = (1.d0 - wij/h3)
                        IF(aws) THEN
                           sij=dtidist2(th(1,i1,i2,i3),
     1                          th(1,jj1,jj2,jj3),bcov)*bii/lambda
                           if(sij.gt.1.d0) CYCLE
                           wij=wij*(1.d0-sij)
                        END IF
                        iwij=iwij+1
                        wwij(iwij)=wij
                        s0ij(iwij)=s0(jj1,jj2,jj3)
                        DO k=1,ngrad
                           siij(iwij,k)=si(jj1,jj2,jj3,k)
                        END DO
                     END DO
                  END DO
               END DO
C
C     now calculate the unbiased estimates of image intensities
C
               call ricefix(s0ij,s0new(i1,i2,i3),sigma2,wwij,iwij,eps)
            call ricevar(s0ij,s0new(i1,i2,i3),sigma2,wwij,iwij,s0var)
               DO k=1,ngrad
                  call ricefix(siij(1,k),sinew(i1,i2,i3,k),
     1                         sigma2,wwij,iwij,eps)
                  call ricevar(siij(1,k),sinew(i1,i2,i3,k),
     1                         sigma2,wwij,iwij,sivar(k))
               END DO
               DO k=1,ngrad
                  sivar(k)=s0var+sivar(k)
                  ls0i(k)=log(s0new(i1,i2,i3))-log(sinew(i1,i2,i3,k))
C   may need to check if the arguments are positive
               END DO
C
C     now we have  log(S_0) - log(S_i)   in   ls0i 
C                  var( log(S_0) - log(S_i) ) in  sivar(k)
C     we now need to estimate the tensor
C

               bi(i1,i2,i3)=sw
               IF(sws0.gt.0.d0) THEN
                  IF(sw.gt.0.d0.and.sw.lt.1.d20) THEN
                     s0new(i1,i2,i3)=sws0/sw
                     sigma2n(i1,i2,i3)=ss2/sw0
                   ELSE
                     s0new(i1,i2,i3)=s0(i1,i2,i3)
                     sigma2n(i1,i2,i3)=sigma2h(i1,i2,i3)
                  END IF
                  sws0=dlog(sws0)
                  DO k=1,ngrad
                     IF(swsi(k).gt.0.d0) THEN
                        swsi(k)=sws0-log(swsi(k))
                     ELSE
                        swsi(k)=0.d0
                     END IF
                  END DO
               ELSE
                  s0new(i1,i2,i3)=s0(i1,i2,i3)
                  DO k=1,ngrad
                     swsi(k)=0.d0
                  END DO 
                  sigma2n(i1,i2,i3)=sigma2h(i1,i2,i3)
               END IF
               DO k=1,6
                  z=0.d0
                  DO l=1,ngrad
                     z=z+solvebtb(k,l)*swsi(l)
                  END DO
                  thnew(k,i1,i2,i3)=z
               END DO
               call eigen3(thnew(1,i1,i2,i3),ew,ev,ierr)
               IF(ierr.ne.0) THEN
C
C       error in eigenvalue decomposition of tensor, set voxel inactive
C
                  thnew(1,i1,i2,i3)=eps
                  thnew(2,i1,i2,i3)=0.d0
                  thnew(3,i1,i2,i3)=0.d0
                  thnew(4,i1,i2,i3)=eps
                  thnew(5,i1,i2,i3)=0.d0
                  thnew(6,i1,i2,i3)=eps
               ELSE IF(dmin1(ew(1),ew(2),ew(3)).lt.eps) THEN
C
C       negative eigenvalue of tensor, project to space of positive definite tensors
C
                  ew(1)=max(ew(1),eps)
                  ew(2)=max(ew(2),eps)
                  ew(3)=max(ew(3),eps)
                  thnew(1,i1,i2,i3)=ew(1)*ev(1,1)*ev(1,1)+
     1                     ew(2)*ev(1,2)*ev(1,2)+ew(3)*ev(1,3)*ev(1,3)
                  thnew(2,i1,i2,i3)=ew(1)*ev(1,1)*ev(2,1)+
     1                     ew(2)*ev(1,2)*ev(2,2)+ew(3)*ev(1,3)*ev(2,3)
                  thnew(3,i1,i2,i3)=ew(1)*ev(1,1)*ev(3,1)+
     1                     ew(2)*ev(1,2)*ev(3,2)+ew(3)*ev(1,3)*ev(3,3)
                  thnew(4,i1,i2,i3)=ew(1)*ev(2,1)*ev(2,1)+
     1                     ew(2)*ev(2,2)*ev(2,2)+ew(3)*ev(2,3)*ev(2,3)
                  thnew(5,i1,i2,i3)=ew(1)*ev(2,1)*ev(3,1)+
     1                     ew(2)*ev(2,2)*ev(3,2)+ew(3)*ev(2,3)*ev(3,3)
                  thnew(6,i1,i2,i3)=ew(1)*ev(3,1)*ev(3,1)+
     1                     ew(2)*ev(3,2)*ev(3,2)+ew(3)*ev(3,3)*ev(3,3)
               ELSE
C
C       valid tensor tensor, copy data
C
C                  DO k=1,6
C                     thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
C                  END DO
               END IF
               mew=(ew(1)+ew(2)+ew(3))/3.d0
               z1=ew(1)-mew
               z2=ew(2)-mew
               z3=ew(3)-mew
               z=3.d0*(z1*z1+z2*z2+z3*z3)
               z1=ew(1)
               z2=ew(2)
               z3=ew(3)
               mew=2.d0*(z1*z1+z2*z2+z3*z3)
               if(mew.le.eps*eps) mew=1.d0
               ani(i1,i2,i3)=sqrt(z/mew)
               det(i1,i2,i3)=ew(1)*ew(2)*ew(3)
               DO k=1,3
                  andir(k,i1,i2,i3)=ev(k,3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
