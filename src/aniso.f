CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsdti(y,th,bi,ani,dir,n1,n2,n3,h,rho,lambda,thnew)
C
C   y        -  observed diffusion tensor data
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
      integer n1,n2,n3
      real*8 y(6,n1,n2,n3),th(6,n1,n2,n3),thnew(6,n1,n2,n3),rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),dir(3,n1,n2,n3)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3a,j3e,jj3,ierr,k
      real*8 wij,adist,sw,swy(6),h3,thi(6),bii,ewert(3),evect(3,3),
     1       mew,z1,z2,z3,dtidist,sij,diri(3),anii
      external adist,dtidist
      logical aws 
      aws=lambda.lt.1e20
      h3=h*h*h
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               sw=0.d0
               bii=bi(i1,i2,i3)
               DO k=1,6
                  swy(k)=0.d0
                  thi(k)=th(k,i1,i2,i3)
               END DO
               anii=ani(i1,i2,i3)
               DO k=1,3
                  diri(k)=dir(k,i1,i2,i3)
               END DO
               call rangex(thi,h,rho,bii,j1a,j1e)
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey(thi,j1,h,rho,bii,j2a,j2e)
                  DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez(thi,j1,j2,h,rho,bii,j3a,j3e)
                     DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        wij=adist(thi,j1,j2,j3,rho,bii)
C     triangular location kernel
                        if(wij.gt.h3) THEN
                           call dblepr("outside",7,wij,1)
                           CYCLE
                        END IF
                        wij = 1.d0 - wij/h3
                        IF(aws) THEN
                        sij=dtidist(diri,dir(1,jj1,jj2,jj3),anii)/bii
                           if(sij.gt.lambda) CYCLE
                           wij=wij*(1.d0-sij/lambda)
                        END IF
                        sw=sw+wij
                        DO k=1,6
                           swy(k)=swy(k)+wij*y(k,jj1,jj2,jj3)
                        END DO
                     END DO
                  END DO
               END DO
               bi(i1,i2,i3)=sw
               DO k=1,6
                  thnew(k,i1,i2,i3)=swy(k)/sw
               END DO
               call eigen3(thnew(1,i1,i2,i3),ewert,evect,ierr)
               mew=(ewert(1)+ewert(2)+ewert(3))/3.d0
               z1=ewert(1)-mew
               z2=ewert(2)-mew
               z3=ewert(3)-mew
               ani(i1,i2,i3)=z1*z1+z2*z2+z3*z3
               DO k=1,3
                  dir(k,i1,i2,i3)=evect(k,1)
               END DO
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
      subroutine initdti(th,n1,n2,n3,ani,dir)
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
      integer n1,n2,n3
      real*8 th(6,n1,n2,n3),ani(n1,n2,n3),dir(3,n1,n2,n3)
      integer i1,i2,i3,ierr,k
      real*8 ewert(3),evect(3,3),mew,z1,z2,z3
C  compute anisotropy index and direction of main anisotropy (nneded in statistical penalty)
         DO i1=1,n1
            DO i2=1,n2
               DO i3=1,n3
                  call eigen3(th(1,i1,i2,i3),ewert,evect,ierr)
                  mew=(ewert(1)+ewert(2)+ewert(3))/3.d0
                  z1=ewert(1)-mew
                  z2=ewert(2)-mew
                  z3=ewert(3)-mew
                  ani(i1,i2,i3)=z1*z1+z2*z2+z3*z3
                  DO k=1,3
                     dir(k,i1,i2,i3)=evect(k,1)
                  END DO
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
      integer i,j,l,ISUPPZ(2),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-20
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
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute anisotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function adist(a,x,y,z,rho,ni)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      integer x,y,z
      real*8 a(6),rho,ni
      real*8 p1,p2,p3,s
      s=rho/ni
      p1=a(1)+s
      p2=a(4)+s
      p3=a(6)+s
      adist=p1*x*x+p2*y*y+p3*z*z+2.d0*(a(2)*x*y+a(3)*x*z+a(5)*y*z)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute statistical penalty for dti
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function dtidist(diri,dirj,ani)
      real*8 diri(3),dirj(3),ani
      dtidist=ani*(diri(1)*dirj(1)+diri(2)*dirj(2)+diri(3)*dirj(3))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangex(a,h,rho,ni,ia,ie)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      integer ia,ie
      real*8 a(6),h,rho,ni
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44
      z=rho/ni
      p1=a(1)+z
      p2=a(4)+z
      p3=a(6)+z
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      t=p5*p4/p2/p3
      z=s/((p1-p44/p2-p55/p3+2.d0*p6*t)*s-p66*p55/p3-p66*p44/p2+
     1                                             2.d0*p66*p6*t)
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
      subroutine rangey(a,ix,h,rho,ni,ja,je)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      integer ja,je,ix
      real*8 a(6),h,rho,ni
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44,x
      x=ix/h
      z=rho/ni
      p1=a(1)+z
      p2=a(4)+z
      p3=a(6)+z
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
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
      subroutine rangez(a,ix,jy,h,rho,ni,ka,ke)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      integer ka,ke,ix,jy
      real*8 a(6),h,rho,ni
      real*8 p1,p2,p3,p4,p5,p6,z,t,x,y
      x=ix/h
      y=jy/h
      z=rho/ni
      p1=a(1)+z
      p2=a(4)+z
      p3=a(6)+z
      p4=a(2)
      p5=a(3)
      p6=a(5)
      t= - (p6*y + p5*x)/p3
      z=(p6*y + p5*x)
      z=(z*z/p3 - p2*y*y - p1*x*x - 2.d0*p4*x*y + 1.d0)/p3
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ka=(t-z)*h
      ke=(t+z)*h
      RETURN
      END
