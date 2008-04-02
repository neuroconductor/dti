CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize anisotropy index and direction of main anisotropy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine projdt(th,n1,n2,n3,thnew,ani,dir,det,mask)
C
C   th       -  observed diffusion tensor data
C   thnew    -  projected tensor data
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 th(6,n1,n2,n3),thnew(6,n1,n2,n3),ani(n1,n2,n3),
     1       dir(3,n1,n2,n3),det(n1,n2,n3),mew
      integer i1,i2,i3,ierr,k
      real*8 ew(3),ev(3,3),z1,z2,z3,z
      logical mask(n1,n2,n3)
C  compute anisotropy index and direction of main anisotropy (nneded in statistical penalty)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF(mask(i1,i2,i3)) THEN
                  call eigen3(th(1,i1,i2,i3),ew,ev,ierr)
                  IF(ierr.ne.0) THEN
C
C       error in eigenvalue decomposition of tensor, set voxel inactive
C
                     DO k=1,6
                        thnew(k,i1,i2,i3)=0.d0
                     END DO
                     mask(i1,i2,i3)=.FALSE.
                  ELSE IF(min(ew(1),ew(2)).lt.1d-5*ew(3)) THEN
C
C       negative eigenvalue of tensor, project to space of positive definite tensors
C
                     ew(1)=max(ew(1),1d-5*ew(3))
                     ew(2)=max(ew(2),1d-5*ew(3))
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
                     DO k=1,6
                        thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
                     END DO
                  END IF
                  IF(mask(i1,i2,i3)) THEN
C
C       compute anisotropy index, direction corresponding to the first eigenvalue
C       and determinant
C
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
                     ani(i1,i2,i3)=sqrt(z/mew)
                     DO k=1,3
                        dir(k,i1,i2,i3)=ev(k,3)
                     END DO
                     z=ew(1)*ew(2)*ew(3)
                     IF(z.le.1d-30) THEN
                        det(i1,i2,i3)=0.d0
                        mask(i1,i2,i3)=.FALSE.
                     ELSE
                        det(i1,i2,i3)=z
                     END IF
                  END IF
               ELSE
C
C    inactive voxel, just copy the data
C
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
      if(abs(z).le.1e-40) THEN
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
