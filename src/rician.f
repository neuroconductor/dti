      subroutine ricefix(y,theta,sigma2,w,n,eps)
      implicit logical (a-z)
      integer n
      real*8 y(n),theta,sigma2,w(n),eps
      integer i
      real*8 besseli,swi,z,z0,z1,theta0,swbi
      external besseli
      swi=0.d0
      theta0=theta+1
      Do i=1,n
         swi=swi+w(i)
      END DO
      DO WHILE (abs(theta0-theta).gt.eps)
         theta0=theta
         swbi=0.d0
         DO i=1,n
            z=y(i)*theta/sigma2
            if(z.lt.27) THEN
               z0=besseli(z,0.d0,1.d0)
               z1=besseli(z,1.d0,1.d0)
               z=z1/z0
            ELSE
               z0=.375/z
               z=z/(z+.5d0+z0+z0/z)
            END IF
            swbi=swbi+w(i)*z*y(i)
         END DO
         theta=swbi/swi
      END DO
      return
      end
      subroutine ricevar(y,theta,sigma2,w,n,vtheta)
C
C  this computes an estimate of D log(\hat{theta})
C
      implicit logical (a-z)
      integer n
      real*8 y(n),theta,sigma2,w(n),vtheta
      integer i
      real*8 besseli,swi,z,z0,z1,z2,swvi,wi,wi2,zv,theta2,pi
      PARAMETER (pi=3.1415927)
      swi=0.d0
      swvi=0.d0
      theta2=theta*theta
      z=theta2/sigma2/2.d0
      z2=z/2.d0
      if(z.lt.2.7d1) THEN
      z=exp(-z2)*((1+z)*besseli(z2,0.d0,1.d0)+z*besseli(z2,1.d0,1.d0))
      z=z*z
      ELSE
         z1=(1.d0+2.d0*z/(1.d0+0.25d0/z+3d-2/z/z))
         z=z1*z1/pi/z
      END IF
      zv=2.d0*sigma2+theta2-0.5D0*pi*sigma2*z
      Do i=1,n
         wi=w(i)
         wi2=wi*wi
         swi=swi+wi
         z=y(i)*theta/sigma2
         if(z.lt.27) THEN
            z0=besseli(z,0.d0,1.d0)
            z1=besseli(z,1.d0,1.d0)
            z=z1/z0
         ELSE
            z0=.375/z
            z=z/(z+.5d0+z0+z0/z)
         END IF
         swvi=swvi+wi2*z*z
      END DO
      vtheta = zv*swvi/theta2/swi/swi
      RETURN
      END
      