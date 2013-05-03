      subroutine ftensor(par,s,nb,b,vinv,gv,fv)
C
C  compute f(par) 
C
      implicit logical (a-z)
      integer nb
      real*8 par(7),s(nb),b(6,nb),vinv(nb),gv(nb),fv
      integer i
      real*8 D(6),th0,rss,res
      th0=par(1)
      call rho2D(par(2),D)
      call sihat(th0,D,b,gv,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in gv
C
      rss=0.d0
      DO i=1,nb
         res=s(i)-gv(i)
         rss=rss+res*res*vinv(i)
      END DO
      fv=rss
C
C   now we have value of the criterion in fv
C
      RETURN
      END
      subroutine gtensor(par,s,nb,b,vinv,gv,fv,grad)
C
C  compute gradient of f(par) 
C
      implicit logical (a-z)
      integer nb
      real*8 par(7),s(nb),b(6,nb),vinv(nb),gv(nb),fv(nb),grad(7)
      integer i
      real*8 D(6),th0,z,z1,res
      th0=par(1)
      call rho2D(par(2),D)
      call sihat(th0,D,b,gv,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in gv
C
      DO i=1,nb
         res=s(i)-gv(i)
         fv(i)=res*vinv(i)
      END DO
C
C     derivative with respect to theta0
C
      z=0.d0
      DO i=1,nb
         z=z+fv(i)*gv(i)
      END DO
      grad(1)=-2.d0/th0*z
C
C     derivatives with respect to r
C
      DO i=2,7
         grad(i)=0.d0
      END DO
      DO i=1,nb
         z=fv(i)*gv(i)
         z1=2.d0*b(1,i)*par(2)+b(2,i)*par(3)+b(3,i)*par(4)
         grad(2)=grad(2)-z*z1
         z1=2.d0*b(4,i)*par(3)+b(2,i)*par(2)+b(5,i)*par(4)
         grad(3)=grad(3)-z*z1
         z1=2.d0*b(6,i)*par(4)+b(3,i)*par(2)+b(5,i)*par(3)
         grad(4)=grad(4)-z*z1
         z1=2.d0*b(4,i)*par(5)+b(5,i)*par(6)
         grad(5)=grad(5)-z*z1
         z1=2.d0*b(6,i)*par(6)+b(5,i)*par(5)
         grad(6)=grad(6)-z*z1
         z1=2.d0*b(6,i)*par(7)
         grad(7)=grad(7)-z*z1         
      END DO       
C
C     We now have the gradient in grad
C      
      RETURN
      END
