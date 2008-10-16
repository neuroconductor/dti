      subroutine besselq(x,n,fw)
      implicit logical (a-z) 
      integer n
      real*8 x(n),fw(n)
      external besseli
      real*8 besseli
      integer i
      DO i=1,n
         x(i)=i*1.d-2
         fw(i)=besseli(x(i),1.d0,2.d0)/besseli(x(i),0.d0,2.d0)
      END DO
      RETURN
      END
      subroutine ricecor1(si,w,n,sw,th,sigma2,fw)
      implicit logical (a-z) 
      integer n,si(n)
      real*8 w(n),th,sigma2,sw,fw(10000)
      integer i,j
      real*8 z,sth
      sth=0.d0
      z=th/sigma2/1.d-2
      DO i=1,n
         j=si(i)*z+1
         if(j.gt.10000) THEN
            sth=sth+si(i)*w(i)
         ELSE
            sth=sth+fw(j)*si(i)*w(i)
         END IF
      END DO
      th=sth/sw
      RETURN
      END
