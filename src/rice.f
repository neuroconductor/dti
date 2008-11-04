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
      subroutine ricecorr(si,w,n,nb,niter,sw,th,s2,sigma2,fw)
      implicit logical (a-z) 
      integer n,nb,si(nb,n),niter(nb)
      real*8 w(n),th(nb),s2(nb),sigma2,sw,fw(10000)
      integer i,j,k,l,iter,sii
      real*8 z,sth,th2,ss2,thk
      iter=1
      DO k=1,nb
         iter=max(iter,niter(k))
      END DO
      DO l=1,iter
         Do k=1,nb
            if(l.le.niter(k)) THEN
               sth=0.d0
               ss2=0.d0
               thk=th(k)
               z=thk/sigma2/1.d-2
               th2=thk*thk
               DO i=1,n
                  sii=si(k,i)
                  j=sii*z+1
                  if(j.gt.10000) THEN
                     sth=sth+sii*w(i)
                     ss2=ss2+w(i)*((sii*sii+th2)*0.5d0-sii*thk)
                  ELSE
                     sth=sth+fw(j)*sii*w(i)
                     ss2=ss2+w(i)*((sii*sii+th2)*0.5d0-fw(j)*sii*thk)
                  END IF
               END DO
               th(k)=sth/sw
               s2(k)=ss2/sw
            END IF
         END DO
         ss2=0.d0
         DO k=1,nb
            ss2=ss2+s2(k)
         END DO
         sigma2=ss2/nb
      END DO
      RETURN
      END
