      subroutine mcorrlag(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      real*8 scorr,res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      real*8 z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0.d0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if (.not.(mask(i1,i2,i3).and.mask(i1+l1,i2+l2,i3+l3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+l1,i2+l2,i3+l3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/zk
               vrmp1=y2/zk
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/zk/dsqrt(vrm)
                  k=k+1
               end if
            enddo
         enddo
      enddo
      scorr=z/k
      return
      end

      subroutine mcorr(res,mask,n1,n2,n3,nv,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,scorr(i1,i2,i3),lag)
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end
