      subroutine initdata(si,n1,n2,n3,nb,maxvalue)
C
C   set si()==0 for all voxel that have si-values <=0 or > maxvalue
C
      integer n1,n2,n3,nb,maxvalue,si(n1,n2,n3,nb)
      logical mask
      integer i1,i2,i3,k,sii
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               mask=.FALSE.
               DO k=1,nb
                  sii=si(i1,i2,i3,k)
                  if(sii.le.0.or.sii.gt.maxvalue) THEN
                     mask=.TRUE.
                     CYCLE
                  END IF
               END DO
               IF(mask) THEN
                  DO k=1,nb
                     si(i1,i2,i3,k)=0.d0
                  END DO
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      real*8 scorr,res(nv,n1,n2,n3),sigma(n1,n2,n3),mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      real*8 vrm,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k,j1,j2,j3
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         j1=i1+l1
         do i2=1,n2-l2
            j2=i2+l2
            do i3=1,n3-l3
               j3=i3+l3
               if (.not.(mask(i1,i2,i3).and.mask(j1,j2,j3))) CYCLE
               zcorr=0.d0
               do i4=1,nv
                  zcorr=zcorr+(res(i4,i1,i2,i3)-mean(i1,i2,i3))*
     1                        (res(i4,j1,j2,j3)-mean(j1,j2,j3))
               enddo
               vrm=sigma(i1,i2,i3)*sigma(j1,j2,j3)
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/vrm
                  k=k+1
               end if
            enddo
         enddo
      enddo
      if( k.gt.0 ) then
         scorr=z/k/nv
      ELSE
         scorr=0
      END IF
      return
      end

      subroutine msd(res,mask,n1,n2,n3,nv,sigma,mean)
      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 sigma(n1,n2,n3),res(nv,n1,n2,n3),mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,iv
      real*8 z,resi,zm
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
                  z=0.d0
                  zm=0.d0
                  DO iv=1,nv
                     resi=res(iv,i1,i2,i3)
                     zm=zm+resi
                     z=z+resi*resi
                  END DO
                  zm=zm/nv
                  z=z/nv-zm*zm
                  sigma(i1,i2,i3)=sqrt(z)
                  mean(i1,i2,i3)=zm
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
      
      subroutine mcorr(res,mask,n1,n2,n3,nv,sigma,mean,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),res(nv,n1,n2,n3),sigma(n1,n2,n3),
     1       mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      call msd(res,mask,n1,n2,n3,nv,sigma,mean)
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,
     1                       scorr(i1,i2,i3),lag)
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end


