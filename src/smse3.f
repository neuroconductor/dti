      subroutine paramw3(h,vext,ind,w,n)
C  compute a description of local weights
C  h    - bandwidth
C  vext - vector (length 2) of relative voxel extensions
C  ind  - integer array dim (3,n) containing relative indices in xyz
C  w    - vector of corresponding weights
C  n    - number of positive weights (initial value
C         (2 int(h)+1)*(2 int(h/vext(1))+1)*(2 int(h/vext(2))+1)
      integer n,ind(3,n)
      double precision h,vext(2),w(n)
      integer i,i1,i2,i3,ih1,ih2,ih3
      double precision hsq,z1,z2,z3
      hsq=h*h
      ih1 = int(h)
      ih2 = int(h/vext(1))
      ih3 = int(h/vext(2))
      i=1
      DO i1=-ih1,ih1
         z1=i1*i1
         DO i2=-ih2,ih2
            z2=i2*vext(1)
            z2=z1+z2*z2
            IF(z2.ge.hsq) CYCLE
            DO i3=-ih3,ih3
               z3=i3*vext(2)
               z3=z2+z3*z3
               IF(z3.ge.hsq) CYCLE
               ind(1,i)=i1
               ind(2,i)=i2
               ind(3,i)=i3
               w(i)=1.d0-z3/hsq
               i=i+1
            END DO
         END DO
      END DO
      n=i-1
      RETURN
      END
