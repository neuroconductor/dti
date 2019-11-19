      subroutine exppm6(p,ex)
      implicit none
      double precision p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=dsin(p)
      ex(1,3)=0.d0
      ex(2,1)=-dsin(p)
      ex(2,2)=dcos(p)
      ex(2,3)=0.d0
      ex(3,1)=0.d0
      ex(3,2)=0.d0
      ex(3,3)=1.d0
      RETURN
      END
      subroutine exppm5(p,ex)
      implicit none
      double precision p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=0.d0
      ex(1,3)=-dsin(p)
      ex(2,1)=0.d0
      ex(2,2)=1.d0
      ex(2,3)=0.d0
      ex(3,1)=dsin(p)
      ex(3,2)=0.d0
      ex(3,3)=dcos(p)
      RETURN
      END
      subroutine exppm4(p,b,ex)
      implicit none
      double precision b,p,ex(3,3)
      double precision D,sqr2,sb,pDs,D2,cpds,spds,spdsh
      sqr2 = dsqrt(2.d0)
      sb = dsin(b)
      D = dsqrt(3.d0-dcos(2.d0*b))
      pDs = p*D/sqr2
      cpds = dcos(pDs)
      spds = dsin(pDs)
      spdsh = dsin(pDs/2.d0)
      D2 = D*D
      ex(1,1)=2.d0*(1.d0+cpDs*sb*sb)/D2
      ex(1,2)=-sqr2/D*sb*spDs
      ex(1,3)=-4.d0/D2*sb*spdsh*spdsh
      ex(2,1)=-ex(1,2)
      ex(2,2)=cpDs
      ex(2,3)=sqr2/D*spDs
      ex(3,1)=ex(1,3)
      ex(3,2)=-ex(2,3)
      ex(3,3)=1.d0-2.d0/D2*(1.d0-cpDs)
      RETURN
      END
      subroutine k456krb(par,b,matm,erg)
C
C  Solve exponential equation for dicrepance parameters
C  compute ||\prod_{i=4}^6 exp(par[i] m_i) - matm||^2
C
      implicit none
      double precision par(3),b,matm(3,3),erg
      integer i1,i2
      double precision s,z,em4(3,3),em5(3,3),em6(3,3),am4(3,3),am5(3,3)
      call exppm4(par(1),b,em4)
      call exppm5(par(2),em5)
      call exppm6(par(3),em6)
      DO i1=1,3
         DO i2=1,3
            am5(i1,i2)=em5(i1,1)*em6(1,i2)+em5(i1,2)*em6(2,i2)+
     1                                     em5(i1,3)*em6(3,i2)
            END DO
         END DO
      DO i1=1,3
         DO i2=1,3
            am4(i1,i2)=em4(i1,1)*am5(1,i2)+em4(i1,2)*am5(2,i2)+
     1                                     em4(i1,3)*am5(3,i2)
         END DO
      END DO
      s=0.d0
      DO i1=1,3
         DO i2=1,3
            z=matm(i1,i2)-am4(i1,i2)
            s=s+z*z
         END DO
      END DO
      erg=s
      RETURN
      END
      subroutine abofg(g,n,bg)
C
C  compute spherical coordinates for gradient vectors
C
      implicit none
      integer n
      double precision g(3,n),bg(2,n)
      integer i
      double precision z,g1,g2,g3,beta,gamma,onemeps
      onemeps=1.d0-1d-10
      DO i=1,n
         g1=g(1,i)
         g2=g(2,i)
         g3=g(3,i)
C standardize to length 1
         z=g1*g1+g2*g2+g3*g3
         z=sqrt(z)
         g1=g1/z
         g2=g2/z
         g3=g3/z
         beta=asin(g1)
         if(abs(g1).lt.onemeps) THEN
            z=g3/cos(beta)
            if(abs(z).lt.onemeps) THEN
               gamma=acos(z)
            ELSE
               gamma=1.570796327d0-sign(1.570796327d0,z)
            END IF
         ELSE
            gamma=0.d0
         END IF
         if(g2.gt.0.d0) gamma = -gamma
         bg(1,i)=beta
         bg(2,i)=gamma
      END DO
      RETURN
      END
      subroutine bgstats(g,n,bg,bghat)
      implicit none
      integer n
      double precision g(3,n),bg(2,n),bghat(2,n,n)
      integer i1,i2
      double precision dgamma,cb1,sb1,cb2,sb2,betah,cbh,z,gammah,cdg
C   first get sperical coordinates in bg
      call abofg(g,n,bg)
      DO i1=1,n
         sb1=sin(bg(1,i1))
         cb1=cos(bg(1,i1))
C    die normalen-vektoren der Gradienten
         DO i2=1,n
            dgamma=bg(2,i1)-bg(2,i2)
            cdg=cos(dgamma)
            if(abs(cdg).gt.1-1d-8) THEN
               bghat(1,i1,i2)=asin(sin(bg(1,i1)-
     1                             sign(1.d0,cdg)*bg(1,i2)))
               bghat(2,i1,i2)=0.d0
               CYCLE
            END IF
            sb2=sin(bg(1,i2))
            cb2=cos(bg(1,i2))
            z=sb1*cb2-cb1*sb2*cdg
            betah=asin(z)
            cbh=cos(betah)
            IF(abs(cbh).gt.1d-8) THEN
               z=cb1*sin(dgamma)/cbh
               z=sign(min(1.d0,abs(z)),z)
               gammah=asin(z)
            ELSE
               if(abs(cb1).gt.1d-6) THEN
C   this should not happen
                  call dblepr("cb1",3,cb1,1)
                  call dblepr("cbh",3,cbh,1)
               END IF
               gammah = dgamma*sign(1d0,cb1*cbh)
            END IF
            bghat(1,i1,i2)=betah
            bghat(2,i1,i2)=gammah
C    die spherischen Koordinaten der Gradientenpaare (Parameter der Rotationsmatix)
         END DO
      END DO
      RETURN
      END
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
