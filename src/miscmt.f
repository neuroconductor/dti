      subroutine inhmtfun(par,si,grad,bv,npar,nc,n,rho,evc,w,krit)
C
C   need  npar=3*(nc+1)
C
      implicit logical (a-z)
      integer n,npar,nc
      real*8 par(npar),si(n),grad(3,n),bv(n),rho,evc(3,nc),w(nc),krit
      integer i,k,l
      real*8 w0,l1,l12,l2,ksi,eta,sksi,sw,z,z1,bvi
      l12=par(1)
      l2=par(2)
      l1=l12+l2
      w0=par(3)
      DO k=1,nc
         ksi = par(k*3+1)
         eta = par(k*3+2)
         w(k) = par(k*3+3)
         sksi = sin(ksi)
         evc(1,k) = sksi*cos(eta)
         evc(2,k) = sksi*sin(eta)
         evc(3,k) = cos(ksi)
      END DO
      sw = 0
      IF(nc.gt.1.and.rho.gt.0) THEN
         DO k=2,nc
            DO l=1,k-1
               z=evc(1,k)*evc(1,l)+evc(2,k)*evc(2,l)+evc(3,k)*evc(3,l)
               sw = sw + z*z
            END DO
         END DO
         sw = sw*rho
      END IF
      DO i=1,n
         bvi=bv(i)
         z = si(i)-w0*exp(-bvi*l1)
         IF(nc.gt.0) THEN
         DO k=1,nc
           z1=evc(1,k)*grad(1,i)+evc(2,k)*grad(2,i)+evc(3,k)*grad(3,i)
           z = z-w(k)*exp(-bvi*(l2+l12*z1*z1))
         END DO
         END IF
         sw = sw+z*z
      END DO
      krit=sw
      RETURN
      END
      subroutine inhmtgrd(par,si,grad,bv,npar,nc,n,rho,evc,w,
     1                    f,evg,el1,el2k,grd)
C
C   need  npar=3*(nc+1)
C
      implicit logical (a-z)
      integer n,npar,nc
      real*8 par(npar),si(n),grad(3,n),bv(n),rho,evc(3,nc),w(nc),
     1       grd(npar),f(n),el1(n),el2k(nc,n),evg(nc,n)
      integer i,k,l
      real*8 w0,l1,l12,l2,ksi,eta,sksi,sw,z,z1,z2,z3,z1i,z2i,bvi,
     1       e1k,e2k,e3k,e1e,e2e,dz2,dz3
      l12=par(1)
      l2=par(2)
      l1=l12+l2
      w0=par(3)
      DO k=1,nc
         ksi = par(k*3+1)
         eta = par(k*3+2)
         w(k) = par(k*3+3)
         sksi = sin(ksi)
         evc(1,k) = sksi*cos(eta)
         evc(2,k) = sksi*sin(eta)
         evc(3,k) = cos(ksi)
      END DO
      sw = 0
      IF(nc.gt.1.and.rho.gt.0) THEN
         DO k=2,nc
            DO l=1,k-1
               z=evc(1,k)*evc(1,l)+evc(2,k)*evc(2,l)+evc(3,k)*evc(3,l)
               sw = sw + z*z
            END DO
         END DO
         sw = sw*rho
      END IF
      DO i=1,n
         bvi=bv(i)
         el1(i)=exp(-bvi*l1)
         z = si(i)-w0*el1(i)
         IF(nc.gt.0) THEN
         DO k=1,nc
           z1=evc(1,k)*grad(1,i)+evc(2,k)*grad(2,i)+evc(3,k)*grad(3,i)
           evg(k,i)=z1
           el2k(k,i)=exp(-bvi*(l2+l12*z1*z1))
           z = z-w(k)*el2k(k,i)
         END DO
         END IF
         f(i) = z
      END DO
C
C  gradient calculations
C  df/dl12, df/dl2, df/dw0
C
      z1=0.d0
      z2=0.d0
      z3=0.d0
      DO i=1,n
         z1i=w0*el1(i)
         z2i=z1i
         z3=z3-f(i)*el1(i)
         IF(nc.gt.0) THEN
            DO k=1,nc
               z1i=z1i+w(k)*evg(k,i)*evg(k,i)*el2k(k,i)
               z2i=z2i+w(k)*el2k(k,i)
            END DO
         END IF
         z1=z1+bv(i)*f(i)*z1i
         z2=z2+bv(i)*f(i)*z2i
      END DO
      grd(1)=2.d0*z1
      grd(2)=2.d0*z2
      grd(3)=2.d0*z3
C
C  gradient calculations
C  df/wk, df/ksik, df/etak
C
      IF(nc.gt.0) THEN
         DO k=1,nc
            z1=0.d0
            z2=0.d0
            z3=0.d0
            ksi = par(k*3+1)
            eta = par(k*3+2)
            e1k=cos(ksi)*cos(eta)
            e2k=cos(ksi)*sin(eta)
            e3k=-sin(ksi)
            e1e=-sin(ksi)*sin(eta)
            e2e=sin(ksi)*cos(eta)
            DO i=1,n
               z1=z1-f(i)*el2k(k,i)
               dz2=e1k*grad(1,i)+e2k*grad(2,i)+e3k*grad(3,i)
               dz3=e1e*grad(1,i)+e2e*grad(2,i)
               z2=z2+f(i)*bv(i)*dz2*evg(k,i)*el2k(k,i)
               z3=z3+f(i)*bv(i)*dz3*evg(k,i)*el2k(k,i)
            END DO
            grd(3*k+3) = 2.d0*z1
            grd(3*k+1) = 4.d0*z2*w(k)*l12
            grd(3*k+2) = 4.d0*z3*w(k)*l12
         END DO
      END IF
C
C  still derivatives for penalty term missing
C
      IF(nc.gt.1) THEN
         DO k=1,nc
            ksi = par(k*3+1)
            eta = par(k*3+2)
            e1k=cos(ksi)*cos(eta)
            e2k=cos(ksi)*sin(eta)
            e3k=-sin(ksi)
            e1e=-sin(ksi)*sin(eta)
            e2e=sin(ksi)*cos(eta)
            DO l=1,nc
               IF(k.eq.l) CYCLE
               z=evc(1,k)*evc(1,l)+evc(2,k)*evc(2,l)+evc(3,k)*evc(3,l)
               grd(3*k+1) = grd(3*k+1)+
     1              2.0*rho*z*(e1k*evc(1,l)+e2k*evc(2,l)+e3k*evc(3,l))
               grd(3*k+2) = grd(3*k+2)+
     1              2.0*rho*z*(e1e*evc(1,l)+e2e*evc(2,l))
            END DO
         END DO
      END IF
      RETURN
      END


      subroutine mthomog(orient,mix,order,n1,n2,n3,mo,mask,minw,
     1                   maxangle,vext,norient)
      implicit logical(a-z)
      integer n1,n2,n3,mo,order(n1,n2,n3)
      real*8 mix(mo,n1,n2,n3),orient(3,mo,n1,n2,n3),minw,maxangle,
     1       vext(3),norient(3,mo,n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j,k,l,oi,ik1(26),ik2(26),ik3(26),j1,j2,j3
      real*8 sw,d(3,26),ddot,dj(3),sd(3),mcos
      external ddot
      mcos = cos(maxangle)
C first remove spurious effects
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               oi=order(i1,i2,i3)
               if(oi.le.1) CYCLE
               DO k=oi,2,-1
                  if(mix(k,i1,i2,i3).lt.minw) THEN
                     order(i1,i2,i3)=k-1
                     sw=1.d0-mix(k,i1,i2,i3)
                     DO l=1,k-1
                        mix(l,i1,i2,i3)=mix(l,i1,i2,i3)/sw
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
C we now only have weights that are larger than minw
C  get orientations of neighboring voxel
      sw = sqrt(vext(1)*vext(1)+vext(2)*vext(2)+vext(3)*vext(3))
      k=1
      DO i1=-1,1
         DO i2=-1,1
            DO i3=-1,1
               if(i1.eq.0.and.i2.eq.0.and.i3.eq.0) CYCLE
               d(1,k) = i1*vext(1)
               d(2,k) = i2*vext(2)
               d(3,k) = i3*vext(3)
               ik1(k) = i1
               ik2(k) = i2
               ik3(k) = i3
               k=k+1
            END DO
         END DO
      END DO
      DO k=1,26
         sw = sqrt(ddot(3,d(1,k),1,d(1,k),1))
         d(1,k) = d(1,k)/sw
         d(2,k) = d(2,k)/sw
         d(3,k) = d(3,k)/sw
      END DO
C now homogeneize directions
      DO i1=2,n1-1
         DO i2=2,n2-1
            DO i3=2,n3-1
               if(.not.mask(i1,i2,i3)) CYCLE
               oi=order(i1,i2,i3)
               DO j=1,oi
                  call dcopy(3,orient(1,j,i1,i2,i3),1,dj,1)
                  call dcopy(3,dj,1,sd,1)
                  DO k=1,26
                     if(abs(ddot(3,dj,1,d(1,k),1)).lt.mcos) CYCLE
                     j1 = i1+ik1(k)
                     j2 = i2+ik2(k)
                     j3 = i3+ik3(k)
                     if(.not.mask(j1,j2,j3)) CYCLE
                     DO l=1,order(j1,j2,j3)
                        sw=ddot(3,dj,1,orient(1,l,j1,j2,j3),1)
                        if(abs(sw).lt.mcos) CYCLE
                        IF(sw.gt.0) THEN
                           sw = sqrt((abs(sw)-mcos)/(1.d0-mcos))
                        ELSE
                            sw = -sqrt((abs(sw)-mcos)/(1.d0-mcos))
                        END IF
                        call daxpy(3,sw,orient(1,l,j1,j2,j3),1,sd,1)
                     END DO
                  END DO
                  sw=sqrt(ddot(3,sd,1,sd,1))
                  norient(1,j,i1,i2,i3)=sd(1)/sw
                  norient(2,j,i1,i2,i3)=sd(2)/sw
                  norient(3,j,i1,i2,i3)=sd(3)/sw
               END DO
C    filtered directions in orient
            END DO
         END DO
      END DO
      RETURN
      END