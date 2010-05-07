      subroutine corrw(z,ngrad,m,siq,erg,th1,w)
C
C  reduce weights if sum(w) > exp(-eps*th1) 
C
      implicit logical (a-z)
      integer m,ngrad
      real*8 z(ngrad,m),siq(ngrad),erg,th1,w(m)
      integer i,j
      real*8 z1,z2,sw
      sw=w(1)
      DO i=2,m
         sw=sw+w(i)
      END DO
      z1 = exp(-1.d-2*th1)
      IF(sw.gt.z1) THEN
         DO i=1,m
            w(i)=w(i)*z1/sw
         END DO
         z2=0.d0
         DO j=1,ngrad
            z1=siq(j)
            DO i=1,m
               z1=z1-w(i)*z(j,i)
            END DO
            z2=z2+z1*z1
         END DO
         erg=sqrt(erg)
      END IF
      RETURN
      END
      
      subroutine mfunpl0(par,siq,grad,m,lpar,ngrad,z,w,work1,erg)
C  
C   only call for m>1  otherwise we have a trivial solution  c1 = 0, c2 = -log(mean(siq))
C   corresponding to an isotropic tensor and weight 1
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg
      integer i,j,i3,ind(10),mode
      real*8 c1,w(m),sw,sth,z1,p0,p1,d1,d2,d3,emc1,
     1       work1(ngrad),work2(10)
      c1 = par(1)
      emc1 = exp(-c1)
      sw = 0
      DO j=1,ngrad
         z(j,1) = emc1
      END DO
      DO i = 2,m
         i3=2*i
         p0=par(i3-2)
         p1=par(i3-1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-c1*z1*z1)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call nnls(z,ngrad,ngrad,m,siq,w,erg,work2,work1,ind,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
      call corrw(z,ngrad,m,siq,erg,c1,w)
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfunpl0g(par,siq,grad,m,lpar,ngrad,z,w,work1,erg,
     1                    fv,dfv,qv,dqv,mm1,dh)
C  
C   only call for m>1  otherwise we have a trivial solution  c1 = 0, c2 = -log(mean(siq))
C   corresponding to an isotropic tensor and weight 1
C
      implicit logical (a-z)
      integer m,mm1,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),fv(ngrad),
     1   erg,dfv(lpar,ngrad),qv(mm1,ngrad),dqv(2,mm1,ngrad),dh(lpar)
      integer i,j,i3,ind(10),mode
      real*8 c1,w(m),sw,sth,cth,spsi,cpsi,z1,p0,p1,d1,d2,d3,emc1,
     1       work1(ngrad),work2(10)
      c1 = par(1)
      emc1 = exp(-c1)
      sw = 0
      DO j=1,ngrad
         z(j,1) = emc1
      END DO
      DO i = 2,m
         i3=2*i
         p0=par(i3-2)
         p1=par(i3-1)
         sth = sin(p0)
         cth = cos(p0)
         spsi = sin(p1)
         cpsi = cos(p1)
         d1 = sth*cpsi
         d2 = sth*spsi
         d3 = cth
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-c1*z1*z1)
            qv(i-1,j) = z1*z1
            dqv(1,i-1,j) = 2.d0*(cth*(cpsi*grad(1,j)+spsi*grad(2,j))-
     1                         sth*grad(3,j))*qv(i-1,j)
            dqv(2,i-1,j) = 2.d0*sth*(cpsi*grad(2,j)-spsi*grad(1,j))*
     1                         qv(i-1,j)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call dcopy(ngrad,siq,1,fv,1)
      call nnls(z,ngrad,ngrad,m,fv,w,erg,work2,work1,ind,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
      call corrw(z,ngrad,m,siq,erg,c1,w)
C
C   now compute what's needed for the gradient
C
      DO j = 1,ngrad
C        call intpr("j1",2,j,1)
         fv(j) = 0.d0 
         fv(j) = w(1)*emc1
         dfv(1,j) = -w(i)*emc1
         DO i = 2,m
            z1 = exp(-c1*qv(i-1,j))
            fv(j) = fv(j)+w(i)*z1
            dfv(1,j) = dfv(1,j)-w(i)*qv(i-1,j)*z1
            dfv(2*i-2,j) = -w(i)*c1*dqv(1,i-1,j)*z1
            dfv(2*i-1,j) = -w(i)*c1*dqv(2,i-1,j)*z1
         END DO
      END DO
C
C   now compute the gradient in dh
C
      DO i=1,lpar
         dh(i) = 0.d0
         DO j = 1,ngrad
           dh(i) = dh(i)-dfv(i,j)*(siq(j)-fv(j))
         END DO
         dh(i) = 2.d0*dh(i)
      END DO
      RETURN
      END 
      subroutine mfunpl(par,siq,grad,m,lpar,ngrad,z,w,work1,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg
      integer i,j,i3,ind(10),mode
      real*8 c1,w(m),sw,sth,z1,p0,p1,d1,d2,d3,
     1       work1(ngrad),work2(10)
C      c1 = exp(par(1))
      c1 = par(1)
      sw = 0
      DO i = 1,m
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-c1*z1*z1)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call nnls(z,ngrad,ngrad,m,siq,w,erg,work2,work1,ind,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfunplgr(par,siq,grad,m,lpar,ngrad,z,w,work1,erg,
     1                    fv,dfv,qv,dqv,dh)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg,
     1   fv(ngrad),dfv(lpar,ngrad),qv(m,ngrad),dqv(2,m,ngrad),dh(lpar)
      integer i,j,i3,ind(10),mode
      real*8 c1,w(m),sw,sth,cth,spsi,cpsi,z1,p0,p1,d1,d2,d3,
     1       work1(ngrad),work2(10)
      c1 = par(1)
      sw = 0
      DO i = 1,m
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         cth = cos(p0)
         spsi = sin(p1)
         cpsi = cos(p1)
         d1 = sth*cpsi
         d2 = sth*spsi
         d3 = cth
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-c1*z1*z1)
            qv(i,j) = z1*z1
            dqv(1,i,j) = 2.d0*(cth*(cpsi*grad(1,j)+spsi*grad(2,j))-
     1                         sth*grad(3,j))*qv(i,j)
            dqv(2,i,j) = 2.d0*sth*(cpsi*grad(2,j)-spsi*grad(1,j))*
     1                         qv(i,j)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call dcopy(ngrad,siq,1,fv,1)
      call nnls(z,ngrad,ngrad,m,fv,w,erg,work2,work1,ind,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
C
C   now compute what's needed for the gradient
C
      DO j = 1,ngrad
C        call intpr("j1",2,j,1)
         fv(j) = 0.d0
         dfv(1,j) = 0.d0 
         DO i = 1,m
            z1 = exp(-c1*qv(i,j))
            fv(j) = fv(j)+w(i)*z1
            dfv(1,j) = dfv(1,j)-w(i)*qv(i,j)*z1
            dfv(2*i,j) = -w(i)*c1*dqv(1,i,j)*z1
            dfv(2*i+1,j) = -w(i)*c1*dqv(2,i,j)*z1
         END DO
      END DO
C
C   now compute the gradient in dh
C
      DO i=1,lpar
         dh(i) = 0.d0
         DO j = 1,ngrad
           dh(i) = dh(i)-dfv(i,j)*(siq(j)-fv(j))
         END DO
         dh(i) = 2.d0*dh(i)
      END DO
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfunpl2(par,siq,grad,m,p,lpar,ngrad,z,w,work1,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg,p
      integer i,j,i3,ind(10),mode
      real*8 c1,c2,w(m),sw,sth,z1,p0,p1,d1,d2,d3,
     1       work1(ngrad),work2(10)

      c1 = exp(par(1))
      c2 = exp(par(2))
      DO i = 1,m
         i3=2*i+1
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-p*log(1.d0+(c2+c1*z1*z1)/p))
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call nnls(z,ngrad,ngrad,m,siq,w,erg,work2,work1,ind,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
      sw = -1.d0
      DO i=1,m
         sw=sw+w(i)
      END DO
C add penalty for sum of weight .neq. 1
      erg=erg+1.d2*sw*sw
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfun(par,siq,grad,m,lpar,ngrad,z,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad),erg
      integer i,j,i3
      real*8 c1,c2,w,sw,ew,sth,z1,p0,p1,d1,d2,d3,dnrm2
      external dnrm2
      c1 = exp(par(1))
      c2 = exp(par(2))
      call dcopy(ngrad,siq,1,z,1)
      sw = 0
      DO i = 1,m
         i3=3*i
         p0=par(i3)
         p1=par(i3+1)
         IF(i.eq.m) THEN
            w = 1.d0 - sw
         ELSE
            ew = exp(par(i3+2))
            w = (1.d0 - sw)*ew/(1.d0+ew)
            sw = sw + w
         END IF
         sth = sin(p0)
C         dir(1) = sth*cos(p1)
C         dir(2) = sth*sin(p1)
C         dir(3) = cos(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j) = z(j) - w*exp(-c2 - c1*z1*z1)
         END DO
      END DO
C compute  ||z||^2
      erg = dnrm2(ngrad,z,1)
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfun0(par,siq,grad,m,lpar,ngrad,z,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad),erg
      integer i,j
      real*8 dotprod3,c1,c2,w,sw,ew,dir(3),sth,z1,z2
      c1 = exp(par(1))
      c2 = exp(par(2))
      sw = 0
      DO i = 1,m
         IF(i.eq.m) THEN
            w = 1.d0 - sw
         ELSE
            ew = exp(par(3*i+2))
            w = (1.d0 - sw)*ew/(1.d0+ew)
            sw = sw + w
         END IF
         sth = sin(par(i*3))
         dir(1) = sth*cos(par(i*3+1))
         dir(2) = sth*sin(par(i*3+1))
         dir(3) = cos(par(i*3))
         DO j = 1,ngrad
            z1 = dotprod3(dir,grad(1,j))
            z(j) = z(j) + w*exp(-c2 - c1*z1*z1)
         END DO
      END DO
C compute  ||sig-z||^2
      z1 = siq(1) - z(1)
      z2 = z1*z1
      DO j = 2,ngrad
         z1 = siq(j) - z(j)
         z2 = z2 + z1*z1
      END DO
      erg = z2
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfun2(par,siq,grad,m,p,lpar,ngrad,z,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad),p,erg
      integer i,j
      real*8 dotprod3,c1,c2,c12,w,sw,ew,dir(3),sth,z1,z2
      c1 = exp(par(1))
      c2 = exp(par(2))
      c12 = c1/(c1+c2)
      sw = 0
      DO i = 1,m
         IF(i.eq.m) THEN
            w = 1.d0 - sw
         ELSE
            ew = exp(par(3*i+2))
            w = (1.d0 - sw)*ew/(1.d0+ew)
            sw = sw + w
         END IF
         sth = sin(par(i*3))
         dir(1) = sth*cos(par(i*3+1))
         dir(2) = sth*sin(par(i*3+1))
         dir(3) = cos(par(i*3))
         DO j = 1,ngrad
            z1 = dotprod3(dir,grad(1,j))
            z(j) = z(j) + w*exp(-p*log(1 + ( c2 + c1*z1*z1 )/p))
         END DO
      END DO
      z1 = siq(1) - z(1)
      z2 = z1*z1
      DO j = 2,ngrad
         z1 = siq(j) - z(j)
         z2 = z2 + z1*z1
      END DO
      erg = z2
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine smsi(si,n,snind,nsn,w,sms)
      implicit logical (a-z)
      integer n,nsn,snind(nsn,n)
      real*8 si(n),sms(n),w(nsn)
      integer i,j
      real*8 z,sw
      sw=0.d0
      DO j=1,nsn
         sw=sw+w(j)
      END DO
      DO i=1,n
         z=0.d0
         DO j=1,nsn
            z=z+w(j)*si(snind(j,i))
         END DO
         sms(i)=z/sw
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine smsiw(si,n,snind,nsn,dgrad,sms)
      implicit logical (a-z)
      integer n,nsn,snind(nsn,n)
      real*8 si(n),sms(n),dgrad(n,n)
      integer i,j,ji
      real*8 z,sw,w
      DO i=1,n
         z=0.d0
         sw=0.d0
         DO j=1,nsn
            ji=snind(j,i)
            w=dgrad(ji,i)
            z=z+w*si(ji)
            sw=sw+w
         END DO
         sms(i)=z/sw
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsiind(si,ngrad,n1,n2,n3,snind,nsn,m,maxc,dgrad,
     1                    sms,siind,ns)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  snind  - array of indices of nearest neighbors directions
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  dgrad  - matrix of pairwise cos of gradient directions
C  sms    - vector for smoothed si
C  siind  - array of indices (output)
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,nsn,ns,snind(nsn,ngrad),
     1        siind(ns,n1,n2,n3),m
      real*8 si(ngrad,n1,n2,n3),sms(ngrad),dgrad(ngrad,ngrad),maxc
      integer i,i1,i2,i3,imin,k
      real*8 zmin,w
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
C  first smooth the si's
               call smsiw(si(1,i1,i2,i3),ngrad,snind,nsn,dgrad,sms)
C  now search for minima of sms (or weighted sms
               siind(1,i1,i2,i3)=m
               DO k=1,m
                  imin=1
                  zmin=sms(1)
                  DO i=2,ngrad
                     if(sms(i).lt.zmin) THEN
                        zmin=sms(i)
                        imin=i
                     END IF
                  END DO
                  siind(k+1,i1,i2,i3)=imin
                  IF(k.lt.m) THEN
                     DO i=1,ngrad
                        w=dgrad(i,imin)
                        if(w.gt.maxc) THEN
                           sms(i)=1d40
                        ELSE
                           sms(i)=sms(i)/(1.d0-w*w)
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getev0(si,ngrad,n1,n2,n3,lev)
      implicit logical (a-z)
      integer n1,n2,n3,ngrad
      real*8 si(ngrad,n1,n2,n3),lev(2,n1,n2,n3)
      integer i1,i2,i3,j
      real*8 z,z1
C      real*8 z,DASUM
C      external DASUM
      z=ngrad
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               lev(1,i1,i2,i3)=0
               z1=si(1,i1,i2,i3)
               DO j=2,ngrad
                  z1=z1+si(j,i1,i2,i3)
               END DO
               lev(2,i1,i2,i3)=-log(z1/z)
C               lev(2,i1,i2,i3)=-log(DASUM(ngrad,si(1,i1,i2,i3),1)/z)
            END DO
         END DO
      END DO
      RETURN
      END
 