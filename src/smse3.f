      subroutine ghfse3i(i4,kstar,k456,nbg,nbghat,ng,kexp,
     1                    kappa0,vext,h,kappa,varred,varred0,n)
      implicit logical (a-z)
      integer ng,i4,kstar,n
      real*8 k456(3,ng,ng),nbg(3,3,ng),nbghat(3,3,ng,ng),vext(2),
     1       kexp,h(kstar),kappa(kstar),varred(kstar),varred0(kstar)
      logical getkappa
      integer k,n0
      real*8 kappa0,hakt,hakt0,vr,vr0,ch,lch,chk,vred,v0r0,v0r
      ch=1.25d0
      lch=-log(ch)
      chk=ch
      getkappa=.FALSE.
      hakt=1.d0
C   initialize kappa
C   loop over steps
      DO k=1,kstar
         call lkfse3i0(hakt,kappa0,i4,k456,nbg,nbghat,ng,
     1                   vext,vr,vr0,n)
         if(k.gt.1.and.log(vr0/vr).lt.lch*kexp*k) getkappa=.TRUE.
C  search for new kappa0 if needed
         DO WHILE (getkappa)
            kappa0=kappa0*0.98d0
            call lkfse3i0(hakt,kappa0,i4,k456,nbg,nbghat,ng,
     1                   vext,vr,vr0,n)
            if(log(vr0/vr).ge.lch*kexp*k) getkappa=.FALSE.
         END DO
C  search for new hakt   
         vred=vr/chk
         DO WHILE (vred.lt.1) 
            hakt=hakt*1.05
            call lkfse3i0(hakt,kappa0,i4,k456,nbg,nbghat,ng,
     1                   vext,vr,vr0,n)
            if(log(vr0/vr).lt.lch*kexp*k) getkappa=.TRUE.
C  search for new kappa0 if needed
            DO WHILE (getkappa)
               kappa0=kappa0*0.98d0
            call lkfse3i0(hakt,kappa0,i4,k456,nbg,nbghat,ng,
     1                   vext,vr,vr0,n)
               if(log(vr0/vr).ge.lch*kexp*k) getkappa=.FALSE.
            END DO
            vred=vr/chk
         END DO
         DO WHILE (vred.gt.1.01) 
            hakt0=hakt
            v0r=vr
            v0r0=vr0
            n0=n
            hakt=hakt/1.005
            call lkfse3i0(hakt,kappa0,i4,k456,nbg,nbghat,ng,
     1                   vext,vr,vr0,n)
            vred=vr/chk
            If (vred.lt.1) THEN
               hakt=hakt0
               vr=v0r
               vr0=v0r0
               n=n0
            END IF
         END DO
         h(k) = hakt
         kappa(k) = kappa0
         varred(k) = vr
         varred0(k) = vr0
         chk=chk*ch
C  number of positive weights for last step in n
      END DO
      RETURN
      END
      subroutine gethse3i(i4,kstar,grad,gr2,gr3,ngrad,kexp,vext,
     1                    h,kappa,varred,varred0,n)
      implicit logical (a-z)
      integer ngrad,i4,kstar,n
      real*8 grad(3,ngrad),gr2(3,ngrad),gr3(2,ngrad),vext(2),
     1       kexp,h(kstar),kappa(kstar),varred(kstar),varred0(kstar)
      logical getkappa
      integer k,n0
      real*8 kappa0,hakt,hakt0,vr,vr0,ch,lch,chk,vred,v0r0,v0r
      ch=1.25d0
      lch=-log(ch)
      chk=ch
      getkappa=.FALSE.
      hakt=1.d0
C   initialize kappa
      call minang(i4,grad,ngrad,kappa0)
      kappa0 = max(12.57d0/ngrad,kappa0)
      call grad23(grad,ngrad,gr2,gr3)
C   loop over steps
      DO k=1,kstar
         call lkse3i0(hakt,kappa0,i4,grad,gr2,gr3,ngrad,vext,vr,vr0,n)
         if(k.gt.1.and.log(vr0/vr).lt.lch*kexp*k) getkappa=.TRUE.
C  search for new kappa0 if needed
         DO WHILE (getkappa)
            kappa0=kappa0*0.98d0
         call lkse3i0(hakt,kappa0,i4,grad,gr2,gr3,ngrad,vext,vr,vr0,n)
            if(log(vr0/vr).ge.lch*kexp*k) getkappa=.FALSE.
         END DO
C  search for new hakt   
         vred=vr/chk
         DO WHILE (vred.lt.1) 
            hakt=hakt*1.05
         call lkse3i0(hakt,kappa0,i4,grad,gr2,gr3,ngrad,vext,vr,vr0,n)
            if(log(vr0/vr).lt.lch*kexp*k) getkappa=.TRUE.
C  search for new kappa0 if needed
            DO WHILE (getkappa)
               kappa0=kappa0*0.98d0
         call lkse3i0(hakt,kappa0,i4,grad,gr2,gr3,ngrad,vext,vr,vr0,n)
               if(log(vr0/vr).ge.lch*kexp*k) getkappa=.FALSE.
            END DO
            vred=vr/chk
         END DO
         DO WHILE (vred.gt.1.01) 
            hakt0=hakt
            v0r=vr
            v0r0=vr0
            n0=n
            hakt=hakt/1.005
         call lkse3i0(hakt,kappa0,i4,grad,gr2,gr3,ngrad,vext,vr,vr0,n)
            vred=vr/chk
            If (vred.lt.1) THEN
               hakt=hakt0
               vr=v0r
               vr0=v0r0
               n=n0
            END IF
         END DO
         h(k) = hakt
         kappa(k) = kappa0
         varred(k) = vr
         varred0(k) = vr0
         chk=chk*ch
C  number of positive weights for last step in n
      END DO
      RETURN
      END
      subroutine minang(i4,grad,ngrad,angle)
C  compute nontrivial min of first discrepancy term
      implicit logical (a-z)
      integer ngrad,i4
      real*8 grad(3,ngrad),angle
      integer j4
      real*8 a12,a13,g1,g2,g3,f1,f2,f3,e2,e3
      angle = 2.d0
      DO j4=1,ngrad
         g1=grad(1,j4)
         g2=grad(2,j4)
         g3=grad(3,j4)
         if(grad(1,i4)*g1+grad(2,i4)*g2+grad(3,i4)*g3.gt..999e0) CYCLE
         if(g1.lt.1.d0-1d-8) THEN
            f1=sqrt(1-g1*g1)
            f2=-g1*g2/f1
            f3=-g1*g3/f1
            e3=-g2/f1
            e2=g3/f1
         ELSE
            f1=0.d0
            f2=-1.d0
            f3=0.d0
            e2=0.d0
            e3=1.d0
C  f2 should be proportional to g2, f3 to g3 
         END IF
         a12 = grad(1,i4)*f1+grad(2,i4)*f2+grad(3,i4)*f3
         a13 = grad(2,i4)*e2+grad(3,i4)*e3
         angle = min(angle,a12*a12*a12*a12+a13*a13*a13*a13)
      END DO
      angle=sqrt(sqrt(angle))
      RETURN
      END
      subroutine grad23(grad,ngrad,grad2,grad3)
      implicit logical (a-z)
      integer ngrad
      real*8 grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad)
      integer j4
      real*8 g1,g2,g3,f1,f2,f3,e2,e3
      DO j4=1,ngrad
         g1=grad(1,j4)
         g2=grad(2,j4)
         g3=grad(3,j4)
         if(g1.lt.1.d0-1d-8) THEN
            f1=sqrt(1-g1*g1)
            f2=-g1*g2/f1
            f3=-g1*g3/f1
            e3=-g2/f1
            e2=g3/f1
         ELSE
            f1=0.d0
            f2=-1.d0
            f3=0.d0
            e2=0.d0
            e3=1.d0
C  f2 should be proportional to g2, f3 to g3 
         END IF
         grad2(1,j4)=f1
         grad2(2,j4)=f2
         grad2(3,j4)=f3
C grad contains only non-zero components, first component is 0
         grad3(1,j4)=e2
         grad3(2,j4)=e3
      END DO
      RETURN
      END
      subroutine ng123(beta,gamma,g123)
      implicit logical (a-z)
      real*8 beta,gamma,g123(3,3)
      real*8 cb,sb,cg,sg
      sb = sin(beta)
      sg = sin(gamma)
      cb = cos(beta)
      cg = cos(gamma)
      g123(1,1) = sb
      g123(2,1) = -cb*sg
      g123(3,1) = cb*cg
      g123(1,2) = cb
      g123(2,2) = sb*sg
      g123(3,2) = -sb*cg
      g123(1,3) = 0.d0
      g123(2,3) = cg
      g123(3,3) = sg
      RETURN
      END
      subroutine expm3(m,ex)
C
C   Compute matrix exponential of 
C
      implicit logical (a-z)
      real*8 m(3,3),ex(3,3)
      integer i1,i2,j
      real*8 mpot(3,3),mpotn(3,3),maxmpot,jfac,z
      DO i1=1,3
         mpot(i1,i1) = 1.d0
         ex(i1,i1)   = 1.d0
         IF(i1.eq.3) CYCLE
         DO i2=i1+1,3
            mpot(i1,i2) = 0.d0
            ex(i1,i2)   = 0.d0
            mpot(i2,i1) = 0.d0
            ex(i2,i1)   = 0.d0
         END DO
      END DO
      maxmpot = 1.d0
      j = 1
      jfac = 1.d0
      DO While (maxmpot.gt.1.d-12)
         jfac = jfac*j
         DO i1=1,3
            DO i2=1,3
               mpotn(i1,i2)=mpot(i1,1)*m(1,i2)+mpot(i1,2)*m(2,i2)+
     1                      mpot(i1,3)*m(3,i2)
            END DO
         END DO
         maxmpot=0.d0
         DO i1=1,3
            DO i2=1,3
               mpot(i1,i2)=mpotn(i1,i2)
               z = mpot(i1,i2)/jfac
               ex(i1,i2)=ex(i1,i2)+z
               z = abs(z)
               IF(z.gt.1d10) EXIT
               maxmpot = max(maxmpot,z)
            END DO
         END DO
         j = j+1
         IF(j.gt.50) EXIT
      END DO
      RETURN
      END
      subroutine k456krit(par,matm,m4,m5,m6,erg)
      implicit logical (a-z)
      real*8 par(3),matm(3,3),m4(3,3),m5(3,3),m6(3,3),erg
      integer i1,i2
      real*8 s,z,am4(3,3),am5(3,3),am6(3,3),em4(3,3),em5(3,3),em6(3,3)
      DO i1=1,3
         DO i2=1,3
            am4(i1,i2)=par(1)*m4(i1,i2)
            am5(i1,i2)=par(2)*m5(i1,i2)
            am6(i1,i2)=par(3)*m6(i1,i2)
         END DO
      END DO
      call expm3(am4,em4)
      call expm3(am5,em5)
      call expm3(am6,em6)
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
      implicit logical (a-z)
      integer n
      real*8 g(3,n),bg(2,n)
      integer i
      real*8 z,g1,g2,g3,beta,gamma,onemeps
      onemeps=1.d0-1d-10
      DO i=1,n
         g1=g(1,i)
         g2=g(2,i)
         g3=g(3,i)
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
      subroutine bgstats(g,n,bg,bghat,nbg,nbghat)
      implicit logical (a-z)
      integer n
      real*8 g(3,n),bg(2,n),bghat(2,n,n),nbg(3,3,n),nbghat(3,3,n,n)
      integer i1,i2
      real*8 dgamma,cb1,sb1,cb2,sb2,betah,cbh,z,gammah
C   first get sperical coordinates in bg
      call abofg(g,n,bg)
      DO i1=1,n
         sb1=sin(bg(1,i1))
         cb1=cos(bg(1,i1))
         call ng123(bg(1,i1),bg(2,i1),nbg(1,1,i1))
         DO i2=1,n
            sb2=sin(bg(1,i2))
            cb2=cos(bg(1,i2))
            dgamma=bg(2,i1)-bg(2,i2)
            z=sb1*cb2-cb1*sb2*cos(dgamma)
            betah=asin(z)
            cbh=cos(betah)
            IF(abs(cbh).gt.1d-8) THEN
               z=cb1*sin(dgamma)/cbh
               z=sign(min(1.d0,abs(z)),z)
               gammah=asin(z)
            ELSE
               gammah=1.570796327d0
            END IF
            bghat(1,i1,i2)=betah
            bghat(2,i1,i2)=gammah
            call ng123(betah,gammah,nbghat(1,1,i1,i2))
         END DO
      END DO  
      RETURN
      END
      subroutine lkse3i0(h,kappa,i4,grad,grad2,grad3,ngrad,vext,
     1                   vred,vred0,n)
      implicit logical (a-z)
      integer ngrad,i4
      real*8 h,kappa,grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad),
     1       vext(2),vred,vred0
      integer ih1,ih2,ih3,i,j1,j2,j3,j4,mj1,mj2,mj3,n
      real*8 v1a,v2a,v3a,w1a,w2a,w3a,d2a,d3a,
     1       d1,d2,pD4,pD4a,kap4,kap2,h2,vd2,vd3,x1,x2,x3,n12,n13,nx1,
     2       h4,gi1,gi2,gi3,gj1,gj2,gj3,sw,sw2,sw0,sw20,wght,anz,rad
      ih1 = 2.0d0*h
      ih2 = 2.0d0*h/vext(1)
      ih3 = 2.0d0*h/vext(2)
      h2 = h*h
      h4 = h2*h2
      kap2 = kappa*kappa
      kap4 = kap2*kap2
      vd2 = vext(1)
      vd3 = vext(2)
      sw=0.d0
      sw0=0.d0
      sw2=0.d0
      sw20=0.d0
      i = 1
      mj1=0
      mj2=0
      mj3=0
      gi1 = grad(1,i4)
      gi2 = grad(2,i4)
      gi3 = grad(3,i4)
      n = 0
      DO j4 = 1,ngrad
         gj1 = grad(1,j4)
         gj2 = grad(2,j4)
         gj3 = grad(3,j4)
         n13 = gi2*grad3(1,j4)+gi3*grad3(2,j4)
         d1 = n13*n13
         pD4 = d1*d1/kap4
         if(pD4.ge.h4) CYCLE  
         n12 = gi1*grad2(1,j4)+gi2*grad2(2,j4)+gi3*grad2(3,j4)
         d1 = n12*n12
         pD4 = pD4 + d1*d1/kap4
         if(pD4.ge.h4) CYCLE  
         w1a=0.5d0*n13*gj1
         w2a=grad3(1,j4)+0.5d0*n13*gj2
         w3a=grad3(2,j4)+0.5d0*n13*gj3
         v1a=grad2(1,j4)+0.5d0*n12*gj1
         v2a=grad2(2,j4)+0.5d0*n12*gj2
         v3a=grad2(3,j4)+0.5d0*n12*gj3
         DO j1 = 0,ih1
            x1 = j1
            DO j2 = -ih2,ih2
               x2 = vd2*j2
               DO j3 = -ih3,ih3
                  x3 = vd3*j3
                  nx1 = x1*gj1+x2*gj2+x3*gj3
                  d2 = nx1*nx1
                  pD4a = pD4 + d2*d2
                  if(pD4a.ge.h4) CYCLE  
C   second term  <\delta_x,n_1(2)> to large
C
C   now third term if something is left
C
                  d3a=x1*w1a+x2*w2a+x3*w3a
                  pD4a = pD4a + d3a*d3a/kap2
                  if(pD4a.ge.h4) CYCLE  
C   third term  <\delta_x,n_3(2)+1/2 n12 n_1(2)> to large for remainder                     
                  d2a=x1*v1a+x2*v2a+x3*v3a
                  pD4a = pD4a + d2a*d2a/kap2
C   third term  <\delta_x,n_2(2)+1/2 n13 n_1(2)> to large for remainder                     
                  if(pD4a.ge.h4) CYCLE  
                  wght= 1.d0-sqrt(pD4a/h4)
C   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice
                  if(j1.eq.0) THEN
                     anz=1.d0
                  ELSE
                     anz=2.d0
                  ENDIF
                  sw=sw+anz*wght
                  if(j4.eq.i4) sw0=sw0+anz*wght
                  wght=wght*wght
                  sw2=sw2+anz*wght
                  if(j4.eq.i4) sw20=sw20+anz*wght
                  n=n+1
                  mj1=max(j1,mj1)
                  mj2=max(abs(j2),mj2)
                  mj3=max(abs(j3),mj3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      vred = sw*sw/sw2
      vred0 = sw0*sw0/sw20
      rad = max(mj1,max(mj2,mj3))
      if(rad.gt.h) THEN
      call dblepr("h",1,h,1)
      call intpr("radius",6,rad,1)
      END IF
      RETURN
      END
      subroutine lkse3i(h,kappa,i4,grad,grad2,grad3,ngrad,vext,ind,
     1                  wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n),i4
      real*8 h,kappa,grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad),
     1       vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4
      real*8 v1a,v2a,v3a,w1a,w2a,w3a,d2a,d3a,
     1       d1,d2,pD4,pD4a,kap4,kap2,h2,vd2,vd3,x1,x2,x3,
     2       n12,n13,nx1,h4,gi1,gi2,gi3,gj1,gj2,gj3
      ih1 = 2.0d0*h
      ih2 = 2.0d0*h/vext(1)
      ih3 = 2.0d0*h/vext(2)
      h2 = h*h
      h4 = h2*h2
      kap2 = kappa*kappa
      kap4 = kap2*kap2
      vd2 = vext(1)
      vd3 = vext(2)
      i = 1
      gi1 = grad(1,i4)
      gi2 = grad(2,i4)
      gi3 = grad(3,i4)
      DO j4 = 1,ngrad
         gj1 = grad(1,j4)
         gj2 = grad(2,j4)
         gj3 = grad(3,j4)
         n13 = gi2*grad3(1,j4)+gi3*grad3(2,j4)
         d1 = n13*n13
         pD4 = d1*d1/kap4
         if(pD4.ge.h4) CYCLE  
         n12 = gi1*grad2(1,j4)+gi2*grad2(2,j4)+gi3*grad2(3,j4)
         d1 = n12*n12
         pD4 = pD4 + d1*d1/kap4
         if(pD4.ge.h4) CYCLE  
         w1a=0.5d0*n13*gj1
         w2a=grad3(1,j4)+0.5d0*n13*gj2
         w3a=grad3(2,j4)+0.5d0*n13*gj3
         v1a=grad2(1,j4)+0.5d0*n12*gj1
         v2a=grad2(2,j4)+0.5d0*n12*gj2
         v3a=grad2(3,j4)+0.5d0*n12*gj3
         DO j1 = 0,ih1
            x1 = j1
            DO j2 = -ih2,ih2
               x2 = vd2*j2
               DO j3 = -ih3,ih3
                  x3 = vd3*j3
                  nx1 = x1*gj1+x2*gj2+x3*gj3
                  d2 = nx1*nx1
                  if(d2.ge.h2) CYCLE  
C   second term  <\delta_x,n_1(2)> to large
                  pD4a = pD4 + d2*d2
                  if(pD4a.ge.h4) CYCLE  
C   first term  <n_1(1),n_2(2)> to large for remainder
C
C   now third term if something is left
C
                  d3a=x1*w1a+x2*w2a+x3*w3a
                  pD4a = pD4a + d3a*d3a/kap2
                  if(pD4a.ge.h4) CYCLE  
C   third term  <\delta_x,n_3(2)+1/2 n12 n_1(2)> to large for remainder                     
                  d2a=x1*v1a+x2*v2a+x3*v3a
                  pD4a = pD4a + d2a*d2a/kap2
C   third term  <\delta_x,n_2(2)+1/2 n13 n_1(2)> to large for remainder                     
                  if(pD4a.ge.h4) CYCLE  
                  if(i.gt.n) THEN
                     call intpr("Exceeded max i",14,i,1)
                     call intpr("for i4",6,i4,1)
                     n = i-1
                     return
                  END IF
                  wght(i)= 1.d0-sqrt(pD4a/h4)
                  ind(1,i) = j1
                  ind(2,i) = j2
                  ind(3,i) = j3
                  ind(4,i) = i4
                  ind(5,i) = j4
                  i = i+1
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine lkfse3i(h,kappa,i4,k456,nbg,nbghat,ng,
     1                   vext,ind,wght,n)
      implicit logical (a-z)
      integer ng,n,ind(5,n),i4
      real*8 h,kappa,k456(3,ng,ng),
     1       nbg(3,3,ng),nbghat(3,3,ng,ng),vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4
      real*8 h2,kap2,x1,x2,x3,xh1,xh2,xh3,k1,k2,k3,k4,k5,k6,z,z1,
     1       vd2,vd3,gi1,gi2,gi3,cb
      ih1 = max(1.d0,5.0d0*h)
      ih2 = max(1.d0,5.0d0*h/vext(1))
      ih3 = max(1.d0,5.0d0*h/vext(2))
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)
      vd3 = vext(2)
      i = 1
      gi1 = nbg(1,1,i4)
      gi2 = nbg(2,1,i4)
      gi3 = nbg(3,1,i4)
      DO j4 = 1,ng
         cb = abs(gi1*nbg(1,1,j4)+gi2*nbg(2,1,j4)+gi3*nbg(3,1,j4))
         k4 = k456(1,i4,j4)
         k5 = k456(2,i4,j4)
         k6 = k456(3,i4,j4)
         z = (k4*k4+k5*k5+k6*k6)/kap2
         if(z.gt.h2) CYCLE
C   last three komponents already to large
         DO j1 = 0,ih1
            x1 = j1
            DO j2 = -ih2,ih2
               x2 = vd2*j2
               DO j3 = -ih3,ih3
                  x3 = vd3*j3
                  xh1=x1*nbg(1,2,j4)+x2*nbg(2,2,j4)+x3*nbg(3,2,j4)
                  xh2=x1*nbg(1,3,j4)+x2*nbg(2,3,j4)+x3*nbg(3,3,j4)
                  xh3=x1*nbg(1,1,j4)+x2*nbg(2,1,j4)+x3*nbg(3,1,j4)
C thats xhat
C now k1       
                  k1=xh1*nbghat(1,2,i4,j4)+xh2*nbghat(2,2,i4,j4)+
     1               xh3*nbghat(3,2,i4,j4)
                  z1=z+k1*k1
                  if(z1.gt.h2) CYCLE
C   last three komponents + first already to large
                  k2=xh1*nbghat(1,3,i4,j4)+xh2*nbghat(2,3,i4,j4)+
     1               xh3*nbghat(3,3,i4,j4)
                  z1=z1+k2*k2
                  if(z1.gt.h2) CYCLE
C   last three komponents + first two already to large
                  k3=xh1*nbghat(1,1,i4,j4)+xh2*nbghat(2,1,i4,j4)+
     1               xh3*nbghat(3,1,i4,j4)
                  z1=z1+k3*k3
                  if(z1.gt.h2) CYCLE
                  if(i.gt.n) THEN
                     call intpr("Exceeded max i",14,i,1)
                     call intpr("for i4",6,i4,1)
                     n = i-1
                     return
                  END IF
                  wght(i)= (1.d0-z1/h2)
C   *cb
                  ind(1,i) = j1
                  ind(2,i) = j2
                  ind(3,i) = j3
                  ind(4,i) = i4
                  ind(5,i) = j4
                  i = i+1
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine lkfse3i0(h,kappa,i4,k456,nbg,nbghat,ng,
     1                   vext,vred,vred0,n)
      implicit logical (a-z)
      integer ng,n,i4
      real*8 h,kappa,k456(3,ng,ng),
     1       nbg(3,3,ng),nbghat(3,3,ng,ng),vext(2),vred,vred0
      integer ih1,ih2,ih3,j1,j2,j3,j4,mj1,mj2,mj3,rad
      real*8 x1,x2,x3,xh1,xh2,xh3,k1,k2,k3,k4,k5,k6,z,z1,
     1       sw,sw2,sw0,sw20,wght,anz,h2,kap2,vd2,vd3,
     2       gi1,gi2,gi3,cb
      ih1 = 5.0d0*h
      ih2 = 5.0d0*h/vext(1)
      ih3 = 5.0d0*h/vext(2)
      sw=0.d0
      sw0=0.d0
      sw2=0.d0
      sw20=0.d0
      mj1=0
      mj2=0
      mj3=0
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)
      vd3 = vext(2)
      n = 1
      gi1 = nbg(1,1,i4)
      gi2 = nbg(2,1,i4)
      gi3 = nbg(3,1,i4)
      DO j4 = 1,ng
         cb = abs(gi1*nbg(1,1,j4)+gi2*nbg(2,1,j4)+gi3*nbg(3,1,j4))
         k4 = k456(1,i4,j4)
         k5 = k456(2,i4,j4)
         k6 = k456(3,i4,j4)
         z = (k4*k4+k5*k5+k6*k6)/kap2
         if(z.gt.h2) CYCLE
C   last three komponents already to large
         DO j1 = 0,ih1
            x1 = j1
            DO j2 = -ih2,ih2
               x2 = vd2*j2
               DO j3 = -ih3,ih3
                  x3 = vd3*j3
                  xh1=x1*nbg(1,2,j4)+x2*nbg(2,2,j4)+x3*nbg(3,2,j4)
                  xh2=x1*nbg(1,3,j4)+x2*nbg(2,3,j4)+x3*nbg(3,3,j4)
                  xh3=x1*nbg(1,1,j4)+x2*nbg(2,1,j4)+x3*nbg(3,1,j4)
C thats xhat
C now k1       
                  k1=xh1*nbghat(1,2,i4,j4)+xh2*nbghat(2,2,i4,j4)+
     1               xh3*nbghat(3,2,i4,j4)
                  z1=z+k1*k1
                  if(z1.gt.h2) CYCLE
C   last three komponents + first already to large
                  k2=xh1*nbghat(1,3,i4,j4)+xh2*nbghat(2,3,i4,j4)+
     1               xh3*nbghat(3,3,i4,j4)
                  z1=z1+k2*k2
                  if(z1.gt.h2) CYCLE
C   last three komponents + first two already to large
                  k3=xh1*nbghat(1,1,i4,j4)+xh2*nbghat(2,1,i4,j4)+
     1               xh3*nbghat(3,1,i4,j4)
                  z1=z1+k3*k3
                  if(z1.gt.h2) CYCLE
                  wght= (1.d0-z1/h2)
C        *cb
C   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice
                  if(j1.eq.0) THEN
                     anz=1.d0
                  ELSE
                     anz=2.d0
                  ENDIF
                  sw=sw+anz*wght
                  if(j4.eq.i4) sw0=sw0+anz*wght
                  wght=wght*wght
                  sw2=sw2+anz*wght
                  if(j4.eq.i4) sw20=sw20+anz*wght
                  n=n+1
                  mj1=max(j1,mj1)
                  mj2=max(abs(j2),mj2)
                  mj3=max(abs(j3),mj3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      vred = sw*sw/sw2
      vred0 = sw0*sw0/sw20
      rad = max(mj1,max(mj2,mj3))
      if(rad.gt.h) THEN
         call dblepr("h",1,h,1)
         call intpr("radius",6,rad,1)
      END IF
      RETURN
      END
      subroutine lkfulse3(h,kappa,k456,nbg,nbghat,ng,
     1                   vext,ind,wght,n)
      implicit logical (a-z)
      integer ng,n,ind(5,n)
      real*8 h(ng),kappa(ng),k456(3,ng,ng),
     1       nbg(3,3,ng),nbghat(3,3,ng,ng),vext(2),wght(n)
      integer ns,ni,i
      ns = 0
      DO i = 1,ng
         ni = n-ns
         call lkfse3i(h(i),kappa(i),i,k456,nbg,nbghat,ng,
     1                vext,ind(1,ns+1),wght(ns+1),ni)
         ns = ns+ni
      END DO
      n = ns
      RETURN
      END      
      subroutine lkse3(h,kappa,grad,grad2,grad3,ngrad,vext,ind,wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n)
      real*8 h(ngrad),kappa(ngrad),grad(3,ngrad),grad2(3,ngrad),
     1       grad3(2,ngrad),vext(2),wght(n)
      integer ns,ni,i
      ns = 0
      call grad23(grad,ngrad,grad2,grad3)
      DO i = 1,ngrad
         ni = n-ns
         call lkse3i(h(i),kappa(i),i,grad,grad2,grad3,ngrad,vext,
     1               ind(1,ns+1),wght(ns+1),ni)
         ns = ns+ni
      END DO
      n = ns
      RETURN
      END      
      subroutine adasmse3(y,th,ni,mask,n1,n2,n3,ngrad,lambda,ind,w,n,
     1                    thn,r,s2inv,sw,swy,swy2)
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C   s2inv containes the inverses of estimated variances
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,n,ind(5,n)
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),th(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1       lambda,w(n),thn(n1,n2,n3,ngrad),sw(ngrad),swy(ngrad),
     2       swy2(ngrad),s2inv(n1,n2,n3,ngrad),r(n1,n2,n3,ngrad)
      integer i,i1,i2,i3,i4,j1,j2,j3,j4
      real*8 z,yj,swyi,nii,thi
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               DO i4=1,ngrad
                  sw(i4)=0.d0
                  swy(i4)=0.d0
                  swy2(i4)=0.d0
               END DO
               i4=0
               DO i=1,n
                  if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
                     i4 = ind(4,i)
                     thi = th(i1,i2,i3,i4)
                     nii = ni(i1,i2,i3,i4)/lambda
                  END IF
                  j1=i1+ind(1,i)
                  j2=i2+ind(2,i)
                  j3=i3+ind(3,i)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  if(.not.mask(j1,j2,j3)) CYCLE          
C                  i4=ind(4,i)
                  j4=ind(5,i)
                  z=thi-th(j1,j2,j3,j4)
                  z=z*z*nii
C                  z=th(i1,i2,i3,i4)-th(j1,j2,j3,j4)
C                  z=z*z*ni(i1,i2,i3,i4)/lambda
C                  call spenalty(th(i1,i2,i3,i4),th(j1,j2,j3,j4),
C     1                          s2inv(i1,i2,i3,i4),ncoil,z)
C                  z=z*ni(i1,i2,i3,i4)/lambda
                  if(z.ge.1.d0) CYCLE
                  z=w(i)*min(1.d0,2.d0-2.d0*z)
                  z=z*s2inv(j1,j2,j3,j4)
                  sw(i4)=sw(i4)+z
                  yj=y(j1,j2,j3,j4)
                  swy(i4)=swy(i4)+z*yj
                  swy2(i4)=swy2(i4)+z*yj*yj
               END DO
               DO i=1,n
                  if(ind(1,i).eq.0) CYCLE
                  if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
                     i4 = ind(4,i)
                     thi = th(i1,i2,i3,i4)
                     nii = ni(i1,i2,i3,i4)/lambda
                  END IF
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
                  j1=i1-ind(1,i)
                  j2=i2-ind(2,i)
                  j3=i3-ind(3,i)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  if(.not.mask(j1,j2,j3)) CYCLE          
C                  i4=ind(4,i)
                  j4=ind(5,i)
                  z=thi-th(j1,j2,j3,j4)
                  z=z*z*nii
C                  z=th(i1,i2,i3,i4)-th(j1,j2,j3,j4)
C                  z=z*z*ni(i1,i2,i3,i4)/lambda
C                  call spenalty(th(i1,i2,i3,i4),th(j1,j2,j3,j4),
C     1                          s2inv(i1,i2,i3,i4),ncoil,z)
C                  z=z*ni(i1,i2,i3,i4)/lambda
C    ni(i1,i2,i3,i4) contains normalization by siinv
                  if(z.ge.1.d0) CYCLE
                  z=w(i)*min(1.d0,2.d0-2.d0*z)
                  z=z*s2inv(j1,j2,j3,j4)
                  sw(i4)=sw(i4)+z
                  yj=y(j1,j2,j3,j4)
                  swy(i4)=swy(i4)+z*yj
                  swy2(i4)=swy2(i4)+z*yj*yj
               END DO
               DO i4=1,ngrad
                  swyi = swy(i4)
                  nii=sw(i4)
                  thn(i1,i2,i3,i4) = swyi/nii
                  ni(i1,i2,i3,i4) = nii
                  if(nii.gt.s2inv(i1,i2,i3,i4)) THEN
                     r(i1,i2,i3,i4)=swyi/sqrt(nii*swy2(i4)-swyi*swyi)
                  ELSE
                     r(i1,i2,i3,i4)=0.d0
                  END IF
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine spenalty(thi,thj,s2inv,ncoil,pen)
      implicit logical (a-z)
C currently not used 
C will be insertet later to handle "Rician" errors
      integer ncoil
      real*8 thi,thj,s2inv,pen
      real*8 z
      if(ncoil.eq.0) THEN
         z=thi-thj
         pen=z*z*s2inv
         RETURN
      END IF
      RETURN
      END
