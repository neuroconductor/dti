C
C __________________________________________________________________
C
      subroutine mfunpl0(par,siq,grad,m,lpar,ngrad,pen,z,w,erg)
C
C   model without isotropic compartment 
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg,pen
      integer i,j,i3,mode,r
      real*8 c1,w(ngrad),sw,sth,z1,p0,p1,d1,d2,d3,work(1000),s(6)
      c1 = par(1)
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
      call dcopy(ngrad,siq,1,w,1)
      call dgelss(ngrad,m,1,z,ngrad,w,ngrad,s,-1.d0,r,work,1000,mode)
      IF(mode.ne.0) THEN
         call intpr("mode",4,mode,1)
         erg = 1d20
      ELSE
C         IF(work(1).gt.1000) THEN
C            call intpr("ngrad",5,ngrad,1)
C            call intpr("m",1,m,1)
C            call dblepr("optimal LWORK",13,work,1)
C         END IF
         sw=0.d0
C penalize for extreme c1 values
         if(c1.gt.8.d0) sw=sw+c1-8.d0
C penalize for negative weights
         if(c1.lt.1.d-2) sw=sw-1.d2*c1+1.d0
         DO i=1,m
            if(w(i).lt.0.d0) sw=sw-pen*w(i)
         END DO
         DO i=m+1,ngrad
            sw=sw+w(i)**2
         END DO
         erg=sw
      END IF
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine getsiin2(si,ngrad,n1,n2,n3,m,
     1         egrad,isample,ntry,sms,z,siind,mval,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  exgrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,ns,siind(ns,n1,n2,n3),m,ntry,
     1       isample(m,ntry)
      real*8 si(ngrad,n1,n2,n3),sms(ngrad),
     1       egrad(ngrad,ngrad),z(ngrad,ns),mval(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,k,ibest,mode,ind(10),l
      real*8 w(1000),krit,work1(1000),work2(10),erg
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
C  now search for minima of sms (or weighted sms
                  ibest=1
                  krit=1e10
                  DO k=1,ntry
                     call dcopy(ngrad,si(1,i1,i2,i3),1,sms,1)
                     DO l=1,m
                  call dcopy(ngrad,egrad(1,isample(l,k)),1,z(1,l),1)
                     END DO
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
                     IF(mode.gt.1) THEN
                        call intpr("mode",4,mode,1)
                        call intpr("isample",7,isample(1,k),m)
                     ELSE 
                        IF(erg.lt.krit) THEN
                           krit=erg
                           ibest=k
                        END IF  
                     END IF
                  END DO
                  siind(1,i1,i2,i3)=m
                  DO l=1,m
                     siind(l+1,i1,i2,i3)=isample(l,ibest)
                  END DO
                  mval(i1,i2,i3)=krit
               ELSE
                  siind(1,i1,i2,i3)=-1
                  mval(i1,i2,i3)=0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsiind(si,ngrad,n1,n2,n3,m,maxc,dgrad,
     1                    sms,siind,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  dgrad  - matrix of pairwise cos of gradient directions
C  sms    - copies of si
C  siind  - array of indices (output)
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,ns,siind(ns,n1,n2,n3),m
      real*8 si(ngrad,n1,n2,n3),sms(ngrad),dgrad(ngrad,ngrad),maxc
      logical mask(n1,n2,n3)
      integer i,i1,i2,i3,imin,k
      real*8 zmin,w
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
C               if(.not.mask(i1,i2,i3)) CYCLE
               if(mask(i1,i2,i3)) THEN
C  now search for minima of sms (or weighted sms
               call dcopy(ngrad,si(1,i1,i2,i3),1,sms,1)
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
               ELSE
                  siind(1,i1,i2,i3)=-1
               END IF
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
 