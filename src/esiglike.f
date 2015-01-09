      subroutine lncchi(sigma,ni,ksi,wj,sj,L,n,work,erg)
C
C  compute local weighted noncentral chi^2 log-likelihood * (-1)
C  
C  sigma - parameter to estimate
C  ni    - sum(wj)
C  ksi   - sum(wj*Sj^2)/ni
C  wj    - local weights
C  sj    - observed values
C  L     - df/2
C  n     - number of local weights/observations
C
      implicit logical (a-z)
      integer n
      real*8 sigma,ni,ksi,wj(n),sj(n),L,work(*),erg
      integer j
      real*8 eta,z,sig2,zs,sl,lm1,za
      real*8 bessliex
      external bessliex
      lm1=L-1
      sig2=sigma*sigma
      eta=0.d0
      sl=ksi-2.d0*L*sig2
         z=sqrt(sl)
         zs=z/sig2
         DO j=1,n
            if(wj(j).gt.0.d0) THEN
               za=sj(j)*zs
               if(za.le.1d2) THEN
                  za=log(bessliex(za,lm1,1.d0,work))
               ELSE
                  za=za-log(za*6.283185d0)/2.d0
C  large value approximation
               END IF
C  avoid Inf in besseli (unlikely in optimum, does not change convexity)
               eta=eta+wj(j)*za
            END IF
         END DO
         erg=ksi/sig2+log(sig2)+lm1/2.d0*log(sl)-eta/ni
      RETURN
      END
      real*8 function lncchi2(sigma,ni,ksi,wj,sj,L,clws,n,work)
C
C  compute local weighted noncentral chi^2 log-likelihood * (-1)
C  
C  sigma - parameter to estimate
C  ni    - sum(wj)
C  ksi   - sum(wj*Sj^2)/ni
C  wj    - local weights
C  sj    - observed values
C  L     - df/2
C  n     - number of local weights/observations
C
      implicit logical (a-z)
      integer n
      real*8 sigma,ni,ksi,wj(n),sj(n),L,work(*)
      integer j
      real*8 eta,z,sig2,zs,pen,sl,lm1,za,clws
      real*8 bessliex
      external bessliex
      lm1=L-1
      sig2=sigma*sigma
      eta=0.d0
      sl=ksi-2.d0*L*sig2
      if(sl.lt.1d-6) THEN
         pen=sl-1d-6
C  \theta estimated as zero: central case
         lncchi2=L*log(sig2)+ksi/sig2+max(sl,0.d0)/2.d0+clws-pen
C clws contains (L-1)\sum_j wj log(sj) + ni (L-1) log2 + lgamma(L)
      ELSE
         z=sqrt(sl)
         zs=z/sig2
         DO j=1,n
C            if(wj(j).gt.0.d0) THEN
               za=sj(j)*zs
               if(za.le.1d2) THEN
                  za=log(bessliex(za,lm1,1.d0,work))
               ELSE
                  za=za-log(za*6.283185d0)/2.d0
C  large value approximation
               END IF
C  avoid Inf in besseli (unlikely in optimum, does not change convexity)
               eta=eta+wj(j)*za
C            END IF
         END DO
         lncchi2=ksi/sig2+log(sig2)+lm1/2.d0*log(sl)-eta/ni
      END IF
      RETURN
      END
C    
C    Minimization of noncentral chi2-likelihood
C    Adapted from procedure localmin in Richard Brent, Algorithms for
C    Minimization without Derivatives, Prentice-Hall, Inc. (1973)
C
      subroutine localmin(low,up,wj,sj,L,n,tol,maxit,work,xmin,fmin)
      implicit logical (a-z)
      integer n,maxit
      real*8 low,up,wj(n),sj(n),L,tol,xmin,fmin,work(*)
      real*8 goldc,a,b,d,e,eps,xm,p,q,r,eps1,eps2,u,v,w,fu,fv,fw,fx,x
      real*8 ni,ksi,sjj,clws,Lm1
      integer it,j
      logical gsect
      real*8 lncchi2,lgammaf
      external lncchi2,lgammaf
      eps=1d-8
      goldc=0.381966
C  Initialize
C  ni    - sum(wj)
C  ksi   - sum(wj*Sj^2)/ni
      Lm1=L-1
      ni=0.d0
      ksi=0.d0
      clws=0.d0
      DO j=1,n
C         if(wj(j).gt.0.d0) THEN
            ni=ni+wj(j)
            sjj=sj(j)
            ksi=ksi+wj(j)*sjj*sjj
            clws=clws+wj(j)*log(sjj)
C         END IF
      END DO
      clws=-Lm1*clws/ni+Lm1*log(2.d0)+lgammaf(L)
      ksi=ksi/ni
      a=low 
      b=up
      v=a+goldc*(b-a)
      w=v
      x=v
      e=0.d0
      d=0.d0
      fx=lncchi2(x,ni,ksi,wj,sj,L,clws,n,work)
      fv=fx
      fw=fx
C  Search for minimum
      DO it=1,maxit   
         xm=0.5d0*(a+b)
         eps1=eps*abs(x)+tol/3.d0
         eps2=2.d0*eps1
         IF (abs(x-xm).le.(eps2-0.5d0*(b-a)) ) EXIT 
         gsect = .TRUE.
         IF (abs(e) .gt. eps1) THEN  
            gsect = .FALSE.
            r=(x-w)*(fx-fv) 
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            IF (q .gt. 0.d0) p=-p
            q=abs(q)
            r=e
            e=d
            IF (abs(p).ge.abs(0.5d0*q*r).or.
     1                  p.le.q*(a-x).or.p.ge.q*(b-x)) THEN
                gsect = .TRUE.
            ELSE
               d=p/q   
               u=x+d
               IF ((u-a).lt.eps2.or.(b-u).lt.eps2) d=sign(eps1,xm-x)
            END IF
         END IF
         IF (gsect) THEN 
            IF (x .ge. xm) THEN  
               e=a-x
            ELSE
               e=b-x
            END IF
            d=goldc*e
         END IF
         IF (abs(d) .ge. eps1) THEN  
            u=x+d
         ELSE
            u=x+sign(eps1,d)
         END IF
         fu=lncchi2(u,ni,ksi,wj,sj,L,clws,n,work)
         IF (fu .le. fx) THEN
            IF (u .ge. x) THEN
               a=x
            ELSE
               b=x
            END IF
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         ELSE
            IF (u.lt.x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu.le.fw.or.w.eq.x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF (fu.le.fv.or.v.eq.x.or.v.eq.w) THEN
               v=u
               fv=fu
            END IF
         END IF
      END DO
      xmin=sqrt(ni/(ni-1.d0))*x
C this seems to correct the bias
      fmin=fx
      RETURN
      END 
 