      subroutine nlrdti(s,nb,n1,n2,n3,mask,b,th0,D,niter,eps,Varth,
     1                  res,rss)
      implicit logical (a-z)
      integer nb,n1,n2,n3,s(nb,n1,n2,n3),niter
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),b(6,nb),Varth(28,n1,n2,n3),res(nb,n1,n2,n3),
     1    th0(n1,n2,n3),eps,rss(n1,n2,n3)
      integer i1,i2,i3,j
      DO i3=1,n3
         call intpr("Nonlinear regression for slice No:",34,i3,1)
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
               call solvedti(s(1,i1,i2,i3),nb,b,th0(i1,i2,i3),
     1                       D(1,i1,i2,i3),Varth(1,i1,i2,i3),
     2                       res(1,i1,i2,i3),niter,eps,
     3                       rss(i1,i2,i3))
               ELSE
                  DO j=1,6
                     D(j,i1,i2,i3)=0.d0
                  END DO
                  rss(i1,i2,i3)=0.d0
               END IF
               if((i1*i2*i3).eq.4913) THEN
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine solvedti(s,nb,b,th0,D,Varth,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(6,nb),th0,Varth(28),F(nb),eps
      integer i,j,k,info,iter,indvar
      real*8 z,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,Dn(6),res,X(7),th0n
      external indvar
      alpha=0.5D0
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      DO i=1,nb
         z=0.d0
         DO j=1,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res
         F(i)=res
      END DO
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO            
         DO i=1,nb
            z=0.d0
            DO j=1,6
               z=z+b(j,i)*D(j)
            END DO
            z=exp(-z)
            X(7)= -z
            z=z*th0
            DO j=1,6
               X(j)=b(j,i)*z
            END DO
            DO j=1,7
               dg(j)=dg(j)+X(j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(j)*X(k)
               END DO
            END DO 
         END DO
         maxabsdg=0.d0
         DO j=1,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-6) THEN
C  prepare things for return if gradient is close to 0
            DO j=1,7
               DO k=j,7
                  Varth(indvar(j,k))=ak(j,k)
               END DO
            END DO
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         DO WHILE (notacc) 
            IF(gamma.lt.1.d0) THEN
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=gamma*ak(j,k)
                  END DO
                  ck(j,j)=ck(j,j)+1.d0-gamma
               END DO
            ELSE
C   we may still need ak and dg so copy them to pk and ck
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=ak(j,k)
                  END DO
               END DO
            END IF
            DO j=1,7
               pk(j)=dg(j)
            END DO
C   Now solve  ak%*%dtheta= dg
	    call dposv("U",7,1,ck,7,pk,7,info)
C  Step 4 we have pk 
            IF(info.ne.0) THEN
               gamma=alpha*gamma
C  thats step 6
            ELSE
C  comute things needed for decision in step 5 
C  if successful F, nrss, and theta will be reused in the  
C  next iteration
               DO j=1,6
                  Dn(j)=D(j)-gamma*pk(j)
               END DO
               th0n=th0-gamma*pk(7)
               nrss=0.d0
               DO i=1,nb
                  z=0.d0
                  DO j=1,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=s(i)-th0n*exp(-z)
                  nrss=nrss+res*res
                  F(i)=res
               END DO
               crss=0.d0
               DO j=1,7
                  crss=crss+dg(j)*pk(j)
               END DO
               crss=rss-delta*gamma*crss
               IF(nrss.le.crss) THEN
                  notacc=.FALSE.
C  accept new estimate, prepare for next iteration
               ELSE
                  gamma=alpha*gamma
C  decrease gamma and try new regularization
               END IF
            END IF
         END DO
         th0=th0n
         DO j=1,6
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      DO j=1,7
         DO k=j,7
            Varth(indvar(j,k))=ak(j,k)
         END DO
      END DO
      RETURN
      END
      subroutine replvar(x,ngrad,n1,n2,n3,rind,tind,lind,nind,sigma2,
     1                   mean)
C
C  estimate variances from replicated gradient images
C  returns  df*variance estimate which should be distributed 
C  chisq(df) * variance
C
      implicit logical (a-z)
      integer ngrad,n1,n2,n3,nind
      integer x(ngrad,n1,n2,n3),rind(ngrad),tind(nind),lind(nind)
      real*8 sigma2(n1,n2,n3),mean(n1,n2,n3)
      integer i,j,k,l,m,r,s
      real*8 z
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               sigma2(i,j,k)=0.d0
            END DO
         END DO
      END DO
      DO l=1,nind
         if(tind(l).le.1) CYCLE
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  mean(i,j,k)=0.d0
               END DO
            END DO
         END DO
         m=lind(l)
         s=0
         DO r=1,ngrad
            if(rind(r).ne.m) CYCLE
            s=s+1
            DO i=1,n1
               DO j=1,n2
                  DO k=1,n3
                     z=x(r,i,j,k)
                     mean(i,j,k)=mean(i,j,k)+z
                     sigma2(i,j,k)=sigma2(i,j,k)+z*z
                  END DO
               END DO
            END DO
         END DO
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  z=mean(i,j,k)
                  sigma2(i,j,k)=sigma2(i,j,k)-z*z/tind(l)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
