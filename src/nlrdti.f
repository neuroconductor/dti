      subroutine nlrdti(s,nb,n1,n2,n3,b,th0,D,niter,eps,Varth,res,X,
     1                  rss)
      implicit logical (a-z)
      integer nb,n1,n2,n3,s(nb,n1,n2,n3),niter
      real*8 D(6,n1,n2,n3),b(6,nb),Varth(28,n1,n2,n3),res(nb,n1,n2,n3),
     1    th0(n1,n2,n3),X(nb,7),F(nb),eps,rss(n1,n2,n3)
      integer i1,i2,i3,j
      real*8 theta(7)
      DO i1=1,n1
         call intpr("i1",2,i1,1)
         DO i2=1,n2
            DO i3=1,n3
               call solvedti(s(1,i1,i2,i3),nb,b,th0(i1,i2,i3),
     1                       D(1,i1,i2,i3),Varth(1,i1,i2,i3),
     2                       res(1,i1,i2,i3),niter,eps,X,
     3                       rss(n1,n2,n3))
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine solvedti(s,nb,b,th0,D,Varth,F,niter,eps,X,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(6,nb),th0,Varth(28),X(nb,7),F(nb),eps,
     1       theta(7)
      integer i,j,k,info,iter
      real*8 ntheta(7),Vtheta(7,7),z,thcorr,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss
      alpha=0.5D0
      delta=0.25D0
      theta(1)=th0
      theta(2)=D(1)
      theta(3)=D(2)
      theta(4)=D(3)
      theta(5)=D(4)
      theta(6)=D(5)
      theta(7)=D(6)
C      call dblepr("theta",5,theta,7)
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      DO i=1,nb
         z=b(1,i)*theta(2)+b(2,i)*theta(3)+b(3,i)*theta(4)+
     1     b(4,i)*theta(5)+b(5,i)*theta(6)+b(6,i)*theta(7)
         z=exp(-z)
         F(i)=s(i)-theta(1)*z
         rss=rss+F(i)*F(i)
      END DO
C      call dblepr("rss",3,rss,1)
      DO iter=1,niter
         DO i=1,nb
            z=b(1,i)*theta(2)+b(2,i)*theta(3)+b(3,i)*theta(4)+
     1        b(4,i)*theta(5)+b(5,i)*theta(6)+b(6,i)*theta(7)
            z=exp(-z)
            X(i,1)= -z
            z=z*theta(1)
            X(i,2)=b(1,i)*z
            X(i,3)=b(2,i)*z
            X(i,4)=b(3,i)*z
            X(i,5)=b(4,i)*z
            X(i,6)=b(5,i)*z
            X(i,7)=b(6,i)*z
         END DO
         DO j=1,7
            dg(j)=X(1,j)*F(1)
            DO k=j,7
               ak(j,k)=X(1,j)*X(1,k)
            END DO
         END DO
         DO i=2,nb
            DO j=1,7
               dg(j)=dg(j)+X(i,j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(i,j)*X(i,k)
               END DO
            END DO 
         END DO
         maxabsdg=abs(dg(1))
         DO j=2,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-5) THEN
C  prepare things for return if gradient is close to 0
            th0=theta(1)
            DO j=1,6
               D(j)=theta(j+1)
            END DO
            i=1
            DO j=1,7
               DO k=j,7
                  Varth(i)=ak(j,k)
                  i=i+1
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
               DO j=1,7
                  ntheta(j)=theta(j)-gamma*pk(j)
               END DO
               nrss=0.d0
               DO i=1,nb
                 z=b(i,1)*ntheta(2)+b(i,2)*ntheta(3)+b(i,3)*ntheta(4)+
     1             b(i,4)*ntheta(5)+b(i,5)*ntheta(6)+b(i,6)*ntheta(7)
                  F(i)=s(i)-ntheta(1)*exp(-z)
                  nrss=nrss+F(i)*F(i)
               END DO
               crss=dg(1)*pk(1)
               DO j=2,7
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
         DO j=1,7
            theta(j)=ntheta(j)
         END DO
C         call dblepr("theta",5,theta,7)
         oldrss=rss
         rss=nrss
C      call dblepr("rss",3,rss,1)
         call rchkusr()
C      call intpr("iter",4,iter,1)
      END DO
      th0=theta(1)
      DO j=1,6
         D(j)=theta(j+1)
      END DO
      i=1
      DO j=1,7
         DO k=j,7
            Varth(i)=ak(j,k)
            i=i+1
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
      integer i,j,k,l,m,r
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
         DO r=1,ngrad
            if(rind(r).ne.m) CYCLE
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
