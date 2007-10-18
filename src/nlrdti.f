      subroutine nlrdti(s,nb,n1,n2,n3,b,th0,D,niter,eps,VarD,res,X,F,
     1                  rss)
      implicit logical (a-z)
      integer nb,n1,n2,n3,s(nb,n1,n2,n3),niter
      real*8 D(6,n1,n2,n3),b(nb,6),VarD(21,n1,n2,n3),res(nb,n1,n2,n3),
     1    th0(n1,n2,n3),X(nb,7),F(nb),eps,rss(n1,n2,n3)
      integer i1,i2,i3,j
      real*8 theta(7),Varth(7,7)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call solvedti(s(1,i1,i2,i3),nb,b,th0(i1,i2,i3),
     1                       D(1,i1,i2,i3),VarD(1,i1,i2,i3),
     2                       res(1,i1,i2,i3),niter,eps,X,F,
     3                       rss(n1,n2,n3))
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine solvedti(s,nb,b,th0,D,VarD,res,niter,eps,X,F,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(nb,6),th0,VarD(21),res(nb),X(nb,7),F(nb),eps,
     1       theta(7)
      integer i,j,k,info,iter
      real*8 ntheta(7),Vtheta(7,7),z,thcorr,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg
      alpha=0.5D0
      delta=0.25D0
      theta(1)=th0
      theta(2)=D(1)
      theta(3)=D(2)
      theta(4)=D(3)
      theta(5)=D(4)
      theta(6)=D(5)
      theta(7)=D(6)
      gamma=1.d0
      alpha=0.7d0
      rss=0.d0
      DO i=1,nb
         z=b(i,1)*theta(2)+b(i,2)*theta(3)+b(i,3)*theta(4)+
     1     b(i,4)*theta(5)+b(i,5)*theta(6)+b(i,6)*theta(7)
         z=exp(-z)
         F(i)=s(i)-theta(1)*z
         rss=rss+F(i)*F(i)
      END DO
      DO iter=1,niter
         DO i=1,nb
            z=b(i,1)*theta(2)+b(i,2)*theta(3)+b(i,3)*theta(4)+
     1        b(i,4)*theta(5)+b(i,5)*theta(6)+b(i,6)*theta(7)
            z=exp(-z)
            X(i,1)= -z
            z=z*theta(1)
            X(i,2)=b(i,1)*z
            X(i,3)=b(i,2)*z
            X(i,4)=b(i,3)*z
            X(i,5)=b(i,4)*z
            X(i,6)=b(i,5)*z
            X(i,7)=b(i,6)*z
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
         IF(maxabsdg.lt.eps) THEN
C  prepare things for return if gradient is close to 0
            th0=theta(1)
            i=1
            DO j=1,6
               D(j)=theta(j+1)
               DO k=j,6
                  VarD(i)=ak(j,k)
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
                  ck(j,j)=ck(j,j)+1-gamma
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
C         DO j=2,7
C            DO k=1,j-1
C               ak(j,k)=ak(k,j)
C            END DO
C         END DO         
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
                  z=exp(-z)
                  F(i)=s(i)-ntheta(1)*z
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
         rss=nrss
         call rchkusr()
      END DO
      th0=theta(1)
      i=1
      DO j=1,6
         D(j)=theta(j+1)
         DO k=j,6
            VarD(i)=ak(j,k)
            i=i+1
         END DO
      END DO
      RETURN
      END