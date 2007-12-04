CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsnrdti(si,nb,n1,n2,n3,mask,btb,sigma2,wlse,
     1                    th0,th0n,D,Dn,Varth,rss,bi,ani,andir,det,
     2                    sigma2h,sigma2n,sigma2r,h,niter,zext,rho,
     3                    lambda,swsi,swsi2,swsi4,F,eps,rician)
C
C   si       -  observed diffusion weighted images
C   nb       -  number of gradients (including zero gradients)
C   n1,n2,n3 -  dimensions of data cube
C   mask     -  logical: 1 for voxel inside region of interest
C   btb      -  b^T%*%b obtained from matrix of gradients b
C   sigma2   -  voxelwise variance estimates (data)
C   th0      -  current estimate of mean s_0 values
C   th0n     -  new  estimate of mean s_0 values (output)
C   D        -  current estimate of tensor
C   Dn       -  new estimate of tensor (output)
C   Varth    -  covariance of (D,th0) (input/output)
C   rss      -  residual sum of squares (output)
C   bi       -  sum of weights
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(D)
C   sigma2h  -  
C   sigma2n  -
C   h        -  actual bandwidth
C   niter    -  number of iterations for imbedded nl-regression
C   zext     -  
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   swsi     -  auxiliary for sum_j w_j si(j)
C   F        -  auxiliary for function values in LSE
C   eps      -  something small and positive
      implicit logical (a-z)
      integer n1,n2,n3,nb,si(nb,n1,n2,n3),niter
      real*8 btb(6,nb),sigma2(n1,n2,n3),swsi2(nb),swsi4(nb),
     1       th0(n1,n2,n3),th0n(n1,n2,n3),sigma2r(n1,n2,n3),
     2       D(6,n1,n2,n3),Dn(6,n1,n2,n3),Varth(28,n1,n2,n3),
     3       bi(n1,n2,n3),ani(n1,n2,n3),andir(3,n1,n2,n3),
     4       det(n1,n2,n3),sigma2h(n1,n2,n3),sigma2n(n1,n2,n3),h,rho,
     5       zext,lambda,swsi(nb),F(nb),eps,
     6       rss(n1,n2,n3)
      logical mask(n1,n2,n3),rician,wlse
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,iz
      real*8 wij,adist,sw,sws0,h3,thi(7),bii,sqrbii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,sij,deti,z,sew,eps3,ss2,sw0,Di(6),dtidisnr,
     2       vth(28),th0i,abessel,az,mswsi2,mswsi2q,mswsi4,s2hat,
     3       rhosw0,crhosw0
      external adist,dtidisnr
      logical aws,lprint
      aws=lambda.lt.1e20
      h3=h*h*h
      eps3=eps*eps*eps
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               lprint=(i1*i2*i3).eq.4913
               if(.not.mask(i1,i2,i3)) CYCLE
               sw=0.d0
               sw0=0.d0
               sws0=0.d0
               ss2=0.d0
               DO k=1,nb
                  swsi(k)=0.d0
                  swsi2(k)=0.d0
                  swsi4(k)=0.d0
              END DO
               DO k=1,28
                  vth(k)=Varth(k,i1,i2,i3)
               END DO
               s2hat = sigma2h(i1,i2,i3)
               deti=exp(log(det(i1,i2,i3))/3)
               bii=bi(i1,i2,i3)
               if(wlse) THEN
                  sqrbii=sigma2h(i1,i2,i3)/sqrt(bii)
               ELSE
                  bii=bii/s2hat
                  sqrbii=sqrt(1.d0/bii)
               END IF
               th0i=th0(i1,i2,i3)
               th0n(i1,i2,i3)=th0i
C    used as initial values
               DO k=1,6
                  Di(k)=D(k,i1,i2,i3)
                  Dn(k,i1,i2,i3)=Di(k)
C    used as initial values
                  thi(k)=Di(k)/deti
               END DO
               thi(1)=thi(1)+rho*sqrbii
               thi(4)=thi(4)+rho*sqrbii
               thi(6)=thi(6)+rho*sqrbii
C  this is scale invariant sice sqrbii scales with sqrt(sigma2) (standard deviation)
               call eigen3(thi,ew,ev,ierr)
               if(ierr.ne.0.or.ew(1).le.eps) THEN
C                  call intpr("ierr",4,ierr,1)
                  thi(1)=1
                  thi(2)=0
                  thi(3)=0
                  thi(4)=1
                  thi(5)=0
                  thi(6)=1
               ELSE
                  sew=ew(1)*ew(2)*ew(3)
                  sew=exp(log(sew)/3.d0)
                  ew(1)=ew(1)/sew
                  ew(2)=ew(2)/sew
                  ew(3)=ew(3)/sew
                  thi(1)=ev(1,1)*ev(1,1)/ew(1)+ev(1,2)*ev(1,2)/ew(2)+
     1                   ev(1,3)*ev(1,3)/ew(3)
                  thi(2)=ev(1,1)*ev(2,1)/ew(1)+ev(1,2)*ev(2,2)/ew(2)+
     1                   ev(1,3)*ev(2,3)/ew(3)
                  thi(3)=ev(1,1)*ev(3,1)/ew(1)+ev(1,2)*ev(3,2)/ew(2)+
     1                   ev(1,3)*ev(3,3)/ew(3)
                  thi(4)=ev(2,1)*ev(2,1)/ew(1)+ev(2,2)*ev(2,2)/ew(2)+
     1                   ev(2,3)*ev(2,3)/ew(3)
                  thi(5)=ev(2,1)*ev(3,1)/ew(1)+ev(2,2)*ev(3,2)/ew(2)+
     1                   ev(2,3)*ev(3,3)/ew(3)
                  thi(6)=ev(3,1)*ev(3,1)/ew(1)+ev(3,2)*ev(3,2)/ew(2)+
     1                   ev(3,3)*ev(3,3)/ew(3)
               END IF
C   create needed estimates of s_i
               call rangex(thi,h,j1a,j1e)
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey(thi,j1,h,j2a,j2e)
                  DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez(thi,j1,j2,h,j3a,j3e,zext)
                      DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE
                        wij=adist(thi,j1,j2,j3,zext)
C     triangular location kernel
                        if(wij.gt.h3) CYCLE
                        wij = (1.d0 - wij/h3)
                        IF(aws) THEN
                           sij = bii*dtidisnr(th0i,Di,
     1               th0(jj1,jj2,jj3),D(1,jj1,jj2,jj3),Vth)/lambda
                           if(sij.gt.1.d0) CYCLE
                           wij=wij*(1.d0-sij)
                        END IF
                        ss2=ss2+wij*sigma2(jj1,jj2,jj3)
                        sw0=sw0+wij
                        if(wlse) wij=wij/sigma2h(jj1,jj2,jj3)
                        if(wij.lt.0.d0) call dblepr("wij",3,wij,1)
                        sw=sw+wij
                        if(rician) THEN
                           DO k=1,nb
                              z=si(k,jj1,jj2,jj3)
                              swsi(k)=swsi(k)+wij*z
                              z=z*z
                              swsi2(k)=swsi2(k)+wij*z
                              swsi4(k)=swsi4(k)+wij*z*z
                           END DO
                        ELSE
                           DO k=1,nb
                              swsi(k)=swsi(k)+wij*si(k,jj1,jj2,jj3)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
               bi(i1,i2,i3)=sw
               if(sw0.lt.1.d0) call dblepr("sw0",3,sw0,1)
               if(rician.and.sw0.gt.1.d0) THEN
                  mswsi2=0.d0
                  mswsi2q=0.d0
                  mswsi4=0.d0
                  DO k=1,nb
                     z = swsi2(k)/sw
                     mswsi2=mswsi2+z
                     mswsi2q=mswsi2q+z*z
                     mswsi4=mswsi4+swsi4(k)
                  END DO
                  mswsi2=mswsi2/nb
                  mswsi2q=mswsi2q/nb
                  mswsi4=mswsi4/nb/sw
                  s2hat = mswsi2q+mswsi2*mswsi2-mswsi4
                  if(s2hat.lt.0.d0) THEN
                     call dblepr("mswsi2",6,mswsi2,1)
                     call dblepr("mswsi2q",7,mswsi2q,1)
                     call dblepr("mswsi4",6,mswsi4,1)
                     call dblepr("s2arg",5,s2hat,1)
                     s2hat=0.d0
                  END IF
                  s2hat = 0.5d0*(mswsi2-sqrt(s2hat))
                  sigma2r(i1,i2,i3)=s2hat
C  thats the joint moment estimate of the Rice variance based on the 2nd and 4th moment
                  rhosw0=sqrt((sw0-1)/sw0)
                  crhosw0=1.d0-rhosw0
                  DO k=1,nb
                  swsi(k)=rhosw0*sqrt(max(0.d0,swsi2(k)/sw-2*s2hat))+
     1                       crhosw0*swsi(k)/sw
                  END DO
               ELSE
                  Do k=1,nb
                     swsi(k)=swsi(k)/sw
                  END DO                  
               END IF
               IF(sw.gt.0.d0.and.sw.lt.1.d20) THEN
                  sigma2n(i1,i2,i3)=ss2/sw0
                  call dsolvdti(swsi,nb,btb,th0n(i1,i2,i3),
     1                          Dn(1,i1,i2,i3),Varth(1,i1,i2,i3),F,
     2                          niter,eps,rss(i1,i2,i3))
               ELSE
                  sigma2n(i1,i2,i3)=sigma2h(i1,i2,i3)
               END IF
               DO k=1,6
                  Di(k)=Dn(k,i1,i2,i3)
               END DO
               call eigen3(Di(1),ew,ev,ierr)
               IF(ierr.ne.0) THEN
C
C       error in eigenvalue decomposition of tensor, set voxel inactive
C
                  Di(1)=eps
                  Di(2)=0.d0
                  Di(3)=0.d0
                  Di(4)=eps
                  Di(5)=0.d0
                  Di(6)=eps
               ELSE IF(dmin1(ew(1),ew(2),ew(3)).lt.eps) THEN
C
C       negative eigenvalue of tensor, project to space of positive definite tensors
C
                  ew(1)=max(ew(1),eps)
                  ew(2)=max(ew(2),eps)
                  ew(3)=max(ew(3),eps)
                  Di(1)=ew(1)*ev(1,1)*ev(1,1)+
     1                     ew(2)*ev(1,2)*ev(1,2)+ew(3)*ev(1,3)*ev(1,3)
                  Di(2)=ew(1)*ev(1,1)*ev(2,1)+
     1                     ew(2)*ev(1,2)*ev(2,2)+ew(3)*ev(1,3)*ev(2,3)
                  Di(3)=ew(1)*ev(1,1)*ev(3,1)+
     1                     ew(2)*ev(1,2)*ev(3,2)+ew(3)*ev(1,3)*ev(3,3)
                  Di(4)=ew(1)*ev(2,1)*ev(2,1)+
     1                     ew(2)*ev(2,2)*ev(2,2)+ew(3)*ev(2,3)*ev(2,3)
                  Di(5)=ew(1)*ev(2,1)*ev(3,1)+
     1                     ew(2)*ev(2,2)*ev(3,2)+ew(3)*ev(2,3)*ev(3,3)
                  Di(6)=ew(1)*ev(3,1)*ev(3,1)+
     1                     ew(2)*ev(3,2)*ev(3,2)+ew(3)*ev(3,3)*ev(3,3)
               END IF
               mew=(ew(1)+ew(2)+ew(3))/3.d0
               z1=ew(1)-mew
               z2=ew(2)-mew
               z3=ew(3)-mew
               z=3.d0*(z1*z1+z2*z2+z3*z3)
               z1=ew(1)
               z2=ew(2)
               z3=ew(3)
               mew=2.d0*(z1*z1+z2*z2+z3*z3)
               if(mew.le.eps*eps) mew=1.d0
               ani(i1,i2,i3)=sqrt(z/mew)
               det(i1,i2,i3)=ew(1)*ew(2)*ew(3)
               DO k=1,3
                  andir(k,i1,i2,i3)=ev(k,3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      integer function indvar(i,j)
      implicit logical (a-z)
      integer i,j
      indvar=i+(j*(j-1))/2
      RETURN 
      END
      real*8 function dtidisnr(th0i,Di,th0j,Dj,Vth)
      implicit logical (a-z)
      real*8 th0i,Di(6),th0j,Dj(6),Vth(28)
      real*8 z,dth(7)
      integer i,j,indvar
      z=0.d0
      DO i = 1,6
         dth(i)=Di(i)-Dj(i)
      END DO
      dth(7)=th0i-th0j
      DO i=1,7
         z=z+dth(i)*dth(i)*Vth(indvar(i,i))
      END DO
      DO i=1,6
         DO j=i+1,7
            z=z+2.D0*dth(i)*dth(j)*Vth(indvar(i,j))
         END DO
      END DO
      dtidisnr=z
      RETURN
      END
      subroutine dsolvdti(s,nb,b,th0,D,Varth,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C  same as solvedti except that s is double
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb),D(6),b(6,nb),th0,Varth(28),F(nb),eps
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
      subroutine sihat(th0i,Di,btb,F,nb)
      implicit logical (a-z)
      integer nb
      real*8 th0i,Di(6),btb(6,nb),F(nb)
      integer j,k
      real*8 z
      DO k=1,nb
         z=0.d0
         DO j=1,6
            z=z+btb(j,k)*Di(j)
         END DO
         F(k)=th0i*exp(-z)
      END DO
      RETURN
      END
      real*8 function besselq(z,bq)
      real*8 z,bq(1001)
      integer iz
      real*8 az,z1
      z1=10*z
C   this depends on the discretization of besselI(.,1)/besselI(.,0)
      iz=z1                  
      az=z1-iz
      besselq=(1.d0-az)*bq(iz)+az*bq(iz+1)
      RETURN
      END
