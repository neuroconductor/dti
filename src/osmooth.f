CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine osmdti(si,siest,varinv,sipred,nb,n1,n2,n3,mask,
     1                  btb,sdcoef,th0,th0n,D,Dn,ani,andir,det,vol,
     2                  niter,vext,lambda,tau,swsi,swsi2,sw,sw2,
     3                  oibd,F,vari,sig2bi,s0ind,eps)
C
C   si       -  observed diffusion weighted images
C   nb       -  number of gradients (including zero gradients)
C   n1,n2,n3 -  dimensions of data cube
C   mask     -  logical: 1 for voxel inside region of interest
C   btb      -  b^T%*%b obtained from matrix of gradients b
C   sigma2   -  voxelwise varinviance estimates (data)
C   th0      -  current estimate of mean s_0 values
C   th0n     -  new  estimate of mean s_0 values (output)
C   D        -  current estimate of tensor
C   Dn       -  new estimate of tensor (output)
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(D)
C   h        -  actual bandwidth
C   niter    -  number of iterations for imbedded nl-regression
C   zext     -   
C   lambda   -  scale factor in the statistical penalty
C   swsi     -  auxiliary for sum_j w_j si(j)
C   F        -  auxiliary for function values in LSE
C   eps      -  something small and positive
      implicit logical (a-z)
      integer n1,n2,n3,nb,si(nb,n1,n2,n3),siest(nb,n1,n2,n3),niter
      real*8 btb(6,nb),swsi2(nb),th0(n1,n2,n3),
     1       th0n(n1,n2,n3),D(6,n1,n2,n3),
     2       Dn(6,n1,n2,n3),sipred(nb,n1,n2,n3),varinv(nb,n1,n2,n3),
     3       ani(n1,n2,n3),andir(3,n1,n2,n3),det(n1,n2,n3),h,
     4       vext(3),lambda,tau,swsi(nb),F(nb),vari(nb),eps,rss,
     5       vol,sw(nb),sw2(nb),sdcoef(4),sig2bi(nb)
      logical mask(n1,n2,n3),s0ind(nb)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,l,center
      real*8 wij,adist0,h2,ew(3),ev(3,3),rssj,la12,la13,oibd(nb),
     1       mew,z1,z2,z3,sij,z,Di(6),dtidisrg,minoibd,wijk,
     2       th0i,rssi,h0,h1,zsd,up,low,zsd2,meanni,taui
      integer nselect
      external adist0,dtidisrg
      logical aws
      aws=lambda.lt.1e20
      h1=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)
      h0=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)/1.4
C  first fill predicted 
      low=sdcoef(1)+sdcoef(3)*sdcoef(2)
      up=sdcoef(1)+sdcoef(4)*sdcoef(2)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call sihat(th0(i1,i2,i3),D(1,i1,i2,i3),btb,
     1                    sipred(1,i1,i2,i3),nb)
            END DO
         END DO
      END DO
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE

               DO k=1,nb
                  z=siest(k,i1,i2,i3)
                  zsd=sdcoef(1)+z*sdcoef(2)
                  if(z.lt.sdcoef(3)) zsd=low
                  if(z.gt.sdcoef(4)) zsd=up
                  sig2bi(k)=zsd*zsd
                  vari(k)=varinv(k,i1,i2,i3)
               END DO
               rssi = dtidisrg(siest(1,i1,i2,i3),
     1                         sipred(1,i1,i2,i3),vari,nb)
               DO k=1,nb
                  sw(k)=0.d0
                  sw2(k)=0.d0
                  swsi(k)=0.d0
                  swsi2(k)=0.d0
               END DO
               th0i=th0(i1,i2,i3)
               th0n(i1,i2,i3)=th0i
C    used as initial values
               IF(ew(1).lt.0.d0) call dblepr("C0",2,ew,3)
               DO k=1,6
                  Di(k)=D(k,i1,i2,i3)
                  Dn(k,i1,i2,i3)=Di(k)
C    used as initial values
               END DO
               call eigen3(Di,ew,ev,ierr)
C produces largest EW in last position
C               if(i1.eq.15.and.i2.eq.15.and.i3.eq.21) THEN
C                  call dblepr("ew",2,ew,3)
C                  call dblepr("ev",2,ev,9)
C               END IF
               
               la12 = (ew(2)-ew(3))/ew(3)
               la13 = (ew(1)-ew(3))/ew(3)
C   now get bandwidth such that the ellopsoid has specified volume
               call gethani0(h0,h1,vol,vext,1.d-2,h)
               h2=h*h
C   create needed estimates of s_i
               nselect=0
               call rangex0(h,j1a,j1e,vext)
               call oriscor2(la12,la13,ev,nb,btb,oibd)
C               call oriscore(la12,la13,ev,nb,btb,oibd)
               meanni=0.d0
               minoibd=-1.d0
               l=0
               DO k=1,nb
                  if(s0ind(k)) CYCLE
                  meanni=meanni+vari(k)*sig2bi(k)
                  minoibd=max(minoibd,oibd(k))
                  l=l+1
               END DO
C               meanni=meanni/l
C               meanni=sqrt(meanni/l)
               meanni=exp(log(meanni/l)/3.d0)
C               taui=ani(i1,i2,i3)/tau
               taui=1.d0/tau
               DO k=1,nb
                  if(s0ind(k)) THEN
                     oibd(k)=1.d0
                  ELSE
C                     oibd(k)=exp((oibd(k)-minoibd)*meanni*taui)
                     oibd(k)=exp((oibd(k)-minoibd)*taui)
C                     oibd(k)=exp(oibd(k)/tau)
                  END IF
               END DO
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey0(j1,h,j2a,j2e,vext)
                  DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez0(j1,j2,h,j3a,j3e,vext)
                     DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE
                        wij=adist0(j1,j2,j3,vext)
C     triangular location kernel
                        if(wij.ge.h2) CYCLE
                        center=abs(i1-jj1)+abs(i2-jj2)+abs(i3-jj3)
                        wij = (1.d0 - wij/h2)
                        IF(aws) THEN
                           rssj = dtidisrg(siest(1,i1,i2,i3),
     1                               sipred(1,jj1,jj2,jj3),vari,nb)
                           sij = max(0.d0,rssj-rssi)/lambda
                           if(center.eq.0) sij=0.d0
                           if(sij.gt.1.d0) CYCLE
C  use Plateau kernel
                           if(sij.gt.0.5d0) THEN
                              wij=wij*2.d0*(1.d0-sij)
                           END IF
                        END IF
                        DO k=1,nb
                           if(center.eq.0) THEN
                              wijk=wij
                           ELSE
                              wijk=wij*oibd(k)
                           END IF
                           z=si(k,jj1,jj2,jj3)
                           zsd=sdcoef(1)+z*sdcoef(2)
                           if(z.lt.sdcoef(3)) zsd=low
                           if(z.gt.sdcoef(4)) zsd=up
                           zsd2=zsd*zsd
                           wijk=wijk/zsd2
                           sw(k)=sw(k)+wijk
                           sw2(k)=sw2(k)+wijk*wijk*zsd2
                           swsi(k)=swsi(k)+wijk*si(k,jj1,jj2,jj3)
                        END DO
                    END DO
                  END DO
               END DO
               Do k=1,nb
                  swsi(k)=swsi(k)/sw(k)
               END DO
               DO k=1,nb
                  siest(k,i1,i2,i3)=swsi(k)
                  vari(k)=sw(k)*sw(k)/sw2(k)
                  varinv(k,i1,i2,i3)=vari(k)
               END DO
C               if(i1.eq.20.and.i2.eq.20.and.i3.eq.21) THEN
C                     call dblepr("oibd",4,oibd,nb)
C                     call dblepr("sw(k)",5,sw,nb)
C                     call dblepr("swsi(k)",7,swsi,nb)
C               END IF
               call dslvdti0(swsi,nb,btb,vari,th0n(i1,i2,i3),
     1                          Dn(1,i1,i2,i3),F,niter,eps,rss)
               call eigen3(Dn(1,i1,i2,i3),ew,ev,ierr)
               IF(ew(1).lt.0.d0) call dblepr("C0",2,ew,3)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C   This version uses a modified orientation dependent statistical penalty
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine osmdti2(si,siest,varinv,sipred,nb,n1,n2,n3,mask,
     1                  btb,sdcoef,th0,th0n,D,Dn,ani,andir,det,vol,
     2                  niter,vext,lambda,tau,swsi,swsi2,sw,sw2,
     3                  oibd,ojbd,F,vari,sig2bi,s0ind,eps)
C
C   si       -  observed diffusion weighted images
C   nb       -  number of gradients (including zero gradients)
C   n1,n2,n3 -  dimensions of data cube
C   mask     -  logical: 1 for voxel inside region of interest
C   btb      -  b^T%*%b obtained from matrix of gradients b
C   sigma2   -  voxelwise varinviance estimates (data)
C   th0      -  current estimate of mean s_0 values
C   th0n     -  new  estimate of mean s_0 values (output)
C   D        -  current estimate of tensor
C   Dn       -  new estimate of tensor (output)
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(D)
C   h        -  actual bandwidth
C   niter    -  number of iterations for imbedded nl-regression
C   zext     -   
C   lambda   -  scale factor in the statistical penalty
C   swsi     -  auxiliary for sum_j w_j si(j)
C   F        -  auxiliary for function values in LSE
C   eps      -  something small and positive
      implicit logical (a-z)
      integer n1,n2,n3,nb,si(nb,n1,n2,n3),siest(nb,n1,n2,n3),niter
      real*8 btb(6,nb),swsi2(nb),th0(n1,n2,n3),
     1       th0n(n1,n2,n3),D(6,n1,n2,n3),oibd(nb),ojbd(nb),
     2       Dn(6,n1,n2,n3),sipred(nb,n1,n2,n3),varinv(nb,n1,n2,n3),
     3       ani(n1,n2,n3),andir(3,n1,n2,n3),det(n1,n2,n3),h,
     4       vext(3),lambda,tau,swsi(nb),F(nb),vari(nb),eps,rss,
     5       vol,sw(nb),sw2(nb),sdcoef(4),sig2bi(nb)
      logical mask(n1,n2,n3),s0ind(nb)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,l,center
      real*8 wij,adist0,h2,ew(3),ev(3,3),rssj,la12,la13,laj12,laj13,
     1       mew,z1,z2,z3,sij,z,Di(6),dtidisrg,minoibd,wijk,
     2       th0i,rssi,h0,h1,zsd,up,low,zsd2,meanni,taui
      integer nselect
      external adist0,dtidisrg
      logical aws
      aws=lambda.lt.1e20
      h1=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)
      h0=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)/1.4
C  first fill predicted 
      low=sdcoef(1)+sdcoef(3)*sdcoef(2)
      up=sdcoef(1)+sdcoef(4)*sdcoef(2)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call sihat(th0(i1,i2,i3),D(1,i1,i2,i3),btb,
     1                    sipred(1,i1,i2,i3),nb)
            END DO
         END DO
      END DO
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE

               DO k=1,nb
                  z=siest(k,i1,i2,i3)
                  zsd=sdcoef(1)+z*sdcoef(2)
                  if(z.lt.sdcoef(3)) zsd=low
                  if(z.gt.sdcoef(4)) zsd=up
                  sig2bi(k)=zsd*zsd
                  vari(k)=varinv(k,i1,i2,i3)
               END DO
               rssi = dtidisrg(siest(1,i1,i2,i3),
     1                         sipred(1,i1,i2,i3),vari,nb)
               DO k=1,nb
                  sw(k)=0.d0
                  sw2(k)=0.d0
                  swsi(k)=0.d0
                  swsi2(k)=0.d0
               END DO
               th0i=th0(i1,i2,i3)
               th0n(i1,i2,i3)=th0i
C    used as initial values
               IF(ew(1).lt.0.d0) call dblepr("C0",2,ew,3)
               DO k=1,6
                  Di(k)=D(k,i1,i2,i3)
                  Dn(k,i1,i2,i3)=Di(k)
C    used as initial values
               END DO
               call eigen3(Di,ew,ev,ierr)
C produces largest EW in last position
               la12 = (ew(2)-ew(3))/ew(3)
               la13 = (ew(1)-ew(3))/ew(3)
C   now get bandwidth such that the ellopsoid has specified volume
               call gethani0(h0,h1,vol,vext,1.d-2,h)
               h2=h*h
C   create needed estimates of s_i
               nselect=0
               call rangex0(h,j1a,j1e,vext)
               call oriscor2(la12,la13,ev,nb,btb,oibd)
C               call oriscore(la12,la13,ev,nb,btb,oibd)
               meanni=0.d0
               minoibd=-1.d0
               l=0
               DO k=1,nb
                  if(s0ind(k)) CYCLE
                  meanni=meanni+vari(k)*sig2bi(k)
                  minoibd=max(minoibd,oibd(k))
                  l=l+1
               END DO
               meanni=meanni/l
C               meanni=sqrt(meanni/l)
C               meanni=exp(log(meanni/l)/3.d0)
C               taui=ani(i1,i2,i3)/tau
               taui=1.d0/tau
               DO k=1,nb
                  if(s0ind(k)) THEN
                     oibd(k)=1.d0
                  ELSE
C                     oibd(k)=exp((oibd(k)-minoibd)*meanni*taui)
C                     oibd(k)=exp((oibd(k)-minoibd)*taui)
                     oibd(k)=exp(oibd(k)*meanni*taui)
                  END IF
               END DO
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey0(j1,h,j2a,j2e,vext)
                  DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez0(j1,j2,h,j3a,j3e,vext)
                     DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE
                        wij=adist0(j1,j2,j3,vext)
C     triangular location kernel
                        if(wij.ge.h2) CYCLE
                        center=abs(i1-jj1)+abs(i2-jj2)+abs(i3-jj3)
                        wij = (1.d0 - wij/h2)
                        IF(aws) THEN
                           rssj = dtidisrg(siest(1,i1,i2,i3),
     1                               sipred(1,jj1,jj2,jj3),vari,nb)
                           sij = max(0.d0,rssj-rssi)/lambda
                           if(center.eq.0) sij=0.d0
                           if(sij.gt.1.d0) CYCLE
C  use Plateau kernel
                           if(sij.gt.0.5d0) THEN
                              wij=wij*2.d0*(1.d0-sij)
                           END IF
                           call eigen3(D(1,jj1,jj2,jj3),ew,ev,ierr)
C produces largest EW in last position
                           laj12 = (ew(2)-ew(3))/ew(3)
                           laj13 = (ew(1)-ew(3))/ew(3)
                           call oriscor2(laj12,laj13,ev,nb,btb,ojbd)
                        END IF
                        DO k=1,nb
                           if(center.eq.0) THEN
                              wijk=1.d0
                           ELSE 
                             IF(s0ind(k)) THEN
                                wijk=0.d0
                             ELSE
                             wijk=wij*oibd(k)*exp(ojbd(k)*meanni*taui)
C                             wijk=wijk*min(1.d0,laj12*laj12/la12*la12)
C  only use information from voxel that have "more" directional information
                             END IF
                           END IF
                           IF(wijk.gt.0.d0) THEN
                              z=si(k,jj1,jj2,jj3)
                              zsd=sdcoef(1)+z*sdcoef(2)
                              if(z.lt.sdcoef(3)) zsd=low
                              if(z.gt.sdcoef(4)) zsd=up
                              zsd2=zsd*zsd
                              wijk=wijk/zsd2
                              sw(k)=sw(k)+wijk
                              sw2(k)=sw2(k)+wijk*wijk*zsd2
                              swsi(k)=swsi(k)+wijk*si(k,jj1,jj2,jj3)
                           END IF
                        END DO
                    END DO
                  END DO
               END DO
               Do k=1,nb
                  swsi(k)=swsi(k)/sw(k)
               END DO
               DO k=1,nb
                  siest(k,i1,i2,i3)=swsi(k)
                  vari(k)=sw(k)*sw(k)/sw2(k)
                  varinv(k,i1,i2,i3)=vari(k)
               END DO
C               if(i1.eq.20.and.i2.eq.20.and.i3.eq.21) THEN
C                     call dblepr("oibd",4,oibd,nb)
C                     call dblepr("sw(k)",5,sw,nb)
C                     call dblepr("swsi(k)",7,swsi,nb)
C               END IF
               call dslvdti0(swsi,nb,btb,vari,th0n(i1,i2,i3),
     1                          Dn(1,i1,i2,i3),F,niter,eps,rss)
               call eigen3(Dn(1,i1,i2,i3),ew,ev,ierr)
               IF(ew(1).lt.0.d0) call dblepr("C0",2,ew,3)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   compute orientation score
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine oriscore(la12,la13,ev,nb,btb,oibd)
      implicit logical (a-z)
      integer nb
      real*8 btb(6,nb),oibd(nb),la12,la13,ev(3,2)
      integer k
      real*8 e1,e2,e3
      e1=ev(1,2)
      e2=ev(2,2)
      e3=ev(3,2)
      DO k=1,nb
         oibd(k)=la12*(e1*e1*btb(1,k)+e1*e2*btb(2,k)+
     1                 e1*e3*btb(3,k)+e2*e2*btb(4,k)+
     2                 e2*e3*btb(5,k)+e3*e3*btb(6,k))
      END DO
      e1=ev(1,1)
      e2=ev(2,1)
      e3=ev(3,1)
      DO k=1,nb
         oibd(k)=oibd(k)+la13*(e1*e1*btb(1,k)+e1*e2*btb(2,k)+
     1                         e1*e3*btb(3,k)+e2*e2*btb(4,k)+
     2                         e2*e3*btb(5,k)+e3*e3*btb(6,k))
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   compute orientation score
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine oriscor2(la12,la13,ev,nb,btb,oibd)
      implicit logical (a-z)
      integer nb
      real*8 btb(6,nb),oibd(nb),la12,la13,ev(3,2)
      integer k
      real*8 e1,e2,e3
      e1=ev(1,2)
      e2=ev(2,2)
      e3=ev(3,2)
      DO k=1,nb
         oibd(k)=la12*(e1*e1*btb(1,k)+e1*e2*btb(2,k)+
     1                 e1*e3*btb(3,k)+e2*e2*btb(4,k)+
     2                 e2*e3*btb(5,k)+e3*e3*btb(6,k))
      END DO
      e1=ev(1,1)
      e2=ev(2,1)
      e3=ev(3,1)
      DO k=1,nb
         oibd(k)=oibd(k)+la13*(e1*e1*btb(1,k)+e1*e2*btb(2,k)+
     1                         e1*e3*btb(3,k)+e2*e2*btb(4,k)+
     2                         e2*e3*btb(5,k)+e3*e3*btb(6,k))
         oibd(k)=oibd(k)*oibd(k)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for isotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangex0(h,ia,ie,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ia,ie
      real*8 h,vext(3)
      real*8 h1
      h1=h/vext(1)
      ia=-h1
      ie=h1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for isotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangey0(ix,h,ja,je,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ja,je,ix
      real*8 h,vext(3)
      real*8 z,x
      x=ix/h*vext(1)
      z=sqrt(1.d0-x*x)*h/vext(2)
      ja=-z
      je=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for isotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangez0(ix,jy,h,ka,ke,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ka,ke,ix,jy
      real*8 h,vext(3)
      real*8 x,y,z,h3
      x=ix/h*vext(1)
      y=jy/h*vext(2)
      h3=h/vext(3)
      z=sqrt(1.d0 - y*y - x*x)*h/vext(3)
      ka=-z
      ke=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for isometric geometry  and given 
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani0(x,y,value,vext,eps,bw)
      implicit logical(a-z)
      real*8 x,y,value,vext(3),eps,bw
      real*8 fw1,fw2,fw3,z
      real*8 sofw3D0
      external sofw3D0
      if(x.ge.y) RETURN
      fw1=sofw3D0(x,vext)
      fw2=sofw3D0(y,vext)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw3D0(x,vext)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw3D0(y,vext)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw3D0(z,vext)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
         call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END  
C  compute sum of weights for isotropic smoothing      
      real*8 function sofw3D0(bw,vext)
      implicit logical(a-z)
      real*8 bw,vext(3)
      integer ia1,ie1,ia2,ie2,ia3,ie3,i1,i2,i3
      real*8 wij,sw,h2,adist0
      external adist0
      h2=bw*bw
      sw=1.d0
      call rangex0(bw,ia1,ie1,vext)
      DO i1=1,ie1
         call rangey0(i1,bw,ia2,ie2,vext)
         DO i2=ia2,ie2
            call rangez0(i1,i2,bw,ia3,ie3,vext)
            DO i3=ia3,ie3
               wij=max(0.d0,1.d0-adist0(i1,i2,i3,vext)/h2)
               sw=sw+2.d0*wij
            END DO
         END DO
      END DO
C  now case i1=0
      call rangey0(0,bw,ia2,ie2,vext)
      DO i2=1,ie2
         call rangez0(0,i2,bw,ia3,ie3,vext)
         DO i3=ia3,ie3
            wij=max(0.d0,1.d0-adist0(0,i2,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      END DO
C  now case i1=i2=0
      call rangez0(0,0,bw,ia3,ie3,vext)
         DO i3=1,ie3
            wij=max(0.d0,1.d0-adist0(0,0,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      sofw3D0=sw
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute isotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function adist0(x,y,z,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer x,y,z
      real*8 xx,yy,zz,vext(3)
      xx=x*vext(1)
      yy=y*vext(2)
      zz=z*vext(3)
      adist0=xx*xx+yy*yy+zz*zz
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Nonlinear heteroskedastic regression
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dslvdti0(s,nb,b,varinv,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C  same as solvedti except that s is double
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb),D(6),b(6,nb),varinv(nb),th0,F(nb),eps
      integer i,j,k,info,iter,icount
      logical negdefin
      real*8 z,gamma,alpha,delta,xzvarinv,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,Dn(6),res,X(7),th0n
C  first check if D defines a positive definite densor
      call regularD(D,negdefin)
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      DO i=1,nb
         z=b(1,i)*D(1)
         DO j=2,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res*varinv(i)
         F(i)=res
      END DO
      th0n = th0
      nrss = rss
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO
         DO i=1,nb
            z=b(1,i)*D(1)
            DO j=2,6
               z=z+b(j,i)*D(j)
            END DO
            z=exp(-z)
            X(7)= -z
            z=z*th0
            DO j=1,6
               X(j)=b(j,i)*z
            END DO
            DO j=1,7
               xzvarinv=X(j)*varinv(i)
               dg(j)=dg(j)+xzvarinv*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+xzvarinv*X(k)
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
            call regularD(D,negdefin)
            if(negdefin) THEN
C  estimate using reparametrization
               call dslvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
            ENDIF
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         icount = 10
         DO WHILE (notacc.and.icount.gt.0)
            icount = icount-1 
            call rchkusr()
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
                  z=b(1,i)*Dn(1)
                  DO j=2,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=s(i)-th0n*exp(-z)
                  nrss=nrss+res*res*varinv(i)
                  F(i)=res
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
         th0=th0n
C  check if tensor is positive definite
         DO j=1,6
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      call regularD(D,negdefin)
      if(negdefin) THEN
C  estimate using reparametrization
         call dslvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
C      else
C         call testreg(D,4)
      ENDIF
      RETURN
      END
