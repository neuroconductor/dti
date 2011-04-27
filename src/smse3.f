      subroutine lkse3n(h,kappa,grad,grad2,grad3,ngrad,vext,ind,wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n)
      real*8 h,kappa,grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad),
     1       vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4,i4,mj1,mj2,mj3
      real*8 e2,e3,g1,g2,g3,f1,f2,f3,v1a,v2a,v3a,w1a,w2a,w3a,d2a,d3a,
     1    d1,d2,pD4,pD4a,kap4,kap2,h2,vd2,vd3,x1,x2,x3,n12,n13,nx1,h4
      ih1 = 1.0d1*h
      ih2 = 1.0d1*h/vext(1)
      ih3 = 1.0d1*h/vext(2)
      h2 = h*h
      h4 = h2*h2
      kap2 = kappa*kappa
      kap4 = kap2*kap2
      vd2 = vext(1)
      vd3 = vext(2)
      DO i4=1,ngrad
         g1=grad(1,i4)
         g2=grad(2,i4)
         g3=grad(3,i4)
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
         grad2(1,i4)=f1
         grad2(2,i4)=f2
         grad2(3,i4)=f3
C grad contains only non-zero components, first component is 0
         grad3(1,i4)=e2
         grad3(2,i4)=e3
      END DO
      i = 1
      mj1=0
      mj2=0
      mj3=0
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = vd2*j2
            DO j3 = -ih3,ih3
               x3 = vd3*j3
               DO j4 = 1,ngrad
                  nx1 = x1*grad(1,j4)+x2*grad(2,j4)+x3*grad(3,j4)
                  d2 = nx1*nx1
                  if(d2.ge.h2) CYCLE  
C   second term  <\delta_x,n_1(2)> to large
                  pD4 = d2*d2
                  DO i4 = 1,ngrad
                     n13 = grad(2,i4)*grad3(1,j4)+grad(3,i4)*grad3(2,j4)
                     d1 = n13*n13
                     pD4a = pD4 + d1*d1/kap4
                     if(pD4a.ge.h4) CYCLE  
C   first term  <n_1(1),n_3(2)> to large for remainder
                     n12 = grad(1,i4)*grad2(1,j4)+
     1                     grad(2,i4)*grad2(2,j4)+grad(3,i4)*grad2(3,j4)
                     d1 = n12*n12
                     pD4a = pD4a + d1*d1/kap4
                     if(pD4a.ge.h4) CYCLE  
C   first term  <n_1(1),n_2(2)> to large for remainder
C
C   now third term if something is left
C
                     w1a=0.5d0*n12*grad(1,j4)
                     w2a=grad3(1,j4)+0.5d0*n12*grad(2,j4)
                     w3a=grad3(2,j4)+0.5d0*n12*grad(3,j4)
                     d3a=x1*w1a+x2*w2a+x3*w3a
                     pD4a = pD4a + d3a*d3a/kap2
                     if(pD4a.ge.h4) CYCLE  
C   third term  <\delta_x,n_3(2)+1/2 n12 n_1(2)> to large for remainder                     
                     v1a=grad2(1,j4)+0.5d0*n13*grad(1,j4)
                     v2a=grad2(2,j4)+0.5d0*n13*grad(2,j4)
                     v3a=grad2(3,j4)+0.5d0*n13*grad(3,j4)
                     d2a=x1*v1a+x2*v2a+x3*v3a
                     pD4a = pD4a + d2a*d2a/kap2
C   third term  <\delta_x,n_2(2)+1/2 n13 n_1(2)> to large for remainder                     
                     if(pD4a.ge.h4) CYCLE  
                     wght(i)= 1.d0-pD4a/h4
                     ind(1,i) = j1
                     ind(2,i) = j2
                     ind(3,i) = j3
                     ind(4,i) = i4
                     ind(5,i) = j4
               if(iabs(j3).gt.mj3) mj3=iabs(j3)
               if(iabs(j2).gt.mj2) mj2=iabs(j2)
               if(j1.gt.mj1) mj1=j1
                     if(i.ge.n) THEN
                        call intpr("Reached max i",13,i,1)
                        return
                     END IF
                     i = i+1
                  END DO
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      call intpr("ih1",3,ih1,1)
      call intpr("mj1",3,mj1,1)
      call intpr("ih2",3,ih2,1)
      call intpr("mj2",3,mj2,1)
      call intpr("ih3",3,ih3,1)
      call intpr("mj3",3,mj3,1)
      RETURN
      END
      subroutine lkse3(h,kappa,grad,grad2,grad3,ngrad,vext,ind,wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n)
      real*8 h,kappa,grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad),
     1       vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4,i4
      real*8 e2,e3,g1,g2,g3,f1,f2,f3,hk,hk2,v1a,v2a,v3a,w1a,w2a,w3a,
     1       v1b,v2b,v3b,w1b,w2b,w3b,d2a,d2b,d3a,d3b,d1,d2,d3
      real*8 h2,vd2,vd3,x1,x2,x3,kofd1,kofd2,kofd3,n12,n13,nx1,twoh2
      ih1 = h
      ih2 = h/vext(1)
      ih3 = h/vext(2)
      h2 = h*h
      twoh2 = 2.d0*h2
      hk = h2*kappa
      hk2 = hk*kappa
      vd2 = vext(1)
      vd3 = vext(2)
      DO i4=1,ngrad
         g1=grad(1,i4)
         g2=grad(2,i4)
         g3=grad(3,i4)
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
         grad2(1,i4)=f1
         grad2(2,i4)=f2
         grad2(3,i4)=f3
C grad contains only non-zero components, first component is 0
         grad3(1,i4)=e2
         grad3(2,i4)=e3
      END DO
      i = 1
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = vd2*j2
            DO j3 = -ih3,ih3
               x3 = vd3*j3
               DO j4 = 1,ngrad
                  nx1 = x1*grad(1,j4)+x2*grad(2,j4)+x3*grad(3,j4)
                  d2 = nx1*nx1
                  if(d2.ge.twoh2) CYCLE  
                  kofd2 = 1.d0-d2/twoh2
                  DO i4 = 1,ngrad
                     n13 = grad(2,i4)*grad3(1,j4)+grad(3,i4)*grad3(2,j4)
                     d1 = n13*n13
                     if(d1.ge.hk2) CYCLE  
                     n12 = grad(1,i4)*grad2(1,j4)+
     1                     grad(2,i4)*grad2(2,j4)+grad(3,i4)*grad2(3,j4)
                     d1 = d1+n12*n12
                     if(d1.ge.hk2) CYCLE                       
                     kofd1 = 1.d0-d1/hk2
                     if(kofd1.le.0.d0) CYCLE  
                     w1a=-0.25d0*grad(1,j4)
                     w2a=grad3(1,j4)-0.25d0*grad(2,j4)
                     w3a=grad3(2,j4)-0.25d0*grad(3,j4)
                     w1b=0.25d0*grad(1,j4)
                     w2b=grad3(1,j4)+0.25d0*grad(2,j4)
                     w3b=grad3(2,j4)+0.25d0*grad(3,j4)
                     d3a=abs(x1*w1a+x2*w2a+x3*w3a)
                     d3b=abs(x1*w1b+x2*w2b+x3*w3b)
                     d3=min(d3a,d3b)
                     if(d3.ge.hk) CYCLE                       
                     v1a=grad2(1,j4)-0.25d0*grad(1,j4)
                     v2a=grad2(2,j4)-0.25d0*grad(2,j4)
                     v3a=grad2(3,j4)-0.25d0*grad(3,j4)
                     v1b=grad2(1,j4)+0.25d0*grad(1,j4)
                     v2b=grad2(2,j4)+0.25d0*grad(2,j4)
                     v3b=grad2(3,j4)+0.25d0*grad(3,j4)
                     d2a=abs(x1*v1a+x2*v2a+x3*v3a)
                     d2b=abs(x1*v1b+x2*v2b+x3*v3b)
                     d3=d3+min(d2a,d2b)
C  this differs from definition 2.3     
                     if(d3.ge.hk) CYCLE  
                     kofd3 = 1.d0-d3/hk
                     wght(i)=kofd1*kofd2*kofd3
                     ind(1,i) = j1
                     ind(2,i) = j2
                     ind(3,i) = j3
                     ind(4,i) = i4
                     ind(5,i) = j4
                     if(i.ge.n) THEN
                        call intpr("Reached max i",13,i,1)
                        return
                     END IF
                     i = i+1
                  END DO
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine lockse3(h,kappa,grad,ngrad,vext,ind,wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n)
      real*8 h,kappa,grad(3,ngrad),vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4,j5
      real*8 h2,d2,d3,z,x1,x2,x3,distx,k1ofx,k2ofxphi,k3ofphi,kap
      kap=kappa/h
      ih1 = h
      ih2 = h/vext(1)
      ih3 = h/vext(2)
      h2 = h*h
      d2 = vext(1)
      d3 = vext(2)
      i = 1
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = d2*j2
            DO j3 = -ih3,ih3
               x3 = d3*j3
               distx = j1*j1+x2*x2+x3*x3
               if(distx.ge.h2) CYCLE
               k1ofx = max(0.d0,1.d0-distx/h2)
               DO j4 = 1,ngrad
                  z = x1*grad(1,j4)+x2*grad(2,j4)+x3*grad(3,j4)
                  if(distx.gt.0.d0) THEN
                  k2ofxphi = max(0.d0,1.d0-kap+kap*z*z/distx)
                  ELSE
                  k2ofxphi = 1.d0
                  END IF
                  if(k2ofxphi.le.0.d0) CYCLE
                  DO j5 = 1,ngrad
                     z = grad(1,j4)*grad(1,j5)+grad(2,j4)*grad(2,j5)+
     1                   grad(3,j4)*grad(3,j5)
                     k3ofphi = max(0.d0,1.d0-kap+kap*z*z)
                     if(k3ofphi.le.0.d0) CYCLE 
                     wght(i)=k1ofx*k2ofxphi*k3ofphi
                     ind(1,i) = j1
                     ind(2,i) = j2
                     ind(3,i) = j3
                     ind(4,i) = j4
                     ind(5,i) = j5
                     if(i.ge.n) THEN
                        call intpr("Reached max i",13,i,1)
                        return
                     END IF
                     i = i+1
                  END DO
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine adasmse3(y,th,ni,mask,n1,n2,n3,ngrad,lambda,ind,w,n,
     1                    thn,r,s2inv,sw,swy,swy2,ncoil)
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C   s2inv containes the inverses of estimated variances
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,n,ind(5,n),ncoil
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),th(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1       lambda,w(n),thn(n1,n2,n3,ngrad),sw(ngrad),swy(ngrad),
     2       swy2(ngrad),s2inv(n1,n2,n3,ngrad),r(n1,n2,n3,ngrad)
      integer i,i1,i2,i3,i4,j1,j2,j3,j4
      real*8 z,yj,swyi,nii
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               DO i4=1,ngrad
                  sw(i4)=0.d0
                  swy(i4)=0.d0
                  swy2(i4)=0.d0
               END DO
               DO i=1,n
                  j1=i1+ind(1,i)
                  j2=i2+ind(2,i)
                  j3=i3+ind(3,i)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  if(.not.mask(j1,j2,j3)) CYCLE          
                  i4=ind(4,i)
                  j4=ind(5,i)
C                  z=th(i1,i2,i3,i4)-th(j1,j2,j3,j4)
C                  z=z*z*ni(i1,i2,i3,i4)/lambda
                  call spenalty(th(i1,i2,i3,i4),th(j1,j2,j3,j4),
     1                          s2inv(i1,i2,i3,i4),ncoil,z)
                  z=z*ni(i1,i2,i3,i4)/lambda
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
                  i4=ind(4,i)
                  j4=ind(5,i)
                  z=th(i1,i2,i3,i4)-th(j1,j2,j3,j4)
                  z=z*z*ni(i1,i2,i3,i4)/lambda
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
      integer ncoil
      real*8 thi,thj,s2inv,pen
      real*8 z
      if(ncoil.eq.0) THEN
         z=thi-thj
         pen=z*z
         RETURN
      END IF
      RETURN
      END
      