      subroutine lkse3(h,kappa,grad,grad2,grad3,ngrad,vext,ind,wght,n)
      implicit logical (a-z)
      integer ngrad,n,ind(5,n)
      real*8 h,kappa,grad(3,ngrad),grad2(3,ngrad),grad3(2,ngrad),
     1       vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4,i4
      real*8 e2,e3,g1,g2,g3,f1,f2,f3,hk,hk2,v1a,v2a,v3a,w1a,w2a,w3a,
     1       v1b,v2b,v3b,w1b,w2b,w3b,w12,w13,d2a,d2b,d3a,d3b
      real*8 h2,d2,d3,x1,x2,x3,dist1,dist2,dist3,n12,n13,nx1
      ih1 = h
      ih2 = h/vext(1)
      ih3 = h/vext(2)
      h2 = h*h
      hk = h2*kappa
      hk2 = hk*kappa
      d2 = vext(1)
      d3 = vext(2)
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
         n12=grad(1,i4)*grad2(1,i4)+
     1       grad(2,i4)*grad2(2,i4)+grad(3,i4)*grad2(3,i4)
         n13 = grad(2,i4)*grad3(1,i4)+grad(3,i4)*grad3(2,i4)
      END DO
      i = 1
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = d2*j2
            DO j3 = -ih3,ih3
               x3 = d3*j3
               DO j4 = 1,ngrad
                  nx1 = x1*grad(1,j4)+x2*grad(2,j4)+x3*grad(3,j4)
                  dist2 = max(0.d0,1.d0-nx1*nx1/2.d0/h2)
                  if(dist2.le.0.d0) CYCLE  
                  DO i4 = 1,ngrad
                     n12 = grad(1,i4)*grad2(1,j4)+
     1                   grad(2,i4)*grad2(2,j4)+grad(3,i4)*grad2(3,j4)
                     n13 = grad(2,i4)*grad3(1,j4)+grad(3,i4)*grad3(2,j4)
                     dist1 = max(0.d0,1.d0-(n12*n12+n13*n13)/hk2)
                     if(dist1.le.0.d0) CYCLE  
                     w12 = 0.25d0*w12
                     v1a=grad2(1,j4)-w12*grad(1,j4)
                     v2a=grad2(2,j4)-w12*grad(2,j4)
                     v3a=grad2(3,j4)-w12*grad(3,j4)
                     v1b=grad2(1,j4)+w12*grad(1,j4)
                     v2b=grad2(2,j4)+w12*grad(2,j4)
                     v3b=grad2(3,j4)+w12*grad(3,j4)
                     w13=0.25d0*w13
                     w1a=-w13*grad(1,j4)
                     w2a=grad3(1,j4)-w13*grad(2,j4)
                     w3a=grad3(2,j4)-w13*grad(3,j4)
                     w1b=w13*grad(1,j4)
                     w2b=grad3(1,j4)+w13*grad(2,j4)
                     w3b=grad3(2,j4)+w13*grad(3,j4)
                     d2a=abs(x1*v1a+x2*v2a+x3*v3a)
                     d2b=abs(x1*v1b+x2*v2b+x3*v3b)
                     d3a=abs(x1*w1a+x2*w2a+x3*w3a)
                     d3b=abs(x1*w1b+x2*w2b+x3*w3b)
                     dist3=min(d2a,d2b)+min(d3a,d3b)
C  this differs from definition 2.3     
                     dist3 = max(0.d0,1.d0-dist3/hk)
                     if(dist3.le.0.d0) CYCLE  
                     wght(i)=dist1*dist2*dist3
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
                  z=th(i1,i2,i3,i4)-th(j1,j2,j3,j4)
                  z=z*z*ni(i1,i2,i3,i4)/lambda
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