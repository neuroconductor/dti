C
C  compute M(c) in Jiang Cheng (2012)
C
      subroutine Mofc(c,kern,nk,ng,e,l,fvmofc)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),fvmofc
C diagonal matrix in l
      integer i,j,k
      real*8 z,ci,ci2,zz
C  first the penalty term
      DO i=1,ng
         ci=c(i)
         z=z+ci*ci*l(i)
      END DO
      DO k=1,ng
         zz=c(1)*c(1)*kern(1,1,k)
         DO i=2,ng
            ci=c(i)
            zz=zz+ci*ci*kern(i,i,k)
            ci2=2.d0*ci
            Do j=1,i-1
               zz=zz+ci2*c(j)*kern(i,j,k)
            END DO
         END DO
         zz=zz-e(k)
         z=z+zz*zz
      END DO
      fvmofc=z/2.d0
      RETURN
      END
C
C  compute dM(c)/dc in Jiang Cheng (2012)
C
      subroutine Dmofcdc(c,kern,nk,ng,e,l,w,mofcdc)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),mofcdc(nk)
C diagonal matrix in l
      integer k
      real*8 z,ddot
      external ddot
C  first the penalty term  Lambda%*%c
      DO k=1,ng
         mofcdc(k) = c(k)* l(k)
      END DO
C  now the sum
      DO k=1,ng
C  compute kern(,,i)%*%c  in w
         call DGEMV("N",nk,nk,1.d0,kern(1,1,k),1,c,1,0.d0,w,1)
C  compute   2*(c^T%*%kern(,,i)%*%c-e(i)) in z
         z= 2.d0*(DDOT(nk,c,1,w,1)-e(k))
C  add 2*(c^T%*%kern(,,i)%*%c-e(i))*kern(,,i)%*%c to mofcdc
         call DAXPY(nk,z,w,1,mofcdc,1)
      END DO
      RETURN
      END
C
C  compute nabla M(c) in Jiang Cheng (2012)
C
      subroutine nablmofc(c,kern,nk,ng,e,l,w,nablam)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),nablam(nk)
      real*8 z,ddot
      external ddot
C  first get  dM(c)/dc in   nablam
      call Dmofcdc(c,kern,nk,ng,e,l,w,nablam)
C  -c^T%*%dM(c)/dc in z
      z = -DDOT(nk,c,1,nablam,1)
      call DAXPY(nk,z,c,1,nablam,1)
      RETURN
      END
C
C  compute Exp_c(dt*v)  in Jiang Cheng (2012)
C
      subroutine expcv(c,v,nk,dt,r) 
      implicit logical (a-z)
      integer nk
      real*8 v(nk),c(nk),r(nk),dt
      integer i
      real*8 sv,cv,vn,dnrm2
      external dnrm2
      vn = dnrm2(nk,v,1)
      cv = dcos(dt*vn)
      sv = dsin(dt*vn)/vn/dt
      DO i=1,nk
         r(i) = c(i)*cv+dt*v(i)*sv
      END DO
      RETURN
      END
C
C  compute c^(k+!) from c^k in Jiang Cheng (2012)
C
      subroutine getnewck(c,kern,nk,ng,e,l,w,nablam,ck1)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),nablam(nk),
     1       ck1(nk)
      real*8 dt0,dt1,dt2,dnrm2,m0,mdt0,mdt1,mdt2,dtn,mdtn
      external dnrm2
C first get nabla M(c) in nablam
      call nablmofc(c,kern,nk,ng,e,l,w,nablam)
C get norm of it 
      dt0 = dnrm2(nk,nablam,1)
C perform line search 
C first M(c) in m0
      call Mofc(c,kern,nk,ng,e,l,m0)
C now M(dt*v/||v||) in mdt
      call expcv(c,nablam,nk,dt0,ck1) 
      call Mofc(ck1,kern,nk,ng,e,l,mdt0)
      winit = .TRUE.
      wextend = .TRUE.
      wndone = .TRUE.
C try to find valid start
      DO while (winit) 
         IF(mdt0.ge.m0) THEN
            dt0 = dt0*.6
            call expcv(c,nablam,nk,dt0,ck1) 
            call Mofc(ck1,kern,nk,ng,e,l,mdt0)
         ELSE
            winit = .FALSE.
         END IF
         if(dt0.lt.1d-6) THEN
C no descent found, keep c(nk)
            winit = .FALSE.
         END IF
      END DO
C best value in dt0, mdt0
C try to get larger stepsize
      DO while (wextend)
         dt1 = 1.5*dt0
         call expcv(c,nablam,nk,dt1,ck1) 
         call Mofc(ck1,kern,nk,ng,e,l,mdt1)
         if(mdt1.lt.mdt0) THEN
            dt0 = dt1
            mdt0 = mdt1
         ELSE
            wextend = .FALSE.
         ENDIF
      END DO
      dt2 = dt0*.8
      call expcv(c,nablam,nk,dt2,ck1) 
      call Mofc(ck1,kern,nk,ng,e,l,mdt2)
      DO while (mdt2.lt.mdt0)
         dt1 = dt0
         mdt1 = mdt0
         dt0 = dt2
         mdt0 =mdt2
         dt2 = dt2*.8
         call expcv(c,nablam,nk,dt2,ck1) 
         call Mofc(ck1,kern,nk,ng,e,l,mdt2)          
      END DO
C now we have dt1 > dt0 > dt2 and mdt1 > mdt0 < mdt2
      DO while (dt2/dt1.gt.0.95d0)
         dtn = 0.5*(dt1+dt0)
         call expcv(c,nablam,nk,dtn,ck1) 
         call Mofc(ck1,kern,nk,ng,e,l,mdtn)
         if(mdtn.lt.mdt0) THEN
            dt2=dt0
            mdt2=mdt0
            dt0=dtn
            mdt0=mdtn
         ELSE
            dt1=dtn
            mdt1=mdtn
         END IF
         dtn = 0.5*(dt2+dt0)
         call expcv(c,nablam,nk,dtn,ck1) 
         call Mofc(ck1,kern,nk,ng,e,l,mdtn)
         if(mdtn.lt.mdt0) THEN
            dt1=dt0
            mdt1=mdt0
            dt0=dtn
            mdt0=mdtn
         ELSE
            dt2=dtn
            mdt2=mdtn
         END IF
      END DO
C now get new c in ck1
      call expcv(c,nablam,nk,dt0,ck1) 
      RETURN
      END
C
C   estimate coefficients
C      
      subroutine sqrteap(ei,kern,nk,ng,n1,n2,n3,l,w,nablam,ck,ck1,
     1                   cres)
      implicit logical (a-z)
      integer nk,ng,n1,n2,n3
      real*8 ck(nk),kern(nk,nk,ng),ei(ng,n1,n2,n3),l(nk),w(nk),
     1       nablam(nk),ck1(nk),cres(nk,n1,n2,n3)
      integer i1,i2,i3,k
      real*8 mck,mck1
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               ck(1)=1.d0
               DO k=2,nk
                  ck(k)=0.d0
               END DO
               call Mofc(ck,kern,nk,ng,ei(1,i1,i2,i3),l,mck)
               ndone=.TRUE.
               DO WHILE (ndone)
                  call getnewck(ck,kern,nk,ng,ei(1,i1,i2,i3),l,w,
     1                          nablam,ck1)
                  call Mofc(ck1,kern,nk,ng,ei(1,i1,i2,i3),l,mck1)
                  call DCOPY(nk,ck1,1,ck,1)
                  ndone = (mck-mck1)/mck.le.1e-5
               END DO
               call DCOPY(nk,ck,1,cres(1,i1,i2,i3),1)
            END DO
         END DO
      END DO
      RETURN
      END
      
      