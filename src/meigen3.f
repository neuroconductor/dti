CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute the largest eigenvalue (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3r(y,lambda,theta,ierr)
      integer ierr
      real*8 y(6),a(3,3),lambda,theta(3)
      integer i,j,l,ISUPPZ(2),lwork,iwork(50),liwork
      real*8 w(3),work(104)
      l=1
      DO i=1,3
         DO j=i,3
	    a(i,j)=y(l)
         END DO
	 l=l+1
      END DO
      lwork=104
      liwork=50
      call dsyevr('V','I','U',3,a,3,vl,vu,3,3,1.d-10,l,w,
     1            theta,3,ISUPPZ,work,lwork,iwork,liwork,ierr)
      if(ierr.eq.0) THEN
         lambda=w(1)
      ELSE
         lambda=0.d0
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3(y,lambda,theta,ierr)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(2),lwork,iwork(50),liwork
      real*8 work(104)
      l=1
      DO i=1,3
         DO j=i,3
	    a(i,j)=y(l)
         END DO
	 l=l+1
      END DO
      lwork=104
      liwork=50
      call dsyevr('V','A','U',3,a,3,vl,vu,1,3,1.d-10,l,lambda,
     1            theta,3,ISUPPZ,work,lwork,iwork,liwork,ierr)
      RETURN
      END
      
      subroutine eig3old(ddiff, diff, evectors, evalues, evaonly, ierr)
      integer ddiff, ierr(ddiff)
      real*8 diff(ddiff,6), evectors(9,ddiff), evalues(3,ddiff)
      logical evaonly
      real*8 xmat(9),work1(3),work2(3)
      integer i,n
      n=3
      DO i=1,ddiff
         xmat(1)=diff(i,1)
         xmat(2)=diff(i,4)
         xmat(3)=diff(i,5)
         xmat(4)=diff(i,4)
         xmat(5)=diff(i,2)
         xmat(6)=diff(i,6)
         xmat(7)=diff(i,5)
         xmat(8)=diff(i,6)
         xmat(9)=diff(i,3)
	 call rs(n, n, xmat, evalues(1,i), evaonly, evectors(1,i), 
     1            work1, work2, ierr(i)) 
      END DO 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Estimate local flow intensity
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine estimdt(theta,try,hd,ht,n1,n2,n3,dt)
C
C   uses rotation symmetric Epanechnicov kernel
C
C
C   theta    array of dominant directions
C   try      the trace of the tensor
C   h        bandwidth
C   n1,n2,n3 dimensions of the grid
C   dt       estimated flow intensity
C
      integer n1,n2,n3
      real*8 theta(n1,n2,n3,3),try(n1,n2,n3),hd,ht,dt(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j2a,j3a,j1e,j2e,j3e,ih1
      real*8 sw,swy,thi1,thi2,thi3,p11,p12,p13,p22,p23,p33,zd,zt,z,
     1     d,d1,d2,d3,ht2
      ht2=ht*ht
      ih1=dsqrt(2*hd*hd+ht2)  
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       swy=0.d0
	       thi1=theta(i1,i2,i3,1)
	       thi2=theta(i1,i2,i3,2)
	       thi3=theta(i1,i2,i3,3)
	       p11=(1-thi1*thi1)/hd
	       p22=(1-thi2*thi2)/hd
	       p33=(1-thi3*thi3)/hd
	       p12=-thi1*thi2/hd
	       p13=-thi1*thi3/hd
	       p23=-thi2*thi3/hd
	       j1a=max0(1,i1-ih1)
	       j1e=min0(n1,i1+ih1)
	       DO j1=j1a,j1e
	          d1=i1-j1
	          DO j2=1,n2
		     d2=i2-j2
		     DO j3=1,n3
		        d3=i3-j3
			d=d1*thi1+d2*thi2+d3*thi3
			zt=d*d/ht2
			if(zt.ge.1.d0) CYCLE
			d=d1*p11+d2*p12+d3*p13
			zd=d*d
			if(zd.ge.1.d0) CYCLE
			d=d1*p12+d2*p22+d3*p23
			zd=zd+d*d
			if(zd.ge.1.d0) CYCLE
			d=d1*p13+d2*p23+d3*p33
			zd=zd+d*d
			if(zd.ge.1.d0) CYCLE
			z=(1.d0-zd)*(1.d0-zt)
			sw=sw+z
			swy=swy+z*try(j1,j2,j3)
        	     END DO
	          END DO
               END DO
	       dt(i1,i2,i3)=swy/sw
	    END DO
	 END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Estimate local flow intensity
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine esttheta(y,dt,n1,n2,n3,h,theta)
C
C   uses rotation symmetric Epanechnicov kernel
C
C
C   theta    array of dominant directions
C   try      the trace of the tensor
C   h        bandwidth
C   n1,n2,n3 dimensions of the grid
C   dt       estimated flow intensity
C
      integer n1,n2,n3
      real*8 theta(n1,n2,n3,3),y(n1,n2,n3,6),h,dt(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j2a,j3a,j1e,j2e,j3e,ih1,ih2,ih3
      real*8 h2,swy(6),lambda
      h2=h*h 
      ih1=h
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       DO k=1,6
	          swy(k)=0.d0
	       END DO
	       j1a=max0(1,i1-ih1)
	       j1e=min0(n1,i1+ih1)
	       DO j1=j1a,j1e
	          d1=i1-j1
		  z1=d1*d1
		  ih2=dsqrt(h2-z1)
	          j2a=max0(1,i2-ih2)
	          j2e=min0(n1,i2+ih2)
	          DO j2=j2a,j2e2
		     d2=i2-j2
		     z2=d2*d2
		     ih3=dsqrt(h2-z1-z2)
	             j3a=max0(1,i3-ih3)
	             j3e=min0(n1,i3+ih3)
		     DO j3=j3a,j3e
		        d3=i3-j3
			z=1.d0-(z1+z2+d3*d3)/h2
			z=z*dt(j1,j2,j3)
			sw=sw+z
			DO k=1,6
			   swy(k)=swy(k)+z*y(j1,j2,j3,k)
			END DO
        	     END DO
	          END DO
               END DO
	       DO k=1,6
	          swy(k)=swy(k)/sw
	       END DO
               call eigen3r(swy,lambda,theta(i1,i2,i3,1),ierr)
	    END DO
	 END DO
      END DO
      RETURN
      END