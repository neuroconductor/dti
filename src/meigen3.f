      subroutine eigen3(ddiff, diff, evectors, evalues, evaonly, ierr)
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
      subroutine estim_dt(theta,try,h,n1,n2,n3,dt)
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
      real*8 theta(n1,n2,n3,3),try(n1,n2,n3),h,dt(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j2a,j3a,j1e,j2e,j3e
      real*8 sw,swy,thi1,thi2,thi3,p11,p12,p13,p22,p23,p33,z,
     1     d,d1,d2,d3  
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       swy=0.d0
	       thi1=theta(i1,i2,i3,1)
	       thi2=theta(i1,i2,i3,2)
	       thi3=theta(i1,i2,i3,3)
	       p11=(1-thi1*thi1)/h
	       p22=(1-thi2*thi2)/h
	       p33=(1-thi3*thi3)/h
	       p12=-thi1*thi2/h
	       p13=-thi1*thi3/h
	       p23=-thi2*thi3/h
	       DO j1=1,n1
	          d1=i1-j1
	          DO j2=1,n2
		     d2=i2-j2
		     DO j3=1,n3
		        d3=i3-j3
			d=d1*p11+d2*p12+d3*p13
			z=d*d
			if(z.ge.1.d0) CYCLE
			d=d1*p12+d2*p22+d3*p23
			z=z+d*d
			if(z.ge.1.d0) CYCLE
			d=d1*p13+d2*p23+d3*p33
			z=z+d*d
			if(z.ge.1.d0) CYCLE
			z=1.d0-z
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