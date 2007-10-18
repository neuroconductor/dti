CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute the largest eigenvalue (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3r(y,lambda,theta,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda,theta(3)
      integer i,j,l,ISUPPZ(2),lwork,iwork(100),liwork,n,m
      real*8 w(3),work(204),vl,vu,eps
      n=3
      m=1
      eps=1.d-20
      l=1
      DO i=1,3
         DO j=i,3
	    a(i,j)=y(l)
	    l=l+1
         END DO
      END DO
      lwork=204
      liwork=100
      call dsyevr('V','I','U',n,a,n,vl,vu,n,n,eps,m,w,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
      lambda=w(1)
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
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 theta(3,n1,n2,n3),try(n1,n2,n3),hd,ht,dt(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j1e,ih1
      real*8 sw,swy,thi1,thi2,thi3,p11,p12,p13,p22,p23,p33,zd,zt,z,
     1     d,d1,d2,d3,ht2
      ht2=ht*ht
      ih1=dsqrt(2*hd*hd+ht2)  
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       swy=0.d0
	       thi1=theta(1,i1,i2,i3)
	       thi2=theta(2,i1,i2,i3)
	       thi3=theta(3,i1,i2,i3)
	       p11=(1-thi1*thi1)/hd
	       p22=(1-thi2*thi2)/hd
	       p33=(1-thi3*thi3)/hd
	       p12=-thi1*thi2/hd
	       p13=-thi1*thi3/hd
	       p23=-thi2*thi3/hd
	       j1a=max(1,i1-ih1)
	       j1e=min(n1,i1+ih1)
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
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 theta(3,n1,n2,n3),y(6,n1,n2,n3),h,dt(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j2a,j3a,j1e,j2e,j3e,ih1,ih2,ih3,
     1        k,ierr,ind(3)
      real*8 h2,swy(6),lambda,sw,d1,d2,d3,z1,z2,z
      h2=h*h 
      ih1=h
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       DO k=1,6
	          swy(k)=0.d0
	       END DO
	       j1a=max(1,i1-ih1)
	       j1e=min(n1,i1+ih1)
	       DO j1=j1a,j1e
	          d1=i1-j1
		  z1=d1*d1
		  ih2=sqrt(h2-z1)
	          j2a=max(1,i2-ih2)
	          j2e=min(n2,i2+ih2)
	          DO j2=j2a,j2e
		     d2=i2-j2
		     z2=d2*d2
		     ih3=sqrt(h2-z1-z2)
	             j3a=max(1,i3-ih3)
	             j3e=min(n3,i3+ih3)
		     DO j3=j3a,j3e
		        d3=i3-j3
			z=1.d0-(z1+z2+d3*d3)/h2
			z=z*dt(j1,j2,j3)
			sw=sw+z
			DO k=1,6
			   swy(k)=swy(k)+z*y(k,j1,j2,j3)
			END DO
        	     END DO
	          END DO
               END DO
	       DO k=1,6
	          swy(k)=swy(k)/sw
	       END DO
               call eigen3r(swy,lambda,theta(1,i1,i2,i3),ierr)
	       if(ierr.gt.0) THEN
	          ind(1)=i1
		  ind(2)=i2
		  ind(3)=i3
		  call intpr("ierr>0 in",9,ind,3)
	       END IF
C   largest eigenvalue in lambda, corresponding EV in theta(1,i1,i2,i3)
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
      subroutine estimdt2(theta,try,ai,lm,hd,ht,n1,n2,n3,dt)
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
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 theta(3,n1,n2,n3),try(n1,n2,n3),hd,ht,dt(n1,n2,n3),
     1       ai(n1,n2,n3),lm(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j1e,ih1
      real*8 sw,swy,thi1,thi2,thi3,p11,p12,p13,p22,p23,p33,zd,zt,z,
     1     d,d1,d2,d3,ht2
      ht2=ht*ht
      ih1=sqrt(2*hd*hd+ht2)  
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       swy=0.d0
	       thi1=theta(1,i1,i2,i3)
	       thi2=theta(2,i1,i2,i3)
	       thi3=theta(3,i1,i2,i3)
	       p11=(1-thi1*thi1)/hd
	       p22=(1-thi2*thi2)/hd
	       p33=(1-thi3*thi3)/hd
	       p12=-thi1*thi2/hd
	       p13=-thi1*thi3/hd
	       p23=-thi2*thi3/hd
	       j1a=max(1,i1-ih1)
	       j1e=min(n1,i1+ih1)
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
C		swy=swy+z*try(j1,j2,j3)*ai(j1,j2,j3)/lm(j1,j2,j3)
		swy=swy+z*ai(j1,j2,j3)
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
      subroutine esttens(y,dt,n1,n2,n3,h,theta,yhat,aihat,lmhat)
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
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 theta(3,n1,n2,n3),y(6,n1,n2,n3),h,dt(n1,n2,n3),
     1    aihat(n1,n2,n3),lmhat(n1,n2,n3),yhat(6,n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,j1a,j2a,j3a,j1e,j2e,j3e,ih1,ih2,ih3,
     1        k,ierr,ind(3),j
      real*8 h2,swy(6),lambda(3),sw,d1,d2,d3,z1,z2,z,th(3,3),lmean
      h2=h*h 
      ih1=h
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       sw=0.d0
	       DO k=1,6
	          swy(k)=0.d0
	       END DO
	       j1a=max(1,i1-ih1)
	       j1e=min(n1,i1+ih1)
	       DO j1=j1a,j1e
	          d1=i1-j1
		  z1=d1*d1
		  ih2=sqrt(h2-z1)
	          j2a=max(1,i2-ih2)
	          j2e=min(n2,i2+ih2)
	          DO j2=j2a,j2e
		     d2=i2-j2
		     z2=d2*d2
		     ih3=sqrt(h2-z1-z2)
	             j3a=max(1,i3-ih3)
	             j3e=min(n3,i3+ih3)
		     DO j3=j3a,j3e
		        d3=i3-j3
			z=1.d0-(z1+z2+d3*d3)/h2
			z=z*dt(j1,j2,j3)
			sw=sw+z
			DO k=1,6
			   swy(k)=swy(k)+z*y(k,j1,j2,j3)
			END DO
        	     END DO
	          END DO
               END DO
	       DO k=1,6
	          swy(k)=swy(k)/sw
	       END DO
               call eigen3(swy,lambda,th,ierr)
C               call eigen3(swy,lambda,theta(1,i1,i2,i3),ierr)
	       if(ierr.gt.0) THEN
	          ind(1)=i1
		  ind(2)=i2
		  ind(3)=i3
		  call intpr("ierr>0 in",9,ind,3)
	       END IF
	       DO j=1,3
	          theta(j,i1,i2,i3)=th(j,3)
	       END DO  
	       DO j=1,6
	          yhat(j,i1,i2,i3)=swy(j)
	       END DO  
	       lmhat(i1,i2,i3)=lambda(3)
	       lmean=(lambda(1)+lambda(2)+lambda(3))/3.d0
	       z1=(lambda(1)-lmean)
	       z=z1*z1
	       z1=(lambda(2)-lmean)
	       z=z+z1*z1
	       z1=(lambda(3)-lmean)
	       z=z+z1*z1
	       z=1.5d0*z/(lambda(1)*lambda(1)+lambda(2)*lambda(2)+
     1                    lambda(3)*lambda(3))
	       aihat(i1,i2,i3)=sqrt(z)
C   largest eigenvalue in lambda, corresponding EV in theta(1,i1,i2,i3)
               call rchkusr()
	    END DO
	 END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute trace of y
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine tracey(y,n1,n2,n3,try)
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(6,n1,n2,n3),try(n1,n2,n3)
      integer i1,i2,i3
      DO i1=1,n1
         DO i2=1,n2
	    DO i3=1,n3
	       try(i1,i2,i3)=y(1,i1,i2,i3)+y(4,i1,i2,i3)+y(6,i1,i2,i3)
	    END DO
	 END DO
      END DO
      return
      end
