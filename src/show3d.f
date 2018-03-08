CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine adcradii(vert,nv,tens,ntens,radii)
      implicit none
      integer nv,ntens
      double precision vert(3,nv),tens(6,ntens),radii(nv,ntens)
      integer i,j
      double precision x,y,z,s,xx,xy,xz,yy,yz,zz
      DO i=1,nv
         x=vert(1,i)
         y=vert(2,i)
         z=vert(3,i)
         xx=x*x
         xy=2.d0*x*y
         xz=2.d0*x*z
         yy=y*y
         yz=2.d0*y*z
         zz=z*z
         DO j=1,ntens
            s=tens(1,j)*xx+tens(2,j)*xy+tens(3,j)*xz+
     1        tens(4,j)*yy+tens(5,j)*yz+tens(6,j)*zz
            radii(i,j)=s
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ellradii(vert,nv,tens,ntens,radii)
      implicit none
      integer nv,ntens
      double precision vert(3,nv),tens(6,ntens),radii(nv,ntens)
      integer ierr,i,j
      double precision qform3,ev(3),edir(3,3)
      DO j=1,ntens
         call eigen3(tens(1,j),ev,edir,ierr)
         if(ev(3).gt.1d-6.and.ierr.eq.0) THEN
            DO i=1,nv
               radii(i,j)=1.d0/sqrt(qform3(vert(1,i),edir,ev))
            END DO
         ELSE
            DO i=1,nv
               radii(i,j)=0.d0
            END DO
         END IF
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for odf in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine odfradii(vert,nv,tens,ntens,radii)
      implicit none
      integer nv,ntens
      double precision vert(3,nv),tens(6,ntens),radii(nv,ntens)
      integer ierr,i,j
      double precision qform3,ev(3),edir(3,3),z,z1
      DO j=1,ntens
         call eigen3(tens(1,j),ev,edir,ierr)
         if(ev(3).gt.1d-6.and.ierr.eq.0) THEN
            z1=0.07957747d0/sqrt(ev(1)*ev(2)*ev(3))
            DO i=1,nv
               z=qform3(vert(1,i),edir,ev)
               radii(i,j)=z1/sqrt(z*z*z)
            END DO
         ELSE
            DO i=1,nv
               radii(i,j)=0.d0
            END DO
         END IF
      END DO
      RETURN
      END

      double precision function qform3(y,edir,ev)
      double precision y(3),edir(3,3),ev(3)
      double precision z1,z2,z3
      z1=y(1)*edir(1,1)+y(2)*edir(2,1)+y(3)*edir(3,1)
      z2=y(1)*edir(1,2)+y(2)*edir(2,2)+y(3)*edir(3,2)
      z3=y(1)*edir(1,3)+y(2)*edir(2,3)+y(3)*edir(3,3)
      qform3=z1*z1/ev(1)+z2*z2/ev(2)+z3*z3/ev(3)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mixtradi(vert,nv,ev,ori,mix,ord,mo,nobj,radii)
      implicit none
C
C     vert   vertices
C     nv     number of vertices
C     ev     ev[1]+ev[2], ev[2], ev[2] are the eigenvalues of components
C     ori    ori[,m]   contain theta and phi of main direction for component m
C     mix    mix[m] mixture coefficient for component m
C     ord    contains number od components
C     mo     maximum number of components
C     nobj   number of objects
C     scale  scaling factor
C     radii  resulting radii of ODF
C
      integer nv,nobj,mo,ord(nobj)
      double precision vert(3,nv),ev(2,nobj),ori(2,mo,nobj),
     1       mix(mo,nobj),radii(nv,nobj)
      integer i,j,k
      double precision dotprod3,c12,sth,dir(3,5),fourpi,c12fp,z,utd,
     1       zk,e1,e2,sm
C     assumes maximum of 5 micxture components
      if(mo.gt.5) THEN
         call intpr("mo restricted to 5, is",18,mo,1)
         return
      END IF
      fourpi = 12.566371d0
      DO j = 1,nobj
         e1 = ev(1,j)
         e2 = ev(2,j)
C   limit maximal excentricity to mex+1 for visualisation purposes
         c12 = (e1-e2)/e1
         c12fp = sqrt(e2/e1)/fourpi
         sm = 1.d0
C   calculate weight for order 0 component in sm
         DO k = 1,ord(j)
            sm = sm - mix(k,j)
            sth = sin(ori(1,k,j))
            dir(1,k) = sth*cos(ori(2,k,j))
            dir(2,k) = sth*sin(ori(2,k,j))
            dir(3,k) = cos(ori(1,k,j))
         END DO
         DO i = 1,nv
            z = sm
            DO k = 1,ord(j)
               utd = dotprod3(dir(1,k),vert(1,i))
               zk = 1.d0-c12*utd*utd
               z  = z + mix(k,j)/sqrt(zk*zk*zk)
            END DO
            radii(i,j) = z*c12fp
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mixandir(ori,mix,ord,mo,nobj,andir)
      implicit none
C
C     ori    ori[,m]   contain theta and phi of main direction for component m
C     mix    mix[m] mixture coefficient for component m
C     ord    contains number od components
C     mo     maximum number of components
C     nobj   number of objects
C     andir  main diffusion directions
C
      integer nobj,mo,ord(nobj)
      double precision ori(2,mo,nobj),mix(mo,nobj),andir(3,mo,nobj)
      integer j,k
      double precision sth
C     assumes maximum of 5 micxture components
      DO j = 1,nobj
         DO k = 1,mo
            if(k.le.ord(j))THEN
               sth = sin(ori(1,k,j))
               andir(1,k,j) = mix(k,j)*sth*cos(ori(2,k,j))
               andir(2,k,j) = mix(k,j)*sth*sin(ori(2,k,j))
               andir(3,k,j) = mix(k,j)*cos(ori(1,k,j))
            ELSE
               andir(1,k,j) = 0.d0
               andir(2,k,j) = 0.d0
               andir(3,k,j) = 0.d0
            END IF
         END DO
      END DO
      RETURN
      END

      double precision function dotprod3(a,b)
      implicit none
      double precision a(3),b(3)
      dotprod3=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      RETURN
      END
