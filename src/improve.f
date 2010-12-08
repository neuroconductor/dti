      subroutine imprpara(maxcomp,maxorder,order,mix,andir,orient,lev,
     1                    vext,par,npar)
      implicit logical (a-z)
      integer maxcomp,maxorder,order(3,3,3),npar
      real*8 mix(maxorder,3,3,3),andir(3,maxorder,3,3,3),minw,mangle,
     1       orient(2,maxorder,3,3,3),lev(2,3,3,3),par(npar),vext(3)
      integer i1,j1,k1,i2,j2,k2,o,o2,ord,ord2,icomp,besti,bestj,
     1        bestk,besto
      real*8 krit(6,3,3,3),dir(3),bestor(2),th,z,bestkrit
C   restricts maximum order to 6
      DO i1=1,3
         i2=4-i1
         DO j1=1,3
            j2=4-j1
            DO k1=1,3
               k2=4-k1
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
               dir(1)=(i1-2)*vext(1)
               dir(2)=(j1-2)*vext(2)
               dir(3)=(k1-2)*vext(3)
               ord = order(i1,j1,k1)
               DO o=1,ord
                  z = abs(dir(1)*andir(1,o,i1,j1,k1)+
     1                    dir(2)*andir(2,o,i1,j1,k1)+
     1                    dir(3)*andir(3,o,i1,j1,k1))
C  this is large if the component shows in direction of the voxel
                  ord2 = order(i2,j2,k2)
                  DO o2=1,ord2
                     z=z+mix(o,i2,j2,k2)*
     1               abs(andir(1,o,i1,j1,k1)*andir(1,o2,i2,j2,k2)+
     1                   andir(2,o,i1,j1,k1)*andir(2,o2,i2,j2,k2)+
     1                   andir(3,o,i1,j1,k1)*andir(3,o2,i2,j2,k2))
C this adds something if the same information is present on the opposite side
                  END DO
                  krit(o,i1,j1,k1)=mix(o,i1,j1,k1)*z
               END DO
            END DO
         END DO
      END DO
C  now prepare initial parameters
      icomp=maxcomp
      th=lev(1,2,2,2)-lev(2,2,2,2)
      DO while (icomp.gt.0) 
C
C   search for highest ranked component (direction)
C
         bestkrit=krit(1,1,1,1)
         bestor(1) = orient(1,1,1,1,1)
         bestor(2) = orient(2,1,1,1,1)
         besti=1
         bestj=1
         bestk=1
         besto=1
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(krit(o,i1,j1,k1).gt.bestkrit) THEN
                        bestkrit=krit(o,i1,j1,k1)
                        bestor(1) = orient(1,o,i1,j1,k1)
                        bestor(2) = orient(2,o,i1,j1,k1)
                        besti=i1
                        bestj=j1
                        bestk=k1
                        besto=o
                     END IF
                  END DO
               END DO
            END DO
         END DO
         th=th+lev(1,besti,bestj,bestk)-lev(2,besti,bestj,bestk)
         par(2*(maxcomp-icomp+1))=bestor(1)
         par(2*(maxcomp-icomp+1)+1)=bestor(2)
C  mark this entry and all that are to close as uninteresting
         dir(1)=andir(1,besto,besti,bestj,bestk)
         dir(2)=andir(2,besto,besti,bestj,bestk)
         dir(3)=andir(3,besto,besti,bestj,bestk)
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(dir(1)*andir(1,o,i1,j1,k1)+
     1                  dir(2)*andir(2,o,i1,j1,k1)+
     2                  dir(3)*andir(3,o,i1,j1,k1).gt..9) THEN
                        krit(o,i1,j1,k1)=0.d0
                     END IF
                  END DO
               END DO
            END DO
         END DO
         icomp=icomp-1
      END DO
      par(1)=th/(maxcomp+1)
      RETURN
      END
      subroutine imprparb(maxcomp,maxorder,order,mix,andir,orient,lev,
     1                    vext,par,npar)
      implicit logical (a-z)
      integer maxcomp,maxorder,order(3,3,3),npar
      real*8 mix(maxorder,3,3,3),andir(3,maxorder,3,3,3),minw,mangle,
     1       orient(2,maxorder,3,3,3),lev(2,3,3,3),par(npar),vext(3)
      integer i1,j1,k1,i2,j2,k2,o,o2,ord,ord2,icomp,besti,bestj,
     1        bestk,besto,ord0,i
      real*8 krit(6,3,3,3),dir(3),bestor(2),th,z,bestkrit
C   restricts maximum order to 6
      ord0=order(2,2,2)
C  keep directions that are found to be informative
      th=1.d0
      IF(ord0.gt.0) THEN
         th=(lev(1,2,2,2)-lev(2,2,2,2))
         DO i=1,ord0
            par(2*(maxcomp-i+1))=orient(1,i,2,2,2)
            par(2*(maxcomp-i+1)+1)=orient(2,i,2,2,2)            
         END DO
      END IF
      DO i1=1,3
         i2=4-i1
         DO j1=1,3
            j2=4-j1
            DO k1=1,3
               k2=4-k1
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
               dir(1)=(i1-2)*vext(1)
               dir(2)=(j1-2)*vext(2)
               dir(3)=(k1-2)*vext(3)
               ord0 = order(i1,j1,k1)
               DO o=1,ord0
                  z = abs(dir(1)*andir(1,o,i1,j1,k1)+
     1                    dir(2)*andir(2,o,i1,j1,k1)+
     1                    dir(3)*andir(3,o,i1,j1,k1))
C  this is large if the component shows in direction of the voxel
                  ord2 = order(i2,j2,k2)
                  DO o2=1,ord2
                     z=z+mix(o,i2,j2,k2)*
     1               abs(andir(1,o,i1,j1,k1)*andir(1,o2,i2,j2,k2)+
     1                   andir(2,o,i1,j1,k1)*andir(2,o2,i2,j2,k2)+
     1                   andir(3,o,i1,j1,k1)*andir(3,o2,i2,j2,k2))
C this adds something if the same information is present on the opposite side
                  END DO
                  krit(o,i1,j1,k1)=mix(o,i1,j1,k1)*z
               END DO
            END DO
         END DO
      END DO
C  make directions that are close to existing ones as uninteresting 
      DO i1=1,3
         DO j1=1,3
            DO k1=1,3
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     z=0.d0
                     DO i=1,ord0
                     z=max(z,abs(
     1                  andir(1,i,2,2,2)*andir(1,o,i1,j1,k1)+
     1                  andir(2,i,2,2,2)*andir(2,o,i1,j1,k1)+
     2                  andir(3,i,2,2,2)*andir(3,o,i1,j1,k1)))
                     END DO
                     krit(o,i1,j1,k1)=krit(o,i1,j1,k1)*(1-z*z)
                  END DO
               END DO
            END DO
         END DO
C  now prepare initial parameters
      icomp=maxcomp
      DO while (icomp.gt.ord0) 
C
C   search for highest ranked component (direction)
C
         bestkrit=krit(1,1,1,1)
         bestor(1) = orient(1,1,1,1,1)
         bestor(2) = orient(2,1,1,1,1)
         besti=1
         bestj=1
         bestk=1
         besto=1
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(krit(o,i1,j1,k1).gt.bestkrit) THEN
                        bestkrit=krit(o,i1,j1,k1)
                        bestor(1) = orient(1,o,i1,j1,k1)
                        bestor(2) = orient(2,o,i1,j1,k1)
                        besti=i1
                        bestj=j1
                        bestk=k1
                        besto=o
                     END IF
                  END DO
               END DO
            END DO
         END DO
         th=th+lev(1,besti,bestj,bestk)-lev(2,besti,bestj,bestk)
         par(2*(maxcomp-icomp+1))=bestor(1)
         par(2*(maxcomp-icomp+1)+1)=bestor(2)
C  mark this entry and all that are to close as uninteresting
         dir(1)=andir(1,besto,besti,bestj,bestk)
         dir(2)=andir(2,besto,besti,bestj,bestk)
         dir(3)=andir(3,besto,besti,bestj,bestk)
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(dir(1)*andir(1,o,i1,j1,k1)+
     1                  dir(2)*andir(2,o,i1,j1,k1)+
     2                  dir(3)*andir(3,o,i1,j1,k1).gt..9) THEN
                        krit(o,i1,j1,k1)=0.d0
                     END IF
                  END DO
               END DO
            END DO
         END DO
         icomp=icomp-1
      END DO
      par(1)=th/(maxcomp+1)
      RETURN
      END