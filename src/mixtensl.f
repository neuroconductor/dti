C
C  version for use with L-BFGS-B
C
      subroutine rskmixl2(par,npar,siq,g,b,ng,risk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,risk
      integer i
      real*8 resi,fval,alpha2
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      alpha2=alpha*alpha
C need to have positive alpha2 to guarantee lambda1 >= lambda2
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-2,lambda,alpha2,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl1(par,npar,siq,g,b,ng,alpha2,risk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha2,risk
      integer i
      real*8 resi,fval
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-1,lambda,alpha2,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl0(par,npar,siq,g,b,ng,lambda,alpha2,risk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha2,risk
      integer i
      real*8 resi,fval
      ncomp=npar/3
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar,lambda,alpha2,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
C
C    Evaluations for tensor-mixture models
C
      subroutine fmixturl(par,npar,lambda,alpha,g,b,fval)
C
C  compute expected value
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, c
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),fval,lambda,alpha
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),b,wi,sw,swi,bclam,f1,f2
      real*8 ddot3sq
      external ddot3sq
C  extract parameters
      ncomp=(npar-2)/3
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
      bclam=b*alpha*lambda
      swi=1.d0/sw
      f1=exp(-b*lambda)*swi
      f2=1.d0
      DO i=1,ncomp
         f2=f2+w(i)*exp(-bclam*ddot3sq(phi(i),eta(i),g))
      END DO
      fval=f1*f2
      RETURN
      END
      subroutine drskml2(par,npar,siq,g,b,ng,drisk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,drisk(npar)
      integer i,j
      real*8 resi,fval,dval(15),drisk0(17)
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
         call dfml2(par,npar,lambda,alpha2,g(1,i),b(i),fval,dval,
     1                  dlam,dalpha)         
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar-1)=drisk0(npar-1)-resi*dlam
         drisk0(npar)=drisk0(npar-1)-resi*dalpha
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
C
C  gradients
C
      subroutine drskml1(par,npar,siq,g,b,ng,alpha2,drisk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha2,drisk(npar)
      integer i
      real*8 resi,fval,dval(15),drisk0(17)
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
        call dfml1(par,npar,lambda,alpha2,g(1,i),b(i),fval,dval,dlam)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar)=drisk0(npar-1)-resi*dlam
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
      subroutine drskml0(par,npar,siq,g,b,ng,lambda,alpha2,drisk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha2,drisk(npar)
      integer i
      real*8 resi,fval,dval(15),drisk0(17)
      ncomp=npar/3
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
         call dfml0(par,npar,lambda,alpha2,g(1,i),b(i),fval,dval)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
      subroutine dfml2(par,npar,lambda,alpha,g,b,fval,dval,dlam,
     1                  dalpha)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i, lambda and alpha
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, alpha
C  fval contains function value
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b,dlam,dalpha
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,swi,bclam
      real*8 embl,embclgd2(5),dgtd(3,5),fval,f1,f2,f1l,f2l,f1c
C  extract parameters
      ncomp=(npar-2)/3
C need to have positive alpha2 to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      bclam=b*alpha*lambda
      swi=1.d0/sw
      embl=exp(-b*lambda)
      f1=exp(-b*lambda)*swi
      f2=1.d0
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         embclgd2(i)=exp(-bclam*dgtd(1,i))
         f2=f2+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=f1*f2
C  compute derivatives
      f1l=-embl*swi*b*alpha
      f2l=0.d0
      f1c=-embl*swi*b*lambda
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=f1*embclgd2(i)-fval/swi
C  derivative with respect to eta_i
      END DO
      dlam=-b*fval+f1l*f2l
C  derivative with respect to lambda
      dalpha=f1c*f2l
C  derivative with respect to c
      RETURN
      END
      subroutine dfml1(par,npar,lambda,alpha,g,b,fval,dval,dlam)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b,dlam
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,swi,bclam
      real*8 embl,embclgd2(5),dgtd(3,5),fval,f1,f2,f1l,f2l,f1c
C  extract parameters
      ncomp=(npar-2)/3
C need to have positive alpha2 to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      bclam=b*alpha*lambda
      swi=1.d0/sw
      embl=exp(-b*lambda)
      f1=exp(-b*lambda)*swi
      f2=1.d0
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         embclgd2(i)=exp(-bclam*dgtd(1,i))
         f2=f2+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=f1*f2
C  compute derivatives
      f1l=-embl*swi*b*alpha
      f2l=0.d0
      f1c=-embl*swi*b*lambda
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=f1*embclgd2(i)-fval/swi
C  derivative with respect to eta_i
      END DO
      dlam=-b*fval+f1l*f2l
C  derivative with respect to lambda
      RETURN
      END
      subroutine dfml0(par,npar,lambda,alpha,g,b,fval,dval)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,swi,bclam
      real*8 embl,embclgd2(5),dgtd(3,5),fval,f1,f2,f1l,f2l,f1c
C  extract parameters
      ncomp=(npar-2)/3
C need to have positive alpha2 to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      bclam=b*alpha*lambda
      swi=1.d0/sw
      embl=exp(-b*lambda)
      f1=exp(-b*lambda)*swi
      f2=1.d0
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         embclgd2(i)=exp(-bclam*dgtd(1,i))
         f2=f2+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=f1*f2
C  compute derivatives
      f1l=-embl*swi*b*alpha
      f2l=0.d0
      f1c=-embl*swi*b*lambda
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=2.d0*f1l*lambda*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=f1*embclgd2(i)-fval/swi
C  derivative with respect to eta_i
      END DO
      RETURN
      END
      
      real*8 function ddot3sq(phi,eta,g)
C  compute (g^T n(phi,eta))^2 
      implicit logical (a-z)
      real*8 phi,eta,g(3),sphi,z
      sphi=dsin(phi)
      z=g(1)*dcos(eta)*sphi+g(2)*dsin(eta)*sphi+g(3)*dcos(phi)
      ddot3sq=z*z
      RETURN
      END
      subroutine dgtddphi(phi,eta,g, d)
C  compute (g^T n(phi,eta))^2 
C           d(g^T n(phi,eta))^2 / d phi
C           d(g^T n(phi,eta))^2 / d eta
C  as components of d
      implicit logical (a-z)
      real*8 phi,eta,g(3),sphi,cphi,ceta,seta,z,ze,zp,d(3)
      sphi=dsin(phi)
      cphi=dcos(phi)
      seta=dsin(eta)
      ceta=dcos(eta)
      z=g(1)*ceta*sphi+g(2)*seta*sphi+g(3)*cphi
      zp=g(1)*ceta*cphi+g(2)*seta*cphi-g(3)*sphi
      ze=g(2)*ceta*sphi-g(1)*seta*sphi
      d(1)=z*z
      d(2)=2.d0*z*zp
      d(3)=2.d0*z*ze
      RETURN
      END
      