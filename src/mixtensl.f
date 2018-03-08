C
C  version for use with L-BFGS-B
C
      subroutine rskmixl2(par,npar,siq,g,b,ng,risk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       risk
      integer i,ncomp
      double precision resi,fval
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-2,lambda,alpha,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl1(par,npar,siq,g,b,ng,alpha,risk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       risk
      integer i,ncomp
      double precision resi,fval
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-1,lambda,alpha,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl0(par,npar,siq,g,b,ng,lambda,alpha,risk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       risk
      integer i
      double precision resi,fval
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar,lambda,alpha,g(1,i),b(i),fval)
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
      implicit none
      integer npar
      double precision par(npar),g(3),fval,lambda,alpha
      integer ncomp,i
      double precision w(5),phi(5),eta(5),b,wi,sw,blam,f1
      double precision ddot3sq
      external ddot3sq
C  extract parameters
      ncomp=npar/3
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
         f1=f1+w(i)*exp(-blam*(1.d0+alpha*ddot3sq(phi(i),eta(i),g)))
      END DO
      fval=f1/sw
      RETURN
      END
      subroutine drskml2(par,npar,siq,g,b,ng,drisk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j,ncomp
      double precision resi,fval,dval(15),drisk0(17),dlam,dalpha
C  calculate numerical gradients for comparison
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfml2(par,npar,lambda,alpha,g(1,i),b(i),fval,dval,
     1                  dlam,dalpha)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar-1)=drisk0(npar-1)-resi*dlam
         drisk0(npar)=drisk0(npar)-resi*dalpha
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)
      END DO
      RETURN
      END
C
C  gradients
C
      subroutine drskml1(par,npar,siq,g,b,ng,alpha,drisk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j,ncomp
      double precision resi,fval,dval(15),drisk0(17),dlam
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfml1(par,npar,lambda,alpha,g(1,i),b(i),fval,dval,dlam)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar)=drisk0(npar)-resi*dlam
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)
      END DO
      RETURN
      END
      subroutine drskml0(par,npar,siq,g,b,ng,lambda,alpha,drisk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j,ncomp
      double precision resi,fval,dval(15),drisk0(17)
      ncomp=npar/3
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfml0(par,npar,lambda,alpha,g(1,i),b(i),fval,dval)
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
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),alpha,lambda,b,dlam,
     1       dalpha
      integer ncomp,i
      double precision w(5),phi(5),eta(5),wi,sw,blam
      double precision embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,f2c,
     1       adgtd(5)
C  extract parameters
      ncomp=(npar-2)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      f1l=-b*(1.d0+alpha)/sw*f1
      f2l=0.d0
      f2c=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         f2c=f2c+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      dlam=f1l-b/sw*f2l
C  derivative with respect to lambda
      dalpha=-blam/sw*(f1+f2c)
C  derivative with respect to c
      RETURN
      END
      subroutine dfml1(par,npar,lambda,alpha,g,b,fval,dval,dlam)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),alpha,lambda,b,dlam
      integer ncomp,i
      double precision w(5),phi(5),eta(5),wi,sw,blam
      double precision embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,adgtd(5)
C  extract parameters
      ncomp=(npar-1)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      f1l=-b*(1.d0+alpha)/sw*f1
      f2l=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      dlam=f1l-b/sw*f2l
C  derivative with respect to lambda
      RETURN
      END
      subroutine dfml0(par,npar,lambda,alpha,g,b,fval,dval)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),alpha,lambda,b
      integer ncomp,i
      double precision w(5),phi(5),eta(5),wi,sw,blam
      double precision embclgd2(5),dgtd(3,5),fval,adgtd(5)
C  extract parameters
      ncomp=npar/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      fval = exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      DO i=1,ncomp
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      RETURN
      END
