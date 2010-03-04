#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

typedef struct {
  int ngrad;
  double *siq;
  double *grad;
} optimex;

// double F77_NAME(dotprod3)(double *, double *);

double fn1(int n, double *par, void *ex) {
  optimex ext = *((optimex*)ex);

  int m = (n-1)/3; // order of mix tensor model

  double z[ext.ngrad]; // value of signal fct. for each grad according to current model
  for (int k = 0; k < ext.ngrad; k++) z[k] = 0;
  double dir[3]; // used for direction vector calculated from angles in par
  double c1 = exp(par[0]);
  double c2 = exp(par[1]);
  double sw = 0;
  double w;
  double ew;
  double sth;
  double z1 = 0;
  int i;
  int j;

  double erg = 0; // result

  for (i = 0; i < m; i++) {
    if (i == m-1 ) {
      w = 1.0 - sw;
    } else {
      ew = exp(par[3*i+4]);
      w = (1.0 - sw)*ew/(1.0+ew);
      sw += w;
    }
    sth = sin(par[i*3+2]);
    dir[0] = sth*cos(par[i*3+3]);
    dir[1] = sth*sin(par[i*3+3]);
    dir[2] = cos(par[i*3+2]);
    for (j = 0; j < ext.ngrad; j++) {
      z1 = dir[0]*ext.grad[0+j] + dir[1]*ext.grad[ext.ngrad+j] + dir[2]*ext.grad[2*ext.ngrad+j];
//      z1 = F77_CALL(dotprod3)(dir, test);
      z[j] += w*exp(-c2 - c1*z1*z1);
    }
  }

  // calculate residual variance
  for (j = 0; j < ext.ngrad; j++) {
    z1 = ext.siq[j] - z[j];
    erg += z1*z1;
  }

  // finished
  return erg;
}

void mixtensor(int *n, double *par, double *x, int *ngrad, double *siq, double *grad, int *maxit, double *reltol, double *Fmin){

  int fail; // failure code: zero is OK
  int fncount; // number of calls to fct.

  optimex myoptimpar;
  myoptimpar.ngrad = *ngrad;
  myoptimpar.siq = siq;
  myoptimpar.grad = grad;

  nmmin(*n, par, x, Fmin, fn1,
        &fail, R_NegInf, *reltol, &myoptimpar,
        1.0, 0.5, 2.0, 0,
        &fncount, *maxit);

}
/*
void mixture(int *n1, int *n2, int *n3, double *mask, double *si, int *ngrad, double *grad, int *sneighbors, int *maxneighb, double *sweights, int *maxcomp, int *maxit, double *reltol,
             double *order, double *lev, double *mix, double *orient, double *p, double *sigma2){

  int i1;
  int i2;
  int i3;
  int mc0;
  // z;
  double *par;
  int lpar;
  double rss = R_PosInf;
  double krit = R_NegInf;
  double *siq;
  double *smsi;

  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      for (int i3 = 0; i3 < n3; i3++) {
        if mask[i1][i2][i3] break; // formuliere als eindim. array! und wir wollen kein break!
        smsi = F77_CALL(smsi)(si[i1][i2][i3][], *ngrad, *sneighbors, *maxneighb, *sweights); // definiere in C und als function!

        z <- locmin(smsi,sneighbors,grad) // wie loesen?
        mc0 = min(*maxcomp, dim(z$grad)[1]) // ueber if

        for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(z$grad[j,])
        par <- numeric(3*mc0+1)
        lpar = 3*mc0+1
        lev[,i1,i2,i3] <- par[1:2] <- log(-log(siq[i1,i2,i3,z$ind[1]])*c(.8,.2))
        par[rep(3*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3] 
        rss <- Inf
        krit <- Inf
        for(k in mc0:1){
          if(k>1) par[3*(2:k)-1] <- -log((k-1):1)
          lpar <- length(par);
        z <-switch(method,
                  "mixtensor" = .C("mixtensor",
                                   as.integer(lpar),
                                   as.double(par),
                                   par = double(lpar),
                                   as.integer(ngrad0),
                                   as.double(siq[i1,i2,i3,]),
                                   as.double(as.matrix(grad)),
                                   as.integer(maxit),
                                   as.double(reltol),
                                   value = as.double(1))[c("par", "value")],
                 "Jian"       = optim(par[1:(3*k+1)],mfun2,siq=siq[i1,i2,i3,],grad=grad,ep=p,control=list(maxit=maxit,reltol=reltol)),
                 "Jian2"      = optim(c(par[1:(3*k+1)],25),mfun3,siq=siq[i1,i2,i3,],grad=grad,ep=p,control=list(maxit=maxit,reltol=reltol)))
        rss <- min(z$value,rss)
        ttt <- z$value+2*(3*k+1)/(ngrad-3*maxcomp-1)*rss
#        cat("risk",z$value,ttt,"\n")
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- as.integer(k)
           lev[,i1,i2,i3] <- z$par[1:2]
           mix[1:k,i1,i2,i3] <- if(k>1) getw(z$par[3*(2:k)-1]) else 1
           or <- matrix(z$par[rep(3*(1:k),rep(2,k))+c(0,1)],2,k)
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           orient[,1:k,i1,i2,i3] <- or
           if(method=="Jian2") p[i1,i2,i3] <- z$par[length(z$par)]
       }
     }
   sigma2[i1,i2,i3] <- rss/(ngrad-3*mc0-1)
     }
   }
 }

}*/