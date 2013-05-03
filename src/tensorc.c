//  L-BFGS-B  version
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

typedef struct
{
  int ngrad;
  double* sig;
  double* btb;
  double* vinv;
  int fnscale;
} optimtens;

optimtens myoptimtens;

extern void F77_NAME(ftensor)(double* param, double* sig, int* ng,
             double* btb, double* vinv, double* gv, double* fv);

extern void F77_NAME(gtensor)(double* param, double* sig, int* ng,
             double* btb, double* vinv, double* gv, 
             double* fv, double* gradient);

double ftens(int param_length, double *param, void* ex){

  double result = 0;
  double* gv = Calloc(&myoptimtens.ngrad, double);
// calculate risk for tensor model
  F77_CALL(ftensor)(param, myoptimtens.sig, &myoptimtens.ngrad, myoptimtens.btb, myoptimtens.vinv, gv, &result);

  Free(gv);
  return result;
}


void gtens(int param_length, double* param, double* result, void* ex){

  double* gv = Calloc(&myoptimtens.ngrad, double);
  double* fv = Calloc(&myoptimtens.ngrad, double);
// calculate gradient of risk for tensor model

  F77_CALL(gtensor)(param, myoptimtens.sig, &myoptimtens.ngrad, myoptimtens.btb, myoptimtens.vinv, gv, fv, result);
  
  Free(fv);
  Free(gv);
}

void dtens( int* n1, double* param, double* sig_in, int* ngrad, 
	    double* btb_in, double* sdcoef, int* maxit, double* reltol){
   
   // tensor optmethod: BFGS
   //
   // i1 and n1 dimx,  current/total respectively
   // ngrad : # gradients
   // maxcomp : maximal order of mix ???
   
   int *mask;
   int trace = 0, nREPORT = 1;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int i, i1, j, l = 0, param_length, ngradc, dimx;
   double low, up; 
   double ttt = 0, abstol = R_NegInf;
   double *param_work = 0, *sig = 0, *vinv = 0;
   
   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim 
   
   //Setting global variables
   param_length=7;
   dimx = *n1;
   ngradc = *ngrad;
   low = sdcoef[0]+sdcoef[2]*sdcoef[1];
   up = sdcoef[0]+sdcoef[3]*sdcoef[1];
   // not yet initialized
   vinv = (double*) R_alloc(ngradc, sizeof(double));
   
   param_work = (double *) R_alloc(param_length, sizeof(double));
   mask = (int *) R_alloc(param_length, sizeof(int));
   for (i = 0; i < param_length; i++) mask[i] = 1;

   myoptimtens.ngrad = ngradc;
   myoptimtens.btb = btb_in;
   myoptimtens.fnscale = 1;
   sig = (double *) R_alloc(ngradc, sizeof(double));
   
   for(i1 = 0; i1 < dimx; i1++){ 
      for(j = 0; j < param_length; j++){
         param_work[l] = param[j+i1*7];
      }
      for(l=0;l<ngradc;l++){
          sig[l] = sig_in[l+i1*ngradc];
	  ttt = sdcoef[0]+sig[l]*sdcoef[1];
	  if(sig[l]<sdcoef[2]) ttt=low;
	  if(sig[l]>sdcoef[3]) ttt=up;
	  vinv[l] = 1./ttt/ttt;
      }
      myoptimtens.sig = sig;
      myoptimtens.vinv = vinv;
      vmmin(param_length, param_work, &Fmin, ftens, gtens,
            *maxit, trace, mask, abstol, *reltol, nREPORT,
            &myoptimtens, &fncount, &grcount, &fail);
      for(j = 0; j < param_length; j++){
          param[j+i1*7] = param_work[j];
      }
   } // end i1
}


