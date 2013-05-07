//  L-BFGS-B  version
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

int ngradd = 0;
double* varinv;
double* signal, *btb;
//double *signal, *varinv, *btb;

typedef struct
{
  int ngrad;
  double* sig;
  double* btb;
  double* varinv;
  int fnscale;
} optimtens;



extern void F77_NAME(ftensor)(double* param, double* sig, int* ng,
             double* btb, double* varinv, double* gv, double* fv);

extern void F77_NAME(gtensor)(double* param, double* sig, int* ng,
             double* btb, double* varinv, double* gv, 
             double* fv, double* gradient);

double ftens(int param_length, double *param, void* ex){

  double result = 0;
  
  double* gv = Calloc(ngradd, double);
//  double* sigi = Calloc(ngradd, double);
//  double* varinvi = Calloc(ngradd, double);
//  int i;
//  Rprintf("vari and sigi %f%f",varinv[100],signal[100]);
//  Rprintf("ngradd %i",ngradd);
/*  for(i=0;i<ngradd;i++){
    varinvi[i] = varinv[i];
    sigi[i] = signal[i];
   Rprintf("vari and sigi %f%f",varinvi[i],sigi[i]);
  }*/
//  Rprintf("\n");
// calculate risk for tensor model
  F77_CALL(ftensor)(param, signal, &ngradd, btb, varinv, gv, &result);
  
//  Free(varinvi);
//  Free(sigi);
  Free(gv);
  return result;
}


void gtens(int param_length, double* param, double* result, void* ex){
 
  double* gv = Calloc(ngradd, double);
  double* fv = Calloc(ngradd, double);
//  double* sigi = Calloc(ngradd, double);
//  double* varinvi = Calloc(ngradd, double);
// calculate gradient of risk for tensor model

//  int i;
/*  for(i=0;i<ngradd;i++){
    sigi[i] = signal[i];
    varinvi[i] = varinv[i];
  }*/
  F77_CALL(gtensor)(param, signal, &ngradd, btb, varinv, gv, fv, result);
  
//  Free(varinvi);
//  Free(sigi);
  Free(fv);
  Free(gv);
}

void dtens( int* n1, double* param, double* sig_in, int* ngrad, double* btb_in, 
	    double* sdcoef, double* sig_tmp, double* vinv_tmp, int* maxit, double* reltol){
   
   // tensor optmethod: BFGS
   //
   // i2 and n1 dimxx,  current/total respectively
   // ngrad : # gradients
   // maxcomp : maximal order of mix ???
   
   int *mask;
   int trace = 0, nREPORT = 1;
   int fncount = 5, grcount=2;    // number of calls to obj fct in optim
   int i, i2, j, l = 0, param_length, dimxx;
   double low, up, signalg; 
   double ttt = 0, abstol = R_NegInf;
   double *param_work = 0; 
//   double* varinv;
   
   int fail;              // failure code for optim: zero is OK
   double Fmin = 0.;          // minimal value of obj fct in optim 
   
   //Setting global variables
   signal = sig_tmp; varinv = vinv_tmp;  ngradd = *ngrad;
   Rprintf("ngrad %i \n",ngradd);
   param_length=7;
   dimxx = *n1;
   btb = btb_in;
   low = sdcoef[0]+sdcoef[2]*sdcoef[1];
   up = sdcoef[0]+sdcoef[3]*sdcoef[1];
   // not yet initialized
//   Rprintf("ngrad %i \n", ngradd);
   
//   Rprintf("param_length %i \n", param_length);
   param_work = (double *) R_alloc(param_length, sizeof(double));
   mask = (int *) R_alloc(param_length, sizeof(int));
   for (i = 0; i < param_length; i++) mask[i] = 1;

   optimtens myoptimtens;
   myoptimtens.ngrad = ngradd;
   myoptimtens.btb = btb_in;
   myoptimtens.fnscale = 1;
//   varinv = varinv;
   Rprintf("ngrad1 %i \n",ngradd);
//   signal = sig_tmp;
   Rprintf("ngrad2 %i \n",ngradd);
//   varinv = vinv_tmp;
   Rprintf("ngrad3 %i \n",ngradd);
   for(i2 = 0; i2 < dimxx; i2++){ 
      for(j = 0; j < param_length; j++){
         param_work[l] = param[j+i2*7];
      }
   Rprintf("ngrad3 %i \n",ngradd);
      for(l=0;l<ngradd;l++){
	  signalg = sig_in[l+i2*ngradd];
   Rprintf("signalg %f \n",signalg);
          sig_tmp[l] = signalg;
   Rprintf("sigl %f \n",sig_tmp[l]);
	  ttt = sdcoef[0]+signalg*sdcoef[1];
	  if(signal[l]<sdcoef[2]) ttt=low;
	  if(signal[l]>sdcoef[3]) ttt=up;
	  vinv_tmp[l] = 1./ttt/ttt;
   Rprintf("vinvl %f \n",vinv_tmp[l]);
      }
   Rprintf("vari and sigi %f%f",varinv[100],signal[100]);
      vmmin(param_length, param_work, &Fmin, ftens, gtens,
            *maxit, trace, mask, abstol, *reltol, nREPORT,
            &myoptimtens, &fncount, &grcount, &fail);
      for(j = 0; j < param_length; j++){
          param[j+i2*7] = param_work[j];
      }
   } // end i2
}


