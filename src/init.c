#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_SUB(d2rall)();
void F77_SUB(projdt2)();
void F77_SUB(r2dall)();
void F77_SUB(rho2d0)();
void F77_SUB(smsigma)();
void F77_SUB(tensres)();
void dtens( int * n1, double * param, double * sig_in, int * ngrad, double * btb_in, 
	    double * sdcoef, double * sig_tmp, double * vinv_tmp, int * maxit, double * reltol);
static R_CMethodDef cEntries[] = {
   {"dtens", (DL_FUNC) &dtens, 10},
   {NULL, NULL, 0, NULL}
};
static R_FortranMethodDef fortranEntries[] = {
   {"D2Rall", (DL_FUNC) &F77_SUB(d2rall), 3},
   {"projdt2", (DL_FUNC) &F77_SUB(projdt2), 9},
   {"R2Dall", (DL_FUNC) &F77_SUB(r2dall), 3},
   {"rho2D0", (DL_FUNC) &F77_SUB(rho2d0), 2},
   {"smsigma", (DL_FUNC) &F77_SUB(smsigma), 7},
   {"tensres", (DL_FUNC) &F77_SUB(tensres), 8},
   {NULL, NULL, 0}
};

void
R_init_dti(DllInfo *info)
{
   R_registerRoutines(info, cEntries, NULL, fortranEntries, NULL);
}
