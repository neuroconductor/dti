#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_SUB(awsadchi)();
void F77_SUB(awslchi2)();
void F77_SUB(awslgaus)();
void F77_SUB(awsvchi)();
void F77_SUB(d2rall)();
void F77_SUB(mediansm)();
void F77_SUB(paramw3)();
void F77_SUB(projdt2)();
void F77_SUB(r2dall)();
void F77_SUB(rho2d0)();
void F77_SUB(smsigma)();
void F77_SUB(tensres)();
void F77_SUB(thcorr)();
void F77_SUB(k456krb)();
void F77_SUB(hg1f1)();
void F77_SUB(inhmtfun)();
void F77_SUB(inhmtgrd)();
void dtens( int * n1, double * param, double * sig_in, int * ngrad, double * btb_in, 
	    double * sdcoef, double * sig_tmp, double * vinv_tmp, int * maxit, double * reltol);
static R_CMethodDef cEntries[] = {
   {"dtens", (DL_FUNC) &dtens, 10},
   {NULL, NULL, 0, NULL}
};
static R_FortranMethodDef fortranEntries[] = {
   {"awsadchi", (DL_FUNC) &F77_SUB(awsadchi), 17},
   {"awsvchi", (DL_FUNC) &F77_SUB(awsvchi), 15},
   {"awslchi2", (DL_FUNC) &F77_SUB(awslchi2), 23},
   {"awslgaus", (DL_FUNC) &F77_SUB(awslgaus), 15},
   {"D2Rall", (DL_FUNC) &F77_SUB(d2rall), 3},
   {"mediansm", (DL_FUNC) &F77_SUB(mediansm), 10},
   {"paramw3", (DL_FUNC) &F77_SUB(paramw3), 5},
   {"projdt2", (DL_FUNC) &F77_SUB(projdt2), 9},
   {"R2Dall", (DL_FUNC) &F77_SUB(r2dall), 3},
   {"rho2D0", (DL_FUNC) &F77_SUB(rho2d0), 2},
   {"smsigma", (DL_FUNC) &F77_SUB(smsigma), 7},
   {"tensres", (DL_FUNC) &F77_SUB(tensres), 8},
   {"thcorr", (DL_FUNC) &F77_SUB(thcorr), 8},
   {"k456krb", (DL_FUNC) &F77_SUB(k456krb), 4},
   {"hg1f1", (DL_FUNC) &F77_SUB(hg1f1), 5},
   {"inhmtfun", (DL_FUNC) &F77_SUB(inhmtfun), 11},
   {"inhmtgrd", (DL_FUNC) &F77_SUB(inhmtfun), 15},
   {NULL, NULL, 0}
};

void
R_init_dti(DllInfo *info)
{
   R_registerRoutines(info, cEntries, NULL, fortranEntries, NULL);
}
