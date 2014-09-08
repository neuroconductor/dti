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
void F77_SUB(imtfunbv)();
void F77_SUB(imtgrdbv)();
void F77_SUB(imtfunb1)();
void F77_SUB(imtgrdb1)();
void F77_SUB(imtfunb0)();
void F77_SUB(imtgrdb0)();
void F77_SUB(mtfunbv)();
void F77_SUB(mtgrdbv)();
void F77_SUB(mtfunb1)();
void F77_SUB(mtgrdb1)();
void F77_SUB(mtfunb0)();
void F77_SUB(mtgrdb0)();
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
   {"imtfunbv", (DL_FUNC) &F77_SUB(imtfunbv), 11},
   {"imtgrdbv", (DL_FUNC) &F77_SUB(imtgrdbv), 15},
   {"imtfunb1", (DL_FUNC) &F77_SUB(imtfunb1), 12},
   {"imtgrdb1", (DL_FUNC) &F77_SUB(imtgrdb1), 16},
   {"imtfunb0", (DL_FUNC) &F77_SUB(imtfunb0), 13},
   {"imtgrdb0", (DL_FUNC) &F77_SUB(imtgrdb0), 17},
   {"mtfunbv", (DL_FUNC) &F77_SUB(mtfunbv), 11},
   {"mtgrdbv", (DL_FUNC) &F77_SUB(mtgrdbv), 14},
   {"mtfunb1", (DL_FUNC) &F77_SUB(mtfunb1), 12},
   {"mtgrdb1", (DL_FUNC) &F77_SUB(mtgrdb1), 15},
   {"mtfunb0", (DL_FUNC) &F77_SUB(mtfunb0), 13},
   {"mtgrdb0", (DL_FUNC) &F77_SUB(mtgrdb0), 16},
   {NULL, NULL, 0}
};

void
R_init_dti(DllInfo *info)
{
   R_registerRoutines(info, cEntries, NULL, fortranEntries, NULL);
}
