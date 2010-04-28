#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>


#include <sys/time.h>
#include <omp.h>

typedef struct
{
  int ngrad;
  double *siq;
  double *grad;
  int p;
} optimex;

typedef struct
{
  int ngrad;
  double *siq;
  double *grad;
  double *w;
} optimexpl;

// double F77_NAME(dotprod3)(double *, double *);
// double z1 = F77_CALL(dotprod3)(dir, test);
//void F77_NAME(nnls)(double *, int, int, int, double *, double *, double *, double *, double *, int * , int *);
void F77_NAME(nnls)(double *a, int *mda, int *m, int *n, double *b, double *X, double *rnorm, double *w, double *zz, int *index, int *mode);
// quicksort

void swapi(int *a, int *b)
{
  int t=*a; *a=*b; *b=t;
}
void swapd(double *a, double *b)
{
  double t=*a; *a=*b; *b=t;
}
void sort(double arr[], int ind[], int beg, int end)
{
  if (end > beg + 1)
  {
    double piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swapd(&arr[l], &arr[--r]);
        swapi(&ind[l], &ind[--r]);
    }
    swapd(&arr[--l], &arr[beg]);
    swapi(&ind[--l], &ind[beg]);
    sort(arr, ind, beg, l);
    sort(arr, ind, r, end);
  }
}

// end quicksort

void paroforient(double *dir, double *angles)
{
  angles[0] = acos(dir[2]);
  double sth = sin(angles[0]);
  angles[1] = 0.;

  if (sth < 1e-8)
  {
    angles[0] = 0.;
  }
  else
  {
    double z = dir[0]/sth;
    angles[1] = (abs(z) >= 1.) ? ((z<0.) ? 0. : M_PI) : acos(z)*sign(dir[1]);
    if (angles[1] < 0.) angles[1] += M_2PI;
  }
  return;
}

inline int get_ind2d(int i, int j, int n1)
{
  return i + j*n1;
}

inline int get_ind3d(int i, int j, int k, int n1, int n2)
{
  return i + j*n1 + k*n1*n2;
}

inline int get_ind4d(int i, int j, int k, int l, int n1, int n2, int n3)
{
  return i + j*n1 + k*n1*n2 + l*n1*n2*n3;
}

inline int get_ind5d(int i, int j, int k, int l, int m, int n1, int n2, int n3, int n4)
{
  return i + j*n1 + k*n1*n2 + l*n1*n2*n3 + m*n1*n2*n3*n4;
}

void smoothsi(int *r1, int *r2, int *r3, double *si, int *n, int *sneighbors, int *nneighbors, double *smoothedsi)
{
  int n1 = *r1, n2 = *r2, n3 = *r3, ngrad = *n, nn = *nneighbors;
  double w[nn];
  int i1, i2, i3, i, j;
  double sw = 0.;
  double z;

  for (j = 0; j < nn; j++)
  {
    w[j] = 1./(j+1.);
    sw += w[j];
  }

  for (i1 = 0; i1 < n1; i1++)
  {
    for (i2 = 0; i2 < n2; i2++)
    {
      for (i3 = 0; i3 < n3; i3++)
      {
        for (i = 0; i < ngrad; i++)
        {
          z = 0.;
          for (j = 0; j < nn; j++) 
            z += w[j]*si[get_ind4d(i1, i2, i3, sneighbors[get_ind2d(j, i, nn)]-1, n1, n2, n3)];
          smoothedsi[get_ind4d(i1, i2, i3, i, n1, n2, n3)] = z/sw;
        }
      }
    }
  }
  return;
}

void getweights(int m, double *cw, double *z)
{
  double ew;
  double sw = 0.0;

  for (int i = 0; i < m; i++)
  {
    // use log and exp to have free optim problem
    ew = exp(cw[i]);
    z[i] = (1. - sw) * ew/(1.+ew);
    sw += z[i];
  }
  z[m] = 1. - sw;
  return;
}

double fn1(int n, double *par, void *ex)
{
  optimex ext = *((optimex*)ex);

  int m = (n-1)/3;     // order of mix tensor model
  double z[ext.ngrad]; // value of signal fct. for each grad according to current model
  for (int k = 0; k < ext.ngrad; k++)
    z[k] = 0;
  double dir[3];       // used for direction vector calculated from angles in par
  double c1 = exp(par[0]), c2 = exp(par[1]), sw = 0, w, ew, sth, z1 = 0;
  int i, j;

  double erg = 0; // result

  for (i = 0; i < m; i++)
  {
    if (i == m-1 )
    {
      w = 1.0 - sw;
    }
    else
    {
      // use log and exp to have free optim problem
      ew = exp(par[3*i+4]);
      w = (1.0 - sw)*ew/(1.0+ew);
      sw += w;
    }
    sth = sin(par[i*3+2]);
    dir[0] = sth*cos(par[i*3+3]);
    dir[1] = sth*sin(par[i*3+3]);
    dir[2] = cos(par[i*3+2]);
    for (j = 0; j < ext.ngrad; j++)
    {
      z1 = dir[0]*ext.grad[j] + dir[1]*ext.grad[ext.ngrad+j] + dir[2]*ext.grad[2*ext.ngrad+j]; 
      z[j] += (ext.p == 0) ? w*exp(-c2 - c1*z1*z1) : w*exp(-ext.p*log(1 + ( c2 + c1*z1*z1 )/ext.p));
    }
  }

  // calculate residual variance
  for (j = 0; j < ext.ngrad; j++)
  {
    z1 = ext.siq[j] - z[j];
    erg += z1*z1;
  }

  // finished
  return erg;
}

double fnpl(int n, double *par, void *ex)
{
  optimexpl ext = *((optimexpl*)ex);

  int m = (n-1)/2;     // order of mix tensor model
  double z[ext.ngrad*m]; // value of signal fct. for each grad according to current model
  double siq[ext.ngrad];
  for (int k = 0; k < ext.ngrad; k++)
  {
    for (int kk = 0; kk < m; kk++)
      z[kk*ext.ngrad + k] = 0;
    siq[k] = ext.siq[k];
  }
  double dir[3];       // used for direction vector calculated from angles in par
  double c1 = exp(par[0]), sth, z1;
  int i, i2, j;
  int ind[10], mode = 0;
  double work2[10];
  double work1[ext.ngrad];

  double erg = 0; // result

  for (i = 0; i < m; i++)
  {
    i2 = 2*i;
    sth = sin(par[i2+1]);
    dir[0] = sth*cos(par[i2+2]);
    dir[1] = sth*sin(par[i2+2]);
    dir[2] = cos(par[i2+1]);
    for (j = 0; j < ext.ngrad; j++)
    {
      z1 = dir[0]*ext.grad[j] + dir[1]*ext.grad[ext.ngrad+j] + dir[2]*ext.grad[2*ext.ngrad+j]; 
      z[j + i*ext.ngrad] += exp(-c1*z1*z1);
    }
  }
 // siq will be replaced, need to copy it if C-version of optim is used
  F77_CALL(nnls)(z, &ext.ngrad, &ext.ngrad, &m, siq, ext.w, &erg, work2, work1, ind, &mode);

  // finished
  return erg;
}

/*
double fn2(int n, double *par, void *ex) {
  optimex ext = *((optimex*)ex);

  int m = (n-1)/3; // order of mix tensor model
  double z[ext.ngrad]; // value of signal fct. for each grad according to current model
  for (int k = 0; k < ext.ngrad; k++) z[k] = 0;
  double dir[3]; // used for direction vector calculated from angles in par
  double c1 = exp(par[0]), c2 = exp(par[1]), sw = 0, w, ew, sth, z1 = 0;
  int i, j;

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
      z[j] += w*exp(-par[n-1]*log(1 + ( c2 + c1*z1*z1 )/par[n-1]));
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
*/

void mixtensor(int *n, double *par, double *x, int *ngrad, double *siq, double *grad, int *maxit, double *reltol, double *Fmin){

  int fail, fncount; // failure code: zero is OK, number of calls to fct.
  optimex myoptimpar;

  myoptimpar.ngrad = *ngrad;
  myoptimpar.siq = siq;
  myoptimpar.grad = grad;
  myoptimpar.p = 0;

  nmmin(*n, par, x, Fmin, fn1,
        &fail, R_NegInf, *reltol, &myoptimpar,
        1.0, 0.5, 2.0, 0,
        &fncount, *maxit);
  return;
}

void mixture(int *method, int *r, int *mask, double *siq, int *siind, int *n, double *grad, int *sneighbors, int *maxcomp, int *ep, int *maxit, double *reltol,
             double *order, double *lev, double *mix, double *orient, double *p, double *sigma2){

  int nv = *r, ngrad = *n, mc = *maxcomp;
  int iv, mc0, lpar;
  double rss, krit, ttt;
  double siiq[ngrad];
  int fail;              // failure code for optim: zero is OK
  int fncount;           // number of calls to obj fct in optim
  double Fmin = 0.;          // minimal value of obj fct in optim
  int i, k;
  int gradind;
  double angles[2];
  double dir[3];
  double par[11];
  double cpar[11];
  double x[11];
  double tmp;
//  struct timeval tp1, tp2;
//  struct timezone tzp;
//  double timm = 0;

// #pragma omp parallel for private(i, siiq, mc0, gradind, tmp, par, dir, angles, rss, krit, k, lpar, cpar, x, fail, Fmin, ttt, fncount, reltol, maxit)
  for (iv = 0; iv < nv; iv++)
  {
    if (mask[iv] != 1)
    {
      order[iv] = 0;
      lev[get_ind2d(0, iv, 2)] = 0;
      lev[get_ind2d(1, iv, 2)] = 0;
      for (i = 0; i < mc; i++)
      {
        mix[get_ind2d(i, iv, mc)] = 0;
        orient[get_ind3d(0, i, iv, 2, mc)] = 0;
        orient[get_ind3d(1, i, iv, 2, mc)] = 0;
      }
//          if (*method == 3) p[get_ind3d(i1, i2, i3, n1, n2)] = 0;
      sigma2[iv] = 0;
    }
    else
    {
      // prepare signal
      for (i=0; i < ngrad; i++)
      {
        siiq[i] = siq[get_ind2d(iv, i, nv)];
      }
      // determine useful prime estimates
      mc0 = siind[get_ind2d(0, iv, mc+1)]; // order of modell 
      gradind = siind[get_ind2d(1, iv, mc+1)];
      tmp = -log(siq[get_ind2d(iv, gradind-1, nv)]);
      par[0] = log(tmp*.8);
      par[1] = log(tmp*.2);
      for (i = 0; i < mc0; i++)
      {
        gradind = siind[get_ind2d(i+1, iv, mc+1)];
        dir[0] = grad[get_ind2d(gradind-1, 0, ngrad)];
        dir[1] = grad[get_ind2d(gradind-1, 1, ngrad)];
        dir[2] = grad[get_ind2d(gradind-1, 2, ngrad)];
        paroforient(dir, angles);    // INDEX????

        par[3*i + 2] = angles[0];
        par[3*i + 3] = angles[1]; 
//            Rprintf("par[] %i, %f %f\n", 3*i+2, par[3*i+2], par[3*i+3]); 
      }
//         if (*method == 3) par[lpar] = 0;

      // estimate models and select
      rss = R_PosInf;
      krit = R_PosInf;
      for (k = mc0; k > 0; k--)
      {
        for (i = 0; i < 3*k+1; i++) cpar[i] = par[i];
       //if (*method == 3) cpar[3*k+1] = ??;
        lpar = (*method == 3) ? 3*k+2 : 3*k+1;
        // use log and exp to have free optim problem
        if (k>1) for (i = 0; i < k-1; i++) cpar[3*i+4] = -log(k-i-1.); // par[3*(2:k)-1] <- -log((k-1):1)
        optimex myoptimpar;
        myoptimpar.ngrad = ngrad;
        myoptimpar.siq = siiq;
        myoptimpar.grad = grad;
        switch (*method)
        { // R code guarantees method is 1, 2, 3
          case 1:
            myoptimpar.p = 0; // unused here
            nmmin(lpar, cpar, x, &Fmin, fn1,
                  &fail, R_NegInf, *reltol, &myoptimpar,
                  1.0, 0.5, 2.0, 0,
                  &fncount, *maxit);
            break;
          case 2:
            myoptimpar.p = *ep; // exp for Jian model
            nmmin(lpar, cpar, x, &Fmin, fn1,
                  &fail, R_NegInf, *reltol, &myoptimpar,
                  1.0, 0.5, 2.0, 0,
                  &fncount, *maxit);
            break;
//          case 3:
//            myoptimpar.p = 0; // unused here
//            nmmin(lpar, par, x, Fmin, fn3,
//                  &fail, R_NegInf, *reltol, &myoptimpar,
//                  1.0, 0.5, 2.0, 0,
//                  &fncount, *maxit);
        }

        if (Fmin < rss) rss = Fmin;
        ttt = Fmin + (6.*k+2.)/(ngrad - 3.*mc - 1.) * rss;
        if (ttt < krit)
        {
          krit = ttt;
          order[iv] = k;
          lev[get_ind2d(0, iv, 2)] = x[0];
          lev[get_ind2d(1, iv, 2)] = x[1];
          if (k == 1)
          {
            mix[get_ind2d(0, iv, mc)] = 1;
          }
          else
          {
            double zmm[k];  // declare maximum?
            double zm[k-1]; // declare maximum?
            for (i = 0; i < k-1; i++) zm[i] = x[3*i+4];
            getweights(k-1, zm, zmm);
            for (i = 0; i < k; i++) mix[get_ind2d(i, iv, mc)] = zmm[i];
          }

          for (i = 0; i < k; i++)
          {
            double theta = x[3*i+2];
            while (theta < 0) 
            {
            //  Rprintf("theta %f %f", theta, M_PI); 
              theta += M_PI;
            //  Rprintf("theta %f\n", theta); 
            }
            while (theta > M_PI) theta -= M_PI;
            double phi = x[3*i+3];
            while (phi < 0) phi += M_2PI;
            while (phi > M_2PI) phi -= M_2PI;
            orient[get_ind3d(0, i, iv, 2, mc)] = theta;
            orient[get_ind3d(1, i, iv, 2, mc)] = phi;
          }
//         if (*method == 3) p[i1*n2*n3 + i2*n3 + i3] = x[lpar-1];
        }
      }
      sigma2[iv] = rss/(ngrad-3.*mc0-1.);
    }
    R_CheckUserInterrupt();
  }
//  gettimeofday(&tp1, &tzp);
//  gettimeofday(&tp2, &tzp);
//  if (tp1.tv_usec > tp2.tv_usec) tp2.tv_usec += 1000000;
//  timm += tp2.tv_usec - tp1.tv_usec;
//  Rprintf("zeit %f\n", timm/1000000.);
  return;
}

void mixturepl(int *method, int *r, int *mask, double *siq, int *siind, int *n, double *grad, int *sneighbors, int *maxcomp, int *ep, int *maxit, double *reltol,
               double *order, double *lev, double *mix, double *orient, double *sigma2){

  int nv = *r, ngrad = *n, mc = *maxcomp;
  int iv, mc0, lpar, ord;
  double rss, krit, ttt, sw = 0;
  double siiq[ngrad];
  int fail;                  // failure code for optim: zero is OK
  int fncount;               // number of calls to obj fct in optim
  double Fmin = 0.;          // minimal value of obj fct in optim
  int i, k;
  int gradind;
  double angles[2];
  double dir[3];
  double par[9];
  double cpar[9];
  double x[9];
  double tmp;

  for (iv = 0; iv < nv; iv++)
  {
    if (mask[iv] != 1)
    {
      order[iv] = 0;
      lev[get_ind2d(0, iv, 2)] = 0;
      lev[get_ind2d(1, iv, 2)] = 0;
      for (i = 0; i < mc; i++)
      {
        mix[get_ind2d(i, iv, mc)] = 0;
        orient[get_ind3d(0, i, iv, 2, mc)] = 0;
        orient[get_ind3d(1, i, iv, 2, mc)] = 0;
      }
      sigma2[iv] = 0;
    }
    else
    {
      // prepare signal
      for (i=0; i < ngrad; i++) siiq[i] = siq[get_ind2d(iv, i, nv)];

      // determine useful prime estimates
      mc0 = siind[get_ind2d(0, iv, mc+1)]; // order of model: possibly constant == mc since intial estimates guarantee mc0 == mc!!
      ord = mc0 + 1;
      gradind = siind[get_ind2d(1, iv, mc+1)];
      tmp = -log(siq[get_ind2d(iv, gradind-1, nv)]);
      par[0] = log(tmp*.8);
      for (i = 0; i < mc0; i++)
      {
        if (i>0) gradind = siind[get_ind2d(i+1, iv, mc+1)];
        dir[0] = grad[get_ind2d(gradind-1, 0, ngrad)];
        dir[1] = grad[get_ind2d(gradind-1, 1, ngrad)];
        dir[2] = grad[get_ind2d(gradind-1, 2, ngrad)];
        paroforient(dir, angles);    // INDEX????

        par[2*i + 1] = angles[0];
        par[2*i + 2] = angles[1];
      }

      // estimate models and select
      rss = R_PosInf;
      krit = R_PosInf;
      for (k = mc0; k > 0; k--)
      {
        if (k < ord)
        {
          lpar = 2*k+1;
          for (i = 0; i < lpar; i++) cpar[i] = par[i];
          optimexpl myoptimpar;
          myoptimpar.ngrad = ngrad;
          myoptimpar.siq = siiq;
          myoptimpar.grad = grad;
          double w[k];
          for (i = 0; i < k; i++) w[i] = 0;
          myoptimpar.w = w;
          switch (*method)
          { // R code guarantees method is 1, 2
            case 1:
              nmmin(lpar, cpar, x, &Fmin, fnpl,
                    &fail, R_NegInf, *reltol, &myoptimpar,
                    1.0, 0.5, 2.0, 0,
                    &fncount, *maxit);
              break;
//            case 2:
//              myoptimpar.p = *ep; // exp for Jian model
//              nmmin(lpar, cpar, x, &Fmin, fnpl,
//                    &fail, R_NegInf, *reltol, &myoptimpar,
//                    1.0, 0.5, 2.0, 0,
//                    &fncount, *maxit);
//              break;
          }

          if (Fmin < rss) rss = Fmin;
          Fmin = fnpl(lpar, x, &myoptimpar); // should return the same value
          int ind[k];
          ord = 0;
          for (i = 0; i < k; i++)
          {
            ind[i] = i;
            if (myoptimpar.w[i] > 0)
            {
              sw = sw + myoptimpar.w[i];
              ord++;
            }
          }
          ttt = Fmin + (6.*ord+2.)/(ngrad - 3.*mc - 1.) * rss;
          revsort(myoptimpar.w, ind, k);
          par[0] = x[0];
          for (i = 0; i < ord; i++)
          {
            par[2*i+1] = x[2*ind[i]+1];
            par[2*i+2] = x[2*ind[i]+2];
          }

          if (ttt < krit)
          {
            krit = ttt;
            order[iv] = ord;
            lev[get_ind2d(0, iv, 2)] = exp(x[0]);
            lev[get_ind2d(1, iv, 2)] = -log(sw);
            if (ord == 1)
            {
              mix[get_ind2d(0, iv, mc)] = 1;
            }
            else
            {
              for (i = 0; i < k; i++) mix[get_ind2d(i, iv, mc)] = myoptimpar.w[i]/sw; // order?
            }

            for (i = 0; i < k; i++)
            {
              double theta = x[2*ind[i]+1];
              while (theta < 0) 
              {
                theta += M_PI;
              }
              while (theta > M_PI) theta -= M_PI;
              double phi = x[2*ind[i]+2];
              while (phi < 0) phi += M_2PI;
              while (phi > M_2PI) phi -= M_2PI;
              orient[get_ind3d(0, i, iv, 2, mc)] = theta;
              orient[get_ind3d(1, i, iv, 2, mc)] = phi;
            }
          }
        }
      }
      sigma2[iv] = rss/(ngrad-3.*mc0-1.);
    }
    R_CheckUserInterrupt();
  }
  return;
}
