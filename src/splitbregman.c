/* 
 * (c) 2013 Tobias Abenius
 *
 * References:
 *   1. Ye, G. & Xie, X. Split Bregman method for large scale fused Lasso. (2010).
 *
 * Speed-up compared to Matlab is 2.8
 *
 * FIXME: 
 *   OpenMP support: as below does not add any speed,
 *          try changing BLAS axpy,scal,dgemv etc to parallell for-loops
 *   Speed: [ ]   Change from .C to .External to avoid data copying time
 *          [x]   starting values a = beta 4.3/4.7
 *          [x]   warmboot with LASSO  2.1/5.9 regardless of lamda1
 *          [x]   randomize starting values 4.7s/5.6s
 *          [x]   support several lambda1 (4.1s), warmboot (3.8s if decreasing lambda values)  5.7/8.3
 *
 * Compilation:
 *   R CMD SHLIB splitbregman.c
 *
 * Usage:
 *   dyn.load('splitbregman.so')
 *   .C('fused',as.integer(c(m,N,p,nolam)),as.integer(c(0,1000),as.double(L),as.double(X),as.double(y),as.double(c(lambda1,lambda2)),as.double(array(0.0,dim=p)))
 *
 *  vim: sw=2:fdm=syntax:foldlevel=3
 */
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

#define USEOPENMP 0
#if USEOPENMP
#  include <omp.h>
#endif

#define nrm2(a,p) F77_CALL(dnrm2)(&p, a, &onei)
#define copy(a,b,p) F77_CALL(dcopy)(&p, a, &onei, b, &onei)
#define scal(a,b,p) F77_CALL(dscal)(&p, a, b, &onei)
#define axpy(a,x,y,p) F77_CALL(daxpy)(&p, a, x, &onei, y, &onei)

#define max(a,b) ( (a)>(b)?(a):(b) )
#define min(a,b) ( (a)<(b)?(a):(b) )
#define newdouble(sz) (double*)R_alloc(sz,sizeof(double));

int my_inverse(double* A, int N) {
  int *IPIV = Calloc(N+1,int);
  int LWORK = N*N;
  double *WORK = Calloc(LWORK,double);
  int INFO;
//  R Calloc does not return if process is out of memory
//  if (!WORK || !IPIV) error(__FILE__ ": Out of memory!");

  F77_NAME(dgetrf)(&N,&N,A,&N,IPIV,&INFO);
  F77_NAME(dgetri)(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
  Free(IPIV);
  Free(WORK);
  return INFO;
}
double my_softthresh(double x,double t) {
  double v = fabs(x) - t;
  double ret = 0.;
  if (v  > 0.) {
    if (x >= 0.) ret = v;
    else ret = -v;
  } 
  return ret;
}
double my_norm1(double *x,int p) {
  double ax = 0.0;
  for(int i=0; i < p; ++i) ax += fabs(x[i]);
  return ax;
}
void fused(int *dim,int *param, double *L,double *X,double *y,double *lambdaLASSO,double *lambdaFUSED,double *beta) {
  // arguments
  int N = dim[2];
  int p = dim[3];
  int m = dim[0];
  int K = dim[1];
  int noL1 = dim[4];
  int noL2 = dim[5];
  int trace = param[0];

  if (trace & 1) REprintf("m=%d K=%d N=%d p=%d noL1=%d noL2=%d\n",m,K,N,p,noL1,noL2);
  //p = K*p;

  // parameters
  double mu1 = 1.0;
  double mu2 = 1.0;
  double delta1 = mu1;
  double delta2 = mu2;
  int maxit = param[1];

  double mu1inv = 1.0/mu1;
  double mu2inv = 1.0/mu2;

  int onei = 1;
  double one = 1.0;
  double mone = -1.0;
  double zero = 0.0;
  int i,j;
  int it;

  // temporary variable
  double *tmpN = newdouble(N);
  double *tmpP = newdouble(p);
  double *tmpu = newdouble(p);
  double *tmpv = newdouble(m);

  // for convergence check
  double *epsilon = newdouble(N);
  double conv;
  double theta;
  double thetaold = 1;

  double *a = newdouble(p);
  double *u = newdouble(p);
  for(i=0; i < p; ++i) {
    a[i] = beta[i]; //norm_rand();
    u[i] = 0.0;
  }

  double *b = newdouble(m);
  double *v = newdouble(m);
  for(i=0; i < m; ++i) {
    b[i] = norm_rand();
    v[i] = 0.0;
  }

  double *betap;

  double *XtX = newdouble(p*p);
  F77_CALL(dgemm)("T","N",&p,&p,&N,&one, X, &N, X, &N, &zero, XtX, &p);

  if (trace & 4) REprintf("XtX: %f %f\n",XtX[0],XtX[1]);

  double *Xty = newdouble(p);
  F77_NAME(dgemv)("T", &N, &p, &one, X, &N, y, &onei, &zero, Xty, &onei);
  if (trace & 4) REprintf("Xty: %f %f %f %f\n",Xty[0],Xty[1],Xty[2],Xty[3]);

  double *LtL = newdouble(p*p); // (p x m) x (m x p)
  // DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  F77_CALL(dgemm)("T","N",   &p,&p,&m,&one,  L, &m,  L, &m, &zero, LtL, &p);
  if (trace & 4) REprintf("L: %f %f\n",L[0],L[1]);
  if (trace & 4) REprintf("LtL: %f %f\n",LtL[0],LtL[1]);

  double *Lbeta = newdouble(m); // (m x p) x (p x 1) -> m x 1

  double *Dinv = newdouble(p*p);
  // D = X'X + mu1*I + mu2*L'L;
  memcpy(Dinv, XtX, sizeof(double)*p*p);
  for(i=0; i < p; ++i) {
    Dinv[p*i+i] += mu1;
    for(j=0; j < p; ++j) {
      Dinv[p*i+j] += mu2*LtL[p*i+j];
    }
  }
  int inverr = my_inverse(Dinv, p);
  if(inverr != 0) {
    if (inverr > 0) {
      REprintf("inverse: Singular matrix at %d\n", inverr);
      error(__FILE__ ": inverse: Singular matrix");
    } else {
      REprintf(":inverse: illegal parameter to dgetri at %d\n",-inverr);
      error(__FILE__ ":inverse: illegal parameter to dgetri!");
    }
  }

  if (trace & 4) REprintf("Dinv: %f %f %f %f %f\n",Dinv[0],Dinv[1],Dinv[2],Dinv[3],Dinv[4]);

  int th_id, nthreads;
  int r = 0;
  for (int l2i=0; l2i < noL2; ++l2i) {
  for (int l1i=0; l1i < noL1; ++l1i,++r) {
    if ( r > 0) memcpy(&beta[r*p], &beta[(r-1)*p], p*sizeof(float));
    betap = &(beta[r*p]);
  for (it=0; it < maxit; ++it) {
    // (1)
    // beta = Dinv * (Xty + mu1*(a - u/mu1) + mu2*L'*(b - v/mu2));
    //                      \_____,______/    \_______,_________/
    //                           (i)                 (ii)
    //        \________________________,________________________/
    //                               (iii)
    if (trace & 4) REprintf("it: %d\n",it+1);
#   if USEOPENMP
#     pragma omp parallel private(th_id) shared(nthreads)
      th_id = omp_get_thread_num();
#     pragma omp sections
      {
#       pragma omp section
        {
#endif
	  if (trace & 2) REprintf("[%d]: section A, it %d\n",th_id,it);
	  // (i) -> a
	  copy(u, tmpu, p);
	  scal(&mu1inv, tmpu, p);
	  axpy(&mone, tmpu, a, p);
	  scal(&mu1, a, p);
	  if (trace & 4) REprintf("(i): %.6f %.6f %.6f %.6f\n",a[0],a[1],a[2],a[3]);
#   if USEOPENMP
        }
#       pragma omp section
        {
#endif
	  if (trace & 2) REprintf("[%d]: section B, it %d\n",th_id,it);
	  // (ii) -> b
	  copy(v, tmpv, m);
	  scal(&mu2inv, tmpv, m);
	  axpy(&mone, tmpv, b, m);
	  copy(b, tmpv,m); // here was a daemon, dgemv does not seem to handle if x identical y
	  if (trace & 4) REprintf("(b - v/mu2): %.6f %.6f %.6f %.6f\n",b[0],b[1],b[2],b[3]);
	  if (trace & 4) REprintf("L: %.6f %.6f %.6f %.6f\n",L[0],L[1],L[2],L[3]);
	  //       DGEMV(TRANS,M,  N, ALPHA, A,LDA, X, INCX,  BETA,  Y, INCY) 
	  F77_NAME(dgemv)("T", &m, &p, &mu2, L, &m, tmpv, &onei, &zero, tmpP, &onei);
	  if (trace & 4) REprintf("(ii): %.6f %.6f %.6f %.6f\n",tmpP[0],tmpP[1],tmpP[2],tmpP[3]);
#if USEOPENMP
        }
      } // sections
#endif

    // (i) + (ii) -> a
    axpy(&one, tmpP, a, p);
    if (trace & 4) REprintf("(i) + (ii): %.6f %.6f %.6f %.6f\n",a[0],a[1],a[2],a[3]);
    // (iii) -> beta
    axpy(&one, Xty, a, p);
    if (trace & 4) REprintf("Xty + (i) + (ii): %.6f %.6f %.6f %.6f\n",a[0],a[1],a[2],a[3]);
    //       DGEMV(TRANS,M,  N, ALPHA, A,   LDA, X, INCX,  BETA,  Y, INCY) 
    F77_NAME(dgemv)("N", &p, &p, &one, Dinv, &p, a, &onei, &zero, betap, &onei);
    if (trace & 4) REprintf("beta: %.6f %.6f %.6f %.6f\n",betap[0],betap[1],betap[2],betap[3]);

#if USEOPENMP
#   pragma omp parallel private(th_id) shared(nthreads)
    {
#endif
#if USEOPENMP
      th_id = omp_get_thread_num();
      /*
      if (th_id == 0) {
	nthreads = omp_get_num_threads();
	if (trace & 2) REprintf("n.o. openmp threads: %d\n",nthreads);
      }
      */
#     pragma omp sections
      {
#       pragma omp section
	{
	  if (trace & 2) REprintf("[%d]: section 1, it %d\n",th_id,it);
#endif
	  // (2) a = ST(beta + u/mu1, lambda1/mu1)
	  copy(u,tmpu,p); // XXX: SEE ABOVE tmpu
	  scal(&mu1inv, tmpu, p);
	  axpy(&one, betap, tmpu, p);
	  for(i=0; i < p; ++i) a[i] = my_softthresh(tmpu[i], lambdaLASSO[l1i]/mu1);
	  if (trace & 4) REprintf("a: %.6f %.6f %.6f %.6f\n",a[0],a[1],a[2],a[3]);

	  // (3) u = u + delta1*(beta - a);    
	  copy(betap,tmpP,p);
	  axpy(&mone,a,tmpP,p);
	  axpy(&delta1,tmpP,u,p);
	  if (trace & 4) REprintf("u: %.6f %.6f %.6f %.6f\n",u[0],u[1],u[2],u[3]);
#if USEOPENMP
	} // section
#       pragma omp section
	{
	  if (trace & 2) REprintf("[%d]: section 2, it %d\n",th_id,it);
#endif
	  // L beta
//                  DGEMV(TRANS,M,  N, ALPHA, A, LDA, X, INCX,  BETA,    Y,   INCY) 
	  F77_NAME(dgemv)("N", &m, &p, &one, L, &m, betap, &onei, &zero, Lbeta, &onei);
	  if (trace & 4) REprintf("L beta: %.6f %.6f %.6f %.6f\n",Lbeta[0],Lbeta[1],Lbeta[2],Lbeta[3]);

	  // (4) b = ST(Lbeta + v/mu2, lambda2/mu2)
	  copy(v,tmpv,m);
	  scal(&mu2inv, tmpv, m);   // XXX: SEE ABOVE tmpv, this and row above could be removed if a new temporary is introduced
	  axpy(&one, Lbeta, tmpv, m);
	  for(i=0; i < m; ++i) b[i] = my_softthresh(tmpv[i], lambdaFUSED[l2i]/mu2);
	  if (trace & 4) REprintf("b: %.6f %.6f %.6f %.6f\n",b[0],b[1],b[2],b[3]);

	  // (5) v = v + delta2*(Lbeta - b);
	  copy(Lbeta,tmpv,m);
	  axpy(&mone,b,tmpv,m);
	  axpy(&delta2,tmpv,v,m);
	  if (trace & 4) REprintf("v: %.6f %.6f %.6f %.6f\n",v[0],v[1],v[2],v[3]);
#if USEOPENMP
	} // section
      }
    }
#endif
    ///// C O N V E R G E N C E
    F77_NAME(dgemv)("N", &N, &p, &one, X, &N, betap, &onei, &zero, epsilon, &onei);
    copy(y, tmpN, N);
    axpy(&mone, epsilon, tmpN, N);
    thetaold = theta;
    theta = nrm2(tmpN,N)/2.0 + lambdaLASSO[l1i] * my_norm1(betap,p) + lambdaFUSED[l2i] * my_norm1(Lbeta,m);
    if (isnan(theta)) {
      error("NaN occured in convergence criterion");
    }
    if (it > 0) {
      conv = fabs(thetaold - theta) / thetaold;
      if (trace & 1) REprintf("[%d] convergence: theta=%f, relE=%f\n", it, theta, conv);
      if (conv <= 1e-5) break;
    }
    ////////
  }
  }
  }
}
