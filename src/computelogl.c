#include<R.h>
#include<Rmath.h>
#include<R_ext/BLAS.h>
#include<R_ext/Lapack.h>
#include<math.h>
#include<stdlib.h>
#include"corrmat.h"

#define ELEM(mat,i,j,nrow) mat[i+j*nrow]
#define SQ(x) ((x)*(x))
static const char* uplo = "Upper";

double log_determinant_chol(const double *mat, int n)
{
  int i;
  double logdet=0.0;
  for(i=0; i<n; ++i) logdet += log(ELEM(mat,i,i,n));
  logdet *= 2.0;
  return logdet;
}
double reduce_sum(const double *vec, int n)
{
  int i;
  double sum=0.0;
  for(i=0; i<n; ++i) sum+=vec[i];
  return sum;
}
double symat_sum(const double *mat, int n)
{
  int i, j;
  double sum = 0.0;
  for(i=0; i<n; ++i)
  {
    sum += ELEM(mat,i,i,n);
    for(j=i+1; j<n; ++j)
      sum += 2.0 * ELEM(mat,i,j,n);
  }
  return sum;
}
void vec_diff_scalar(const double *vin, double scalar, int n, double* vout)
{
  int i;
  for(i = 0; i < n; ++i)
    vout[i] = vin[i] - scalar;
}
void theta2rho(const double *theta, int n, double* rho)
{
  int i;
  for(i=0; i<n; ++i) rho[i] = exp(-0.25*exp(theta[i]));
}
void computelogl(const double *rho, const double *design,
		 const double *response, const double *alpha,
		 const double *cond, const int *N, const int *p,
		 double *neglogl)
{
  int nsq, info;
  int incx=1, incy=1;
  double sigma2, betahat, logdetR, dn;
  double alpblas=1.0, betablas=0.0;
  double *corr, *soly, *dev;
  nsq = (*N)*(*N);
  corr = (double*) malloc(sizeof(double)*nsq);
  sigma2 = 1.0;
  corrmat(design,rho,alpha,&sigma2,cond,N,p,corr);
  F77_CALL(dpotrf)(uplo,N,corr,N,&info);
  if(info)
  {
    free(corr);
    error("bad Cholesky decomp (info=%d)", info);
  }
  logdetR = log_determinant_chol(corr,*N);
  F77_CALL(dpotri)(uplo,N,corr,N,&info);
  soly = (double*) malloc(sizeof(double)*(*N));
  F77_CALL(dsymv)(uplo,N,&alpblas,corr,N,response,&incx,&betablas,soly,&incy);
  betahat = reduce_sum(soly,*N)/symat_sum(corr,*N);
  dev = (double*) malloc(sizeof(double)*(*N));
  vec_diff_scalar(response,betahat,*N,dev);
  /* reuse soly and sigma2 */
  F77_CALL(dsymv)(uplo,N,&alpblas,corr,N,dev,&incx,&betablas,soly,&incy);
  sigma2 = F77_CALL(ddot)(N,dev,&incx,soly,&incy);
  dn = (double)(*N);
  sigma2 /= dn;
  *neglogl = dn*log(sigma2)+logdetR;
  free(corr);
  free(soly);
  free(dev);
}

void computeloglMerr(const double* theta, const double *design,
		     const double* response, const double *alpha,
		     const double* cond, const int* N, const int *np,
		     double *neglogl)
{
  int nsq, info, p;
  int incx = 1, incy = 1;
  double sigma2, sigma2eps, nug, betahat, logdetR;
  double *rho, *corr, *soly, *dev;
  double alpblas = 1.0, betablas = 0.0;
  nsq = (*N)*(*N);
  p = *np - 2;
  corr = (double*) malloc(sizeof(double)*nsq);
  rho = (double*) malloc(sizeof(double)*p);
  theta2rho(theta, p, rho);
  sigma2 = exp(theta[*np-2]);
  sigma2eps = exp(theta[*np-1]);
  nug = sigma2eps + (*cond);
  corrmat(design,rho, alpha, &sigma2, &nug, N, &p, corr);
  F77_CALL(dpotrf)(uplo,N,corr,N,&info);
  if(info)
  {
    free(corr);
    free(rho);
    error("bad Cholesky decomp (info=%d)", info);
  }
  logdetR = log_determinant_chol(corr,*N);
  F77_CALL(dpotri)(uplo,N,corr,N,&info);
  soly = (double*) malloc(sizeof(double)*(*N));
  F77_CALL(dsymv)(uplo,N,&alpblas,corr,N,response,&incx,&betablas,soly,&incy);
  betahat = reduce_sum(soly,*N)/symat_sum(corr,*N);
  dev = (double*) malloc(sizeof(double)*(*N));
  vec_diff_scalar(response,betahat,*N,dev);
  /* reuse soly */
  F77_CALL(dsymv)(uplo,N,&alpblas,corr,N,dev,&incx,&betablas,soly,&incy);
  sigma2 = F77_CALL(ddot)(N,dev,&incx,soly,&incy);
  *neglogl = 0.5*sigma2 + 0.5*logdetR;
  free(corr);
  free(rho);
  free(soly);
  free(dev);
}
/* prediction for only one point, achieve better performance for the optim
 of EI */
void calcEI(const double *newx, const double *design,
	    const double *solres, const double *rho,
	    const double *betahat, const double *sigma2,
	    const double *sigma2eps, const double *alpha,
	    const double *Einv, const double *fmin,
	    const int *N, const int *p,
	    double *yhat, double *mse, double *nei)
{
  int n0 = 1;
  int nsq = (*N)*(*N);
  int incx = 1, incy = 1;
  double *Rpred, *solER;
  double cres, quadER, bilER1, quadE1;
  double part1, part2, pmse;
  double alpblas = 1.0, betablas = 0.0;
  double shat, norm, ei;
  Rpred = (double*) malloc(sizeof(double) * (*N));
  crosscorrmat(newx, design, rho, alpha, &n0, N, p, Rpred);
  cres = F77_CALL(ddot)(N, Rpred, &incx, solres, &incy);
  *yhat = *betahat + (*sigma2) * cres;

  solER = (double*) malloc(sizeof(double) * (*N));
  F77_CALL(dsymv)(uplo,N,&alpblas,Einv,N,Rpred,&incx,&betablas,solER,&incy);
  quadER = F77_CALL(ddot)(N,Rpred,&incx,solER,&incy);
  bilER1 = reduce_sum(solER, *N);
  quadE1 = reduce_sum(Einv,nsq);
  part1 = SQ(*sigma2)*quadER;
  part2 = SQ(1.0-(*sigma2)*bilER1)/quadE1;
  pmse = (*sigma2)+(*sigma2eps) - part1 + part2;
  *mse = (pmse>0.0)? pmse: 0.0;
  free(Rpred);
  free(solER);
  /* calculate ei */

  if(*mse <= 0.0)
  {
    *nei = 0.0;
    return;
  }
  shat = sqrt(*mse);
  norm = ((*fmin)-(*yhat))/shat;
  ei = norm*shat*pnorm(norm,0.0,1.0,1,0)+shat*dnorm(norm,0.0,1.0,0);
  *nei = -ei;
}
