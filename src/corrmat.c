#include<math.h>
#define ELEM(mat,i,j,nrow) mat[i+j*nrow]

void corrmat(const double *design, const double *rho,
	     const double *alpha, const double *sigma2,
	     const double *nug, const int *nrow,
	     const int* ncol, double *corr)
{
  int i, j, k, n, p;
  double temp, sig2, dist, nugget;
  n = *nrow;
  p = *ncol;
  sig2 = *sigma2;
  nugget = *nug;
  for(i=0; i<n; ++i)
  {
    ELEM(corr,i,i,n) = sig2 + nugget;
    for(j=0; j<i; ++j)
    {
      temp = 1.0;
      for(k=0; k<p; ++k)
      {
	dist = fabs(ELEM(design,i,k,n) - ELEM(design,j,k,n));
	dist = 4.0*pow(dist,*alpha);
	temp *= pow(rho[k],dist);
      }
      ELEM(corr,i,j,n) = ELEM(corr,j,i,n) = sig2*temp;
    }
  }
}

void crosscorrmat(const double *design1, const double *design2, const double *rho,
		  const double *alpha, const int *nrow1, const int *nrow2,
		  const int* ncol, double *corr)
{
  int i, j, k, n1, n2, p;
  double temp, dist;
  n1 = *nrow1;
  n2 = *nrow2;
  p = *ncol;
  for(i=0; i<n1; ++i)
  {
    for(j=0; j<n2; ++j)
    {
      temp = 1.0;
      for(k = 0; k < p; ++k)
      {
	dist = fabs(ELEM(design1,i,k,n1) - ELEM(design2,j,k,n2));
	dist = 4.0 * pow(dist,*alpha);
	temp *= pow(rho[k],dist);
      }
      ELEM(corr,i,j,n1) = temp;
    }
  }
}
