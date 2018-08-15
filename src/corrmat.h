#ifndef __CORRMAT_H__
#define __CORRMAT_H__
void corrmat(const double *design, const double *rho,
	     const double *alpha, const double *sigma2,
	     const double *nug, const int *nrow,
	     const int* ncol, double *corr);

void crosscorrmat(const double *design1, const double *design2, const double *rho,
		  const double *alpha, const int *nrow1, const int *nrow2,
		  const int* ncol, double *corr);

#endif
