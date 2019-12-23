#ifndef __header_Blas_h__
#define __header_Blas_h__
#include "Config.h"
#ifndef MNS_WITH_INTEL_LAPACK
#error  undefined macro MNS_WITH_INTEL_LAPACK
#endif
     
#if !MNS_WITH_INTEL_LAPACK

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
#else
#define _WIN64 0
#if __MNS_ILP__
#define MKL_ILP64 1
#else
#undef MKL_ILP64
#endif
#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#endif

#if ! MNS_WITH_INTEL_LAPACK

extern integer 		idamax_	(integer *	n, 
				 doublereal *	dx, 
				 integer *	incx);
extern int 		dcopy_	(integer *	n, 
				 doublereal *	dx, 
				 integer *	incx, 
				 doublereal *	dy, 
				 integer *	incy);

extern doublereal 	ddot_	(integer*,
				 doublereal*,
				 integer*,
				 doublereal*,
				 integer*);
extern int 		dger_	(integer *m, 
				 integer *n, 
				 doublereal *alpha, 
				 doublereal *x, 
				 integer *incx, 
				 doublereal *y, 
				 integer *incy, 
				 doublereal *a, 
				 integer *lda);

extern int 		dgemm_(char *transa,
			       char *transb,
			       integer *m,
			       integer *n,
			       integer *k,
			       doublereal *alpha,
			       doublereal *a,
			       integer *lda,			       
			       doublereal *b,
			       integer *ldb,
			       doublereal *beta,
			       doublereal *c__,			       
			       integer *ldc);

extern int 		dscal_(integer *n,
			       doublereal *da,
			       doublereal *dx,			       
			       integer *incx);

extern int 		dgemv_(char *trans,
			       integer *m,
			       integer *n,
			       doublereal * alpha,
			       doublereal *a,
			       integer *lda,
			       doublereal *x,
			       integer *incx,
			       doublereal *beta,
			       doublereal *y,
			       integer *incy);

extern doublereal 	dnrm2_(integer *n,
			       doublereal *x,
			       integer *incx);

extern int 		daxpy_(integer *n,
			       doublereal *da,
			       doublereal *dx,			       
			       integer *incx,
			       doublereal *dy,
			       integer *incy);

extern doublereal 	dasum_(integer *n,
			       doublereal *dx,
			       integer *incx);

#define Blas_idamax(n,dx,incx) 			idamax_((integer*)n, (doublereal*)dx, (integer*)incx)
#define Blas_dcopy(n, dx, incx, dy, incy) 	dcopy_((integer*)n, (doublereal*)dx, (integer*)incx, (doublereal*)dy, (integer*)incy)
#define Blas_ddot(n,dx,incx,dy,incy) 		ddot_((integer*)n,(doublereal*)dx,(integer*)incx,(doublereal*)dy,(integer*)incy)
#define Blas_dnrm2(n, x, incx) 			dnrm2_((integer*)n, (doublereal*)x, (integer*)incx)
#define Blas_daxpy(n, da, dx, incx, dy, incy) 	daxpy_((integer*)n, (doublereal*)da, (doublereal*)dx, (integer*)incx, (doublereal*)dy, (integer*)incy)
#define Blas_dasum(n, dx, incx) 		dasum_((integer*)n, (doublereal*)dx, (integer*)incx)
#define Blas_dscal(n, da, dx,  incx) 		dscal_((integer*)n, (doublereal*)da, (doublereal*)dx,  (integer*)incx)
#define Blas_dger(m, n, alpha, x, incx, y, incy, a, lda)		\
  dger_((integer*)m, (integer*)n, (doublereal*)alpha, (doublereal*)x, (integer*)incx, (doublereal*)y, (integer*)incy, (doublereal*)a, (integer*)lda)
#define Blas_dgemm(transa, transb, m,n, k,alpha, a, lda,  b, ldb, beta, c__,  ldc) \
  dgemm_((char*)transa,(char*)transb, (integer*)m, (integer*)  n, (integer*)k, (doublereal*)alpha, (doublereal*)a, (integer*)lda,  (doublereal*)b, (integer*)ldb, (doublereal*)beta, (doublereal*)c__,  (integer*)ldc)
#define Blas_dgemv(trans, m, n,  alpha, a, lda, x, incx, beta, y, incy) \
  dgemv_((char*)trans, (integer*)m, (integer*)n, (doublereal*) alpha, (doublereal*)a, (integer*)lda, (doublereal*)x, (integer*)incx, (doublereal*)beta, (doublereal*)y, (integer*)incy)

#else

#define Blas_daxpy(n, da, dx, incx, dy, incy) 			daxpy((pI)n, (pR)da, (pR)dx, (pI)incx, (pR)dy, (pI)incy)
#define Blas_dscal(n, da, dx,  incx) 				dscal((pI)n, (pR)da, (pR)dx,  (pI)incx)
#define Blas_dnrm2(n, x, incx) 					dnrm2((pI)n, (pR)x, (pI)incx)
#define Blas_dasum(n, dx, incx) 				dasum((pI)n, (pR)dx, (pI)incx)
#define Blas_idamax(n,dx,incx) 					idamax((pI)n, (pR)dx, (pI)incx)
#define Blas_dcopy(n, dx, incx, dy, incy) 			dcopy((pI)n, (pR)dx, (pI)incx, (pR)dy, (pI)incy)
#define Blas_ddot(n,dx,incx,dy,incy) 				ddot((pI)n,(pR)dx,(pI)incx,(pR)dy,(pI)incy)
#define Blas_dger(m, n, alpha, x, incx, y, incy, a, lda) 	dger((pI)m, (pI)n, (pR)alpha, (pR)x, (pI)incx, (pR)y, (pI)incy, (pR)a, (pI)lda)
#define Blas_dgemm(transa, transb, m,n, k,alpha, a, lda,  b, ldb, beta, c__,  ldc) \
  dgemm((char*)transa,(char*)transb, (pI)m, (pI)  n, (pI)k, (pR)alpha, (pR)a, (pI)lda,  (pR)b, (pI)ldb, (pR)beta, (pR)c__,  (pI)ldc)
#define Blas_dgemv(trans, m, n,  alpha, a, lda, x, incx, beta, y, incy) \
  dgemv((char*)trans, (pI)m, (pI)n, (pR) alpha, (pR)a, (pI)lda, (pR)x, (pI)incx, (pR)beta, (pR)y, (pI)incy)

#endif

#endif
