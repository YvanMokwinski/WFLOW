#ifndef __HEADER_NS_CONFIG_LAPACK_H__
#define __HEADER_NS_CONFIG_LAPACK_H__

#ifndef __CONFIG_LAPACK_NOGLOBAL__
#endif

#ifndef __MK_WITH_INTEL_LAPACK__
#error  undefined macro __MK_WITH_INTEL_LAPACK__
#endif

     
#if !__MK_WITH_INTEL_LAPACK__

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
#else
#define _WIN64 0


#if 0
#ifndef __MNS_ILP__
#error NS CONFIG LAPACK . H ddddddddddddd
#endif

#ifndef __MK_ILP64__
#error NS CONFIG LAPACK . H ddddddddddddd
#endif
#endif

#ifdef __MNS_ILP__
#define MKL_ILP64 1
#else
#undef MKL_ILP64
#endif

#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#endif



#if NOT __MK_WITH_INTEL_LAPACK__
integer idamax_(integer *n, doublereal *dx, integer *incx);
extern int dcopy_(integer *n, doublereal *dx, integer *incx, 
		  doublereal *dy, integer *incy);

extern doublereal ddot_(integer*,doublereal*,integer*,doublereal*,integer*);
extern int dger_(integer *m, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
			   doublereal *a, integer *lda);
extern /* Subroutine */ int dgesv_(integer *n, integer *nrhs, doublereal *a, integer 
				   *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);

extern int dgemm_(char *transa, char *transb, integer *m, integer *
		  n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
		  doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
		  integer *ldc);
extern int dscal_(integer *n, doublereal *da, doublereal *dx, 
			    integer *incx);
extern int dgemv_(char *trans, integer *m, integer *n, doublereal * alpha, 
		  doublereal *a, integer *lda, doublereal *x, integer *incx, 
		  doublereal *beta, doublereal *y, integer *incy);
extern doublereal dnrm2_(integer *n, doublereal *x, integer *incx);

extern int daxpy_(integer *n, doublereal *da, doublereal *dx, 
		  integer *incx, doublereal *dy, integer *incy);
extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
extern /* Subroutine */ int dgetrs_(char *trans, integer *n, integer *nrhs, 	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *    ldb, integer *info);
extern doublereal dasum_(integer *n, doublereal *dx, integer *incx);

extern /* Subroutine */ int dgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
				    iwork, integer *info);

extern /* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
				    lda, integer *ipiv, integer *info);
extern /* Subroutine */ int dgetrs_(char *trans, integer *n, integer *nrhs, 
				    doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
				    ldb, integer *info);

extern  int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
		    doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
		    ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
		    integer *info);

#define nsblas_idamax(n,dx,incx) idamax_((integer*)n, (doublereal*)dx, (integer*)incx)
#define nsblas_dcopy(n, dx, incx, dy, incy) dcopy_((integer*)n, (doublereal*)dx, (integer*)incx, (doublereal*)dy, (integer*)incy)
#define nsblas_ddot(n,dx,incx,dy,incy) ddot_((integer*)n,(doublereal*)dx,(integer*)incx,(doublereal*)dy,(integer*)incy)
#define nsblas_dger(m, n, alpha, x, incx, y, incy, a, lda) dger_((integer*)m, (integer*)n, (doublereal*)alpha, (doublereal*)x, (integer*)incx, (doublereal*)y, (integer*)incy, (doublereal*)a, (integer*)lda)
#define nslapack_dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) dgesv_((integer*)n, (integer*)nrhs, (doublereal*)a, (integer   *)lda, (integer*)ipiv, (doublereal*)b, (integer*)ldb, (integer*)info)
#define nsblas_dgemm(transa, transb, m,n, k,alpha, a, lda,  b, ldb, beta, c__,  ldc) dgemm_((char*)transa,(char*)transb, (integer*)m, (integer*)  n, (integer*)k, (doublereal*)alpha, (doublereal*)a, (integer*)lda,  (doublereal*)b, (integer*)ldb, (doublereal*)beta, (doublereal*)c__,  (integer*)ldc)
#define nsblas_dscal(n, da, dx,  incx) dscal_((integer*)n, (doublereal*)da, (doublereal*)dx,  (integer*)incx)
#define nsblas_dgemv(trans, m, n,  alpha, a, lda, x, incx, beta, y, incy) dgemv_((char*)trans, (integer*)m, (integer*)n, (doublereal*) alpha, (doublereal*)a, (integer*)lda, (doublereal*)x, (integer*)incx, (doublereal*)beta, (doublereal*)y, (integer*)incy)
#define nsblas_dnrm2(n, x, incx) dnrm2_((integer*)n, (doublereal*)x, (integer*)incx)
#define nsblas_daxpy(n, da, dx, incx, dy, incy) daxpy_((integer*)n, (doublereal*)da, (doublereal*)dx, (integer*)incx, (doublereal*)dy, (integer*)incy)
#define nslapack_dgetrf(m, n, a, lda, ipiv, info) dgetrf_((integer*)m, (integer*)n, (doublereal*)a, (integer*)lda, (integer*)ipiv, (integer*)info)
#define nslapack_dgetrs(trans, n, nrhs,a,lda,ipiv,b,ldb, info) dgetrs_((char*)trans, (integer*)n, (integer*)nrhs, 	(doublereal*)a, (integer*)lda, (integer*)ipiv, (doublereal*)b, (integer*)ldb, (integer*)info)
#define nslapack_dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info) dgesvd_(jobu,jobvt,(integer*)m,(integer*)n,(doublereal*)a,(integer*)lda,(doublereal*)s,(doublereal*)u,(integer*)ldu,(doublereal*)vt,(integer*)ldvt,(doublereal*)work,(integer *)lwork,(integer*)info)

#define nsblas_dasum(n, dx, incx) dasum_((integer*)n, (doublereal*)dx, (integer*)incx)
#define nslapack_dgesvx(fact, trans,n,nrhs, a,lda, af, ldaf,ipiv, equed, r__, c__,b, ldb, x,ldx,rcond, ferr, berr, work, iwork,info) dgesvx_((char*)fact, (char*)trans, (integer*)n, (integer*)	nrhs, a, (integer*)lda, af, (integer*)ldaf,(integer*)ipiv, equed, r__, c__,b, (integer*)ldb, x, (integer*)ldx, 	rcond, ferr, berr, work, (integer*) iwork, (integer*)info)
 

#define nslapack_dgelsy	dgelsy_
#define nslapack_dgeev 	dgeev_

#else

#define nsblas_idamax(n,dx,incx) idamax((I*)n, (R*)dx, (I*)incx)
#define nsblas_dcopy(n, dx, incx, dy, incy) dcopy((I*)n, (R*)dx, (I*)incx, (R*)dy, (I*)incy)
#define nsblas_ddot(n,dx,incx,dy,incy) ddot((I*)n,(R*)dx,(I*)incx,(R*)dy,(I*)incy)
#define nsblas_dger(m, n, alpha, x, incx, y, incy, a, lda) dger((I*)m, (I*)n, (R*)alpha, (R*)x, (I*)incx, (R*)y, (I*)incy, (R*)a, (I*)lda)
#define nslapack_dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) dgesv((I*)n, (I*)nrhs, (R*)a, (I   *)lda, (I*)ipiv, (R*)b, (I*)ldb, (I*)info)
#define nsblas_dgemm(transa, transb, m,n, k,alpha, a, lda,  b, ldb, beta, c__,  ldc) dgemm((char*)transa,(char*)transb, (I*)m, (I*)  n, (I*)k, (R*)alpha, (R*)a, (I*)lda,  (R*)b, (I*)ldb, (R*)beta, (R*)c__,  (I*)ldc)
#define nsblas_dscal(n, da, dx,  incx) dscal((I*)n, (R*)da, (R*)dx,  (I*)incx)
#define nsblas_dgemv(trans, m, n,  alpha, a, lda, x, incx, beta, y, incy) dgemv((char*)trans, (I*)m, (I*)n, (R*) alpha, (R*)a, (I*)lda, (R*)x, (I*)incx, (R*)beta, (R*)y, (I*)incy)
#define nsblas_dnrm2(n, x, incx) dnrm2((I*)n, (R*)x, (I*)incx)
#define nsblas_daxpy(n, da, dx, incx, dy, incy) daxpy((I*)n, (R*)da, (R*)dx, (I*)incx, (R*)dy, (I*)incy)
#define nslapack_dgetrf(m, n, a, lda, ipiv, info) dgetrf((I*)m, (I*)n, (R*)a, (I*)lda, (I*)ipiv, (I*)info)
#define nslapack_dgetrs(trans, n, nrhs,a,lda,ipiv,b,ldb, info) dgetrs((char*)trans, (I*)n, (I*)nrhs, 	(R*)a, (I*)lda, (I*)ipiv, (R*)b, (I*)ldb, (I*)info)
#define nslapack_dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info) dgesvd(jobu,jobvt,(I*)m,(I*)n,(R*)a,(I*)lda,(R*)s,(R*)u,(I*)ldu,(R*)vt,(I*)ldvt,(R*)work,(I *)lwork,(I*)info)
#define nsblas_dasum(n, dx, incx) dasum((I*)n, (R*)dx, (I*)incx)
#define nslapack_dgesvx(fact, trans,n,nrhs, a,lda, af, ldaf,ipiv, equed, r__, c__,b, ldb, x,ldx,rcond, ferr, berr, work, iwork,info) dgesvx((char*)fact, (char*)trans, (I*)n, (I*)	nrhs, a, (I*)lda, af, (I*)ldaf,(I*)ipiv, equed, r__, c__,b, (I*)ldb, x, (I*)ldx, 	rcond, ferr, berr, work, (I*) iwork, (I*)info)
 
#define nslapack_dgelsy	dgelsy
#define nslapack_dgeev 	dgeev

#endif

#if 0
extern const char __vmps_blas_transN[2];
extern const char __vmps_blas_transT[2];
extern const char __vmps_blas_transV[2];
extern const R __vmps_blas_regal1 ;
extern const R __vmps_blas_mregal1;
extern const R __vmps_blas_regal0 ;
extern const I  __vmps_blas_negal1 ;

#define nsblas_dgemv1(_n,_m,_A,_x,_y) 				nsblas_dgemv(__vmps_blas_transN,(_n),(_m),&__vmps_blas_regal1,(_A),_n,(_x),&__vmps_blas_negal1,&__vmps_blas_regal0,(_y),&__vmps_blas_negal1)
#define nsblas_dgemv2(_n,_m,_sa,_A,_x,_sy,_y) 			nsblas_dgemv(__vmps_blas_transN,(_n),(_m),(_sa),(_A),_n,(_x),&__vmps_blas_negal1,(_sy),(_y),&__vmps_blas_negal1)
#define nsblas_dgemtv1(_n,_m,_A,_x,_y) 				nsblas_dgemv(__vmps_blas_transT,(_n),(_m),&__vmps_blas_regal1,(_A),(_n),(_x),&__vmps_blas_negal1,&__vmps_blas_regal0,(_y),&__vmps_blas_negal1)
#define nsblas_dgemtv2(_n,_m,_sa,_A,_x,_sy,_y) 			nsblas_dgemv(__vmps_blas_transT,(_n),(_m),(_sa),(_A),(_n),(_x),&__vmps_blas_negal1,(_sy),(_y),&__vmps_blas_negal1)
#define nsblas_dgemtm1(_m,_n,_p,_a,_aoff,_b,_boff,_c,_coff) 	nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,(_m),(_n),(_p),&__vmps_blas_regal1,(_a),(_aoff),(_b),(_boff),&__vmps_blas_regal0,(_c),(_coff))
#define nsblas_dgemtmt1(_m,_n,_p,_a,_aoff,_b,_boff,_c,_coff) 	nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transT,(_m),(_n),(_p),&__vmps_blas_regal1,(_a),(_aoff),(_b),(_boff),&__vmps_blas_regal0,(_c),(_coff))
#define nsblas_dgemmt1(_m,_n,_p,_a,_aoff,_b,_boff,_c,_coff) 	nsblas_dgemm(__vmps_blas_transN,__vmps_blas_transT,(_m),(_n),(_p),&__vmps_blas_regal1,(_a),(_aoff),(_b),(_boff),&__vmps_blas_regal0,(_c),(_coff))
#define nsblas_dgemm1(_m,_n,_p,_a,_aoff,_b,_boff,_c,_coff) 	nsblas_dgemm(__vmps_blas_transN,__vmps_blas_transN,(_m),(_n),(_p),&__vmps_blas_regal1,(_a),(_aoff),(_b),(_boff),&__vmps_blas_regal0,(_c),(_coff))
#endif

#endif
