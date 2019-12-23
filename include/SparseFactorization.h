#ifndef __header_SparseFactorization_h__
#define __header_SparseFactorization_h__

#include "Sparse.h"
#include "eSparseFactorizationMethod.h"
#include "ExternUmfpack.h"
#include "ExternSparsekit.h"
#include "ExternPardiso.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
  eSparseFactorizationMethod 	factorization_method;    
  pSparse			matrix;
  pExternPardiso		extern_pardiso;
  pExternUmfpack 		extern_umfpack;
  pExternSparsekit		extern_sparsekit;
} SparseFactorization,*RESTRICT pSparseFactorization;


typedef const SparseFactorization * cst_pSparseFactorization;

/**
   \brief Kill the SparseFactorization object
   @param self_ pointer to the sparse linear solver algorithm
   @return NULL
   @note if factorization_ is NULL, nothing is done
*/
pSparseFactorization 	SparseFactorization_kill	(pSparseFactorization		self_);

/**
   \brief Create a SparseFactorization object
   @param method_  factorization method
   @param  err_  error
   @return new SparseFactorization object
*/
pSparseFactorization 	SparseFactorization_new		(const eSparseFactorizationMethod	method_,
							 pErr					err_);

/**
   \brief Compute the factorization 
   @param self_ pointer to the factorization algorithm
   @param S_ pointer to the sparse matrix 
*/
void  			SparseFactorization_precompute	(pSparseFactorization 		self_,
							 pSparse  			S_);

/**
   \brief Solve \f$A x = b\f$
   @param self_ pointer to the factorization algorithm
   @param tr_ transpose
   @param x_ solution
   @param y_ right hand side
   @note
   - solve (LU) y_ = x_
   - this routine HAVE TO BE compatible with lusol or lutsol routine of nsITERATIVE,
   that's why a void pointer is used
*/
void			SparseFactorization_compute	(void   * 			self_,
							 const char * 			tr_,
							 pR 				x_,
							 cst_pR 			y_);


#ifdef __cplusplus
}
#endif

  
#endif
