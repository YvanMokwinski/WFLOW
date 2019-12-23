#ifndef __header_ExternSparsekit_h__
#define __header_ExternSparsekit_h__

#include "Sparse.h"

typedef struct
{

  /**
     \brief SPARSEKIT incomplete factorization
   */
  pSparse 		preconditioner;

  /**
     \brief SPARSEKIT size of the permutation array
  */
  I 			sparsekit_iperm_n;

  /**
     \brief SPARSEKIT permutation array
  */
  I * 			sparsekit_iperm;

  /**
     \brief SPARSEKIT error
  */
  I   			sparsekit_ierr;

  /**
     \brief SPARSEKIT size of the integer working array
  */
  I 			sparsekit_iwork_n;

  /**
     \brief SPARSEKIT integer working array
  */
  I * 			sparsekit_iwork;

  /**
     \brief SPARSEKIT size of the real working array
  */
  I 			sparsekit_rwork_n;

  /**
     \brief SPARSEKIT real working array
  */
  R * 			sparsekit_rwork;

  /**
     \brief SPARSEKIT tolerance for incomplete factorization
  */
  R 			sparsekit_tol;

  /**
     \brief SPARSEKIT alph parameter for incomplete factorization
  */
  R 			sparsekit_alph;

  /**
     \brief SPARSEKIT level of fil for incomplete factorization
  */
  I 			sparsekit_lfil;

  /**
     \brief SPARSEKIT permutation tolerance for incomplete factorization with permutation
  */
  R 			sparsekit_permtol;

  /**
     \brief SPARSEKIT levels
  */
  I  * 			sparsekit_levs;

  /**
     \brief SPARSEKIT float work
  */
  R * 			sparsekit_fpar;

} ExternSparsekit,*RESTRICT pExternSparsekit;

pExternSparsekit 	ExternSparsekit_new		();
pExternSparsekit 	ExternSparsekit_kill		(pExternSparsekit 	const self_);

void 			ExternSparsekit_precompute	(pExternSparsekit 	const self_,
							 pSparse		S_);

void  			ExternSparsekit_compute		(void   * 		self_,
							 const char * 		tr_,
							 pR 			sol_,
							 cst_pR 		rhs_);


#endif
