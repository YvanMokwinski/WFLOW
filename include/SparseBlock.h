#ifndef __header_SparseBlock_h__
#define __header_SparseBlock_h__

#include "Sparse.h"
#ifdef __cplusplus
extern "C"
{
#endif

/**\brief struct to define a sub block of a sparse matrix */
typedef struct
{
  /**\brief pointer to the parent matrix */
  pSparse	self_parent;
  /**\brief number of rows */
  I 		n_i;
  /**\brief number of columns */
  I 		n_j;
  /**\brief starting row */
  I 		start_j;
  /**\brief starting column */
  I 		start_i;
  /**\brief begin */
  pI   		b;
  /**\brief number of coeff in sequence */
  pI   		bn;
} SparseBlock,*RESTRICT pSparseBlock;

typedef const SparseBlock * RESTRICT cst_pSparseBlock;
 
#ifndef NDEBUG

cst_pSparse 		cst_SparseBlock_get_parent	(cst_pSparseBlock 	const B_);
pSparse 		SparseBlock_get_parent		(pSparseBlock		const B_);
I 			SparseBlock_get_ni		(cst_pSparseBlock 	const B_);
I 			SparseBlock_get_nj		(cst_pSparseBlock 	const B_);
I 			SparseBlock_get_starti		(cst_pSparseBlock 	const B_);
I 			SparseBlock_get_startj		(cst_pSparseBlock 	const B_);
pI 			SparseBlock_get_bn		(pSparseBlock		const B_);
pI 			SparseBlock_get_b		(pSparseBlock		const B_);
cst_pI 			cst_SparseBlock_get_bn		(cst_pSparseBlock 	const B_);
cst_pI 			cst_SparseBlock_get_b		(cst_pSparseBlock 	const B_);

#else

#define 		cst_SparseBlock_get_parent(_self) 	(_self)->self_parent
#define 		SparseBlock_get_parent(_self)		(_self)->self_parent
#define 		SparseBlock_get_ni(_self) 		(_self)->n_i
#define 		SparseBlock_get_nj(_self) 		(_self)->n_j
#define 		SparseBlock_get_starti(_self) 		(_self)->start_i
#define 		SparseBlock_get_startj(_self) 		(_self)->start_i
#define 		cst_SparseBlock_get_bn(_self) 		(_self)->bn
#define 		cst_SparseBlock_get_b(_self) 		(_self)->b
#define 		SparseBlock_get_bn(_self) 		(_self)->bn
#define 		SparseBlock_get_b(_self) 		(_self)->b

#endif

void 			SparseBlock_gemv_sub	(cst_pR 			a_,
						 cst_pSparseBlock const 	block_,
						 cst_pR 			rhs_,
						 cst_pR 			b_,
						 pR 				sol_);

void  			SparseBlock_ass_blank	(pSparseBlock const		B_,
						 cst_pI 			ni_,
						 cst_pI 			ddli_,
						 cst_pI 			nj_,
						 cst_pI 			ddlj_,
						 cst_pR 			locmat_,
						 cst_pI 			locmatoff_,
						 pI				blank_);

void 			SparseBlock_clr		(pSparseBlock const		B_);
void 			SparseBlock_free	(pSparseBlock const		B_);
void 			SparseBlock_clear	(pSparseBlock const		B_);

void 			SparseBlock_def		(pSparseBlock const		B_,
						 pSparse const			A_,
						 cst_pI 			start_i_,
						 cst_pI 			ni_,
						 cst_pI 			start_j_,
						 cst_pI 			nj_);

void 			SparseBlock_gemv	(cst_pR 			a_,
						 cst_pSparseBlock const 	B_,
						 cst_pR 			global_rhs_,
						 cst_pR 			b_,
						 pR 				global_sol_);
#ifdef __cplusplus
}
#endif

#endif
