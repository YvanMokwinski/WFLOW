#ifndef __header_ExternUmfpack_h__
#define __header_ExternUmfpack_h__

#include "Sparse.h"

typedef struct
{
} ExternUmfpack,*RESTRICT pExternUmfpack;

pExternUmfpack 	ExternUmfpack_new		();
pExternUmfpack 	ExternUmfpack_kill		(pExternUmfpack 	const self_);
void 		ExternUmfpack_precompute	(pExternUmfpack 	const self_,
						 pSparse			S_);
void  		ExternUmfpack_compute		(void   * 			self_,
						 const char * 			tr_,
						 pR 				sol_,
						 cst_pR 			rhs_);

#endif
