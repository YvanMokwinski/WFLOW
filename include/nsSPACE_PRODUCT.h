#ifndef __header_nsSPACE_PRODUCT_h__
#define __header_nsSPACE_PRODUCT_h__

#include "nsSPACE.h"

typedef struct
{
  I 			nspace;
  I 			nddl;
  I 			nddlelm;
  pSpaceReadOnly	spaces[16];
  pI 			dep;
} nsSPACE_PRODUCT_ST,*RESTRICT nsSPACE_PRODUCT;

typedef const nsSPACE_PRODUCT_ST 	cst_nsSPACE_PRODUCT_ST;
typedef cst_nsSPACE_PRODUCT_ST*RESTRICT	cst_nsSPACE_PRODUCT;

I 		nsSPACE_PRODUCT_get_nddl	(cst_nsSPACE_PRODUCT const 	self_);
I 		nsSPACE_PRODUCT_get_nddlelm	(cst_nsSPACE_PRODUCT const 	self_);
I 		nsSPACE_PRODUCT_get_nspace	(cst_nsSPACE_PRODUCT const 	self_);
pSpaceReadOnly 	nsSPACE_PRODUCT_get_space	(cst_nsSPACE_PRODUCT const 	self_,
						 const I 			i_);
void 		nsSPACE_PRODUCT_def		(nsSPACE_PRODUCT 		self_,
						 cst_pI				dep_,
						 const I 			nspace_,
						 ...);
void  		nsSPACE_PRODUCT_free		(nsSPACE_PRODUCT 		self_);
pSparse		nsSPACE_PRODUCT_galerkin_def	(cst_nsSPACE_PRODUCT		self_);

#endif
