#ifndef __header_ExternPardiso_h__
#define __header_ExternPardiso_h__

#include "Sparse.h"
#ifdef __cplusplus
extern "C"
{
#endif
typedef struct
{
#if __MNS_WITH_PARDISO__  
  pSparse			matrix;
  I				N;
  void * 			pardiso_pt[64];
  I 				pardiso_iparm[64];
  I      			pardiso_maxfct;
  I 				pardiso_mnum;
  I  				pardiso_phase;
  I  				pardiso_error;
  I  				pardiso_msglvl;
  double 			pardiso_ddum[64];
  I  				pardiso_idum[64];
  I  				pardiso_mtype;
  I				pardiso_nrhs;
#endif
} ExternPardiso,*RESTRICT pExternPardiso;

pExternPardiso 	ExternPardiso_kill		(pExternPardiso 	const 	self_);
pExternPardiso 	ExternPardiso_new		();

void 		ExternPardiso_precompute	(pExternPardiso 	const 	self_,
						 pSparse			S_);

void  		ExternPardiso_compute		(void   * 			self_,
						 const char * 			tr_,
						 pR 				sol_,
						 cst_pR 			rhs_);


#ifdef __cplusplus
}
#endif

#endif
