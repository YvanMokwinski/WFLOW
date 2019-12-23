#ifndef __header_SparseStokes_h__
#define __header_SparseStokes_h__

#include "SparseBlock.h"
#ifdef __cplusplus
extern "C"
{
#endif

  //
  // UU UB US
  // BU BB
  // SU    SS
  //
typedef struct
{
  /** \brief global matrix */
  pSparse 		A; 

  I 			nddlu_dirichlet;
  I 			nddlv_dirichlet;
  I 			nddlw_dirichlet;
  I 			nddlp_dirichlet;
  I 			nddls_dirichlet;
  I 			dec_ddlu;
  I 			dec_ddlv;
  I 			dec_ddlw;
  I 			dec_ddlp;
  I 			dec_ddls;

  /** \brief flag for linear  */
  I 			linear[16];
  /** \brief pointers matrix to blocks */
  pSparseBlock		block[16];
  /** \brief blocks */
  SparseBlock		block_BU;
  SparseBlock		block_UB;
  SparseBlock		block_BB;
  SparseBlock		block_UU;
  SparseBlock		block_SS;
  SparseBlock		block_SU;
  SparseBlock		block_US;

  pI			ddlu_dirichlet;
  pI			ddlv_dirichlet;
  pI			ddlw_dirichlet;
  pI			ddlp_dirichlet;
  pI			ddls_dirichlet;
  pI			slip_iperm;
} SparseStokes;

typedef SparseStokes*RESTRICT pSparseStokes;
typedef const SparseStokes*RESTRICT cst_pSparseStokes;

  void 	SparseStokes_residu		(pSparseStokes		sparse_stokes_,
					 pR 				global_rhs_,
					 cst_pR 			x_,
					 cst_pR 			xi_,
					 cst_pR 			xii_,
					 void * 			usrptr_);
  
  void 	SparseStokes_free		(pSparseStokes const		S_);
  
  void 	SparseStokes_clr		(pSparseStokes const		S_);
  
  void  	SparseStokes_dirichlet_rhs	(cst_pSparseStokes const	S_,
						 pR 				global_rhs_);
  
  void  	SparseStokes_pressure_rhs	(cst_pSparseStokes const	S_,
						 cst_pR 			x_,
						 pR 				global_rhs_);
  
  void  SparseStokes_slip_rhs		(cst_pSparseStokes const	S_,
					 cst_pR 			x_,
					 pR 				global_rhs_);
  
  
  void SparseStokes_dirichlet_init(pSparseStokes 	const 	S_,
				   const I 	nedge_boundary_,
				   const I 	nshape_face_u_,
				   const I 	nshape_face_p_,
				   const I 	dec_ddlu_,
				   const I 	dec_ddlv_,
				   const I 	dec_ddlp_,
				   const I 	dec_ddls_);
  
  cst_pSparse		cst_SparseStokes_get_A		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_UU		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_UB		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_UF		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_FF		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_FU		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_BB		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_BU		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_SS		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_SU		(cst_pSparseStokes 	const S_);
  cst_pSparseBlock 	cst_SparseStokes_get_US		(cst_pSparseStokes 	const S_);

pSparse			SparseStokes_get_A		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_UU		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_UB		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_UF		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_FF		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_FU		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_BB		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_BU		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_SS		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_SU		(pSparseStokes 		const S_);
pSparseBlock 		SparseStokes_get_US		(pSparseStokes 		const S_);

#if 0
#ifndef NDEBUG

#else

#define cst_SparseStokes_get_A(_S)  (_S)->A
#define cst_SparseStokes_get_UU(_S) &(_S)->block_UU
#define cst_SparseStokes_get_UB(_S) &(_S)->block_UB
#define cst_SparseStokes_get_UF(_S) &(_S)->block_UF
#define cst_SparseStokes_get_FU(_S) &(_S)->block_FU
#define cst_SparseStokes_get_FF(_S) &(_S)->block_FF
#define cst_SparseStokes_get_BU(_S) &(_S)->block_BU
#define cst_SparseStokes_get_BB(_S) &(_S)->block_BB
#define cst_SparseStokes_get_SS(_S) &(_S)->block_SS
#define cst_SparseStokes_get_SU(_S) &(_S)->block_SU
#define cst_SparseStokes_get_US(_S) &(_S)->block_US

#define SparseStokes_get_A(_S)  (_S)->A
#define SparseStokes_get_UU(_S) &(_S)->block_UU
#define SparseStokes_get_UB(_S) &(_S)->block_UB
#define SparseStokes_get_UF(_S) &(_S)->block_UF
#define SparseStokes_get_FU(_S) &(_S)->block_FU
#define SparseStokes_get_FF(_S) &(_S)->block_FF
#define SparseStokes_get_BU(_S) &(_S)->block_BU
#define SparseStokes_get_BB(_S) &(_S)->block_BB
#define SparseStokes_get_SS(_S) &(_S)->block_SS
#define SparseStokes_get_SU(_S) &(_S)->block_SU
#define SparseStokes_get_US(_S) &(_S)->block_US

#endif
#endif

#ifdef __cplusplus
}
#endif


#endif
