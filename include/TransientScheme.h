#ifndef __header_TransientScheme_h__
#define __header_TransientScheme_h__

#include "eTransientMethod.h"
#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
  eTransientMethod 	method;
  I 			nstep;
  R 			scheme[8];
} TransientScheme,* RESTRICT pTransientScheme;

Err 	TransientScheme_def		(pTransientScheme 	const	self_,
					 cst_eTransientMethod 		method_);

Err 	TransientScheme_precompute	(pTransientScheme 	const	scheme_,
					 const I 			itime_,
					 const R 			dt_,
					 const R 			idt_,
					 const R 			dti_,
					 const R 			ratio_dti_,
					 const R 			iratio_dti_);

void 	TransientScheme_compute		(pTransientScheme	const	self_,
					 cst_pI			const	fn_,
					 cst_pR			const	f_,
					 cst_pI			const	foff_,
					 pR			const	r_);

#ifdef __cplusplus
}
#endif

#endif
