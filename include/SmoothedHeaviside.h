#ifndef __header_SmoothedHeaviside_h__
#define __header_SmoothedHeaviside_h__

#include "eHeaviside.h"
#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif
enum __ensHEAVISIDE_iinfo 	{   __ensHEAVISIDE_iinfo_error=0,
				    __ensHEAVISIDE_iinfo_lev,
				    __ensHEAVISIDE_iinfo_n };

enum __ensHEAVISIDE_rinfo 	{   __ensHEAVISIDE_rinfo_error=0,
				    __ensHEAVISIDE_rinfo_eps,
				    __ensHEAVISIDE_rinfo_epsmin,
				    __ensHEAVISIDE_rinfo_epsmax,
				    __ensHEAVISIDE_rinfo_n };

cst_pS 				__ensHEAVISIDE_rinfo_string	(const enum __ensHEAVISIDE_rinfo flag_);
enum __ensHEAVISIDE_rinfo 	__ensHEAVISIDE_rinfo_enum	(cst_pS name);

cst_pS 				__ensHEAVISIDE_iinfo_string	(const enum __ensHEAVISIDE_iinfo flag_);
enum __ensHEAVISIDE_iinfo 	__ensHEAVISIDE_iinfo_enum	(cst_pS name);

typedef struct
{
  eHeaviside 	m_heaviside;
  cst_pR  	m_eps;
} SmoothedHeaviside,* RESTRICT pSmoothedHeaviside;

typedef const SmoothedHeaviside * RESTRICT pSmoothedHeavisideReadOnly;

void SmoothedHeaviside_def	(pSmoothedHeaviside 	const 		heaviside_info_,
				 cst_eHeaviside 			heaviside_,
				 cst_pR 		const 		eps_,
				 pErr			const 		err_);

void SmoothedHeaviside_eval	(pSmoothedHeavisideReadOnly	const	self_, 
				 const I 				rn_,
				 pR 				const	r_,
				 pErr					err_);
#ifdef __cplusplus
}
#endif

#endif
