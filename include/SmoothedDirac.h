#ifndef __header_SmoothedDirac_h__
#define __header_SmoothedDirac_h__
#include "eDirac.h"
#include "Err.h"

enum __ensDIRAC_iinfo 	{   __ensDIRAC_iinfo_error=0,
			    __ensDIRAC_iinfo_lev,
			    __ensDIRAC_iinfo_n };

enum __ensDIRAC_rinfo 	{   __ensDIRAC_rinfo_error=0,
			    __ensDIRAC_rinfo_eps,
			    __ensDIRAC_rinfo_epsmin,
			    __ensDIRAC_rinfo_epsmax,
			    __ensDIRAC_rinfo_n };

cst_pS 			__ensDIRAC_rinfo_string	(const enum __ensDIRAC_rinfo flag_);
enum __ensDIRAC_rinfo 	__ensDIRAC_rinfo_enum	(cst_pS name);
cst_pS 			__ensDIRAC_iinfo_string	(const enum __ensDIRAC_iinfo flag_);
enum __ensDIRAC_iinfo 	__ensDIRAC_iinfo_enum	(cst_pS name);

typedef struct
{
  eDirac dirac;
  cst_pR  eps;
} SmoothedDirac,* RESTRICT pSmoothedDirac;

typedef const SmoothedDirac * RESTRICT pSmoothedDiracReadOnly;

void SmoothedDirac_def	(pSmoothedDirac 	const 	self_,
			 cst_eDirac 			dirac_,
			 cst_pR 		const 	eps_,
			 pErr				err_);

void SmoothedDirac_eval	(pSmoothedDiracReadOnly	const	self_, 
			 const I 			rn_,
			 pR 			const 	r_,
			 pErr				err_);

#endif
