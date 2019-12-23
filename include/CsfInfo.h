#ifndef __header_CsfInfo_h__
#define __header_CsfInfo_h__

#include "nsOP_INFO.h"
#include "nsISO_INFO.h"
#include "SmoothedHeaviside.h"
#include "SmoothedDirac.h"
#include "eOperator.h"
#include "eSmoothingMethod.h"
#include "Err.h"

enum __ensCSF_iinfo{__ensCSF_iinfo_error=0,
		    __ensCSF_iinfo_kappa_filtering,
		    __ensCSF_iinfo_n };

enum __ensCSF_rinfo{__ensCSF_rinfo_error=0,
		    __ensCSF_rinfo_kappa_tol,
		    __ensCSF_rinfo_normal_tol,
		    __ensCSF_rinfo_distance_iso,
		    __ensCSF_rinfo_n };

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
  L			color;
  eSmoothingMethod 	smoothing_method;

  nsOP_INFO_ST		op_f;
  nsOP_INFO_ST		op_kappa;
  nsOP_INFO_ST		op_normal;
  nsOP_INFO_ST		op_tension;
  nsOP_INFO_ST		op_nabla;
  nsOP_INFO_ST		op_distance;

  nsISO_INFO_ST 	iso_info;
  
  SmoothedDirac		m_smoothedDirac;
  SmoothedHeaviside	m_smoothedHeaviside;

  I 			iinfo[__ensCSF_iinfo_n];
  R 			rinfo[__ensCSF_rinfo_n];  
} CsfInfo,* RESTRICT pCsfInfo;

typedef const CsfInfo *RESTRICT pCsfInfoReadOnly;

void CsfInfo_def	(pCsfInfo 			csf_,
			 const L			color_,
			 cst_eOperator			f_op_,
			 cst_eOperator			nabla_op_,
			 cst_eOperator			normal_op_,
			 const R 			normal_tol_,
			 cst_eOperator			kappa_op_,
			 const R 			kappa_tol_,
			 cst_eOperator			tension_op_,
			 cst_eOperator			distance_op_,
			 const R 			distance_iso_,
			 cst_ens_method_capturing	kind_curve_,
			 const I 			distance_level_,
			 cst_eSmoothingMethod 		smoothing_method_,
			 const I 			qh_degree_,
			 const eHeaviside		heaviside_,
			 const eDirac			dirac_,
			 cst_pR				eps_,
			 pErr 				err_);

#ifdef __cplusplus
}
#endif

#endif
