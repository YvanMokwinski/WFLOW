#ifndef __header_eSparseIterativeMethod_h__
#define __header_eSparseIterativeMethod_h__

#include "ConfigEnum.h"

ConfigEnum(eSparseIterativeMethod,
	   __eSparseIterativeMethod_GMRES,/*!< value for gmres method */
	   __eSparseIterativeMethod_FGMRES,/*!< value for fgmres method */
	   __eSparseIterativeMethod_TFQMR,/*!< value for tfqmr method */
	   __eSparseIterativeMethod_BCGSTAB,/*!< value for bcgstab method */
	   __eSparseIterativeMethod_DBCG,/*!< value for double conjugate gradient method */
	   __eSparseIterativeMethod_CG,/*!< value for conjugate gradient method */
	   __eSparseIterativeMethod_CGNR,/*!< value for conjugate gradient NR method  */
	   __eSparseIterativeMethod_BCG,/*!< value for bi-conjugate gradient method */
	   __eSparseIterativeMethod_FOM,/*!< value for full orthogonalization method */
	   __eSparseIterativeMethod_DQGMRES/*!< value for dqgmres */);


ConfigOption(cst_pS,eSparseIterativeMethod);


#endif
