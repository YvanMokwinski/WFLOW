#ifndef __header_eSparseFactorizationMethod_h__
#define __header_eSparseFactorizationMethod_h__

#include "ConfigEnum.h"

ConfigEnum(eSparseFactorizationMethod,
	   __eSparseFactorizationMethod_UMFPACKLU,/*!< lu with UMFPACK */
	   __eSparseFactorizationMethod_PARDISOLU,/*!< lu with PARDISO */
	   __eSparseFactorizationMethod_CHOL,/*!< cholesky  */
	   __eSparseFactorizationMethod_ILU0,/*!< ilu0 with SPARSEKIT2 */
	   __eSparseFactorizationMethod_MILU0,/*!< milu0 with SPARSEKIT2 */
	   __eSparseFactorizationMethod_ILUT,/*!< ilut with SPARSEKIT2  */
	   __eSparseFactorizationMethod_ILUTP,/*!< ilutp with SPARSEKIT2 */
	   __eSparseFactorizationMethod_ILUK,/*!< iluk with SPARSEKIT2  */
	   __eSparseFactorizationMethod_ILUD/*!<ilud with SPARSEKIT2  */);

#endif



