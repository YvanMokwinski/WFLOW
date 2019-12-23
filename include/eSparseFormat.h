#ifndef __header_eSparseFormat_h__
#define __header_eSparseFormat_h__

#include "ConfigEnum.h"

ConfigEnum(eSparseFormat,
	   __eSparseFormat_CSR,/*!< value for Compressed Sparse Row format */
	   __eSparseFormat_CSC,/*!< value for Compressed Sparse Column format */
	   __eSparseFormat_COO/*!< value for Coordinates format */);

#endif
