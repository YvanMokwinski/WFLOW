#ifndef __header_eSparsePreconditioningMethod_h__
#define __header_eSparsePreconditioningMethod_h__

#include "ConfigEnum.h"

ConfigEnum(eSparsePreconditioningMethod,
	   __eSparsePreconditioningMethod_NO,/*!< value for no preconditioning */
	   __eSparsePreconditioningMethod_LEFT,/*!< value for left preconditioning */
	   __eSparsePreconditioningMethod_RIGHT/*!< value for right preconditioning */);


#endif
