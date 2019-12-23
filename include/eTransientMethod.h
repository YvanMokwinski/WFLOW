#ifndef __header_eTransientMethod_h__
#define __header_eTransientMethod_h__

#include "ConfigEnum.h"

ConfigEnum(eTransientMethod,
	   __eTransientMethod_NO,/*!< value for gmres method */
	   __eTransientMethod_EULER,/*!< value for fgmres method */
	   __eTransientMethod_IMR,/*!< value for tfqmr method */
	   __eTransientMethod_GEAREULER,/*!< value for bcgstab method */
	   __eTransientMethod_GEARIMR,/*!< value for double conjugate gradient method */
	   __eTransientMethod_TRAPEZE /*!< value for conjugate gradient method */ );


ConfigOption(L,eTransientMethod);


#endif
