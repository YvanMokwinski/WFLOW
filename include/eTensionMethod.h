#ifndef __HEADER_eTensionMethod_H__
#define __HEADER_eTensionMethod_H__

#include "ConfigEnum.h"

ConfigEnum(eTensionMethod,
	   __eTensionMethod_NO, /**< value for no surface force */
	   __eTensionMethod_CSF,/**< value for CSF model */
	   __eTensionMethod_SSF	/**< value for SSF model  */);

ConfigOption(L,eTensionMethod);


#endif




