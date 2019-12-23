#ifndef __HEADER_ePrecision_H__
#define __HEADER_ePrecision_H__

#include "ConfigEnum.h"

ConfigEnum(ePrecision,
	   __ePrecision_SINGLE,
	   __ePrecision_DOUBLE,
	   __ePrecision_QUAD);

/**
   \brief Define a macro to get sizeof of defined real value R
*/
#define PrecisionSizeof(_self,_n)			\
  switch(_self)						\
    {							\
    case __ePrecision_SINGLE:				\
      {							\
	(_n) = sizeof(float);				\
	break;						\
      }							\
    case __ePrecision_DOUBLE:				\
      {							\
	(_n) = sizeof(double);				\
	break;						\
      }							\
    case __ePrecision_QUAD:				\
      {							\
	(_n) = sizeof(long double);			\
	break;						\
      }							\
    case __ePrecision_ALL:				\
    case __ePrecision_ERROR:				\
      {							\
	(_n) = (size_t)0;				\
	break;						\
      }							\
    }


#endif




