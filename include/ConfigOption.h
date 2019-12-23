#ifndef __header_ConfigOption_h__
#define __header_ConfigOption_h__

#include "Config.h"
#include "Type.h"
#include <string.h>


#define ConfigOption_cDeclaration(_type,_enum)				\
  									\
  typedef struct {							\
    _enum flag;								\
    _type default_value;						\
    cst_pS description;							\
    cst_pS name_cmdline;						\
    cst_pS name_configfile;} _enum##Options;				\
  									\
  cst_pS 	_enum##_get_description		(const _enum self_);	\
  cst_pS 	_enum##_get_name_cmdline	(const _enum self_);	\
  _type  	_enum##_get_default		(const _enum self_);	\


#define ConfigOption_cImplementation(_type,_enum)			\
  									\
  cst_pS _enum##Option_get (const _enum i_)				\
  {									\
    return __tableoptions_##_enum[i_].name_cmdline;			\
  }									\
  									\
  cst_pS  _enum##Option_get_description (const _enum i_)		\
  {									\
    return __tableoptions_##_enum[i_].description;			\
  }									\
  									\
  _type  _enum##Option_get_default (const _enum i_)			\
  {									\
    return __tableoptions_##_enum[i_].default_value;			\
  }									\
  									\


#if 0
  _enum	Cmdline_get_##_t(pCmdline const self_,const char * name_)  
  _enum Cmdline_get_##_enum 	(pCmdline const cmdline_,const char * name_) \
  {									\
    Str ctmp;								\
    if (Cmdline_get_name(cmdline_,name_,ctmp))				\
      {									\
	return _enum##_get_enum(ctmp);					\
      }									\
    else								\
      {									\
	return __##_enum##_ERROR;					\
      }									\
  }
#endif

#ifndef __cplusplus

#define ConfigOption_DECLARATION(_type,_enum) 		ConfigOption_cDeclaration(_type,_enum)
#define ConfigOption_IMPLEMENTATION(_type,_enum) 	ConfigOption_cImplementation(_type,_enum)

#else

#include <ostream>

#define ConfigOption_DECLARATION(_type,_enum)				\
  extern "C"  {	ConfigOption_cDeclaration(_type,_enum); }		\
  
#define ConfigOption_IMPLEMENTATION(_type,_enum) extern "C" { ConfigOption_cImplementation(_type,_enum); }

#endif

#define ConfigOption(_type,_enum) 	     ConfigOption_DECLARATION(_type,_enum)

#define ConfigOptionTable(_type,_enum,...) static const _enum##Options  __tableoptions_##_enum[__##_enum##_ALL] =   { { __##_enum##_ERROR,((_type)0),NULL,NULL}, __VA_ARGS__ };\
  ConfigOption_IMPLEMENTATION(_type,_enum)

/**/
#endif
