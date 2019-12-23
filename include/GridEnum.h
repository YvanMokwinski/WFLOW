#ifndef __header_GridEnum_h__
#define __header_GridEnum_h__

#include "Type.h"
#include <string.h>
#define GridEnum_cDeclaration(_token)				\
  /**\brief Get enum @param str_ string @return enum */		\
  _token _token##_get_enum	(cst_pS 	str_);		\
  /**\brief Get string @param self_ self enum @return string */	\
  cst_pS _token##_get_string	(const _token 	self_)


#define GridEnum_cImplementation(_token)				\
  _token _token##_get_enum	(cst_pS	str_)				\
  {									\
    _token method = __##_token##_ERROR;					\
    if (str_)								\
      {									\
	_token imethod = __##_token##_ERROR;				\
	for (++imethod;imethod<__##_token##_ALL;++imethod)		\
	  if (!strcmp(str_,__tablenames_##_token[imethod]))		\
	    break;							\
	method = (imethod<__##_token##_ALL) ? imethod : __##_token##_ERROR; \
      }									\
    if (!method)							\
      {									\
	fprintf(stderr,"'%s' is not a "#_token" name\n",(str_)?str_:"(NULL)"); \
	fprintf(stderr,"\tlist of available method(s) :\n");		\
	{ _token imethod = __##_token##_ERROR;				\
	  for (++imethod;imethod<__##_token##_ALL;++imethod)		\
	    {								\
	      fprintf(stderr,"\t\t'%s'\n",__tablenames_##_token[imethod]); \
	    } }								\
      }									\
    return method;							\
  }									\
  									\
  cst_pS   _token##_get_string	(const _token 	self_)			\
  {									\
    return __tablenames_##_token[self_];				\
  }									\





#ifndef __cplusplus

#define GridEnum_DECLARATION(_token) 	  GridEnum_cDeclaration(_token)
#define GridEnum_IMPLEMENTATION(_token) GridEnum_cImplementation(_token)

#else

#include <ostream>



#define GridEnum_DECLARATION(_token)					\
  extern "C"  {	GridEnum_cDeclaration(_token); }			\
  /**\brief Overload of increment operator for C++*/			\
  inline _token & operator ++(_token & self_) { self_=static_cast<_token>(self_+1); return self_; }; \
  /**\brief Overload of output stream operator for C++*/		\
  inline std::ostream& operator<<(std::ostream&s_,const _token d_)

#define GridEnum_IMPLEMENTATION(_token) extern "C" { GridEnum_cImplementation(_token); }


#endif

/* __attribute__((aligned(256))) */

#define GridEnum(_a,...) enum __##_a {__##_a##_ERROR=0, __VA_ARGS__ ,__##_a##_ALL } ; typedef enum __##_a _a; typedef const _a cst_##_a; GridEnum_DECLARATION(_a)
#define GridEnumStrings(_a,...) static cst_pS __tablenames_##_a[__##_a##_ALL] =   { #_a"_ERROR",__VA_ARGS__ };GridEnum_IMPLEMENTATION(_a)

#endif
