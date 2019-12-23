#ifndef __header_nsISO_INFO_h__
#define __header_nsISO_INFO_h__

#include "ens_method_capturing.h"
typedef struct
{
  ens_method_capturing 	method;
  I 			lev;
  R 			iso;
  R 			tol;
} nsISO_INFO_ST,*nsISO_INFO;

typedef const nsISO_INFO_ST cst_nsISO_INFO_ST;
typedef cst_nsISO_INFO_ST *cst_nsISO_INFO;

void nsISO_INFO_def		(nsISO_INFO 			iso_info_,
				 cst_ens_method_capturing 	method_,		    
				 const R 			iso_,
				 const R 			tol_,
				 const I 			lev_);



#endif
