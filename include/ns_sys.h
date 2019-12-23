#ifndef __header_ns_sys_h__
#define __header_ns_sys_h__

#include "Config.h"
#include "Type.h"
#include <stdio.h>
#include <stdlib.h>
#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif

  Err 		ns_fclose		(FILE * 	fich_);

  FILE * 	ns_fopen		(cst_pS  	filename_,
					 cst_pS  	mode_);

  FILE * 	ns_fopen_wsecure	(cst_pS  	filename_);

  pS 		ns_string_basename	(cst_pS  	t_,
					 cst_pS  	suffix);

  void 		ns_string_basenam	(cst_pS  	t_,
					 pS 		o_);

  pS 		ns_string_get_basename	(cst_pS  	t_);

#ifdef __cplusplus
}
#endif

#endif
