#ifndef __HEADER_eFileSuffix_H__
#define __HEADER_eFileSuffix_H__

#include "ConfigEnum.h"

ConfigEnum(eFileSuffix,
	   __eFileSuffix_CmgData,
	   __eFileSuffix_MnsBinary,
	   __eFileSuffix_Latex,
	   __eFileSuffix_Matlab,
	   __eFileSuffix_MeditMesh,
	   __eFileSuffix_MeditSolution,
	   __eFileSuffix_Pirate,
	   __eFileSuffix_Text,
	   __eFileSuffix_VtkAscii,
	   __eFileSuffix_ZippedMnsBinary,
	   __eFileSuffix_ZippedMeditMesh,
	   __eFileSuffix_ZippedMeditSolution);

#ifdef __cplusplus
extern "C"
{
#endif
const char * 	eFile_desc	(cst_eFileSuffix 	format_);
eFileSuffix 	eFile_suffix	(cst_pS		const	filename_);

cst_eFileSuffix eFileSuffix_get_suffix			(cst_pS const 	filename_);
cst_eFileSuffix eFileSuffix_get_suffixAndBasename	(cst_pS const 	filename_,
							 STR 		basename_);
#ifdef __cplusplus
};
#endif

#endif




