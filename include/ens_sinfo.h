#ifndef __header__ens_sinfo_h__
#define __header__ens_sinfo_h__


#include "Type.h"
#include "ConfigEnum.h"
// #include "Cmdline.h"



ConfigEnum(ens_sinfo,
	   __ens_sinfo_name,
	   __ens_sinfo_pblmname,
	   __ens_sinfo_ofilename);


ConfigOption(cst_pS,ens_sinfo);
#if 0
void ens_sinfo_from_Cmdline(STR		sinfo_[__ens_sinfo_ALL],
			    const L	have_configfile_,
			    pCmdline  const 	cmdline_);
#endif
#endif
