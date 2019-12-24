#ifndef __header__ens_iinfo_h__
#define __header__ens_iinfo_h__

#include "ConfigEnum.h"
// #include "Cmdline.h"

ConfigEnum(ens_iinfo,
	   __ens_iinfo_newton_maxiter,
	   __ens_iinfo_ntime,
	   __ens_iinfo_noboundary_vcod,
	   __ens_iinfo_mesh_adaptivity_maxiter,
	   __ens_iinfo_ntime_interval,
	   __ens_iinfo_heaviside,
	   __ens_iinfo_dirac,
	   __ens_iinfo_restart,
	   __ens_iinfo_nproc);


ConfigOption(I,ens_iinfo);
#if 0
void ens_iinfo_from_Cmdline(I 			iinfo_[__ens_iinfo_ALL],
			    const L		have_configfile_,
			    pCmdline  const 	cmdline_);
#endif

#endif
