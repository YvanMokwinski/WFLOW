#ifndef __header__ens_rinfo_h__
#define __header__ens_rinfo_h__

#include "ConfigEnum.h"
// #include "Cmdline.h"


ConfigEnum(ens_rinfo,
	   __ens_rinfo_hmin,
	   __ens_rinfo_hmax,
	   __ens_rinfo_ratio_viscosity,
	   __ens_rinfo_ratio_density,
	   __ens_rinfo_dtmin,
	   __ens_rinfo_dtmax,
	   __ens_rinfo_epsmin,
	   __ens_rinfo_epsmax,
	   __ens_rinfo_newton_tol_residu,
	   __ens_rinfo_newton_tol_correc,
	   __ens_rinfo_reynold,
	   __ens_rinfo_weber,
	   __ens_rinfo_froude);


ConfigOption(R,ens_rinfo);
#if 0
void ens_rinfo_from_Cmdline(R 			rinfo_[__ens_rinfo_ALL],
			    const L	have_configfile_,
			    pCmdline  const 	cmdline_);
#endif


#endif
