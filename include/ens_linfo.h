#ifndef __header__ens_linfo_h__
#define __header__ens_linfo_h__
#include "ConfigEnum.h"
//#include "Cmdline.h"


ConfigEnum(ens_linfo,
	   __ens_linfo_verbose,
	   __ens_linfo_mesh_adaptivity,
	   __ens_linfo_time_adaptivity,
	   __ens_linfo_transport,
	   __ens_linfo_transport_uncoupled,
	   __ens_linfo_transport_galerkin,
	   __ens_linfo_pressure_uncoupled,
	   __ens_linfo_pressure_freematrix,
	   __ens_linfo_slip,
	   __ens_linfo_color,
	   __ens_linfo_tension,
	   __ens_linfo_redistance,
	   __ens_linfo_pspg,
	   __ens_linfo_usupg,
	   __ens_linfo_axisymetric_x,
	   __ens_linfo_axisymetric_y,
	   __ens_linfo_dynamic_boundary_condition,
	   __ens_linfo_vnsdg,
	   __ens_linfo_skip_verif);


ConfigOption(L,ens_linfo);

#if 0
void ens_linfo_from_Cmdline(L 		linfo_[__ens_linfo_ALL],
			    const L  	have_configfile_,
			    pCmdline const 	cmdline_);
#endif

#endif








