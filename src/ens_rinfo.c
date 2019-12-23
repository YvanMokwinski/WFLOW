#include "ns_sys.h"
#include "Cmdline.h"
#include "ens_rinfo.h"

ConfigEnumStrings(ens_rinfo,
		  "hmin",
		  "hmax",
		  "ratio_viscosity",
		  "ratio_density",
		  "dtmin",
		  "dtmax",
		  "epsmin",
		  "epsmax",
		  "newton_tol_residu",
		  "newton_tol_correction",
		  "reynold",
		  "weber",
		  "froude");

ConfigOptionTable(R,
		  ens_rinfo,
		  {__ens_rinfo_hmin,((R)1.0e-4),"hmin","--hmin","HMIN"},
		  {__ens_rinfo_hmax,((R)0.25),"hmax","--hmax","HMAX"},
		  {__ens_rinfo_ratio_viscosity,((R)1.0),"ratio viscosity","--nu","NU"},
		  {__ens_rinfo_ratio_density,((R)1.0),"ratio density","--rho","RHO"},
		  {__ens_rinfo_dtmin,((R)1.0e-4),"dtmin","--dt","DT"},
		  {__ens_rinfo_dtmax,((R)0.5),"dtmax","--dtmax","DTMAX"},
		  {__ens_rinfo_epsmin,((R)0.25),"epsmin","--epsmin","EPSMIN"},
		  {__ens_rinfo_epsmax,((R)0.9),"epsmax","--epsmax","EPSMAX"},
		  {__ens_rinfo_newton_tol_residu,((R)1.0e-8),"stopping criteria based on residu","--ntolr","NEWTON_TOLR"},
		  {__ens_rinfo_newton_tol_correc,((R)1.0e-5),"stopping criteria based on correction","--ntolc","NEWTON_TOLC"},
		  {__ens_rinfo_reynold,((R)0.5),"reynold number","--reynold","REYNOLD"},
		  {__ens_rinfo_weber,((R)0.5),"weber number","--weber","WEBER"},
		  {__ens_rinfo_froude,((R)0.5),"froude number","--froude","FROUDE"});


void ens_rinfo_from_Cmdline(R 			rinfo_[__ens_rinfo_ALL],
			    const L  	have_configfile_,
			    pCmdline const 	cmdline_)

{
  { ens_rinfo i = __ens_rinfo_ERROR;
    for (++i;i<__ens_rinfo_ALL;++i)
      {
	if (__tableoptions_ens_rinfo[i].flag!=i)
	  {
	    fprintf(stderr,"ens_rinfo_argv:wrong dat nbase\n");
	  }
	const L a = Cmdline_get_real(cmdline_,__tableoptions_ens_rinfo[i].name_cmdline,&rinfo_[i]);
	if ( (NOT a) AND (NOT have_configfile_) )
	  {
	    rinfo_[i] = __tableoptions_ens_rinfo[i].default_value;
	  }
      } }
}
