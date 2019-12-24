#include "ns_sys.h"
#include "ens_iinfo.h"
//#include "Cmdline.h"
#include "eHeaviside.h"
#include "eDirac.h"



ConfigEnumStrings(ens_iinfo,
		  "newton_maxiter",
		  "ntime",
		  "noboundary_vcod",
		  "mesh_adaptivity_maxiter",
		  "ntime_interval",
		  "heaviside",
		  "dirac",
		  "restart",
		  "nproc");

ConfigOptionTable(I,
		  ens_iinfo,
		  {__ens_iinfo_newton_maxiter,((I)20),"newton max iter","--newton-niter","NEWTON_MAXITER"},
		  {__ens_iinfo_ntime,((I)7),"n time step","--ntime","NTIME"},
#if 0
		  {__ens_iinfo_transient_scheme_velocity,((I)__ens_scheme_transient_euler),"transient scheme velocity","--ut","TRANSIENT_SCHEME_VELOCITY"},
		  {__ens_iinfo_transient_scheme_transport,((I)__ens_scheme_transient_euler),"transient scheme transport","--ft","TRANSIENT_SCHEME_TRANSPORT"},
#endif
		  {__ens_iinfo_noboundary_vcod,((I)100),"no boundary cod","--noboundary-cod","NOBOUNDARY_COD"},
#if 0
		  {__ens_iinfo_linsys_solver,((I)__ens_linsys_method_pardiso),"choose linsys solver","--linsys","LINSYS"},
		  {__ens_iinfo_tension_method,((I)__ens_tension_method_csf),"choose method (0 ou 1)","--csf","CSF"},
#endif
		  {__ens_iinfo_mesh_adaptivity_maxiter,((I)5),"n time step","--adapt-niter","ADAPT_NITER"},
#if 0
		  {__ens_iinfo_mesh_adaptivity_method,((I)__ens_mesh_adaptivity_method_anisotropic),"n time step","--ns-mesh-adaptivity-method","MESH_ADAPTIVITY_METHOD"},
		  {__ens_iinfo_dg_transient_scheme,((I)__ens_scheme_transient_euler),"dg scheme","--ns-dg-transient-scheme","DG_TRANSIENT_SCHEME"},
#endif
		  {__ens_iinfo_ntime_interval,((I)5),"ntime interval","--ns-ntime-interval","NTIME_INTERVAL"},
		  {__ens_iinfo_heaviside,((I)__eHeaviside_m0p6),"heaviside smooth","--ns-heaviside","HEAVISIDE"},
		  {__ens_iinfo_dirac,((I)__eDirac_m0p6),"dirac smooth","--ns-dirac","DIRAC"},
		  {__ens_iinfo_restart,((I)0),"restart time step","--ns-restart","RESTART"},
		  {__ens_iinfo_nproc,((I)1),"nproc","-n","NPROC"});

#if 0
#include "Cmdline.h"

void ens_iinfo_from_Cmdline(I 			iinfo_[__ens_iinfo_ALL],
			    const L 	have_configfile_,
			    pCmdline const 	cmdline_)
{
  ens_iinfo i=__ens_iinfo_ERROR;
  for (++i;i<__ens_iinfo_ALL;++i)
    {
      if (__tableoptions_ens_iinfo[i].flag!=i)
	{
	  fprintf(stderr,"ens_iinfo_argv:wrong dat nbase\n");
	}
      const L have_option = Cmdline_get_integer(cmdline_,__tableoptions_ens_iinfo[i].name_cmdline,&iinfo_[i]);
      if ( (NOT have_option) AND (NOT have_configfile_) )
	{
	  iinfo_[i] = __tableoptions_ens_iinfo[i].default_value;
	}
    }
}
#endif

