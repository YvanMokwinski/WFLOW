#include "ns_sys.h"
#include "ens_linfo.h"

ConfigEnumStrings(ens_linfo,
		  "verbose",
		  "mesh_adaptivity",
		  "time_adaptivity",
		  "transport",
		  "transport_uncoupled",
		  "transport_galerkin",
		  "pressure_uncoupled",
		  "pressure_freematrix",
		  "slip",
		  "color",
		  "tension",
		  "redistance",
		  "pspg",
		  "usupg",
		  "axisymetric_x",
		  "axisymetric_y",
		  "dynamic_boundary_condition",
		  "vnsdg",
		  "skip_verif"		  );

ConfigOptionTable(L,
		  ens_linfo,
		  {__ens_linfo_verbose,__emnsNO,"activate verbose","-v","VERBOSE"},
		  {__ens_linfo_mesh_adaptivity,__emnsNO,"mesh adaptivity","--mesh-adaptivity","MESH_ADAPTIVITY"},
		  {__ens_linfo_time_adaptivity,__emnsNO,"time adaptivity","--time-adaptivity","TIME_ADAPTIVITY"},
		  {__ens_linfo_transport,__emnsYES,"apply transport","--transport","TRANSPORT"},
		  {__ens_linfo_transport_uncoupled,__emnsNO,"apply uncoupled transport","--transport-uncoupled","TRANSPORT_UNCOUPLED"},
		  {__ens_linfo_transport_galerkin,__emnsYES,"apply dg for transport","--poulou","DG"},
		  {__ens_linfo_pressure_uncoupled,__emnsNO,"uncouple pressure","--pressure-uncoupled","PRESSURE_UNCOUPLED"},
		  {__ens_linfo_pressure_freematrix,__emnsNO,"do not build pressure B matrix ","--pfreem","PFREEM"},
		  {__ens_linfo_slip,__emnsNO,"apply slip boundary condition","--slip","SLIP"},
		  {__ens_linfo_color,__emnsYES,"choose color signed distance","--sd","SD"},
		  {__ens_linfo_tension,__emnsYES,"do not apply tension","--notension","NOTENSION"},
		  {__ens_linfo_redistance,__emnsNO,"apply redistance algorithm","--redistance","REDISTANCE"},
		  {__ens_linfo_pspg,__emnsNO,"apply PSPG to pressure equations","--pspg","PSPG"},
		  {__ens_linfo_usupg,__emnsNO,"apply SUPG to velocity equations","--usupg","USUPG"},
		  {__ens_linfo_axisymetric_x,__emnsNO,"set axisymetric x","--axix","AXIX"},
		  {__ens_linfo_axisymetric_y,__emnsNO,"set axisymetric y","--axiy","AXIY"},
		  {__ens_linfo_dynamic_boundary_condition,__emnsNO,"no no ","--no","DYNAMIC_BOUNDARY_CONDITION"},
		  {__ens_linfo_vnsdg,__emnsNO,"enable vnsdg","--dg","VNSDG"},
		  {__ens_linfo_skip_verif,__emnsNO,"skip verif","-f","SKIP_VERIF"} );


#if 0
#include "Cmdline.h"
void ens_linfo_from_Cmdline(L 		linfo_[__ens_linfo_ALL],
			    const L 	have_configfile_,
			    pCmdline const 	cmdline_)
{
  ens_linfo i = __ens_linfo_ERROR;
  for (++i;i<__ens_linfo_ALL;++i)
    {
      if (__tableoptions_ens_linfo[i].flag!=i)
	{
	  fprintf(stderr,"ens_linfo_argv:wrong dat nbase\n");
	}
      const L a = Cmdline_get_logical(cmdline_,__tableoptions_ens_linfo[i].name_cmdline);
      if (a)
	{
	  linfo_[i] = (__tableoptions_ens_linfo[i].default_value)?__emnsNO:__emnsYES;
	}
      else
	{
	  if ( (NOT a) AND (NOT have_configfile_) )
	    {
	      linfo_[i] = __tableoptions_ens_linfo[i].default_value;
	    }
	}
    }
}
#endif
