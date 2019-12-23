#include "ns_sys.h"
#include "Cmdline.h"
#include "ens_sinfo.h"


ConfigEnumStrings(ens_sinfo,
		  "name",
		  "pblmname",
		  "ofilename");

ConfigOptionTable(cst_pS,
		  ens_sinfo,
		  {__ens_sinfo_name,"widji","name of the process","--name","NAME"},
		  {__ens_sinfo_pblmname,"laplace","name of the problem","--pblm","PBLMNAME"},
		  {__ens_sinfo_ofilename,"ns.out","name of the output file","-o","OFILENAME"} );


void ens_sinfo_from_Cmdline(STR 			sinfo_[__ens_sinfo_ALL],
			    const L  		have_configfile_,
			    pCmdline const 	cmdline_)
  
{
  { ens_sinfo i=__ens_sinfo_ERROR;
    for (++i;i<__ens_sinfo_ALL;++i)
      {
	if (__tableoptions_ens_sinfo[i].flag!=i)
	  {
	    fprintf(stderr,"ens_sinfo_argv:wrong dat nbase\n");

	  }
	const L a = Cmdline_get_string(cmdline_,__tableoptions_ens_sinfo[i].name_cmdline,&sinfo_[i][0]);
	if ( (NOT a) AND (NOT have_configfile_) )
	  {
	    sprintf(&sinfo_[i][0],"%s",__tableoptions_ens_sinfo[i].default_value);
	  }
      } }
}
