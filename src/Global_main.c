#include "nsGLOBAL.h"
#include "ns_config.h"
#include "libns_exact.h"
#include "Monitor.h"


pParameters 		Global_get_Parameters		(pGlobal 		const self_)
{
  return self_->parameters; 
}

pParametersReadOnly 		GlobalReadOnly_get_Parameters		(pGlobalReadOnly 		const self_)
{
  return self_->parameters; 
}


#if 0
static Err Global_default(pGlobal const  self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif

  pParameters const Params 					= Global_get_Parameters(&self_[0]);
  
  Params->tension_method 					= __eTensionMethod_CSF;
  Params->method_transient_scheme[__eKindEquation_TRANSPORT] 	= __eTransientMethod_EULER;
  Params->method_transient_scheme[__eKindEquation_VELOCITY]  	= __eTransientMethod_EULER;
  Params->method_linearsolver[__eKindEquation_VELOCITY]  	= __eLinearSolver_DIRECT;
  Params->linfo[__ens_linfo_color]				= __emnsYES;
  Params->linfo[__ens_linfo_transport]				= __emnsYES;
  Params->rinfo[__ens_rinfo_froude]				= 1e+15;  
  sprintf(Params->sinfo[__ens_sinfo_pblmname],"%s","canal");
  return __eErr_no;
}


void Global_main(pGlobal const 		self_,
		 pCmdline const 	cmdline_,
		 STR 			errmsg_,
		 Err*			err_)
{
  /*
    0/ READ CONFIG FILE ? 
    1/ OVERWRITE CONFIG WITH COMMAND LINE 
    2/ CHECK INVALID ARGUMENT FROM COMMAND LINE 
    3/ RESTART MODE ? 
    4/ DEFAULT
    5/ DEFINE TRANSIENT SCHEME FOR EACH EQUATION    
    6/ INITIALIZE PARAMETERS
    7/ VERIFICATION
  */
  err_[0] 		= __eErr_no;
  L 	have_configfile = __emnsNO;

  //
  // Initialize structure to zero.
  //
  memset(self_,0,sizeof(nsGLOBAL_ST));

  self_->nddlelmu 	= ensBASIS_n(_shape_u,__eFace_TRIANGLE);
  self_->nddlelmp 	= ensBASIS_n(_shape_p,__eFace_TRIANGLE);
  self_->nddlelmdg	= ensBASIS_n(_shape_dg,__eFace_TRIANGLE);
  
  Workelm_init(Global_get_Workelm(self_));
  pParameters Params 	=  Global_get_Parameters(self_);
  
  /*
    0/ READ CONFIG FILE ? 
  */  
  { STR ctmp;
    if ( (Cmdline_get_string(cmdline_,"-c",ctmp)) )
      {
#ifndef NDEBUG
	Monitor_msg(0,"ns_global_main:mkmake_ns_read_configfile(%s) ...\n",ctmp);
#endif
	have_configfile = __emnsNO;
	if (NOT have_configfile)
	  {
	    Monitor_errmsg(self_->iproc,"ns_global_main:no config file '%s'",ctmp);
	    err_[0] = __eErr_user;
	    return;
	  }
	if (err_[0])
	  {
	    Monitor_errmsg(self_->iproc,"ns_global_main:mkmake_ns_read_configfile('%s') failed",ctmp);
	    return;
	  }
      } }

  if (Params->linfo[__ens_linfo_verbose])
    {
      Monitor_msg(0,"nsGLOBAL_main:step 0/5:read cmdline ...");
    }

  //
  // 1/ OVERWRITE CONFIG WITH COMMAND LINE 
  //
  err_[0] = Parameters_def(Params,
			   cmdline_,
			   have_configfile);
  if (err_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:nsPARAMS_def failed");
      return;
    }


  //
  // 2/ CHECK INVALID ARGUMENT FROM COMMAND LINE 
  //
  err_[0] = Cmdline_check_invalid(cmdline_);
  if (err_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:Cmdline_check_invalid failed");
      return;
    }
  const L  cmdline_isempty = Cmdline_isempty(cmdline_);

  
  //
  // 3/ RESTART MODE ? 
  //
  if ( (NOT Params->iinfo[__ens_iinfo_restart])
       AND
       (cmdline_isempty) )
    { 
      Monitor_errmsg(self_->iproc,"missing mesh file");
      err_[0]= __eErr_user;
      return;
    } 
  else if (Params->iinfo[__ens_iinfo_restart])
    {
      pTime const gTimeInfo = Global_get_Time(self_);
      char ctmp[512];
      sprintf(ctmp,"integrator.u.%.5"nsFORMAT_INTEGER".timeinfo.txt",Params->iinfo[__ens_iinfo_restart]);
      FILE*fich = fopen(ctmp,"r");
      int ii = fscanf(fich,""ifmt" "ifmt" %lf",&gTimeInfo->itime_interval,&gTimeInfo->itime,&gTimeInfo->ti);
      if (ii==4321)
	{}
      fclose(fich);
      sprintf(ctmp,"integrator.P1.%.5"nsFORMAT_INTEGER".mesh",gTimeInfo->itime_interval);
      extern char * strdup(const char*);
      self_->filename = strdup(ctmp);
      ns_basename(self_->filename,self_->basename);
    }
  else
    {
      self_->filename = Cmdline_get_arg(cmdline_,1);
      ns_basename(self_->filename,self_->basename);
    }
  


  //
  // 4/ DEFAULT
  //
  err_[0] = Global_default(self_);
  if (err_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:nsGLOBAL_default failed");
      return;
    }


  
  //
  // 5/ DEFINE TRANSIENT SCHEME FOR EACH EQUATION    
  //
  { eKindEquation ieq = __eKindEquation_ERROR;
    for (++ieq;ieq<__eKindEquation_ALL;++ieq)
      {
	if (ieq!=__eKindEquation_PRESSURE)
	  {
	    err_[0] = TransientScheme_def(&(self_->transient_schemes[ieq]),
					 Params->method_transient_scheme[ieq]);
	    if (err_[0])
	      {
		Monitor_errmsg(self_->iproc,"nsGLOBAL_main:TransientScheme_def failed on equation '%s' method : %s",
			  eKindEquation_get_string(ieq),
			  eTransientMethod_get_string(Params->method_transient_scheme[ieq]));
		break;
	      }
	  }
      } }  
  if (err_[0])
    {
      return;
    }
  


  //
  // 6/ INITIALIZE PARAMETERS
  //
  self_->rv[__ens_rv_iweber]		 		= ((R)1.0)/Params->rinfo[__ens_rinfo_weber];
  self_->rv[__ens_rv_ifroude]	 			= ((R)1.0)/Params->rinfo[__ens_rinfo_froude];
  self_->rv[__ens_rv_ireynold] 				= ((R)1.0)/Params->rinfo[__ens_rinfo_reynold];
  self_->rv[__ens_rv_coeff_ratio_viscosity] 		= Params->rinfo[__ens_rinfo_ratio_viscosity];
  self_->rv[__ens_rv_coeff_ratio_density]   		= Params->rinfo[__ens_rinfo_ratio_density];
  self_->rv[__ens_rv_dt] 				= Params->rinfo[__ens_rinfo_dtmin];
  self_->rv[__ens_rv_idt]				= ((R)1.0)/self_->rv[__ens_rv_dt];
  self_->rv[__ens_rv_midt]				= -self_->rv[__ens_rv_idt];
  self_->rv[__ens_rv_eps] 				= Params->rinfo[__ens_rinfo_epsmin];
  self_->rv[__ens_rv_ieps] 				= ((R)1.0)/self_->rv[__ens_rv_eps];
  self_->rv[__ens_rv_vmin] 				= ((R)1.0)/(Params->rinfo[__ens_rinfo_hmax]*Params->rinfo[__ens_rinfo_hmax]);
  self_->rv[__ens_rv_vmax] 				= ((R)1.0)/(Params->rinfo[__ens_rinfo_hmin]*Params->rinfo[__ens_rinfo_hmin]);
  

  //
  // 7/ VERIFICATION
  //
  if (NOT Params->linfo[__ens_linfo_skip_verif])
    {
      Monitor_msg(0,"");
      Monitor_msg(0," NS CONFIGURATION");
      Monitor_msg(0,"");
      Monitor_msg(0," general parameters");
      Monitor_msg(0,"        verbose    : %s",(Params->linfo[__ens_linfo_verbose])?"yes":"no");
      Monitor_msg(0,"        verify     : %s",(Params->linfo[__ens_linfo_skip_verif])?"yes":"no");
      Monitor_msg(0,"        DG         : %s",(Params->linfo[__ens_linfo_vnsdg])?"yes":"no");
      Monitor_msg(0,"        dynamic bc : %s",(Params->linfo[__ens_linfo_dynamic_boundary_condition])?"yes":"no");
      Monitor_msg(0,"        axix       : %s",(Params->linfo[__ens_linfo_axisymetric_x])?"yes":"no");
      Monitor_msg(0,"        axiy       : %s",(Params->linfo[__ens_linfo_axisymetric_y])?"yes":"no");
      Monitor_msg(0,"        nobc  vcod : %s",(Params->iinfo[__ens_iinfo_noboundary_vcod])?"yes":"no");
      Monitor_msg(0," discrete formulation");
      Monitor_msg(0,"    velocity");
      Monitor_msg(0,"        method     : %s",(Params->linfo[__ens_linfo_usupg])?"supg":"galerkin");
      Monitor_msg(0,"        nbasis     : "ifmt"",(I)self_->nddlelmu);
      Monitor_msg(0,"        slip       : %s",(Params->linfo[__ens_linfo_slip])?"yes":"no");
      Monitor_msg(0,"    pressure");
      Monitor_msg(0,"        uncoupled  : %s",(Params->linfo[__ens_linfo_pressure_uncoupled])?"yes":"no");
      Monitor_msg(0,"        freematrix : %s",(Params->linfo[__ens_linfo_pressure_freematrix])?"yes":"no");
      Monitor_msg(0,"        method     : %s",(Params->linfo[__ens_linfo_pspg])?"pspg":"galerkin");
      Monitor_msg(0,"        nbasis     : "ifmt"",(I)self_->nddlelmp);
      Monitor_msg(0,"    transport : %s",(Params->linfo[__ens_linfo_transport])?"yes":"no");
      Monitor_msg(0,"        eps        : %.3e",Params->rinfo[__ens_rinfo_epsmin]);
      Monitor_msg(0,"        marker     : %s",(Params->linfo[__ens_linfo_color])?"pseudo-concentration":"signed-distance");
      Monitor_msg(0,"        uncoupled  : %s",(Params->linfo[__ens_linfo_transport_uncoupled])?"yes":"no");
      Monitor_msg(0,"        method     : %s",(Params->linfo[__ens_linfo_transport_galerkin])?"supg":"dg");
      Monitor_msg(0,"        nbasis     : "ifmt"",(I)self_->nddlelmp);
      Monitor_msg(0,"        nbasis dg  : "ifmt"",(I)self_->nddlelmdg);
      Monitor_msg(0,"        redistance : %s",(Params->linfo[__ens_linfo_redistance])?"yes":"no");
      Monitor_msg(0,"    tension : %s",(Params->linfo[__ens_linfo_tension])?"yes":"no");
      Monitor_msg(0,"        method     : %s",eTensionMethod_get_string(Params->tension_method));
      Monitor_msg(0,"    time : %s",(Params->linfo[__ens_linfo_time_adaptivity])?"yes":"no");
      Monitor_msg(0,"        restart        : "ifmt"",Params->iinfo[__ens_iinfo_restart]);
      Monitor_msg(0,"        ntime          : "ifmt"",Params->iinfo[__ens_iinfo_ntime]);
      Monitor_msg(0,"        ntime_interval : "ifmt"",Params->iinfo[__ens_iinfo_ntime_interval]);
      Monitor_msg(0,"        velocity   : %s",eTransientMethod_get_string(Params->method_transient_scheme[__eKindEquation_VELOCITY]));
      Monitor_msg(0,"        transport  : %s",eTransientMethod_get_string(Params->method_transient_scheme[__eKindEquation_TRANSPORT]));
      Monitor_msg(0," discrete solution");
      Monitor_msg(0,"        maxiter    : "ifmt"",Params->iinfo[__ens_iinfo_newton_maxiter]);
      Monitor_msg(0,"        residu tol : %.3e",Params->rinfo[__ens_rinfo_newton_tol_residu]);
      Monitor_msg(0,"        correc tol : %.3e",Params->rinfo[__ens_rinfo_newton_tol_correc]);
      Monitor_msg(0,"        linsol     : %s",eLinearSolver_get_string(Params->method_linearsolver[__eKindEquation_VELOCITY]));
      Monitor_msg(0," mesh adaptivity : %s",(Params->linfo[__ens_linfo_mesh_adaptivity])?"yes":"no");
      Monitor_msg(0,"        maxiter    : "ifmt"",Params->iinfo[__ens_iinfo_mesh_adaptivity_maxiter]);
      Monitor_msg(0,"        hmin       : %.3e",Params->rinfo[__ens_rinfo_hmin]);
      Monitor_msg(0,"        hmax       : %.3e",Params->rinfo[__ens_rinfo_hmax]);
      Monitor_msg(0," adimensional numbers");
      Monitor_msg(0,"        viscosity ratio  : %.3e",Params->rinfo[__ens_rinfo_ratio_viscosity]);
      Monitor_msg(0,"        density   ratio  : %.3e",Params->rinfo[__ens_rinfo_ratio_density]);
      Monitor_msg(0,"        1/Re             : %.3e",self_->rv[__ens_rv_ireynold]);
      Monitor_msg(0,"        1/We             : %.3e",self_->rv[__ens_rv_iweber]);
      Monitor_msg(0,"        1/Fr             : %.3e",self_->rv[__ens_rv_ifroude]);
      Monitor_msg(0," disk data");
      Monitor_msg(0," file : %s",self_->filename);
      { char titi;
	int please = 0;
	titi = 'a';
	while ( (titi!='y') AND (titi!='n') )
	  {
	    ++please;
	    fprintf(stdout,"are you ok ? ('y' or 'n') ");
	    if (please>3)
	      {
		fprintf(stdout,"please ... are you ok ? ('y' or 'n') ");
	      }
	    if (please>4)
	      {
		fprintf(stderr,"you are stupid, go to sleep and try tomorrow\n");
		exit(1);
	      }
	    int ii = scanf("%s",&titi);
	    if (ii)
	      {

	      }
	  }    
	if (titi=='n')
	  {
	    err_[0]=__eErr_user;
	    return;
	  } }
    }
#if 0
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    ns_printf("nsGLOBAL_main:step 0/5:read cmdline done.");
#endif
#endif
  return;
}



#endif
