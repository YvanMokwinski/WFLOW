
#include "nsGLOBAL.h"
#include "ns_config.h"
#include "libns_exact.h"
#include "Monitor.h"



cst_pParameters cst_Global_get_Parameters	(cst_pGlobal self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->parameters;
}

pParameters Global_get_Parameters	(pGlobal self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->parameters;
}

cst_pTime	cst_Global_get_Time	(cst_pGlobal self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->timeInfo;
}

pTime		Global_get_Time		(pGlobal self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->timeInfo;
}


static void ns_basename(const char * 	t_,
			char		o_[512])
{
  I j = 0;
  while ((j<512)AND(t_[j]!='\0'))++j;
  I n = ((j<512)?j:0);
  j=0;
  { I i;for (i=0;i<n;++i) if (t_[i]=='.') j=i; }
  const char * pp = (!j)? NULL: (&t_[j+1]);
  if (!pp) j = n; 
  { I i;for (i=0;i<j;++i) { o_[i] = t_[i]; } }
  o_[j] = '\0';
}

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
  sprintf(Params->sinfo[__ens_sinfo_pblmname],"%s","canal");
  Params->linfo[__ens_linfo_color]				= __emnsYES;
  Params->linfo[__ens_linfo_transport]				= __emnsYES;
  Params->rinfo[__ens_rinfo_froude]				= 1e+15;
  return __eErr_no;
}


void nsWORKELM_init(nsWORKELM_ST*const workelm_)
{
  workelm_->uelm 	= &workelm_->ddlelm[_ju];
  workelm_->velm 	= &workelm_->ddlelm[_jv];
  workelm_->pelm 	= &workelm_->ddlelm[_jp];

  workelm_->uelmi 	= &workelm_->ddlelmi[_ju];
  workelm_->velmi 	= &workelm_->ddlelmi[_jv];
  workelm_->pelmi 	= &workelm_->ddlelmi[_jp];

  workelm_->uelmii 	= &workelm_->ddlelmii[_ju];
  workelm_->velmii 	= &workelm_->ddlelmii[_jv];
  workelm_->pelmii 	= &workelm_->ddlelmii[_jp];
}

Err Global_main(pGlobal const 	self_,
		pCmdline const 	cmdline_)
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
  Err 	err 		= __eErr_no;
  L 	have_configfile = __emnsNO;
  memset(self_,0,sizeof(nsGLOBAL_ST));

  self_->nddlelmu 	= ensBASIS_n(_shape_u,__eElement_TRIA);
  self_->nddlelmp 	= ensBASIS_n(_shape_p,__eElement_TRIA);
  self_->nddlelmdg	= ensBASIS_n(_shape_dg,__eElement_TRIA);
  
  nsWORKELM_init(&self_->workelm);
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
#if 0
	have_configfile = mkmake_ns_read_configfile	(ctmp,
							 Params->rinfo,
							 Params->iinfo,
							 Params->sinfo,
							 Params->linfo,
							 &err);
#endif
	if (NOT have_configfile)
	  {

	    
	    Monitor_errmsg(self_->iproc,"ns_global_main:no config file '%s'",ctmp);
	    return __eErr_user;
	  }
	if (err)
	  {
	    Monitor_errmsg(self_->iproc,"ns_global_main:mkmake_ns_read_configfile('%s') failed",ctmp);
	    return err;
	  }
      } }

  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(0,"nsGLOBAL_main:step 0/5:read cmdline ...");

  /* 
     1/ OVERWRITE CONFIG WITH COMMAND LINE 
  */
  err = Parameters_def(Params,
		       cmdline_,
		       have_configfile);
  if (err)
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:nsPARAMS_def failed");
      return err;
    }


  /*
    2/ CHECK INVALID ARGUMENT FROM COMMAND LINE 
  */
  err = Cmdline_check_invalid(cmdline_);
  if (err)
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:Cmdline_check_invalid failed");
      return err;
    }
  const L  cmdline_isempty = Cmdline_isempty(cmdline_);

  
  /*
    3/ RESTART MODE ? 
  */
  if ( (NOT Params->iinfo[__ens_iinfo_restart])
       AND
       (cmdline_isempty) )
    { 
      Monitor_errmsg(self_->iproc,"missing mesh file");
      return __eErr_user;
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
  


  /*
    4/ DEFAULT
  */
  err = Global_default(self_);
  if (err)
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_main:nsGLOBAL_default failed");
      return err;
    }


  
  /*
    5/ DEFINE TRANSIENT SCHEME FOR EACH EQUATION    
  */
  { eKindEquation ieq = __eKindEquation_ERROR;
    for (++ieq;ieq<__eKindEquation_ALL;++ieq)
      {
	if (ieq!=__eKindEquation_PRESSURE)
	  {
	    err = TransientScheme_def(&(self_->transient_schemes[ieq]),
					 Params->method_transient_scheme[ieq]);
	    if (err)
	      {
		Monitor_errmsg(self_->iproc,"nsGLOBAL_main:TransientScheme_def failed on equation '%s' method : %s",
			  eKindEquation_get_string(ieq),
			  eTransientMethod_get_string(Params->method_transient_scheme[ieq]));
		break;
	      }
	  }
      } }  
  if (err)
    {
      return err;
    }
  


  /*
    6/ INITIALIZE PARAMETERS
  */
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
  

  /*
    7/ VERIFICATION
  */
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
	    return __eErr_user;
	  } }
    }
#if 0
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    ns_printf("nsGLOBAL_main:step 0/5:read cmdline done.");
#endif
#endif
  return err;
}



Err Global_precompute(pGlobal 		const self_,
		      ns_mesh * 	const mesh_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif

  Err 		err 	= __eErr_no;
  pParameters 	Params 	=  Global_get_Parameters(self_);
  /*1 ###############################################*/
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 1/5:reading mesh from file '%s'",self_->filename);
#endif
  self_->mesh_usrptr = mesh_;
  nsWORKELM_init(&self_->workelm);
  /*2 ###############################################*/
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data ...");
#endif
  err = nsFEM_DATA_def(&self_->fem,
		       mesh_->elm);  
  if (err)
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_precompute:ns_exact_init failed");
      return err;
    }
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data done");
#endif
  /*3 ###############################################*/
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize ...");
#endif
  err = nsGLOBAL_initialize(self_,
			    self_->mesh_usrptr);
  if (err)
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_precompute:nsGLOBAL_initialize failed");
      return err;
    }
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
  Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize done.");
#endif

  /*4 ###############################################*/
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 4/5:nsCSF_INFO_def ...");
#endif
  CsfInfo_def	(&self_->csfInfo,
		 Params->linfo[__ens_linfo_color],
		 __eOperator_ih,
		 __eOperator_qh,
		 __eOperator_ih,
		 ((R)0.125),
		 __eOperator_qh,
		 ((R)0.99),
		 __eOperator_ph,
		 __eOperator_ph,
		 0.5,
		 __ens_method_capturing_linear,
		 1,
		 __eSmoothingMethod_LQ_l,
		 1,
		 (eHeaviside)Params->iinfo[__ens_iinfo_heaviside],
		 (eDirac)Params->iinfo[__ens_iinfo_dirac],
		 &self_->rv[__ens_rv_eps]);  
  
#ifndef NDEBUG
  if (Params->linfo[__ens_linfo_verbose])
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 4/5:CsfInfo_def done.");
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 5/5:completed.");
    }
#endif
  return err;
}
