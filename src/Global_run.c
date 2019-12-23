#include "nsGLOBAL.h"
#include "ns_constantes.h"
#include "ns_config.h"
#include "Monitor.h"


static void Global_solve_time_step(pGlobal const 	self_)
{
  
#if 0
  { R nrmdtf,nrmdtu;
    integrator_compute_derivative(&nrmdtf,&nrmdtu); }
#endif
  /*
    COMPUTE VOLUME
  */
#if 0
  /*
    COMPUTE CFL APPROXIMATION
  */
    nsGLOBAL_compute_cfl	(&self_[0]);
#endif

    { I jj=0;
      for (jj=0;jj<50;++jj)
	{
	  nsGLOBAL_NEWTON_DRIVER_ST newton_driver = {ns_build_residu,
						     ns_build_system};
	  
	  Err err = Global_newton(self_,
				  &newton_driver); 
	  if (err)
	    {
	      Monitor_errmsg(self_->iproc,"nsGLOBAL_run:nsGLOBAL_newton failed");
	      return;
	    }
	  
	  nsGLOBAL_NEWTON_RESULTS_ST*const res  = &self_->newton_results;
	  if ( (res->nrmc[self_->newton_iter]<1.0e-8)AND(res->nrmr[self_->newton_iter]<1.0e-8) )
	    break;
	} }
}


void 	nsGLOBAL_get_ddlcnc_u		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm(self_->space_u,
		 &ielm_,
		 ddl_);
}

void 	nsGLOBAL_get_ddlcnc_v		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm(self_->space_u,
		 &ielm_,
		 ddl_);
  { I i;
    for (i=0;i<_nu;++i)
      {
	ddl_[i] += self_->space_u->nddl;
      } }
}

void 	nsGLOBAL_get_ddlcnc_p		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm	(self_->space_p,
			 &ielm_,
			 ddl_);  
  { I i;
    for (i=0;i<_np;++i)
      {
	ddl_[i] += _dim*self_->space_u->nddl;
      } }
}


void nsGLOBAL_get_ddlcnc(nsGLOBAL_ST*const 	self_,
			 const I 		ielm_,
			 pI 			ddl_)
{
  nsGLOBAL_get_ddlcnc_u	(self_,ielm_,ddl_ + _ju);
  nsGLOBAL_get_ddlcnc_v	(self_,ielm_,ddl_ + _jv);
  nsGLOBAL_get_ddlcnc_p	(self_,ielm_,ddl_ + _jp);
}


void nsGLOBAL_next_time_interval(nsGLOBAL_ST*const self_)
{
  pParametersReadOnly 	const gParameters	= GlobalReadOnly_get_Parameters(self_);
  pTimeReadOnly 	const gTime		= GlobalReadOnly_get_Time(&self_[0]);
  ns_mesh * mesh				= (ns_mesh*)self_->mesh_usrptr;
  const I numNodes	= ns_mesh_get_numNodes(mesh);
  const I numEdges 	= ns_mesh_get_numEdges(mesh);

  nsSPACE_write_medit((pSpaceReadOnly)mesh->spaces[__ensBASIS_LAGRANGE_1],"integrator.P1.%.5"nsFORMAT_INTEGER"",gTime->itime_interval); 
  nsSPACE_write_medit((pSpaceReadOnly)mesh->spaces[__ensBASIS_LAGRANGE_2],"integrator.P2.%.5"nsFORMAT_INTEGER"",gTime->itime_interval); 
  if (gParameters->linfo[__ens_linfo_vnsdg])
    {
      printf("print dg_print_mesh\n");
      nsGLOBAL_dg_print_mesh_interval	(&self_[0],
					 "dg",		 
					 _shape_dg);      
    }      
  else
    {
#ifndef NDEBUG
      Monitor_msg(self_->iproc,"nsGLOBAL_run:does not call dg_print_mesh\n");
#endif
    }
  
  if  ( (gParameters->linfo[__ens_linfo_axisymetric_x]) OR (gParameters->linfo[__ens_linfo_axisymetric_y]) )
    {
      const I axi_axe 	= (gParameters->linfo[__ens_linfo_axisymetric_x])?0:1;
      pI axi_perm  	= (pI)malloc(numNodes*sizeof(I));
      pI axi_rperm 	= (pI)malloc(numNodes*sizeof(I));
      const I axi_nrot 	= 40;
      I axi_nvertex	= 0;
      I noaxi_nvertex	= 0;
      print_mesh_extrude_axi_permutation	(mesh,
						 &axi_axe,
						 axi_perm,
						 axi_rperm,
						 &axi_nvertex,
						 &noaxi_nvertex);	  
      char ctmp[512];
      sprintf(ctmp,"integrator.P1axi.%.5"nsFORMAT_INTEGER"",gTime->itime_interval);	
      print_mesh_extrude_axi		(mesh,
					 ctmp,
					 &axi_nrot,
					 &axi_axe,
					 axi_perm,
					 axi_rperm,			   			   
					 &axi_nvertex,
					 &noaxi_nvertex);   
      free(axi_perm);
      free(axi_rperm); 
    }
	
  if  ( (gParameters->linfo[__ens_linfo_axisymetric_x]) OR (gParameters->linfo[__ens_linfo_axisymetric_y]) )
    {
      const I axi_axe 	= (gParameters->linfo[__ens_linfo_axisymetric_x])?0:1;
      pI axi_perm  	= (pI)malloc((numNodes+numEdges)*sizeof(I));
      pI axi_rperm 	= (pI)malloc((numNodes+numEdges)*sizeof(I));
      const I axi_nrot 	= 40;
      I axi_nvertex	= 0;
      I noaxi_nvertex	= 0;
      print_mesh_extrude_axi_permutationP2	(mesh,
						 &axi_axe,
						 axi_perm,
						 axi_rperm,
						 &axi_nvertex,
						 &noaxi_nvertex);	  
      char ctmp[512];
      sprintf(ctmp,"integrator.P2axi.%.5"nsFORMAT_INTEGER"",gTime->itime_interval);	
      print_mesh_extrude_axiP2		(mesh,
					 ctmp,
					 &axi_nrot,
					 &axi_axe,
					 axi_perm,
					 axi_rperm,			   			   
					 &axi_nvertex,
					 &noaxi_nvertex);   
      free(axi_perm);
      free(axi_rperm); 
    }
#if 0
  /* 
     PRINT METRIC AT THE END 
  */      
  { char ctmp[512];
    sprintf(ctmp,"integrator.metric.%.5"nsFORMAT_INTEGER"",gTime->itime_interval);	
    nsMETRIC_print	(self_->metric_end_time_interval,ctmp); }
#endif
}




void Global_runOld		(pGlobal const 		self_,
				 eStateSolver 		job_[1],
				 cst_pI 		rwork_n_,
				 pR 			rwork_,
				 cst_pI 		iwork_n_,
				 pI 			iwork_)  
{
  if (self_->required_rwork_n<rwork_n_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_run:not enough rwork_");
    }
  pParametersReadOnly const gParameters	= GlobalReadOnly_get_Parameters(self_);
  pTime const gTimeInfo 		= Global_get_Time(self_);

#ifndef NDEBUG
  Monitor_msg(self_->iproc,"nsGLOBAL_run:input job = "efmt"",job_[0]);
#endif

  const L dynamic_boundary_condition 	= ParametersReadOnly_getl(gParameters,
								  __ens_linfo_dynamic_boundary_condition);
  const L slip_condition 		= ParametersReadOnly_getl(gParameters,
								  __ens_linfo_slip);
  const L firstTimeStep 		= (gTimeInfo->itime==0) ? __emnsYES : __emnsNO;


  /*__ WHAT HAS BEEN DONE  ____________________________________*/
  switch (job_[0])
    {
      
    case __eStateSolver_initial_condition:
      {
	goto state___eStateSolver_initial_condition;
      }
      
    case __eStateSolver_initial_condition_dg:
      {
	goto state___eStateSolver_initial_condition_dg;
      }
      
    case __eStateSolver_slip_velocity_symbolic:
      {
	if (gParameters->linfo[__ens_linfo_slip])
	  {	      
	    goto state___eStateSolver_slip_velocity_symbolic;
	  }
	else
	  {
	    fprintf(stderr,"errororor slip\n");
	    exit(1);
	  }
	break;
      }
      
    case __eStateSolver_dirichlet_pressure_symbolic:
      {
	goto state___eStateSolver_dirichlet_pressure_symbolic;
      }
      
    case __eStateSolver_dirichlet_velocity_symbolic:
      {
	goto state___eStateSolver_dirichlet_velocity_symbolic;
      }
      
    case __eStateSolver_dirichlet_marker_symbolic:
      {
	goto state___eStateSolver_dirichlet_marker_symbolic;
      }
      
    case __eStateSolver_dirichlet_pressure:
      {
	goto state___eStateSolver_dirichlet_pressure;
      }
      
    case __eStateSolver_dirichlet_velocity:
      {
	goto state___eStateSolver_dirichlet_velocity;
      }
      
    case __eStateSolver_dirichlet_marker:
      {
	goto state___eStateSolver_dirichlet_marker;
      }
      
    case __eStateSolver_dirichlet_dg:
      {
	goto state___eStateSolver_dirichlet_dg;
      }
      
    case __eStateSolver_ERROR:
      {
	break;
      }
    case __eStateSolver_exit:
      {
	Monitor_errmsg(self_->iproc,"nsGLOBAL_run:wrong input job(=__eStateSolver_exit)");
	return;
      }
      
    case __eStateSolver_ALL:
      {
	fprintf(stderr,"integrator job error "efmt"\n",job_[0]);
	exit(1);
      }
    }


  Time_init(gTimeInfo);
  /*___ NEXT TIME INTERVAL _______________________________________*/
  
 state_next_time_interval:  
  Time_next_time_interval(gTimeInfo,
			    self_->rv[__ens_rv_dt]*((R)gParameters->iinfo[__ens_iinfo_ntime]));
  Monitor_msg(self_->iproc,"time interval [%e,%e], mass %e",gTimeInfo->time_interval_t0,gTimeInfo->time_interval_tf,
	    ((R)0.0));
  /*nsSTATISTICS_get(&Z->stats,Z->adapt_iter,(Time->itime>1)?(Time->itime-1):0,__ens_rres_area)*/
  self_->adapt_iter = negal0;
  


  /*__ RESTART TIME INTERVAL ____________________________________*/
 state_restart_time_interval:      
  Time_restart_time_interval(gTimeInfo);
#if 0
  fprintf(stdout,"restart interval["ifmt"] [%e,%e,%e||||,%e,%e]\n",gTimeInfo->itime_interval,tii,ti,t,time_interval_t0,time_interval_tf);
#endif
  if (gParameters->linfo[__ens_linfo_vnsdg])
    {
      if ( (gParameters->iinfo[__ens_iinfo_ntime]>0) AND (gTimeInfo->itime_interval==0) )
	{
	  self_->ielm = 0;      
	redo_state___eStateSolver_initial_condition_dg:
#ifdef TONTON
	  
	  self_->dgcoo_n	= self_->cubature_triangle.q_n;
	  self_->cooelm[0] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+0]-1)+0];
	  self_->cooelm[3] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+0]-1)+1];
	  self_->cooelm[1] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+1]-1)+0];
	  self_->cooelm[4] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+1]-1)+1];
	  self_->cooelm[2] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+2]-1)+0];
	  self_->cooelm[5] 	= self_->ddlcoo[2*(self_->ddlcnc[6*self_->ielm+2]-1)+1];
	  { I i;
	    for (i=0;i<self_->cubature_triangle.q_n;++i)
	      {
		self_->dgcoo[2*i]   = (((R)1.0)-self_->cubature_triangle.q_p[i]-self_->cubature_triangle.q_p[self_->cubature_triangle.q_n+i])*self_->cooelm[0]+self_->cubature_triangle.q_p[i]*self_->cooelm[1]+self_->cubature_triangle.q_p[self_->cubature_triangle.q_n+i]*self_->cooelm[2];
		self_->dgcoo[2*i+1] = (((R)1.0)-self_->cubature_triangle.q_p[i]-self_->cubature_triangle.q_p[self_->cubature_triangle.q_n+i])*self_->cooelm[3]+self_->cubature_triangle.q_p[i]*self_->cooelm[4]+self_->cubature_triangle.q_p[self_->cubature_triangle.q_n+i]*self_->cooelm[5];
	      } }
#endif
	  job_[0] = __eStateSolver_initial_condition_dg;
	  return;
	  
	state___eStateSolver_initial_condition_dg:
#if 0
	  if (gParameters->linfo[__ens_linfo_color])
	    {	 
     
	      ns_heaviside(gParameters->iinfo[__ens_iinfo_heaviside],self_->cubature_triangle.q_n,rwork_,&self_->rv[__ens_rv_eps]);
	    }

	  nsblas_dgemv(transN,&self_->dg_n,&self_->cubature_triangle.q_n,&regal1,self_->qelm_dg_wbasis,&self_->dg_n,rwork_,&negal1,&regal0,&self_->var_dg.xi[self_->ielm*self_->dg_n],&negal1); 
	  nsblas_dcopy(&self_->dg_n,&self_->var_dg.xi[self_->ielm*self_->dg_n],&negal1,&self_->var_dg.x[self_->ielm*self_->dg_n],&negal1);	 	  
#endif
	  ++self_->ielm;
	  if (self_->ielm<self_->nelm)
	    goto redo_state___eStateSolver_initial_condition_dg;	  
	}

#ifndef NDEBUG
      Monitor_msg(self_->iproc,"nsGLOBAL_run:initial dg condition done");
#endif
    }



  /*___ INITIAL CONDITION  _______________________________________*/
  /* IF NOT STEADY STATE AND FIRST TIME INTERVAL */
  if ( (gParameters->iinfo[__ens_iinfo_ntime]>0) AND (gTimeInfo->itime_interval==0) )
    {
      job_[0] = __eStateSolver_initial_condition;
      return;
    }  

 state___eStateSolver_initial_condition:
#ifndef NDEBUG
  Monitor_msg(self_->iproc,"nsGLOBAL_run:initial condition back");
#endif  
  if ( (gParameters->iinfo[__ens_iinfo_ntime]>0) AND (gTimeInfo->itime_interval==0) )
    {
#if 0
      nsMETRIC_setidentity	(self_->metric);
#endif
      ns_var_zeroi		(&self_->var_u);
      ns_var_zeroi		(&self_->var_p);      
      ns_var_cpy_xi2x		(&self_->var_u);
      ns_var_cpy_xi2x		(&self_->var_p);      
    }  
  /*___ FIN INITIAL CONDITION  _______________________________________*/

  /*
    DANS LE CAS OU ON A INTERPOLE
    ON REPASSE LES VARIABLES xi et xii A LA MOULINETTE IDENTITE
  */
#if 0
  if (ns_linfo[__ens_linfo_mesh_adaptivity])
    {
      if (self_->itime_interval>0)
	{
	  /*
	    REDEFINE INTERPOLATE SOLUTION
	  */  
	  printf("allo\n");
	  build_residu_interpolate	(mesh->linsys_vectors.rhs,mesh->linsys_vectors.x,mesh->linsys_vectors.xi,mesh->linsys_vectors.xii,&self_->var_volumic);    
	  build_system_interpolate	(&negal3,mesh->linsys_vectors.x,&mesh->linsys_vectors.N); 
	  ns_linsys_compute_constraint	(&mesh->linsys_A,&mesh->sparse_stokes.A,mesh->linsys_vectors.rhs,mesh->linsys_vectors.c);
	}	  
    }
#endif


  /*__________________________________________________________ PROCHAIN PAS DE TEMPS _______________________________________*/  
 state_compute_next_time:
  {        
    Time_next_time(gTimeInfo,
		   gParameters->rinfo[__ens_rinfo_dtmin]);
    
    self_->rv[__ens_rv_dti] 			= TimeReadOnly_get_dti(gTimeInfo);
    self_->rv[__ens_rv_dt]  			= TimeReadOnly_get_dt(gTimeInfo);
    self_->rv[__ens_rv_idt] 			= ((R)1.0)/self_->rv[__ens_rv_dt];
    self_->rv[__ens_rv_ratio_dti]   		= self_->rv[__ens_rv_idt]*self_->rv[__ens_rv_dti];
    self_->rv[__ens_rv_iratio_dti]		= ((R)1.0)/self_->rv[__ens_rv_ratio_dti];
    self_->rv[__ens_rv_midt]			= -self_->rv[__ens_rv_idt];
    
    self_->ruser[ruser_t] 			= TimeReadOnly_get_t(gTimeInfo);
    self_->ruser[ruser_ti] 			= TimeReadOnly_get_ti(gTimeInfo);
    self_->ruser[ruser_tii] 			= TimeReadOnly_get_tii(gTimeInfo); 


    /* _____ SYMBOLIC PRESSURE ___________________________________________________ */
    if ( (dynamic_boundary_condition) OR (firstTimeStep) )
      {
	job_[0] = __eStateSolver_dirichlet_pressure_symbolic;
	return;
      }
    
  state___eStateSolver_dirichlet_pressure_symbolic:
    /* ____________________________________________________________________________ */

    
    /* _____ SYMBOLIC EULERIAN MARKER _____________________________________________ */
    if ( (dynamic_boundary_condition) OR (firstTimeStep) )
      {
	job_[0] = __eStateSolver_dirichlet_marker_symbolic;
	return;
      }
    
  state___eStateSolver_dirichlet_marker_symbolic:
    /* ____________________________________________________________________________ */
    
    
    /* _____ SYMBOLIC VELOCITY  ___________________________________________________ */
    if ( (dynamic_boundary_condition) OR (firstTimeStep) )
      {
	job_[0] = __eStateSolver_dirichlet_velocity_symbolic;
	return;    
      } 
  state___eStateSolver_dirichlet_velocity_symbolic:
    /* ____________________________________________________________________________ */

    
    /* _____ SLIP SYMBOLIC VELOCITY _______________________________________________ */
    if ( ( (dynamic_boundary_condition) OR (firstTimeStep) )
	 AND (slip_condition) )
      {
	job_[0] = __eStateSolver_slip_velocity_symbolic;
	return;
      }    
  state___eStateSolver_slip_velocity_symbolic:
    /* ____________________________________________________________________________ */

    
    
    if ( (dynamic_boundary_condition) OR (firstTimeStep) )
      {
#ifndef NDEBUG
	Monitor_msg(self_->iproc,"nsGLOBAL_run:slip velocity sym back");
#endif    
      }

    /* valeur des conditions de dirichlets (peut dependre du temps, donc il faut les recharger ) */    
    /* I L FAUT METTTRE UN FLAG POUR LA DEPENDENCE EN TEMPS DES CONDITIONS DE DIRICHLETS */
    job_[0] = __eStateSolver_dirichlet_pressure;
    return;
  state___eStateSolver_dirichlet_pressure:
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:dirichlet pressure back job "efmt"",job_[0]);
#endif    
    /* ____________________________________________________________________________ */

    
    /* ____________________________________________________________________________ */
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:go to job "efmt"",job_[0]);
#endif 
    job_[0] = __eStateSolver_dirichlet_velocity;
    return;        
  state___eStateSolver_dirichlet_velocity:    
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:dirichlet velocity back");
#endif
    /* ____________________________________________________________________________ */

    
    /* ____________________________________________________________________________ */
    if (gParameters->linfo[__ens_linfo_transport])
      {
	job_[0] = __eStateSolver_dirichlet_marker;
	return;    
      }
  state___eStateSolver_dirichlet_marker:
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:dirichlet marker back");
#endif
    /* ____________________________________________________________________________ */
    
    if (gParameters->linfo[__ens_linfo_transport])
      {
	if (gParameters->linfo[__ens_linfo_vnsdg])
	  {
	    job_[0] = __eStateSolver_dirichlet_dg;
	    return;    
	  }
      }
  state___eStateSolver_dirichlet_dg:
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:dirichlet dg back");
#endif
    
#ifndef NDEBUG
    Monitor_msg(self_->iproc,"nsGLOBAL_run:solve time step ...");
#endif


    Global_solve_time_step(self_);
    
#if 0
    { I ielm;
      for (ielm=0;ielm<self_->pointeur_maillage->nelm;++ielm)
	{
	  const R nx = self_->var_dgnormal.x[ielm*_ndg+0];
	  const R ny = self_->var_dgnormal.y[ielm*_ndg+0];
	  self_->var_dgnormal.x[ielm*_ndg+0]=self_->var_dgnormal.x[ielm*_ndg+0]/nsSQRT(nx*nx+ny*ny)/nsSQRT(2.0);
	  self_->var_dgnormal.y[ielm*_ndg+0]=self_->var_dgnormal.y[ielm*_ndg+0]/nsSQRT(nx*nx+ny*ny)/nsSQRT(2.0);
	}
    }
#endif

#if 0
    {
      R cooelm[12] 	= {0.0,0.1,0.0,0.5,0.5,0.0,
				   0.0,0.0,1.0,0.0,0.8,0.5};    
      I nn= 60;
      R pt[200*2]; 
      I i;
      for (i=0;i<nn;++i)
	{
	  pt[i] 	= cos(2.0*acos(-1.0)* ((R)i)/((R)(nn)));
	  pt[nn+i] 	= sin(2.0*acos(-1.0)* ((R)i)/((R)(nn)));
	}
      FILE*fich = fopen("titi","w");
      for (i=0;i<nn;++i)
	{
	  const R r = pt[i];
	  const R s = pt[nn+i];
	  R b_[6];
	  b_[0] 		= s*(2.0e0*s+4.0e0*r-3.0e0)+r*(2.0e0*r-3.0e0)+1.0e0;
	  b_[1] 		= r*(2.0e0*r-1.0e0);
	  b_[2] 		= s*(2.0e0*s-1.0e0);
	  b_[3] 		= r*(4.0e0-4.0e0*r)-4.0e0*r*s;
	  b_[4] 		= 4.0e0*r*s;
	  b_[5] 		= s*(-4.0e0*s-4.0e0*r+4.0e0);
	  static const I n6=6;
	  R x 		= ddot(&n6,cooelm,&negal1,b_,&negal1);
	  R y 		= ddot(&n6,cooelm+6,&negal1,b_,&negal1);
	  fprintf(fich,"%e %e\n",x,y);
	}
      fclose(fich);
    }
#endif    

    /*    nsGLOBAL_compute_volume	(&self_[0]);*/
    if ( (gTimeInfo->itime>0) AND ((gTimeInfo->itime+1)%(gParameters->iinfo[__ens_iinfo_ntime])==0) )
      {
#if 0
	printf("reinit\n");
	printf("reinit done\n");
#endif
      }

    if (gParameters->iinfo[__ens_iinfo_ntime]>1)
      {	
	
	if ( (gParameters->linfo[__ens_linfo_redistance]) AND (gTimeInfo->itime>0) AND (gTimeInfo->itime%5==5) )
	  {
	    /* 
	       COMPUTE DISCRETIself_ATION
	    */	
	  }
	/* 
	   COMPUTE METRIC 
	*/
#if 0
	integrator_metric();
#endif	

	/* 
	   PRINT   VARIABLES  
	*/
	nsGLOBAL_transient_print(self_);

	/* next */
	++gTimeInfo->itime;

	LinsysVectors_update_previous_x(self_,
					self_->linsysVectors);

	if (gTimeInfo->itime<(gTimeInfo->itime_interval+1)*gParameters->iinfo[__ens_iinfo_ntime])
	  {
	    goto  state_compute_next_time;
	    exit(1);
	  }
      }       
  }     
  
  /* MESH ADAPTIVITY ___________________________________________________________________*/
  if (gParameters->linfo[__ens_linfo_mesh_adaptivity])
    {
      if ( ((self_->adapt_iter<4)AND(gTimeInfo->itime_interval==0)) OR (self_->adapt_iter<gParameters->iinfo[__ens_iinfo_mesh_adaptivity_maxiter]) )
	{	
	  nsGLOBAL_change_mesh(self_,self_->basename);
	  if (gTimeInfo->itime_interval>0)
	    {
#if 0
	      ns_transfert_solution(__ens_method_transfert_static);	  
#endif
	    }	  
	  /* 
	     RESTART TIME INTERVAL
	     NEXT ADAPTIVITY ITERATION
	  */  	  
	  gTimeInfo->itime = gTimeInfo->itime_interval*gParameters->iinfo[__ens_iinfo_ntime];	  
	  ++self_->adapt_iter;	  
	  
	  
	  if (gTimeInfo->itime_interval==0)
	    {
	      /*	      nsMASS_reset_iter(&self_->mass_info);*/
	    }	  
	  
	  
	  fprintf(stdout,"itime set to "ifmt"\n",gTimeInfo->itime);	  
	  goto  state_restart_time_interval;
	}      
    } 

  /*__________________________________________________________________________*/
  if (gTimeInfo->itime_interval<=gParameters->iinfo[__ens_iinfo_ntime_interval])
    {   
      nsGLOBAL_next_time_interval(self_);
#if 0
      if (self_->metric)
	{
	  nsMETRIC_copy	(self_->metric,
			 self_->metric_end_time_interval);
	}
#endif
      
      Time_next_interval(gTimeInfo);
      if (gTimeInfo->itime_interval<=gParameters->iinfo[__ens_iinfo_ntime_interval])
	{
	  goto  state_next_time_interval;
	}
    }
  
  job_[0] = __eStateSolver_exit;
  return;
}

