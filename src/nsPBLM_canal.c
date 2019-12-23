#include "ns_sys.h"
#include "nsGLOBAL.h"


void boundaryf(cst_pI 		n_,
	       cst_pI 		xm_,
	       pR     		x_,
	       cst_pI       	xoff_,
	       cst_pI       	gl_m_,
	       cst_pR       	gl_,
	       cst_pI       	gloff_,
	       cst_pI       	cod_,
	       cst_pR       	t_,
	       void * 		usrptr_)
{
  I i;
  for (i=0;i<n_[0];++i)
    {
      x_[i] = gl_[gloff_[0]+i];
#if 1
      if (x_[i]>0.0)
	x_[i]=1.0;
      else if (x_[i]<0.0)
	x_[i]=0.0;
      else
	x_[i]=0.5;
#endif
#if 0
      ns_heaviside(global[0].params.iinfo[__ens_iinfo_heaviside],1,&x_[i],&ns_rv[__ens_rv_eps]);
#endif
    }
}

void boundaryu(cst_pI 		n_,
	       cst_pI 		xm_,
	       pR     		x_,
	       cst_pI       	xoff_,
	       cst_pI       	gl_m_,
	       cst_pR       	gl_,
	       cst_pI       	gloff_,
	       cst_pI       	cod_,
	       cst_pR       	t_,
	       void * 		usrptr_)
{
  I i;
  for (i=0;i<n_[0];++i)
    {
      const R x = gl_[i];
      const R y = gl_[gloff_[0]+i];
      if (y>0.5)
	{
	  x_[i] = ((((R)0.5)-y)*(y-1.0)*((R)16.0));
	  if (x_[i]<0.0)
	    {
	      printf("AAAAAAAAAAAA %e %e %e\n",x_[i],x,y);
	      exit(1);
	    }
	}
      else
	{
	  x_[i] = ((((R)0.5)-y)*(y)*((R)16.0));
	  if (x_[i]<0.0)
	    {
	      printf("AAAAAAAAAAAA222 %e %e %e\n",x_[i],x,y);
	      exit(1);
	    }
	}
      x_[xoff_[0]+i] = ((R)0.0);
    }
}


void use_integrator_probleme_canal(eStateSolver job_[1],
				   void * 	usrptr_)
{
  //  nsGLOBAL_ST*const Z = &global[iproc_];
    nsGLOBAL_ST*const Z = usrptr_;
  static const I dim = ((I)_dim);
  pParameters const Params	= Global_get_Parameters(Z);
  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(Z);


  pSparseStokes M = &Z->sparseStokes;

  const I 	rwork_n			= Z->required_rwork_n;
  const I 	iwork_n			= Z->required_iwork_n;
  pR 		rwork			= (pR)malloc((rwork_n+1)*sizeof(R));
  pI 		iwork			= (pI)malloc((iwork_n+1)*sizeof(I));

  Params->method_transient_scheme[__eKindEquation_VELOCITY]	= __eTransientMethod_NO;
  Params->method_transient_scheme[__eKindEquation_TRANSPORT]	= __eTransientMethod_NO;
  const ns_mesh*mesh = (const ns_mesh*)Z->mesh_usrptr;
  const I nelm  = ns_mesh_nelm(mesh);
  const I nddlu = Z->nddlu;
#if 0
  const R 	current_t 		= Z->ruser[ruser_t];
  const I 	nddlp 			= Z->nddlp;
#endif
  pVariables variables = Global_get_Variables(Z);
  job_[0] 				= __eStateSolver_ERROR;
L 	we_continue 		= __emnsYES;
  while (we_continue)    
    {
      Global_runOld(//&global[iproc_],
		 Z,
		 job_,
		 &rwork_n,
		 rwork,
		 &iwork_n,
		 iwork);
#ifndef NDEBUG      
      printf("---> job_[0] = "efmt"\n",job_[0]);
#endif            
      switch(job_[0])
	{
	case __eStateSolver_slip_velocity_symbolic:
	  {
	    M->nddls_dirichlet = 0;
	    break;
	  }
	  
	case __eStateSolver_dirichlet_dg:
	  {
#if 0
	    pTimeReadOnly const gTimeInfo	= cst_Global_get_Time(Z);
	    R uvw[128];
	    R loc[128];
	    R pos[128];
	    if (0)
	    memset(Z->dg_rhs,0,sizeof(R)*Z->nddlelmdg*nelm);
#endif
	    { I jelm;
	      for (jelm=0;jelm<nelm;++jelm)
		{
		  { I jadj;
		    for (jadj=0;jadj<3;++jadj)
		      {
			I kelm = mesh->adj[3*jelm+jadj];
			if (kelm==0)
			  {
			    if (Z->ddlcod[Z->ddlcnc[jelm*6+3+jadj]-1]==23)
			      {

#ifdef TONTON
				const I N = ns_dg_fluxelm_boundary_get_coo(&global[iproc_],jelm,jadj,pos);
#if 0
				I i;
				for (i=0;i<N;++i){printf("%e %e\n",pos[i],pos[N+i]);}
#endif
				
				boundaryu(&N,&dim,uvw,&N,&dim,pos,&N,NULL,&gTimeInfo->t,NULL);
				boundaryf(&N,&negal1,loc,&negal1,&dim,pos,&N,NULL,&gTimeInfo->t,NULL);

				ns_dg_fluxelm_boundary(&global[iproc_],&mregal1,jelm,jadj,Z->dg_rhs,loc,uvw);
				
				{
				  I i;
				  for (i=0;i<N;++i)
				    {
				      loc[i] = 1.0-loc[i];
				    }
				}

				ns_dg_fluxelm_boundary(&global[iproc_],&mregal1,jelm,jadj,Z->dg_rhs2,loc,uvw);
#endif
			      }
			  }
		      } }
		} }
	    break;
	  }
	  
	case __eStateSolver_initial_condition_dg:
	  {
#if 0
	    R loc[128];
	    R pos[128];
	    { I j;
	      for (j=0;j<nelm;++j)
		{
		  const I N = ns_dg_projectionelm_get_coo(jelm,pos);
		  { I k;
		    for (k=0;k<N;++k)
		      {
			loc[k] = pos[N+k];
		      } }
		  ns_dg_projectionelm(&mregal1,jelm,rwork,loc);
		} }
#endif
	    break;
	  }
	  
	case __eStateSolver_initial_condition:
	  {
	    /* DO NOTHING */
	    break;
	  }
	  
	case __eStateSolver_dirichlet_marker:
	  {
	    /* DO NOTHING */
	    break;
	  }
	  
	case __eStateSolver_dirichlet_marker_symbolic:
	  {
	    /* DO NOTHING */
	    break;
	  }	    

	case __eStateSolver_dirichlet_velocity_symbolic:
	  {
	    M->nddlu_dirichlet = 0;
	    M->nddlv_dirichlet = 0;
	    { I i;
	      for (i=0;i<nddlu;++i)
		{

		  if ( 
		      (Z->ddlcod[i]==1) 
		      OR (Z->ddlcod[i]==2)
		      OR (Z->ddlcod[i]==3)
		      OR (Z->ddlcod[i]==4)
		      OR (Z->ddlcod[i]==5)
		      OR (Z->ddlcod[i]==12)
		      OR (Z->ddlcod[i]==34)
		      OR (Z->ddlcod[i]==45)
		      OR (Z->ddlcod[i]==23) )
		    {
		      M->ddlu_dirichlet[M->nddlu_dirichlet++] = i;
		    }

		  if ( 
		      (Z->ddlcod[i]==1) 
		      OR (Z->ddlcod[i]==2)
		      OR (Z->ddlcod[i]==3)
		      OR (Z->ddlcod[i]==4)
		      OR (Z->ddlcod[i]==5)
		      OR (Z->ddlcod[i]==51)
		      OR (Z->ddlcod[i]==12)
		      OR (Z->ddlcod[i]==34)
		      OR (Z->ddlcod[i]==45)
		      OR (Z->ddlcod[i]==23) )
		    {
		      M->ddlv_dirichlet[M->nddlv_dirichlet++] = i;
		    }

#if 0
		  if ( (Z->ddlcod[i]==1) OR (Z->ddlcod[i]==2)OR(Z->ddlcod[i]==3)OR(Z->ddlcod[i]==4)OR(Z->ddlcod[i]==12)OR(Z->ddlcod[i]==34)OR(Z->ddlcod[i]==23)OR(Z->ddlcod[i]==41) )
		    {
		      M->ddlv_dirichlet[M->nddlv_dirichlet++] = i;
		    }
#endif
		} } 
	    break;	      
	  }


	case __eStateSolver_dirichlet_pressure_symbolic:
	  {	      
	    /* 
	       DO NOTHING
	       BUT SPECIFY NO DOF IS SPECIFIED DIRICHLET 
	    */
	    M->nddlp_dirichlet = 0;

#if 0
	    { I i;
	      for (i=0;i<nddlp;++i)
		{
		  if ( (ddlcod[i]==1) )
		    {
		      ++M->nddlp_dirichlet;
		      M->ddlp_dirichlet[0] = i;
		      break;
		    }
		} }
#endif	      
	    break;
	  }
	    	    
	case __eStateSolver_dirichlet_velocity:
	  {
	    const eVariable variable_u = __eVariable_U;
	    pR vu = Variables_get(variables,&variable_u,&negal0);
	    pR vv = vu + Z->nddlu;

	    { I i;
	      for (i=0;i<M->nddlu_dirichlet;++i)
		{
		  const I k = M->ddlu_dirichlet[i];

		  if( ( (Z->ddlcod[k]==23)OR(Z->ddlcod[k]==2)OR(Z->ddlcod[k]==3))
		    OR
		      ( (Z->ddlcod[k]==34)OR(Z->ddlcod[k]==3)OR(Z->ddlcod[k]==4)) )
		    {
		      R uu[3];
		      boundaryu(&negal1,&dim,uu,&negal1,&dim,Z->ddlcoo+k*2,&negal1,NULL,&gTimeInfo->t,NULL);
		      /*		      printf("%e\n",uu[0]-*(hdl->u+k));*/
		      vu[k] =uu[0];


		      /*		      *(hdl->v+k)=uu[1];*/
#if 0
		      Z->u[k] = (0.5-Z->ddlcoo[k*2+1])*(Z->ddlcoo[k*2+1]+0.5)*regal4;
#endif
		      /*		      Z->u[k] = (regal1-Z->ddlcoo[k*2+1])*(Z->ddlcoo[k*2+1]+regal1)*regal4;*/


#if 0
		      if (Z->ddlcoo[k*2+1]<0.0)
			hdl->u[k] = -(0.5+Z->ddlcoo[k*2+1])*(Z->ddlcoo[k*2+1])*regal4*regal4;
		      else
			hdl->u[k] = (0.5-Z->ddlcoo[k*2+1])*(Z->ddlcoo[k*2+1])*regal4*regal4*2.0;
#endif
		    }
#if 0
		  else
		    {
		      *(hdl->u+k)=0.0;
		    }
#endif
		} }

	    { I i;
	      for (i=0;i<M->nddlv_dirichlet;++i)
		{
		  const I k = M->ddlv_dirichlet[i];
		  //		  hdl->v[k] = regal0;
		  vv[k] = regal0;
		} }      
	    break;
	  }

	case __eStateSolver_dirichlet_pressure:
	  {
#if 0
	    { I i;
	      for (i=0;i<nddlp_dirichlet;++i)
		{
		  const I k = ddlp_dirichlet[i];
		  p[k] = 1.0;
		} }
#endif
	    break;
	  }
	  
	case __eStateSolver_exit:
	  {
	    we_continue = __emnsNO;
	    break;
	  }
	case __eStateSolver_ERROR:
	  {
	    we_continue = __emnsNO;
#if 0
	    ns_errmsg("nsPBLM_canal:nsGLOBAL_run failed");
#endif
	    break;
	  }
	case __eStateSolver_ALL:
	  {
	    we_continue = __emnsNO;
#if 0
	    ns_errmsg("nsPBLM_canal:nsGLOBAL_run unexpected job");
#endif
	    break;
	  }
	}       	      
    }    
  free(rwork);
}



