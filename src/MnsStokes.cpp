

#include <signal.h>
#include <pthread.h>
#include "Cmdline.h"
#include "ns_config.h"


#include "ns_sys.h"
#include "ns_constantes.h"
#include "ns_config_lapack.h"
#include "Monitor.h"
#include <math.h>
#include <time.h>

#include "Global.hpp"
//#include "nsPBLM.h"
#include "Monitor.h"
#include "ns_sys.h"
#include "MetricTensor.h"
#if 0
static void ns_main_signal(int sigid) 
{
  

  
  fflush(stdout);
  fflush(stderr);
  fprintf(stderr,"\n");
  Monitor_errmsg(0,"ns_main_signal:Unexpected error"); 
  switch(sigid) {
    case SIGABRT:
      Monitor_errmsg(0,"ns_main_signal:<Abnormal stop>");  
      exit(1);
      break;
    case SIGFPE:
      Monitor_errmsg(0,"ns_main_signal:<Floating-point exception>"); 
      exit(1);
      break;
    case SIGILL:
      Monitor_errmsg(0,"ns_main_signal:<Illegal instruction>"); exit(1);break;
    case SIGSEGV:
      Monitor_errmsg(0,"ns_main_signal:<Segmentation fault>");  exit(1);break;
    case SIGTERM:
      {
	Monitor_errmsg(0,"<Program killed>2");  
	exit(1);
	break;
      }
    case SIGINT:
      {
	Monitor_errmsg(0,"<Program killed>");  
	exit(1);
	break;
      }
  default:
    {
      Monitor_errmsg(0,"  ??? unrecognized error");  
      exit(1);
      break;
    }
  }

}


static void ns_main_atexit()
{
#ifndef NDEBUG
  fprintf(stderr,"ns_main_atexit ...\n");
#endif

#ifndef NDEBUG
  fprintf(stderr,"ns_main_atexit done.\n");
#endif
}
#endif

void  Global_checkError(const char * 	filename_,
			const int 	fileline_,
			STR 		errmsg_,
			Err*		err_)
{
  if (err_[0])
    {
      Monitor_errmsg(0,"ERROR DETECTED");
      Monitor_errmsg(0,"ERROR CODE      : '%s' (err=" efmt ")",eErr_get_string(err_[0]),err_[0]);
      Monitor_errmsg(0,"ERROR MESSAGE   : '%s'",errmsg_);
      Monitor_errmsg(0,"ERROR FROM FILE : '%s'",filename_);
      Monitor_errmsg(0,"ERROR AT LINE   : %d",fileline_);
      exit(err_[0]);
    }
}

#define CHECK_ERROR(_errmsg,_err) Global_checkError(__FILE__,__LINE__,_errmsg,_err)

class ProblemDefinition
{
protected:
  ProblemDefinition()
  {
  };
  virtual ~ProblemDefinition()
  {
  };
public:
  virtual void use_integrator(eStateSolver job_[1],
			      Global * 	global_) = 0;
  
};


class SolverCanal : public ProblemDefinition
{
public:

  SolverCanal();
  virtual ~SolverCanal();
  void Run(Global * global_);
  
protected:

  static void boundaryf(cst_pI 		n_,
			cst_pI 		xm_,
			pR     		x_,
			cst_pI       	xoff_,
			cst_pI       	gl_m_,
			cst_pR       	gl_,
			cst_pI       	gloff_,
			cst_pI       	cod_,
			cst_pR       	t_,
			void * 		usrptr_);
  
  static void boundaryu	(cst_pI 	n_,
			 cst_pI 	xm_,
			 pR     	x_,
			 cst_pI       	xoff_,
			 cst_pI       	gl_m_,
			 cst_pR       	gl_,
			 cst_pI       	gloff_,
			 cst_pI       	cod_,
			 cst_pR       	t_,
			 void * 	usrptr_);

  static void ComputeU(const double 	x,
		       const double 	y,
		       const double 	z,
		       double&		u,
		       double&		v,
		       double&		w);
  
  void use_integrator	(eStateSolver job_[1],
			 Global* 	global_);      
};

SolverCanal::SolverCanal() {};
SolverCanal::~SolverCanal() {};

void SolverCanal::Run(Global * global_)
{
  eStateSolver job = __eStateSolver_ERROR;
  use_integrator(&job,
		 global_);
};
  
void SolverCanal::boundaryf(cst_pI 		n_,
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

      if (x_[i]>0.0)
	x_[i]=1.0;
      else if (x_[i]<0.0)
	x_[i]=0.0;
      else
	x_[i]=0.5;

#if 0
      ns_heaviside(global[0].params.iinfo[__ens_iinfo_heaviside],1,&x_[i],&ns_rv[__ens_rv_eps]);
#endif
      
    }
  };


 void SolverCanal::ComputeU(const double 	x,
		       const double 	y,
		       const double 	z,
		       double&		u,
		       double&		v,
		       double&		w)
  {	  
    if (y>0.5)
      {
	u = ((((R)0.5)-y)*(y-1.0)*((R)16.0));
      }
    else
      {
	u = ((((R)0.5)-y)*(y)*((R)16.0));
      }
	  
    v = ((R)0.0);
    w = ((R)0.0);
  };
  
   void SolverCanal::boundaryu(cst_pI 		n_,
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

    const I n = n_[0];

    { I i;
      for (i=0;i<n;++i)
	{
	  //	  const R x = gl_[i];
	  const R y = gl_[gloff_[0]+i];
//	  const R z = 0.0;
//
//	  R u,v,w;
//	  ComputeU(x,
//		   y,
//		   z,
//		   u,
//		   v,
//		   w);
//
//	  x_[i] = u;
//	  x_[xoff_[0]+i] = v;
//	  
	  if (y>0.5)
	    {
	      x_[i] = ((((R)0.5)-y)*(y-1.0)*((R)16.0))*2.0;
	    }
	  else
	    {
	      x_[i] = ((((R)0.5)-y)*(y)*((R)16.0));
	    }
	  	  
	  x_[xoff_[0]+i] = ((R)0.0);
	  
	} }
    
  };
  
  void SolverCanal::use_integrator(eStateSolver job_[1],
		      Global* 	global_)
  {
    static const I dim = ((I)_dim);
    const ns_mesh*	mesh 		= global_->GetMesh();    
    Variables* 		variables 	= global_->GetVariables();
    Parameters* 	parameters	= global_->GetParameters();
    
    pTimeReadOnly 	gTimeInfo	= global_->GetTime();       
    pSparseStokes 	M 		= global_->GetSparseStokes();
    if (!M)
      {
	std::cout << "Sparse stokes matrix does not exist" << std::endl;
	exit(1);
      }
    const I 		rwork_n		= global_->GetRequiredRworkN();
    const I 		iwork_n		= global_->GetRequiredIworkN();
    pR 			rwork		= (pR)malloc((rwork_n+1)*sizeof(R));
    pI 			iwork		= (pI)malloc((iwork_n+1)*sizeof(I));

    //
    // Do not do transient.
    //
    parameters->SetKindTransientMethod(KindEquation::VELOCITY,
				       KindTransientMethod::UNDEFINED);
    parameters->SetKindTransientMethod(KindEquation::TRANSPORT,
				       KindTransientMethod::UNDEFINED);

    const I nelm  = ns_mesh_nelm(mesh);
    const I nddlu = global_->GetNumDofs(__eVariable_U);
    
#if 0
    const R 	current_t 		= Z->ruser[ruser_t];
    const I 	nddlp 			= Z->nddlp;
#endif

    job_[0] = __eStateSolver_ERROR;
    
  run:
    global_->Run(job_,
		 &rwork_n,
		 rwork,
		 &iwork_n,
		 iwork);
    
#ifndef NDEBUG      
    printf("---> job_[0] = " efmt "\n",job_[0]);
#endif
    
    switch(job_[0])
      {
	
      case __eStateSolver_slip_velocity_symbolic:
	{
	  M->nddls_dirichlet = 0;
	  goto run;
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
#ifdef TONTON
			  if (Z->ddlcod[Z->ddlcnc[jelm*6+3+jadj]-1]==23)
			    {

			      
			      
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


			    }
#endif

			  
			}
		    } }
	      } }
	  goto run;
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
	  goto run;
	}
	
      case __eStateSolver_initial_condition:
	{
	  /* DO NOTHING */
	  goto run;
	}
	    
      case __eStateSolver_dirichlet_marker:
	{
	  /* DO NOTHING */
	  goto run;
	}
	
      case __eStateSolver_dirichlet_marker_symbolic:
	{
	  /* DO NOTHING */
	  goto run;
	}	    
	
      case __eStateSolver_dirichlet_velocity_symbolic:
	{
	  
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_velocity_symbolic ... ");
	      } }
	    

	  M->nddlu_dirichlet = 0;
	  M->nddlv_dirichlet = 0;
	  { I i;
	    for (i=0;i<nddlu;++i)
	      {
		const I ddlcod = global_->GetDofTopology(i);
		if ( 
		    (ddlcod==1) 
		    || (ddlcod==2)
		    || (ddlcod==3)
		    || (ddlcod==4)
		    || (ddlcod==5)
		    || (ddlcod==12)
		    || (ddlcod==34)
		    || (ddlcod==45)
		    || (ddlcod==23) )
		  {
		    //		    std::cout << "nddlu[" << M->nddlu_dirichlet << "] = " << i << std::endl;
		    M->ddlu_dirichlet[M->nddlu_dirichlet++] = i;
		  }
		
		if ( 
		    (ddlcod==1) 
		    || (ddlcod==2)
		    || (ddlcod==3)
		    || (ddlcod==4)
		    || (ddlcod==5)
		    || (ddlcod==51)
		    || (ddlcod==12)
		    || (ddlcod==34)
		    || (ddlcod==45)
		    || (ddlcod==23) )
		  {
		    //		    std::cout << "nddlv[" << M->nddlv_dirichlet << "] = " << i << std::endl;
		    M->ddlv_dirichlet[M->nddlv_dirichlet++] = i;
		  }
		
#if 0
		if ( (ddlcod==1) OR (ddlcod==2)OR(ddlcod==3)OR(ddlcod==4)OR(ddlcod==12)OR(ddlcod==34)OR(ddlcod==23)OR(ddlcod==41) )
		  {
			M->ddlv_dirichlet[M->nddlv_dirichlet++] = i;
		  }
#endif
		
	      } }
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_velocity_symbolic done. ");
	      } }

	  goto run;
	}

	
      case __eStateSolver_dirichlet_pressure_symbolic:
	{	      
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_pressure_symbolic ... ");
	      } }
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
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_pressure_symbolic done. ");
	      } }
	  goto run;
	}
	
      case __eStateSolver_dirichlet_velocity:
	{
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_velocity ... ");
	      } }
	  const eVariable variable_u = __eVariable_U;
	  pR vu = variables->Get(variable_u,0);
	  pR vv = vu + nddlu;

	  cst_pR ddlcoo = global_->GetDofCoordinates();
	  
	  { I i;
	    for (i=0;i<M->nddlu_dirichlet;++i)
	      {
		const I k = M->ddlu_dirichlet[i];
		const I ddlcod = global_->GetDofTopology(k);
		if( ( (ddlcod==23)OR(ddlcod==2)OR(ddlcod==3))
		    OR
		    ( (ddlcod==34)OR(ddlcod==3)OR(ddlcod==4)) )
		  {
		    
		    R uu[3];		    
		    boundaryu(&negal1,&dim,uu,&negal1,&dim,ddlcoo+k*2,&negal1,NULL,&gTimeInfo->t,NULL);
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

	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_velocity done. ");
	      } }

	  
	  goto run;
	}
	
      case __eStateSolver_dirichlet_pressure:
	{
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_pressure ... ");
	      } }
#if 0
	  { I i;
	    for (i=0;i<nddlp_dirichlet;++i)
	      {
		const I k = ddlp_dirichlet[i];
		p[k] = 1.0;
	      } }
#endif
	  { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	    if (verbose)
	      {
		Monitor_msg(0,"StateSolver dirichlet_pressure done. ");
	      } }
	  goto run;
	}
	
	
	
      case __eStateSolver_exit:
	{
	  break;
	}
	
      case __eStateSolver_ERROR:
	{
#if 0
	  ns_errmsg("nsPBLM_canal:nsGLOBAL_run failed");
#endif
	  break;
	}

	
	
      case __eStateSolver_ALL:
	{
#if 0
	  ns_errmsg("nsPBLM_canal:nsGLOBAL_run unexpected job");
#endif
	  break;
	}
	
      }       	      
    
    free(rwork);
  };
  







class Solver
{
protected:

  Cmdline 		m_cmdline;
  Parameters* 		m_parameters {};
  Global* 		m_solver {};
  ns_mesh 		m_mesh;
  
public:
  Solver(int 		argc,
	 const char**	argv,
	 STR 		errmsg_,
	 Err*		err_)
  {
    //
    // Initialization of the command line
    //
    Cmdline_def (&m_cmdline,
		 argc,
		 argv);    
    
    I nproc = 1;
    if (false == Cmdline_get_integer(&m_cmdline,
				     "-n",
				     &nproc))
      {
	Monitor_warn(0,"set nproc to 1");
	nproc=1;
      }

    std::cout << "MNSSTOKES::Initialization " << std::endl;
    //
    // MAIN
    //
    {
      this->m_solver = new Global(&m_cmdline,
				  errmsg_,
				  err_);
      
      CHECK_ERROR(errmsg_,err_);
      
    }
    std::cout << "MNSSTOKES::Initialization done. " << std::endl;
  
    this->m_parameters = this->m_solver->GetParameters();

    
    const bool verbose = m_parameters->GetInfoLogical(InfoLogical::verbose);
    //
    // READ MESH
    //
    {
      if (verbose)
	{
	  Monitor_msg(0,"ReadMesh ..?");
	}  
      
      ns_mesh_read(&this->m_mesh,
		   errmsg_,
		   err_,
		   "%s.p%.2lld.mesh",
		   argv[1],
		   0);
      
      CHECK_ERROR(errmsg_,
		  err_);
    
      if (verbose)
	{
	  Monitor_msg(0,"ReadMesh Success.");
	}  
      
#if 0

      pMetricTensor metric = ns_mesh_addmetric(&m_mesh);


      
      MetricTensor_geometric_statistical(metric);

      {
	I i;
	for (i=0;i<metric->n;++i)
	  {

	    metric->x[i] *=100.0;
	  }
      }

      
      ns_mesh_write_medit(&m_mesh,
			  "allo");

      MetricTensor_write_medit(metric,
			       "allo");

      
#endif
      

    }

  };

  virtual ~Solver()
  {
    if (this->m_solver)
      {
	delete this->m_solver;
	this->m_solver = nullptr;
      }
  };

  void precompute(Err*		err_)
  {
    STR errmsg;
    const bool verbose = this->m_parameters->GetInfoLogical(InfoLogical::verbose);
    if (verbose)
      {
	Monitor_msg(0,"precompute ...");
      }
    
    this->m_solver->Precompute(&this->m_mesh,
			       errmsg,
			       err_);
    
    CHECK_ERROR(errmsg,
		err_);
    
    if (verbose)
      {
	Monitor_msg(0,"precompute success.");
      }
  };

  void compute(Err*		err_)
  {

    //
    // COMPUTE
    //
    const bool verbose = m_parameters->GetInfoLogical(InfoLogical::verbose);
    
    if (verbose)
      {
	Monitor_msg(0,"Run ..?");
      }    
    
    SolverCanal solverCanal;
    solverCanal.Run(this->m_solver);
    
    if (verbose)
      {
	Monitor_msg(0,"Run Success.");
      }
    
  };
  
  
};



int main(int 		argc,
	 const char**	argv)
{

  Monitor_def(0,
	      argv[0],	      
	      MonitorMode_STD);

  STR errmsg;
  Err err;
  Solver solver(argc,
		argv,
		errmsg,
		&err);
    
  solver.precompute	(&err);
  if (err)
    {
      std::cerr << "error precompute " << eErr_get_string(err) << std::endl;
      return err;
    }

  solver.compute	(&err);
  if (err)
    {
      std::cerr << "error compute " << eErr_get_string(err) << std::endl;
      return err;
    }
    
  return err;

}
