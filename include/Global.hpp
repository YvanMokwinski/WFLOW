
#pragma once

#include "nsGLOBAL.h"
#include "Parameters.hpp"
#include "Monitor.h"
#include "NonLinearSolverNewton.hpp"
#include "StokesNonLinearSolverDriver.hpp"
#include "LinsysVectors.hpp"
#include "Variables.hpp"

#include "LinearSolver/Direct/MKL/Pardiso.hpp"
#include "LinearSolver/Iterative/MKL/Fgmres.hpp"
#include "LinearSolver/Preconditioner/Jacobi.hpp"
#include "LinearSolver/Preconditioner/MKL/Ilu0.hpp"
#include "LinearSolver/Preconditioner/MKL/Ilut.hpp"

class MyLinearOperator : public LinearSolver::ILinearOperator
{
  SparseStokes*m_a;

public: MyLinearOperator(SparseStokes*a)
  {
    m_a = a;
  };
  
public: void Apply(const char * transpose_,
		   pR 		y_,
		   cst_pR 	x_)
  {
    R r1 = 1.0;
    R r0 = 0.0;    
    Sparse_gemv	(&r1,
		 m_a->A,
		 x_,
		 &r0,
		 y_);    
  };
  
};






class Global
{
  
private:  pGlobal m_self;  
protected: I iproc;
protected: SolverConstants* 			m_solverConstants;
protected: Variables*				m_variables;
protected: NonLinearSolverNewton* 		m_nonLinearSolverNewton;
protected: LinsysVectors* 			m_linsysVectors;
protected: Parameters*				m_parameters;
protected: LinearSolver::IInverseOperator*	m_linsysA;
protected: LinearSolver::IInverseOperator*	m_linsysF;
protected: LinearSolver::IInverseOperator*	m_linsysB;
protected: LinearSolver::IInverseOperator*	m_linsysS;
protected: LinearSolver::IInverseOperator*	m_linsysC;
protected: StokesNonLinearSolverDriver* 	m_stokesNonLinearSolverDriver;
  
public: Err Global_default();
  
private:  void Global_main(pGlobal const 	self_,
			   cmdline* 		cmdline_,
			   STR 			errmsg_,
			   Err*			err_);
  
public:  Global(cmdline* 	 	cmdline_,
		STR 			errmsg_,
		Err*			err_);

public:  virtual ~Global();
  
public: inline const Parameters* 	GetParameters() const noexcept;
public: inline Parameters* 		GetParameters() noexcept;

public: inline const Variables* 	GetVariables() const noexcept;
public: inline Variables*       	GetVariables() noexcept;
  
public: inline pTimeReadOnly  		GetTime() const noexcept;
public: inline pTime 			GetTime() noexcept;
  
public: inline const ns_mesh * 		GetMesh() const noexcept;
public: inline pSparseStokes 		GetSparseStokes() noexcept;
public: inline cst_pR 			GetDofCoordinates() const noexcept;
public: inline I 			GetDofTopology(const I index_) const noexcept;

public: inline I 			GetRequiredRworkN() const noexcept;
public: inline I 			GetRequiredIworkN() const noexcept;

public: inline I 			GetNumDofs(const eVariable evariable_) const noexcept;

public: void 				Global_initializePart1	(pGlobal 	const   self_,
								 ns_mesh *	const  	mesh_,
								 STR 			errmsg_,
								 Err*			err_);
  
public: inline void 			Precompute(ns_mesh * 	const 	mesh_,
						   STR 			errmsg_,
						   Err*			err_);
  
public: inline void SolveTimeStep();

public: inline void NextTimeInterval();
  
public: inline void Run(eStateSolver 	job_[1],
			cst_pI 		rwork_n_,
			pR 			rwork_,
			cst_pI 		iwork_n_,
			pI 			iwork_);
  
public: inline void PrintVariables();
public: inline void UpdatePreviousX();

};

Err Global::Global_default()
{
  Parameters* parameters = this->GetParameters();
    
  parameters->SetKindTensionMethod	(KindTensionMethod::CSF);
  parameters->SetKindTransientMethod	(KindEquation::TRANSPORT,
					 KindTransientMethod::EULER);
  parameters->SetKindTransientMethod	(KindEquation::VELOCITY,
					 KindTransientMethod::EULER);
  parameters->SetInfoLogical		(InfoLogical::color,
					 true);

  parameters->SetInfoLogical		(InfoLogical::transport,
					 true);
    
  parameters->SetInfoReal		(InfoReal::froude,
					 1e+15);
    
  parameters->SetKindLinearSolver	(KindEquation::VELOCITY,
					 KindLinearSolver::Direct);

  parameters->SetInfoString		(InfoString::pblmname,
					 "canal");
    
  return __eErr_no;
};

  
 void Global::Global_main(pGlobal const 	self_,
			  cmdline* 	cmdline_,
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
  
    /*
      0/ READ CONFIG FILE ? 
    */  
    { STR ctmp;
      if ( (cmdline_->get_string("-c",ctmp)) )
	{
#ifndef NDEBUG
	  Monitor_msg(0,"ns_global_main:mkmake_ns_read_configfile(%s) ...\n",ctmp);
#endif
	  have_configfile = __emnsNO;
	  if (NOT have_configfile)
	    {
	      Monitor_errmsg(this->iproc,"ns_global_main:no config file '%s'",ctmp);
	      err_[0] = __eErr_user;
	      return;
	    }
	  if (err_[0])
	    {
	      Monitor_errmsg(this->iproc,"ns_global_main:mkmake_ns_read_configfile('%s') failed",ctmp);
	      return;
	    }
	} }

    //
    // 1/ OVERWRITE CONFIG WITH COMMAND LINE 
    //
    Parameters* parameters= new Parameters(cmdline_,
					   have_configfile);
    this->m_parameters = parameters;

    self_->parameters = (void*)parameters;

  


    //
    // 2/ CHECK INVALID ARGUMENT FROM COMMAND LINE 
    //
    err_[0] = (cmdline_->check_invalid()) ? __eErr_user : __eErr_no;
    if (err_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_main:Cmdline_check_invalid failed");
	return;
      }
    const bool cmdline_isempty = cmdline_->isempty();

    const I restart = parameters->GetInfoInteger(InfoInteger::restart);
    //
    // 3/ RESTART MODE ? 
    //
    if ( (0 == restart)
	 AND
	 (cmdline_isempty) )
      { 
	Monitor_errmsg(this->iproc,"missing mesh file");
	err_[0]= __eErr_user;
	return;
      } 
    else if (restart)
      {
	pTime const gTimeInfo = Global_get_Time(self_);
	char ctmp[512];
	sprintf(ctmp,"integrator.u.%.5" nsFORMAT_INTEGER ".timeinfo.txt",restart);
	FILE*fich = fopen(ctmp,"r");
	int ii = fscanf(fich,"" ifmt " " ifmt " %lf",&gTimeInfo->itime_interval,&gTimeInfo->itime,&gTimeInfo->ti);
	if (ii==4321)
	  {}
	fclose(fich);
	sprintf(ctmp,"integrator.P1.%.5" nsFORMAT_INTEGER ".mesh",gTimeInfo->itime_interval);
	extern char * strdup(const char*);
	self_->filename = strdup(ctmp);
	ns_basename(self_->filename,self_->basename);
      }
    else
      {
	self_->filename = cmdline_->get_arg(1);
	ns_basename(self_->filename,self_->basename);
      }
  


    //
    // 4/ DEFAULT
    //
    err_[0] = Global_default();
    if (err_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_main:nsGLOBAL_default failed");
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
					    (eTransientMethod)parameters->GetKindTransientMethod((KindEquation::EnumType)ieq));
	      if (err_[0])
		{
		  Monitor_errmsg(this->iproc,"nsGLOBAL_main:TransientScheme_def failed on equation '%s' method : %s",
				 eKindEquation_get_string(ieq),
				 eTransientMethod_get_string((eTransientMethod)parameters->GetKindTransientMethod((KindEquation::EnumType)ieq)));
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
    self_->rv[__ens_rv_iweber]		 		= ((R)1.0)/parameters->GetInfoReal(InfoReal::weber);
    self_->rv[__ens_rv_ifroude]	 			= ((R)1.0)/parameters->GetInfoReal(InfoReal::froude);
    self_->rv[__ens_rv_ireynold] 				= ((R)1.0)/parameters->GetInfoReal(InfoReal::reynold);
    self_->rv[__ens_rv_coeff_ratio_viscosity] 		= parameters->GetInfoReal(InfoReal::ratio_viscosity);
    self_->rv[__ens_rv_coeff_ratio_density]   		= parameters->GetInfoReal(InfoReal::ratio_density);
    self_->rv[__ens_rv_dt] 				= parameters->GetInfoReal(InfoReal::dtmin);
    self_->rv[__ens_rv_idt]				= ((R)1.0)/self_->rv[__ens_rv_dt];
    self_->rv[__ens_rv_midt]				= -self_->rv[__ens_rv_idt];
    self_->rv[__ens_rv_eps] 				= parameters->GetInfoReal(InfoReal::epsmin);
    self_->rv[__ens_rv_ieps] 				= ((R)1.0)/self_->rv[__ens_rv_eps];
    self_->rv[__ens_rv_vmin] 				= ((R)1.0)/(parameters->GetInfoReal(InfoReal::hmax)*parameters->GetInfoReal(InfoReal::hmax));
    self_->rv[__ens_rv_vmax] 				= ((R)1.0)/(parameters->GetInfoReal(InfoReal::hmin)*parameters->GetInfoReal(InfoReal::hmin));
  

    //
    // 7/ VERIFICATION
    //
    if (NOT parameters->GetInfoLogical(InfoLogical::skip_verif))
      {
	Monitor_msg(0,"");
	Monitor_msg(0," NS CONFIGURATION");
	Monitor_msg(0,"");
	Monitor_msg(0," general parameters");
	Monitor_msg(0,"        verbose    : %s",(parameters->GetInfoLogical(InfoLogical::verbose))?"yes":"no");
	Monitor_msg(0,"        verify     : %s",(parameters->GetInfoLogical(InfoLogical::skip_verif))?"yes":"no");
	Monitor_msg(0,"        DG         : %s",(parameters->GetInfoLogical(InfoLogical::vnsdg))?"yes":"no");
	Monitor_msg(0,"        dynamic bc : %s",(parameters->GetInfoLogical(InfoLogical::dynamic_boundary_condition))?"yes":"no");
	Monitor_msg(0,"        axix       : %s",(parameters->GetInfoLogical(InfoLogical::axisymetric_x))?"yes":"no");
	Monitor_msg(0,"        axiy       : %s",(parameters->GetInfoLogical(InfoLogical::axisymetric_y))?"yes":"no");
	Monitor_msg(0,"        nobc  vcod : %s",(parameters->GetInfoInteger(InfoInteger::noboundary_vcod))?"yes":"no");
	Monitor_msg(0," discrete formulation");
	Monitor_msg(0,"    velocity");
	Monitor_msg(0,"        method     : %s",(parameters->GetInfoLogical(InfoLogical::usupg))?"supg":"galerkin");
	Monitor_msg(0,"        nbasis     : " ifmt "",(I)self_->nddlelmu);
	Monitor_msg(0,"        slip       : %s",(parameters->GetInfoLogical(InfoLogical::slip))?"yes":"no");
	Monitor_msg(0,"    pressure");
	Monitor_msg(0,"        uncoupled  : %s",(parameters->GetInfoLogical(InfoLogical::pressure_uncoupled))?"yes":"no");
	Monitor_msg(0,"        freematrix : %s",(parameters->GetInfoLogical(InfoLogical::pressure_freematrix))?"yes":"no");
	Monitor_msg(0,"        method     : %s",(parameters->GetInfoLogical(InfoLogical::pspg))?"pspg":"galerkin");
	Monitor_msg(0,"        nbasis     : " ifmt "",(I)self_->nddlelmp);
	Monitor_msg(0,"    transport : %s",(parameters->GetInfoLogical(InfoLogical::transport))?"yes":"no");
	Monitor_msg(0,"        eps        : %.3e",parameters->GetInfoReal(InfoReal::epsmin));
	Monitor_msg(0,"        marker     : %s",(parameters->GetInfoLogical(InfoLogical::color))?"pseudo-concentration":"signed-distance");
	Monitor_msg(0,"        uncoupled  : %s",(parameters->GetInfoLogical(InfoLogical::transport_uncoupled))?"yes":"no");
	Monitor_msg(0,"        method     : %s",(parameters->GetInfoLogical(InfoLogical::transport_galerkin))?"supg":"dg");
	Monitor_msg(0,"        nbasis     : " ifmt "",(I)self_->nddlelmp);
	Monitor_msg(0,"        nbasis dg  : " ifmt "",(I)self_->nddlelmdg);
	Monitor_msg(0,"        redistance : %s",(parameters->GetInfoLogical(InfoLogical::redistance))?"yes":"no");
	Monitor_msg(0,"    tension : %s",(parameters->GetInfoLogical(InfoLogical::tension))?"yes":"no");
	Monitor_msg(0,"        method     : %s",eTensionMethod_get_string((eTensionMethod)parameters->GetKindTensionMethod()));
	Monitor_msg(0,"    time : %s",(parameters->GetInfoLogical(InfoLogical::time_adaptivity))?"yes":"no");
	Monitor_msg(0,"        restart        : " ifmt "",parameters->GetInfoInteger(InfoInteger::restart));
	Monitor_msg(0,"        ntime          : " ifmt "",parameters->GetInfoInteger(InfoInteger::ntime));
	Monitor_msg(0,"        ntime_interval : " ifmt "",parameters->GetInfoInteger(InfoInteger::ntime_interval));
	Monitor_msg(0,"        velocity   : %s",eTransientMethod_get_string((eTransientMethod)parameters->GetKindTransientMethod(KindEquation::VELOCITY)));
	Monitor_msg(0,"        transport  : %s",eTransientMethod_get_string((eTransientMethod)parameters->GetKindTransientMethod(KindEquation::TRANSPORT)));
	Monitor_msg(0," discrete solution");
	Monitor_msg(0,"        maxiter    : " ifmt "",parameters->GetInfoInteger(InfoInteger::newton_maxiter));
	Monitor_msg(0,"        residu tol : %.3e",parameters->GetInfoReal(InfoReal::newton_tol_residu));
	Monitor_msg(0,"        correc tol : %.3e",parameters->GetInfoReal(InfoReal::newton_tol_correc));
	Monitor_msg(0,"        linsol     : %s",eLinearSolver_get_string((eLinearSolver)parameters->GetKindLinearSolver(KindEquation::VELOCITY)));
	Monitor_msg(0," mesh adaptivity : %s",(parameters->GetInfoLogical(InfoLogical::mesh_adaptivity))?"yes":"no");
	Monitor_msg(0,"        maxiter    : " ifmt "",parameters->GetInfoInteger(InfoInteger::mesh_adaptivity_maxiter));
	Monitor_msg(0,"        hmin       : %.3e",parameters->GetInfoReal(InfoReal::hmin));
	Monitor_msg(0,"        hmax       : %.3e",parameters->GetInfoReal(InfoReal::hmax));
	Monitor_msg(0," adimensional numbers");
	Monitor_msg(0,"        viscosity ratio  : %.3e",parameters->GetInfoReal(InfoReal::ratio_viscosity));
	Monitor_msg(0,"        density   ratio  : %.3e",parameters->GetInfoReal(InfoReal::ratio_density));
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
    return;
  };
  
Global::Global(cmdline*		cmdline_,
	 STR 			errmsg_,
	 Err*			err_)
  {
    this->m_solverConstants = new SolverConstants();
    this->m_self = (pGlobal)calloc(1,sizeof(nsGLOBAL_ST));
    Global_main	(this->m_self,
		 cmdline_,
		 errmsg_,
		 err_);    
  };

Global::~Global()
  {
    if (this->m_linsysVectors)
      {
	delete this->m_linsysVectors;
	this->m_linsysVectors = NULL;
      }
    
    Global_free(this->m_self);
    free(this->m_self);
    this->m_self = NULL;
  };
  
  inline const Parameters* Global::GetParameters() const noexcept { return this->m_parameters; };
  inline Parameters* Global::GetParameters() noexcept { return this->m_parameters; };

  inline const Variables* Global::GetVariables() const noexcept  { return this->m_variables; };  
  inline Variables*       Global::GetVariables() noexcept { return this->m_variables; };
  
  inline pTimeReadOnly  Global::GetTime() const noexcept
  {
    return GlobalReadOnly_get_Time(this->m_self);
  };

  inline pTime Global::GetTime() noexcept
  {
    return Global_get_Time(this->m_self);
  };
  
  inline const ns_mesh * Global::GetMesh() const noexcept
  {
    return (const ns_mesh*)this->m_self->mesh_usrptr;
  };

  inline pSparseStokes Global::GetSparseStokes() noexcept
  {
    return this->m_self->sparseStokes;
  };

  inline cst_pR Global::GetDofCoordinates() const noexcept
  {
    return this->m_self->ddlcoo;
  };

  inline I Global::GetDofTopology(const I index_) const noexcept
  {
    return this->m_self->ddlcod[index_];
  };

  inline I Global::GetRequiredRworkN() const noexcept
  {
    return this->m_self->required_rwork_n;
  };
  
  inline I Global::GetRequiredIworkN() const noexcept
  {
    return this->m_self->required_iwork_n;
  };

  inline I Global::GetNumDofs(const eVariable evariable_) const noexcept
  {
    switch(evariable_)
      {
      case __eVariable_U:
	{
	  return this->m_self->nddlu;
	}
      case __eVariable_P:
	{
	  return this->m_self->nddlp;
	}
      case __eVariable_F:
      case __eVariable_ERROR:
      case __eVariable_ALL:
	{
	  break;
	}
      }
    fprintf(stderr,"Global::GetNumDofs failed\n");
    exit(1);
    return 0;
  };

  void Global::Global_initializePart1		(pGlobal 	const   self_,
					 ns_mesh *	const  	mesh_,
					 STR 			errmsg_,
					 Err*			err_)
  {
    
    err_[0] = __eErr_no;
    
    const Parameters *parameters 	= this->GetParameters();
    const bool  slip 			= parameters->GetInfoLogical(InfoLogical::slip);  
    const bool  transport 		= parameters->GetInfoLogical(InfoLogical::transport);
    const bool  transport_uncoupled 	= parameters->GetInfoLogical(InfoLogical::transport_uncoupled);
    const bool  include_dg 		= parameters->GetInfoLogical(InfoLogical::vnsdg);  
    const bool  pressure_uncoupled 	= parameters->GetInfoLogical(InfoLogical::pressure_uncoupled);
    const bool  verbose			= parameters->GetInfoLogical(InfoLogical::verbose);  
    const bool  pressure_freematrix  	= parameters->GetInfoLogical(InfoLogical::pressure_freematrix);
    
    cst_eFace element  = ns_mesh_elm(mesh_);
    
    self_->mesh_usrptr 		= mesh_;
    self_->nelm 		= ns_mesh_nelm(mesh_);
    self_->nvertex		= ns_mesh_get_numNodes(mesh_);
    self_->nedge_boundary	= ns_mesh_nedge_boundary(mesh_);
    self_->nedge		= ns_mesh_get_numEdges(mesh_);
    
    memset		(&self_->var_u,0,sizeof(ns_var));
    memset		(&self_->var_p,0,sizeof(ns_var));  
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize ...");
#endif
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:add space u ...");
#endif
    
    self_->space_u = ns_mesh_addspace(mesh_,_shape_u);
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:add spaces p ...");
#endif
    
    self_->space_p = ns_mesh_addspace(mesh_,_shape_p);
    if (include_dg)
      {
#ifndef NDEBUG
	Monitor_msg(this->iproc,"nsGLOBAL_initialize:add spaces dg ...");
#endif
	self_->space_dg = ns_mesh_addspace(mesh_,_shape_dg);
      }  
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:compute dof ...");
#endif
    
    self_->nddlu		= ns_mesh_nddlspace(mesh_,_shape_u);
    self_->nddlp		= ns_mesh_nddlspace(mesh_,_shape_p);
    self_->dec_ddlu		= 0;
    self_->dec_ddlp		= self_->nddlu*_dim;
    self_->total_nddl 		= self_->nddlu*_dim+self_->nddlp;      
    self_->blank_n		= self_->total_nddl;
    self_->iwork_n		= self_->total_nddl;    
    self_->ddlcoo		= mesh_->coo;
    self_->ddlcnc		= mesh_->cnc;
    self_->ddlcod		= mesh_->cod;
    self_->jacelm     		= mesh_->jacelm;
    self_->required_rwork_n	= MAX(self_->nddlp,self_->nddlu*_dim);
    if (verbose)
      {
	fprintf(stdout,"nedge boundary " ifmt "\n",self_->nedge_boundary);
	fprintf(stdout,"nedge interior " ifmt "\n",self_->nedge-self_->nedge_boundary);
	fprintf(stdout,"nddlelmu       " ifmt "\n",self_->nddlelmu);
	fprintf(stdout,"nddlelmp       " ifmt "\n",self_->nddlelmp);
	fprintf(stdout,"nddlu          " ifmt "\n",self_->nddlu);
	fprintf(stdout,"nddlp          " ifmt "\n",self_->nddlp);
	fprintf(stdout,"nddl           " ifmt "\n",self_->total_nddl);      
      }
    
    self_->blank = (pI)calloc(self_->total_nddl+2*self_->nedge_boundary,sizeof(I));
    self_->iwork = (pI)malloc(self_->total_nddl*sizeof(I));
    
#if 0
    I kind = 0;
    if (transport)
      {
	if (transport_uncoupled)
	  {
	    kind = 1;
	    if (pressure_uncoupled)
	      {
		kind=2;
	      }
	  }  
      }
    else
      {
	kind = 1;
	if (pressure_uncoupled)
	  {
	    kind=2;
	  }
      }
#endif
    
    self_->patch 		= mesh_->patch;
    self_->bpatch		= mesh_->bpatch;
    self_->adj   		= mesh_->adj;
    self_->jacedge       	= mesh_->jacedge;
    self_->normaledge    	= mesh_->normaledge;
    self_->trelm         	= mesh_->trelm;

    //
    // Create the sparse matrix.
    //
    self_->sparseStokes = (SparseStokes*)calloc(1,sizeof(SparseStokes));
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:sparse_A_init ...");
#endif
    sparse_A_init(self_,__eKindSystem_up);    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:sparse_A_init done.");
#endif
    
    Sparse_fortran_indexation(SparseStokes_get_A(self_->sparseStokes));
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_UU ...");
#endif
    
    { const I nddluxdim 		= self_->nddlu*_dim;
      const I start_i_block_UU 	= 0;
      const I start_j_block_UU 	= 0;
      SparseBlock_def	(SparseStokes_get_UU(self_->sparseStokes),
			 SparseStokes_get_A(self_->sparseStokes),
			 &start_i_block_UU,
			 &nddluxdim,
			 &start_j_block_UU,
			 &nddluxdim);  }
    /**/
    if (pressure_uncoupled)
      {
	if (NOT pressure_freematrix)
	  {
	    if (verbose)
	      Monitor_msg(this->iproc,"nsGLOBAL_initialize:sparse_B_init ...");
	    
	    sparse_B_init(self_);
	    
	    Sparse_fortran_indexation(self_->sparseB);
#if 0
	    sparse_C_init(iproc_);
	    nsSPARSE_fortran_indexation(&self_->sparse_C);
#endif
	    
	  }
      }
    else
      {
#ifndef NDEBUG
	Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_UB ... " ifmt "",this->iproc);
#endif
	{
	  const I nddluxdim		= self_->nddlu*_dim;
	  const I start_i_block_BU 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlp : nddluxdim;
	  const I start_j_block_BU 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlu : 0;
	  const I start_i_block_UB 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlu : 0;
	  const I start_j_block_UB 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlp : nddluxdim;
	  const I start_i_block_BB 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlp : nddluxdim;
	  const I start_j_block_BB 	= ( (transport) AND (NOT transport_uncoupled) ) ? self_->dec_ddlp : nddluxdim;      
	  
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_BU ...");
#endif
	  
	  SparseBlock_def	(SparseStokes_get_BU(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_i_block_BU,
				 &self_->nddlp,
				 &start_j_block_BU,
				 &nddluxdim);            
	  
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_UB ...");
#endif
	  SparseBlock_def	(SparseStokes_get_UB(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_i_block_UB,
				 &nddluxdim,
				 &start_j_block_UB,
				 &self_->nddlp);
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_BB ...");
#endif
	  SparseBlock_def	(SparseStokes_get_BB(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_i_block_BB,
				 &self_->nddlp,
				 &start_j_block_BB,
				 &self_->nddlp);      
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"nsGLOBAL_initialize:block_BB done.");
#endif
	  
	  
	}
      }
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:decalage ...");
#endif
    
#if 0
    self_->metric			= ns_mesh_addmetric(mesh_);
    self_->metric_end_time_interval	= ns_mesh_addmetric(mesh_);  
    /* 
       CALCUL DE LAMETRIQUE GEOMETRIQUE (DOIT ETRE FAIT AVANT LE DECALAGE DU PATCH
    */
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:compute statistical metric ...");
#endif
    nsMETRIC_geometric_statistical(self_->metric);    
#endif
    /*    
	  ON DECALAE PATCH POUR QU'IL PUOISSE ETRE UTILISE PAR MK
	  DONC APRES SI ON UTILISE PLUS LE PATCH DANS 'NS'
    */    
    if (transport)
      {
#if 0
	nsDG_init();
	self_->transport_usrptr = nsDG_transport_new(mesh_);
#endif
    }
    
    /**/
    { I i;
    for (i=0;i<self_->nelm*self_->nddlelmu;++i)
      {
	mesh_->own_patch[i]+=1;
      } }

    if (slip)
      {
	{	const I start_i_block_SS 	= self_->total_nddl;
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"define slip block " ifmt " " ifmt "\n",start_i_block_SS,start_i_block_SS);
#endif


	  SparseBlock_def	(SparseStokes_get_SS(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_i_block_SS,
				 &self_->nddl_slip,
				 &start_i_block_SS,
				 &self_->nddl_slip); }      

	{ const I nddluxdim 		= self_->nddlu*_dim;
	  const I start_i_block_SU 	= self_->total_nddl;
	  const I start_j_block_SU 	= 0;
#ifndef NDEBUG
	  Monitor_msg(this->iproc,"define slip block " ifmt " " ifmt " nddl_slip " ifmt "\n",start_i_block_SU,start_j_block_SU,self_->nddl_slip);
#endif
	  SparseBlock_def	(SparseStokes_get_SU(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_i_block_SU,
				 &self_->nddl_slip,
				 &start_j_block_SU,
				 &nddluxdim);    

#ifndef NDEBUG
	  Monitor_msg(this->iproc,"define slip block " ifmt " " ifmt "\n",start_i_block_SU,start_j_block_SU);
#endif
	  SparseBlock_def	(SparseStokes_get_US(self_->sparseStokes),
				 SparseStokes_get_A(self_->sparseStokes),
				 &start_j_block_SU,
				 &nddluxdim,
				 &start_i_block_SU,
				 &self_->nddl_slip);
	


	} 
      }

    /*  I global_rhs_n  		= mesh_->sparse_stokes.A.n;*/

    const I nshape_on_face_u 		= ensBASIS_n_onface(_shape_u,element);
    const I nshape_on_face_p 		= ensBASIS_n_onface(_shape_p,element);

#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_initialize:stokes dirichlet ...");
#endif

    SparseStokes_dirichlet_init	(self_->sparseStokes,
				 mesh_->nedge_boundary,
				 nshape_on_face_u,
				 nshape_on_face_p,
				 self_->dec_ddlu,
				 self_->dec_ddlu+self_->nddlu,
				 self_->dec_ddlp,
				 (slip)
				 ? self_->dec_ddlp+self_->nddlp
				 : 0);
    
  };

  
  
  inline void Global::Precompute(ns_mesh * 	const 	mesh_,
			 STR 			errmsg_,
			 Err* 			err_)
  {

#ifndef NDEBUG
    DebugVerif(mesh_);
#endif
    err_[0] 				= __eErr_no;
    const Parameters* parameters 	= this->GetParameters();
    pWorkelm	workelm			= Global_get_Workelm(this->m_self);
    const bool 	color 			= parameters->GetInfoLogical(InfoLogical::color);
  
    /*1 ###############################################*/
  
#ifndef NDEBUG
    const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 1/5:reading mesh from file '%s'",this->m_self->filename);
      }
#endif

    this->m_self->mesh_usrptr = mesh_;
    Workelm_init(workelm);

    //
    // 2
    //
#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data ...");
      }
#endif
  
    Fem_def(Global_get_Fem(this->m_self),
	    mesh_->elm,
	    errmsg_,
	    err_);
  
    if (err_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_precompute:ns_exact_init failed");
	return;
      }

#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data done");
      }
#endif

    //
    // 3
    //

#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize ...");
      }
#endif

    STR errmsg;
    
    Global_initializePart1(this->m_self,
			   (ns_mesh*)this->m_self->mesh_usrptr,
			   errmsg,
			   err_);
    if (err_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_precompute:nsGLOBAL_initialize failed");
	return;
      }

    {
      const I ntimesteps = 3;
      this->m_linsysVectors = new LinsysVectors(&this->m_self->total_nddl,
						&ntimesteps);  
      
      pR x[4];      
      x[0] = m_linsysVectors->GetX();
      x[1] = m_linsysVectors->GetXi(&negal1);
      x[2] = m_linsysVectors->GetXi(&negal2);
      x[3] = nullptr;
      const I n = m_linsysVectors->GetN();


      {
	err_[0] = __eErr_no;  
#if 0
	STR errmsg;
#endif
	I nddls		[__eVariable_ALL];
	I shifts	[__eVariable_ALL];

	{ eVariable variable = __eVariable_ERROR;
	  for (++variable;variable<__eVariable_ALL;++variable)
	    {
	      nddls[variable] = 0;
	    } }
  
	{ eVariable variable = __eVariable_ERROR;
	  for (++variable;variable<__eVariable_ALL;++variable)
	    {
	      shifts[variable] = 0;
	    } }

	const bool  include_dg 		= parameters->GetInfoLogical(InfoLogical::vnsdg);  
	const bool  transport_galerkin	= parameters->GetInfoLogical(InfoLogical::transport_galerkin); 
	const bool  transport 		= parameters->GetInfoLogical(InfoLogical::transport);  
	const bool  transport_uncoupled	= parameters->GetInfoLogical(InfoLogical::transport_uncoupled);
	const bool  pressure_uncoupled 	= parameters->GetInfoLogical(InfoLogical::pressure_uncoupled);
	const bool  slip 		= parameters->GetInfoLogical(InfoLogical::slip);

	nddls[__eVariable_U] 	= this->m_self->nddlu;
	nddls[__eVariable_P] 	= this->m_self->nddlp;
	shifts[__eVariable_U]	= this->m_self->dec_ddlu;
	shifts[__eVariable_P] 	= this->m_self->dec_ddlp;


	m_variables = new Variables(n,
				    ntimesteps,
				    x,
				    nddls,
				    shifts);
#if 0
	Variables_def(&this->m_self->m_variables,
		      &n,
		      &ntimesteps,
		      x,
		      nddls,
		      shifts,
		      errmsg,
		      err_);
	if (err_[0])
	  {
	    Monitor_errmsg(this->iproc,"nsGLOBAL_initialize:Variables_def failed '%s' ",errmsg);
	    return;
	  }
#endif

#ifndef NDEBUG
	Monitor_msg(this->iproc,"nsGLOBAL_initialize:create variables ...");
#endif
	static const eVariable variable_u = __eVariable_U;
	static const eVariable variable_p = __eVariable_P;
	pTimeReadOnly gTimeInfo = this->GetTime();
	ns_var_def	(&this->m_self->var_u,
			 mesh_,
			 "U",
			 _shape_u,
			 __emnsNO,
			 &negal2,
			 &gTimeInfo->t,
			 &gTimeInfo->ti,
			 &gTimeInfo->tii,
			 this->m_variables->Get(variable_u,0),
			 this->m_variables->Get(variable_u,1),
			 this->m_variables->Get(variable_u,2));
	
	ns_var_def	(&this->m_self->var_p,
			 mesh_,
			 "P",
			 _shape_p,
			 __emnsNO,
			 &negal1,
			 &gTimeInfo->t,
			 &gTimeInfo->ti,
			 &gTimeInfo->tii,
			 this->m_variables->Get(variable_p,0),
			 this->m_variables->Get(variable_p,1),
			 this->m_variables->Get(variable_p,2));
#ifndef NDEBUG
	Monitor_msg(this->iproc,"nsGLOBAL_initialize:divergence matrix ...");
#endif

	{ Sparse_clear(this->m_self->sparseStokes->A);    
	  /* COMPUTE DIVERGENCE MATRIX */
	  ns_build_system_divergence(this->m_self);  
	  /* COMPUTE SLIP MATRIX */
	  if (slip)
	    {
	      ns_build_system_slip(this->m_self,
				   this->m_self->slip_iperm,
				   &this->m_self->nddlelmu);
	    } }
  
#if 0
	pLinsys 	Linsys_new			(cst_eLinearSolver 			kind_,
							 cst_eSparseFactorizationMethod		factorization_method_,
							 cst_eSparseIterativeMethod		iterative_method_,
							 cst_eSparsePreconditioningMethod	preconditioning_method_);
#endif
	//	const eLinearSolver kindLinearSolverVelocity = parameters->GetKindLinearSolver(KindEquation::VELOCITY);

	if (NULL == this->m_self->sparseStokes->A)
	  {
	    fprintf(stderr,"A is NULL!\n");
	    exit(1);
	  }

	this->m_linsysA = new LinearSolver::Direct::MKL::Pardiso(this->m_self->sparseStokes->A);        
	//	this->m_linsysA = new LinearSolver::MKLPardiso(__eSparseFactorizationMethod_PARDISOLU,this->m_self->sparseStokes->A);        
	if (pressure_uncoupled)
	  {
	    //	    this->m_linsysB = new LinearSolver::Direct::Pardiso::Pardiso(this->m_self->sparseStokes->A);        
	    //	    this->m_linsysB  = new LinearSolver::MKLPardiso(__eSparseFactorizationMethod_PARDISOLU,this->m_self->sparseStokes->A);
	  }   
	if (((transport)AND(transport_uncoupled)AND(transport_galerkin)))
	  {
	    //	    this->m_linsysF = new LinearSolver::Direct::Pardiso::Pardiso(this->m_self->sparseStokes->A);        
	    //	    this->m_linsysF  = new LinearSolver::MKLPardiso(__eSparseFactorizationMethod_PARDISOLU,this->m_self->sparseStokes->A);
	  }
	
	if (include_dg)
	  {      
#if 0
	    ns_dg_decoupled_init(this->m_self,&this->m_self->nddlelmdg);
#endif
	  }

	if (pressure_uncoupled)
	  {
#if 0
	    this->m_self->nonlinear_solver = nsNONLINEAR_STOKES_DECOUPLED_GLOBAL_new	(&this->iproc,
											 &this->m_self->total_nddl,
											 &this->m_self->sparse_stokes.A,
											 &this->m_self->sparse_stokes.B,
											 __eFactorizationMethod_PARDISOLU,
											 __ensMETHOD_ITERATIVE_gmres,
											 __ensMETHOD_ITERATIVE_PRECOND_left);
#endif
	  }
	else
	  {
#if 0
	    this->m_self->nonlinear_solver = nsSTOKES_DIRECT_SOLVER_new	(this->iproc,
									 __ensMETHOD_ITERATIVE_gmres,
									 __ensMETHOD_ITERATIVE_PRECOND_left,
									 &this->m_self->sparse_stokes.A,
									 NULL,
									 __eFactorizationMethod_PARDISOLU);
#endif
#if 0
	    this->m_self->nonlinear_solver = nsSTOKES_DIRECT_SOLVER_new	(this->iproc,
									 &this->m_self->sparse_stokes.A,
									 __eFactorizationMethod_PARDISOLU);
#endif
	  }

      }

      
      if (err_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_precompute:nsGLOBAL_initialize failed");
	return;
      }
    }
#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize done.");
      }
#endif

    //
    // 4
    //
  
#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 4/5:nsCSF_INFO_def ...");
      }
#endif

    CsfInfo_def	(&this->m_self->m_csfInfo,
		 color ? __emnsYES : __emnsNO,
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
		 (eHeaviside)parameters->GetInfoInteger(InfoInteger::heaviside),
		 (eDirac)parameters->GetInfoInteger(InfoInteger::dirac),
		 &this->m_self->rv[__ens_rv_eps],
		 err_);  
  
#ifndef NDEBUG
    if (verbose)
      {
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 4/5:CsfInfo_def done.");
	Monitor_msg(this->iproc,"nsGLOBAL_precompute:step 5/5:completed.");
      }
#endif
    
  };

  
  void Global::SolveTimeStep()
  {
    
#if 0
    { R nrmdtf,nrmdtu;
      integrator_compute_derivative(&nrmdtf,&nrmdtu); }
    /*
      COMPUTE VOLUME
    */
    /*
      COMPUTE CFL APPROXIMATION
    */
    nsGLOBAL_compute_cfl	(&self_[0]);
#endif
    MyLinearOperator *h = new MyLinearOperator(this->GetSparseStokes());
    for (I jj=0;jj<1;++jj)
      {
	//	std::cerr << "give the matrix vector product" << std::endl;
	//	exit(1);	
	m_nonLinearSolverNewton->Solve(h,
				       this->m_linsysA,
				       this->m_parameters,
				       this->m_linsysVectors);
	
	auto nonLinearSolverResults = m_nonLinearSolverNewton->Results();
	if ( (nonLinearSolverResults->nrmc[this->m_self->newton_iter]<1.0e-8)
	     && (nonLinearSolverResults->nrmr[this->m_self->newton_iter]<1.0e-8) )
	  break;
      } 
    
  };

  void Global::NextTimeInterval()
  {
    pGlobal self_ = this->m_self;

    const Parameters* parameters 	= this->GetParameters();
    pTimeReadOnly 	gTime		= this->GetTime();
    const ns_mesh * 	mesh		= this->GetMesh();
    
    const I numNodes	= ns_mesh_get_numNodes(mesh);
    const I numEdges 	= ns_mesh_get_numEdges(mesh);

    nsSPACE_write_medit((pSpaceReadOnly)mesh->spaces[__ensBASIS_LAGRANGE_1],"integrator.P1.%.5" nsFORMAT_INTEGER "",gTime->itime_interval); 
    nsSPACE_write_medit((pSpaceReadOnly)mesh->spaces[__ensBASIS_LAGRANGE_2],"integrator.P2.%.5" nsFORMAT_INTEGER "",gTime->itime_interval); 
    if (parameters->GetInfoLogical(InfoLogical::vnsdg))
      {
	printf("print dg_print_mesh\n");
	nsGLOBAL_dg_print_mesh_interval	(&self_[0],
					 "dg",		 
					 _shape_dg);      
      }      
    else
      {
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:does not call dg_print_mesh\n");
#endif
    }

    const bool isAxisymetricX = parameters->GetInfoLogical(InfoLogical::axisymetric_x);
    const bool isAxisymetricY = parameters->GetInfoLogical(InfoLogical::axisymetric_y);
    if  (isAxisymetricX || isAxisymetricY)
    {
      const I axi_axe 	= isAxisymetricX?0:1;
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
      sprintf(ctmp,"integrator.P1axi.%.5" nsFORMAT_INTEGER "",gTime->itime_interval);	
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
	
    if  (isAxisymetricX || isAxisymetricY)
    {
      const I axi_axe 	= isAxisymetricX?0:1;
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
      sprintf(ctmp,"integrator.P2axi.%.5" nsFORMAT_INTEGER "",gTime->itime_interval);	
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
    sprintf(ctmp,"integrator.metric.%.5" nsFORMAT_INTEGER "",gTime->itime_interval);	
    nsMETRIC_print	(self_->metric_end_time_interval,ctmp); }
#endif

  };
  
  void Global::Run(eStateSolver 	job_[1],
		   cst_pI 		rwork_n_,
		   pR 			rwork_,
		   cst_pI 		iwork_n_,
		   pI 			iwork_)
  {
    pGlobal self_ = this->m_self;
    const I required_rwork_n = this->GetRequiredRworkN();
    if (required_rwork_n<rwork_n_[0])
      {
	Monitor_errmsg(this->iproc,"nsGLOBAL_run:not enough rwork_");
      }
    const Parameters* parameters	= this->GetParameters();
    pTime const gTimeInfo 			= this->GetTime();

    
    if (!this->m_stokesNonLinearSolverDriver)
      {
	this->m_stokesNonLinearSolverDriver = new StokesNonLinearSolverDriver(this->m_self);
	this->m_stokesNonLinearSolverDriver->SetParams(this->m_parameters);
	this->m_nonLinearSolverNewton = new NonLinearSolverNewton(m_stokesNonLinearSolverDriver);
      }
    
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_run:input job = " efmt "",job_[0]);
#endif
    
    const bool dynamic_boundary_condition 	= parameters->GetInfoLogical(InfoLogical::dynamic_boundary_condition);
    const bool slip_condition 		= parameters->GetInfoLogical(InfoLogical::slip);
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
	  if (slip_condition)
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
	  Monitor_errmsg(this->iproc,"nsGLOBAL_run:wrong input job(=__eStateSolver_exit)");
	  return;
	}
      
      case __eStateSolver_ALL:
	{
	  fprintf(stderr,"integrator job error " efmt "\n",job_[0]);
	  exit(1);
	}
      }


    Time_init(gTimeInfo);
    /*___ NEXT TIME INTERVAL _______________________________________*/
  
  state_next_time_interval:  

    Time_next_time_interval(gTimeInfo,self_->rv[__ens_rv_dt]*((R)parameters->GetInfoInteger(InfoInteger::ntime)));

    Monitor_msg(this->iproc,"time interval [%e,%e], mass %e",gTimeInfo->time_interval_t0,gTimeInfo->time_interval_tf,
		((R)0.0));

    /*nsSTATISTICS_get(&Z->stats,Z->adapt_iter,(Time->itime>1)?(Time->itime-1):0,__ens_rres_area)*/
    self_->adapt_iter = negal0;

    /*__ RESTART TIME INTERVAL ____________________________________*/
  state_restart_time_interval:
    Time_restart_time_interval(gTimeInfo);
#if 0
    fprintf(stdout,"restart interval[" ifmt "] [%e,%e,%e||||,%e,%e]\n",gTimeInfo->itime_interval,tii,ti,t,time_interval_t0,time_interval_tf);
#endif
    if (parameters->GetInfoLogical(InfoLogical::vnsdg))
      {
	if ( (parameters->GetInfoInteger(InfoInteger::ntime)>0) AND (gTimeInfo->itime_interval==0) )
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
	    if (parameters->GetInfoLogical(InfoLogical::linfo_color))
	      {	 
     
		ns_heaviside(parameters->GetInfoLogical(InfoLogical::heaviside),self_->cubature_triangle.q_n,rwork_,&self_->rv[__ens_rv_eps]);
	      }

	    nsblas_dgemv(transN,&self_->dg_n,&self_->cubature_triangle.q_n,&regal1,self_->qelm_dg_wbasis,&self_->dg_n,rwork_,&negal1,&regal0,&self_->var_dg.xi[self_->ielm*self_->dg_n],&negal1); 
	    nsblas_dcopy(&self_->dg_n,&self_->var_dg.xi[self_->ielm*self_->dg_n],&negal1,&self_->var_dg.x[self_->ielm*self_->dg_n],&negal1);	 	  
#endif
	    ++self_->ielm;
	    if (self_->ielm<self_->nelm)
	      goto redo_state___eStateSolver_initial_condition_dg;	  
	  }

#ifndef NDEBUG
	Monitor_msg(this->iproc,"nsGLOBAL_run:initial dg condition done");
#endif
      }

    /*___ INITIAL CONDITION  _______________________________________*/
    /* IF NOT STEADY STATE AND FIRST TIME INTERVAL */
    {
      const I nt = parameters->GetInfoInteger(InfoInteger::ntime);
      if ( (nt>0) && (gTimeInfo->itime_interval==0) )
	{
	  job_[0] = __eStateSolver_initial_condition;
	  return;
	}
    }
    
  state___eStateSolver_initial_condition:
#ifndef NDEBUG
    Monitor_msg(this->iproc,"nsGLOBAL_run:initial condition back");
#endif
    {
      const I nt = parameters->GetInfoInteger(InfoInteger::ntime);
      if ( (nt>0) AND (gTimeInfo->itime_interval==0) )
	{
#if 0
	  nsMETRIC_setidentity	(self_->metric);
#endif
	  ns_var_zeroi		(&self_->var_u);
	  ns_var_zeroi		(&self_->var_p);      
	  ns_var_cpy_xi2x	(&self_->var_u);
	  ns_var_cpy_xi2x	(&self_->var_p);      
	}
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


    /*_________ PROCHAIN PAS DE TEMPS _______________________________________*/  
  state_compute_next_time:
    {        
      Time_next_time(gTimeInfo,
		     parameters->GetInfoReal(InfoReal::dtmin));

      (*this->m_solverConstants)[SolverConstants::dti] 		= TimeReadOnly_get_dti(gTimeInfo);
      (*this->m_solverConstants)[SolverConstants::dt] 		= TimeReadOnly_get_dt(gTimeInfo);
      (*this->m_solverConstants)[SolverConstants::idt] 		= ((R)1.0)/(*this->m_solverConstants)[SolverConstants::dt];
      (*this->m_solverConstants)[SolverConstants::ratio_dti] 	= (*this->m_solverConstants)[SolverConstants::idt] * (*this->m_solverConstants)[SolverConstants::dti];
      (*this->m_solverConstants)[SolverConstants::iratio_dti] 	= ((R)1.0)/(*this->m_solverConstants)[SolverConstants::ratio_dti];
      (*this->m_solverConstants)[SolverConstants::midt] 	= -(*this->m_solverConstants)[SolverConstants::idt];


      //
      // SELF USE
      //
      self_->rv[__ens_rv_dti] 			= TimeReadOnly_get_dti(gTimeInfo);
      self_->rv[__ens_rv_dt]  			= TimeReadOnly_get_dt(gTimeInfo);
      self_->rv[__ens_rv_idt] 			= ((R)1.0)/self_->rv[__ens_rv_dt];
      self_->rv[__ens_rv_ratio_dti]   		= self_->rv[__ens_rv_idt]*self_->rv[__ens_rv_dti];
      self_->rv[__ens_rv_iratio_dti]		= ((R)1.0)/self_->rv[__ens_rv_ratio_dti];
      self_->rv[__ens_rv_midt]			= -self_->rv[__ens_rv_idt];

      //
      // SELF USE
      //      
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
	  Monitor_msg(this->iproc,"nsGLOBAL_run:slip velocity sym back");
#endif    
	}

      /* valeur des conditions de dirichlets (peut dependre du temps, donc il faut les recharger ) */    
      /* I L FAUT METTTRE UN FLAG POUR LA DEPENDENCE EN TEMPS DES CONDITIONS DE DIRICHLETS */
      job_[0] = __eStateSolver_dirichlet_pressure;
      return;
    state___eStateSolver_dirichlet_pressure:
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:dirichlet pressure back job " efmt "",job_[0]);
#endif    
      /* ____________________________________________________________________________ */

    
      /* ____________________________________________________________________________ */
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:go to job " efmt "",job_[0]);
#endif 
      job_[0] = __eStateSolver_dirichlet_velocity;
      return;        
    state___eStateSolver_dirichlet_velocity:    
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:dirichlet velocity back");
#endif
      /* ____________________________________________________________________________ */

    
      /* ____________________________________________________________________________ */
      if (parameters->GetInfoLogical(InfoLogical::transport))
	{
	  job_[0] = __eStateSolver_dirichlet_marker;
	  return;    
	}
    state___eStateSolver_dirichlet_marker:
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:dirichlet marker back");
#endif
      /* ____________________________________________________________________________ */
    
      if (parameters->GetInfoLogical(InfoLogical::transport))
	{
	  if (parameters->GetInfoLogical(InfoLogical::vnsdg))
	    {
	      job_[0] = __eStateSolver_dirichlet_dg;
	      return;    
	    }
	}
    state___eStateSolver_dirichlet_dg:
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:dirichlet dg back");
#endif
    
#ifndef NDEBUG
      Monitor_msg(this->iproc,"nsGLOBAL_run:solve time step ...");
#endif


      
      { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	if (verbose)
	  {
	    Monitor_msg(0,"StateSolver SOLVE TIME STEP. ");
	  } }
      
      this->SolveTimeStep();

      { const bool verbose = parameters->GetInfoLogical(InfoLogical::verbose);
	if (verbose)
	  {
	    Monitor_msg(0,"StateSolver SOLVE TIME STEP DONE. ");
	  } }

    
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
      if ( (gTimeInfo->itime>0) AND ((gTimeInfo->itime+1)%(parameters->GetInfoInteger(InfoInteger::ntime))==0) )
	{
#if 0
	  printf("reinit\n");
	  printf("reinit done\n");
#endif
	}

      if (parameters->GetInfoInteger(InfoInteger::ntime)>1)
	{	
	
	  if ( (parameters->GetInfoLogical(InfoLogical::redistance)) AND (gTimeInfo->itime>0) AND (gTimeInfo->itime%5==5) )
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
	  this->PrintVariables();

	  /* next */
	  ++gTimeInfo->itime;
	  this->UpdatePreviousX();

#if 0
	  LinsysVectors_update_previous_x(self_,
					  self_->linsysVectors);
#endif
	  
	  if (gTimeInfo->itime<(gTimeInfo->itime_interval+1)*parameters->GetInfoInteger(InfoInteger::ntime))
	    {
	      goto  state_compute_next_time;
	      exit(1);
	    }
	}       
    }     
  
    /* MESH ADAPTIVITY ___________________________________________________________________*/
    if (parameters->GetInfoLogical(InfoLogical::mesh_adaptivity))
      {
	if ( ((self_->adapt_iter<4)AND(gTimeInfo->itime_interval==0)) OR (self_->adapt_iter<parameters->GetInfoInteger(InfoInteger::mesh_adaptivity_maxiter)) )
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
	    gTimeInfo->itime = gTimeInfo->itime_interval*parameters->GetInfoInteger(InfoInteger::ntime);	  
	    ++self_->adapt_iter;	  
	  
	  
	    if (gTimeInfo->itime_interval==0)
	      {
		/*	      nsMASS_reset_iter(&self_->mass_info);*/
	      }	  
	  
	  
	    fprintf(stdout,"itime set to " ifmt "\n",gTimeInfo->itime);	  
	    goto  state_restart_time_interval;
	  }      
      } 

    /*__________________________________________________________________________*/
    if (gTimeInfo->itime_interval<=parameters->GetInfoInteger(InfoInteger::ntime_interval))
      {
	this->NextTimeInterval();
#if 0
	if (self_->metric)
	  {
	    nsMETRIC_copy	(self_->metric,
				 self_->metric_end_time_interval);
	  }
#endif
      
	Time_next_interval(gTimeInfo);
	if (gTimeInfo->itime_interval<=parameters->GetInfoInteger(InfoInteger::ntime_interval))
	  {
	    goto  state_next_time_interval;
	  }
      }
  
    job_[0] = __eStateSolver_exit;

#if 0    
    Global_run(this->m_self,
	       job_,
	       rwork_n_,
	       rwork_,
	       iwork_n_,
	       iwork_);
#endif
  };

  void Global::PrintVariables()
  {
    pTimeReadOnly gTimeInfo = this->GetTime();

    ns_var_print(&this->m_self->var_p,
		 "p",
		 &gTimeInfo->itime_interval,
		 &gTimeInfo->itime);

    ns_var_print(&this->m_self->var_u,
		 "u",
		 &gTimeInfo->itime_interval,
		 &gTimeInfo->itime);
    
#if 0
    const nsPARAMS_ST*const Params = &G_[0].params;
#endif
    
#if 0  
    if (  G_[0].params.linfo[__ens_linfo_vnsdg])
      {
#if 0
	{ mkSTRING aa;
	  sprintf(aa,"integrator.dg.%.5" nsFORMAT_INTEGER ".bin",tim->itime);
	  ns_var_write_bin(&G_[0].var_dg,aa); }
#endif
#if 0
	print_dscalar1	("fdisc",x,_nf);
#endif
	if (0)
	  {
	    
#if 0       
	    { mkSTRING aa;
	      sprintf(aa,"integrator.dg.%.5" nsFORMAT_INTEGER ".bb",tim->itime);
	      FILE * fich = fopen(aa,"w");
	      fprintf(fich,"2 1 " ifmt " 1\n",G_[0].pointeur_maillage->nelm);
	      {I i;for(i=0;i<G_[0].pointeur_maillage->nelm;++i)fprintf(fich,"%e\n",G_[0].var_dg.x[i*G_[0].nddlelmdg]*sqrt(2.0));}
	      fclose(fich);
#endif
	    }
	    
	  }
#if 0
	mkSTRING aa;
	sprintf(aa,"integrator.dg.%.5" nsFORMAT_INTEGER "",tim->itime);
	ns_var_write_medit_with_mesh(&G_[0].var_dg,aa); 
	
#endif
	
#if 0
	if (pele_pele==0)
	  {
	    pele_pele=1;
	    ns_print_dgspace	("dg.mesh",
				 _shape_dg);
	    
	  }
#endif
	
#if 0
	nsGLOBAL_print_dgspace	(&G_[0],
				 "visudg.mesh",
				 _shape_dg);
	nsGLOBAL_dg_print_sol	(&G_[0],
				 "visudg",		 
				 _shape_dg,
				 G_[0].var_dg.x,
				 &G_[0].dg_n);      
#endif
	
      }
#if 0
    if  ( (parameters->GetInfoLogical(InfoLogical::axisymetric_x]) OR (parameters->GetInfoLogical(InfoLogical::axisymetric_y]) )
      {
	const I axi_axe 	= (parameters->GetInfoLogical(InfoLogical::axisymetric_x])?0:1;
	I * axi_perm  	= imalloc(G_[0].pointeur_maillage->nvertex);
	I * axi_rperm 	= imalloc(G_[0].pointeur_maillage->nvertex);
	const I axi_nrot 	= 40;
	I axi_nvertex		= 0;
	I noaxi_nvertex	= 0;
	print_mesh_extrude_axi_permutation	(G_[0].pointeur_maillage,
						 &axi_axe,
						 axi_perm,
						 axi_rperm,
						 &axi_nvertex,
						 &noaxi_nvertex);	  
	char ctmp[512];
	sprintf(ctmp,"integrator.faxi.%.5" nsFORMAT_INTEGER "",tim->itime);
	print_scalar_extrude_axi		(G_[0].pointeur_maillage,
						 ctmp,
						 &axi_nrot,
						 axi_rperm,			   			   
						 &axi_nvertex,
						 &noaxi_nvertex,
						 G_[0].var_f.x);
	ifree(axi_perm);
	ifree(axi_rperm); 
      }	
#endif
#endif
  
  };

  
  void Global::UpdatePreviousX()
  {
    const Parameters* parameters	= this->GetParameters();
    pTimeReadOnly gTimeInfo		= this->GetTime();
    
    const I vN 	= m_linsysVectors->GetN();
    cst_pR x 	= m_linsysVectors->GetX();
    
    switch(parameters->GetKindTransientMethod(KindEquation::VELOCITY))
      {
	
      case KindTransientMethod::UNDEFINED:
	{
	  break;
	}
	
      case KindTransientMethod::EULER:
	{
	  static const I istep = 1;
	  cst_pR xi = m_linsysVectors->GetXi(&istep);
	  nsblas_dcopy(&vN,x,&negal1,xi,&negal1);	
	  break;
	}
	
      case KindTransientMethod::IMR:
      {
	static const I istep = 1;
	cst_pR xi = m_linsysVectors->GetXi(&istep);
	
	const I N  = _dim*this->m_self->nddlu;
	const I N2 = vN-N;
	/* 
	   mise a jour imr
	*/
	nsblas_dscal	(&N,&mregal1,xi,&negal1);
	nsblas_daxpy	(&N,&regal2,x,&negal1,xi,&negal1);	
	nsblas_dcopy	(&N2,&x[N],&negal1,&xi[N],&negal1);
	break;
      }
      
    case KindTransientMethod::GEAREULER:
      {
	static const I istep1 = 1;
	cst_pR xi = m_linsysVectors->GetXi(&istep1);

	static const I istep = 2;
 	cst_pR xii = m_linsysVectors->GetXi(&istep);

	nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	nsblas_dcopy(&vN,x,&negal1,xi,&negal1);		       	
	break;
      }
      
    case KindTransientMethod::GEARIMR:
      {
	static const I istep1 = 1;
	cst_pR xi = m_linsysVectors->GetXi(&istep1);

	static const I istep = 2;
 	cst_pR xii = m_linsysVectors->GetXi(&istep);

	if (gTimeInfo->itime>0)
	  {
	    nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	    nsblas_dcopy(&vN,x,&negal1,xi,&negal1);		       	
	  }
	else
	  {
	    nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	    const I N  =  _dim*this->m_self->nddlu;
	    const I N2 = vN-N;
	    /* 
	       mise a jour imr
	    */
	    nsblas_dscal	(&N,&mregal1,xi,&negal1);
	    nsblas_daxpy	(&N,&regal2,x,&negal1,xi,&negal1);	
	    nsblas_dcopy	(&N2,&x[N],&negal1,&xi[N],&negal1);
	  }
	break;
      }
    case KindTransientMethod::TRAPEZE:
      {
	Monitor_errmsg(this->iproc,"__eTransientMethod_TRAPEZE not yet");
	break;
      }
      
    case KindTransientMethod::ERROR:
    case KindTransientMethod::ALL_VALUES:
      {
	Monitor_errmsg(this->iproc,"nsLINSYS_VECTORS_update_previous_x:switch failed on __eTransientMethod \n");
	break;
      }
    }	      
  };
