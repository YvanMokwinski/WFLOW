#pragma once

#include "INonLinearSolverDriver.hpp"
#include "NonLinearSolverResults.hpp"
#include "LinsysVectors.hpp"

#include "LinearSolver/IInverseOperator.hpp"
#include "LinearSolver/ILinearOperator.hpp"

#include "Performance.hpp"

class NonLinearSolverNewton
{
protected:
  INonLinearSolverDriver * 	m_nonLinearSolverDriver;
  NonLinearSolverResults 	m_nonLinearSolverResults;
  Performance 			m_performance;
protected:
  
  static void DoLineSearch(const I 	iproc_,
			   pI 		job_,
			   cst_pI      	n_,
			   cst_pR 	fold_,
			   cst_pR 	stpmax_,
			   cst_pR 	slope_,
			   cst_pR 	descent_,
			   cst_pR 	xold_,
			   pR 		x_,
			   cst_pR 	r_,
			   R 		rinfo_[])
  {
    R optimize_residu,a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,tmplam,test=((R)1.0e+30);
    I i;  

    static const R ALF 	= ((R)1.0e-4);
    static const R TOLX = ((R)1.0e-7);
    
    if (job_[0]==1)
      {
	alamin 	= rinfo_[0];
	alam 	= rinfo_[1];
	f2	= rinfo_[2];
	fold2   = rinfo_[3];
	alam2	= rinfo_[4];
	goto start_backtrack_back;
      }
    
    for (i=0;i<n_[0];++i) 
      {
	const R bb	= (xold_[i]<regal1)
	  ? regal1
	  : xold_[i];
	const R temp 	= nsFABS(descent_[i])/bb;
	if (temp>test)
	  {
	    test=temp;
	  }
      }
    
    alamin 	= TOLX/test;
    alam   	= regal1;  
    rinfo_[0] 	= alamin;
    rinfo_[1] 	= alam;
    
    /* residu */
  start_backtrack:
    if (alam<((R)1.0e-2))
      {
	
#ifndef NDEBUG
	Monitor_warn	(iproc_,
			 "too small alam %e < 1.0e-2, set alam=1/4",
			 alam);      
#endif
	
#if 0
#ifndef NDEBUG
	ns_warn("too small alam %e < 1.0e-2, set alam=1/4",alam);      
#endif
#endif
	alam 	= ((R)0.25);
	rinfo_[1] = alam;
	nsblas_dcopy(n_,xold_,&negal1,x_,&negal1);
	nsblas_daxpy(n_,&alam,descent_,&negal1,x_,&negal1);      
	job_[0]   = 0;
	return;
      }
    
    nsblas_dcopy(n_,xold_,&negal1,x_,&negal1);
    nsblas_daxpy(n_,&alam,descent_,&negal1,x_,&negal1);    
    job_[0]   = 1;
    return;
  start_backtrack_back:  
    alamin 		= rinfo_[0];
    alam   		= rinfo_[1];
    
    rinfo_[0] = alamin;
    rinfo_[1] = alam;
    
    optimize_residu 	= r_[0];
    if (alam<alamin)
      {
	nsblas_dcopy(n_,xold_,&negal1,x_,&negal1);
	rinfo_[1] = alam;
	job_[0] = 0;
	return;
      }  
    else if (optimize_residu<=fold_[0]+ALF*alam*slope_[0]) 
      {
	rinfo_[1] = alam;
	job_[0] = 0;
	return;
      }
    else 
      {
	if (alam==regal1)
	  {
	    tmplam = -slope_[0] / (regal2*(optimize_residu-fold_[0]-slope_[0]));
	  }
	else
	  {
	    rhs1 	=  optimize_residu-fold_[0] - alam*slope_[0];
	    rhs2 	= f2-fold2-alam2*slope_[0];
	    a 		= (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	    b		= (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	    if (a==regal0) 
	      {
		tmplam = -slope_[0]/(b*regal2);
	      }
	    else
	      {		  
		disc = b*b-regal3*a*slope_[0];
		if (disc<regal0)
		  {
		    Monitor_warn	(iproc_,
					 "linesearch roundoff error halt set tau=1");
		    alam  	= regal1;
		    nsblas_dcopy(n_,xold_,&negal1,x_,&negal1);
		    nsblas_daxpy(n_,&alam,descent_,&negal1,x_,&negal1);      
		    rinfo_[1] 	= alam;
		    job_[0] 	= 0;
		    return;
		  }
		else 
		  {
		    tmplam = (-b+nsSQRT(disc))/(regal3*a);
		  }
	      }
	    if (tmplam>regal1demi*alam)
	      {
		tmplam=regal1demi*alam;
	      }
	  }	  
      }
    alam2	= alam;
    f2		= optimize_residu;
    fold2	= fold_[0];
    alam 	*= regal1demi;
    if (alam<tmplam)
      {
	alam=tmplam;
      }
    rinfo_[0] = alamin;
    rinfo_[1] = alam;
    rinfo_[2] = f2;
    rinfo_[3] = fold2;
    rinfo_[4] = alam2;
    goto  start_backtrack;
  };
  
  
  void ApplyLineSearch(const I 	iproc_,
		       pR 	x_,
		       cst_pI 	n_,
		       pR 	global_grad_,
		       pR 	global_rhs_,
		       pR 	global_c_,
		       pR 	global_x_,
		       pR 	global_xi_,
		       pR 	global_xii_)
  { 

    const R 	stmpmax_ref 		= ((R)100.0);
    const R 	a0 			= nsblas_ddot(n_,global_x_,&negal1,global_x_,&negal1);
    const R 	a1 			= nsSQRT(a0);
    const R 	a2 			= (a1 < ((R)n_[0]))?((R)n_[0]):a1;
    const R 	stpmax  		= stmpmax_ref * a2;    
    R 		rinfo_backtrack[32]{};
    const R 	optimal_residuold 	= nsblas_ddot(n_,global_rhs_,&negal1,global_rhs_,&negal1)*regal1demi;    
    
    /* on inverse la correction */
    nsblas_dscal(n_,
		 &mregal1,
		 global_c_,
		 &negal1);
    /**/
    { const R sum = nsSQRT(nsblas_ddot(n_,
				       global_c_,
				       &negal1,
				       global_c_,
				       &negal1));
      if (sum>stpmax)
	{      
	  const R aa = stpmax / sum;
	  nsblas_dscal(n_,
		       &aa,
		       global_c_,
		       &negal1);
	} }
    
    const R slope = nsblas_ddot(n_,
				global_grad_,
				&negal1,
				global_c_,
				&negal1);
    
    nsblas_dcopy(n_,
		 global_x_,
		 &negal1,
		 global_grad_,
		 &negal1);
    
    { R optimal_residu = regal0;
      I job = 0;
      DoLineSearch(iproc_,
		   &job,
		   n_,
		   &optimal_residuold,
		   &stpmax,
		   &slope,
		   global_c_,
		   global_grad_,
		   global_x_,
		   &optimal_residu,
		   rinfo_backtrack);
      
      while (job>0)
	{

	  for (I i=0;i<n_[0];++i)
	    {
	      global_rhs_[i]=((R)0.0);
	    }
	  
	  m_nonLinearSolverDriver->BuildResidu(global_rhs_,
					       global_x_,
					       global_xi_,
					       global_xii_);
	  
	  optimal_residu = nsblas_ddot(n_,
				       global_rhs_,
				       &negal1,
				       global_rhs_,
				       &negal1)*regal1demi;	

	  
	  DoLineSearch(iproc_,
		       &job,
		       n_,
		       &optimal_residuold,
		       &stpmax,
		       &slope,
		       global_c_,
		       global_grad_,
		       global_x_,
		       &optimal_residu,
		       rinfo_backtrack);
	} }
    /**/
    x_[0] = rinfo_backtrack[1];    
  };


 public:

  const NonLinearSolverResults* Results() const
  {
    return &this->m_nonLinearSolverResults;
  }
  
  NonLinearSolverNewton(INonLinearSolverDriver * nonLinearSolverDriver_)
  {
    m_nonLinearSolverDriver 	= nonLinearSolverDriver_;
  };


  
  void Solve(LinearSolver::ILinearOperator *    matrixVectorOperator_,
	     LinearSolver::IInverseOperator * 	inverseOperator_,
	     Parameters * 			parameters_,
	     LinsysVectors * 			linsysVectors_)
  {

#define ns_btime()  		((R)-1.0)*((R) clock()/CLOCKS_PER_SEC)
#define ns_time(_t) 		((_t)+((R) clock()/CLOCKS_PER_SEC))
    
    const bool  verbose 			= parameters_->GetInfoLogical(InfoLogical::verbose);
    //    Err 	err 			= __eErr_no;
    I maxiter[1]			= {parameters_->GetInfoInteger(InfoInteger::newton_maxiter)};
    R tol_correction[1]			= {parameters_->GetInfoReal(InfoReal::newton_tol_residu)};
    R tol_residu[1]			= {parameters_->GetInfoReal(InfoReal::newton_tol_correc)};

    static const I step2 = 2;
    static const I step1 = 1;

    pR rhs 				= linsysVectors_->GetRhs();
    pR c 				= linsysVectors_->GetCorr();
    pR grad 				= linsysVectors_->GetGrad();
    pR x 				= linsysVectors_->GetX();
    
    pR xi 				= linsysVectors_->GetXi(&step1);
    pR xii 				= linsysVectors_->GetXi(&step2);    
    const I N				= linsysVectors_->GetN();


    m_nonLinearSolverResults.Reset();

    m_performance.Reset();
    m_performance[Performance::Initial] = ns_btime();

    I newton_iter = 0;
    

#define clr_dvect(_n,_x) { I _k; for (_k=0;_k<(_n);++_k) { (_x)[_k] = regal0;  }}
    
    /*0-A CLEAR GLOBAL_RHS ############################################################### */
    clr_dvect(N,rhs);

#undef clr_dvect
    /*0-B COMPUTE RESIDU ################################################################# */
    
    m_performance[Performance::Work] = ns_btime();
    
    {
      R perf_t = ns_btime();

      m_nonLinearSolverDriver->BuildResidu(rhs,
					   x,
					   xi,
					   xii);
      
      m_nonLinearSolverResults.nrmr[newton_iter]  = nsblas_dnrm2(&N,
								 rhs,
								 &negal1);
      
      perf_t = ns_time(perf_t);
      m_performance[Performance::Residu] = perf_t;
      m_performance[Performance::TotalResidu] += perf_t;
    }

    /*1-LOOP ############################################################################  */
    for (newton_iter=0;newton_iter<maxiter[0];++newton_iter)
      {	
	/*  ruser[ruser_newton_iter] = iter; */
	if (newton_iter>0)
	  {
	    /* 1-A CALCUL DU GRADIENT global_grad=sparse_A*global_rhs ######################### */
	    matrixVectorOperator_->Apply("No transpose",grad,rhs);
#if 0
	    Sparse_gemv(&regal1,
			m_self->sparseStokes.A,
			rhs,
			&regal0,
			grad);  
#endif	  
	    /* 1-B MISE A JOUR AVEC LINESEARCH ################################################ */
	    m_performance[Performance::Work] = ns_btime();

	    int iproc = 0;
	    ApplyLineSearch(iproc,
			    &m_nonLinearSolverResults.nrmh[newton_iter],
			    &N,
			    grad,
			    rhs,
			    c,
			    x,
			    xi,
			    xii);
#if 0
	    if (err)
	      {
		Monitor_errmsg(0,"nsGLOBAL_newton:ns_applylinesearch failed");
		break;
	      }
#endif
	    if (verbose)
	      {
		Monitor_msg(0,"nsGLOBAL_newton:newton step " rfmt "",    m_nonLinearSolverResults.nrmh[newton_iter]);
	      }

	    m_performance[Performance::Residu] = ns_time(m_performance[Performance::Work]);
	    m_performance[Performance::TotalResidu] += m_performance[Performance::Residu];

	    /* RESULTS  ############################################ */
	    m_nonLinearSolverResults.nrmr[newton_iter] 	= nsblas_dnrm2(&N,rhs,&negal1);	    
	    m_nonLinearSolverResults.nrmc[newton_iter] 	= nsblas_dnrm2(&N,c,&negal1);
	    if (verbose)
	      {
		Monitor_msg(0,"%.2" nsFORMAT_INTEGER " c:%.2e r:%.2e",
			    newton_iter,
			    m_nonLinearSolverResults.nrmc[newton_iter],
			    m_nonLinearSolverResults.nrmr[newton_iter]);
	      }
#if 0	    
	    if (verbose)
	      {
		Monitor_msg(0,"%.4" nsFORMAT_INTEGER " %.2e %.2" nsFORMAT_INTEGER " c:%.2e r:%.2e",
			    itime_[0],
			    t_[0],
			    newton_iter,
			    m_nonLinearSolverResults.nrmc[newton_iter],
			    m_nonLinearSolverResults.nrmr[newton_iter]);
	      }
#endif
	    /* STOPPING CRITERIA */
	    if ( (    m_nonLinearSolverResults.nrmc[newton_iter]<tol_correction[0])
		 AND
		 (    m_nonLinearSolverResults.nrmr[newton_iter]<tol_residu[0]) )
	      {
		Monitor_msg(0,"%.2" nsFORMAT_INTEGER " c:%.2e r:%.2e",
			    newton_iter,
			        m_nonLinearSolverResults.nrmc[newton_iter],
			        m_nonLinearSolverResults.nrmr[newton_iter]);

#if 0
		Monitor_msg(0,"%.4" nsFORMAT_INTEGER " %.2e %.2" nsFORMAT_INTEGER " c:%.2e r:%.2e",
			    itime_[0],
			    t_[0],
			    newton_iter,
			        m_nonLinearSolverResults.nrmc[newton_iter],
			        m_nonLinearSolverResults.nrmr[newton_iter]);
#endif
		break;	  	 	  
	      }
	  }	

	/* POUR NE PAS CALCULER LE SYSTEME POUR RIEN EN FIN DE BOUCLE SI ON N'A PAS CONVERGE */
	if (newton_iter==maxiter[0]-1)
	  {

	    Monitor_warn	(0,
				 "nsGLOBAL_newton:warning stopping criteria not satisfied\nconvergence history :");
	    { I i;
	      for (i=0;i<=newton_iter;++i)
		{
		  Monitor_warn	(0,
				 "%.2" nsFORMAT_INTEGER " h_c:%.2e c:%.2e r:%.2e",
				 i,
				 m_nonLinearSolverResults.nrmh[i],
				 m_nonLinearSolverResults.nrmc[i],
				 m_nonLinearSolverResults.nrmr[i]);

#if 0
		  Monitor_warn	(0,
				 "%.4" nsFORMAT_INTEGER " %.2" nsFORMAT_INTEGER " h_c:%.2e c:%.2e r:%.2e t:%e",
				 itime_[0],i,
				 m_nonLinearSolverResults.nrmh[i],
				 m_nonLinearSolverResults.nrmc[i],
				 m_nonLinearSolverResults.nrmr[i],
				 t_[0]);
#endif
		} }
	    break;
	  }
            
	/* ANNULATION PARTIELLE DE LA MATRICE ######################### */
	m_performance[Performance::Solve] = ((R)0.0);

	/* BUILD JACOBIAN ############################################# */
	m_performance[Performance::Work] = ns_btime();

	m_nonLinearSolverDriver->BuildSystem(&negal3,
					     x,
					     &N);

	m_performance[Performance::System] = ns_time(m_performance[Performance::Work]);
	m_performance[Performance::TotalSystem] += m_performance[Performance::System];

	/* RESOLUTION DU SYSTEME LINEAIRE ############################# */
	m_performance[Performance::Work] = ns_btime();


	{
	  bool precomputeHasFailed;
	  inverseOperator_->Compute(&precomputeHasFailed);
	  if (precomputeHasFailed)
	    {
	      std::cerr << inverseOperator_->GetErrorMessage() << std::endl;
	      exit(1);
	    }
	}
	
	m_performance[Performance::Solve] = ns_time(m_performance[Performance::Work]);
	/* ############################################################ */      
	m_performance[Performance::Work] = ns_btime();

	{ I i;
	  for (i=0;i<N;++i)
	    {
	      c[i]=0.0;
	    } }
	
	{
	  bool applyHasFailed;
	  const I sizeOfTemporaryVector = inverseOperator_->GetSizeOfTemporaryVector();
	  pR temporaryVector = (sizeOfTemporaryVector>0) ? (pR)malloc(sizeof(R)*sizeOfTemporaryVector) : NULL;
	  inverseOperator_->Apply("No transpose",
				  c,
				  rhs,
				  sizeOfTemporaryVector,
				  temporaryVector,
				  &applyHasFailed);
	  if (applyHasFailed)
	    {
	      std::cerr << inverseOperator_->GetErrorMessage() << std::endl;
	      exit(1);
	    }

	  if (NULL != temporaryVector)
	    {
	      free(temporaryVector);
	    }
	}
	
	m_performance[Performance::Solve] = ns_time(m_performance[Performance::Work]);
	/* ############################################################ */       
	m_performance[Performance::TotalSolve] += m_performance[Performance::Solve];

	m_performance[Performance::Iter] = m_performance[Performance::Solve]
	  + m_performance[Performance::System]
	  + m_performance[Performance::Residu];
	
	m_performance[Performance::TotalIter] += m_performance[Performance::Iter];
#if 0
	if (verbose)
	  {
	    Monitor_msg(0,
			"newton time percent(residu,sys,solve) %.0f %.0f %.0f",
			(perf->p[perf_t_residu]/perf->p[perf_t_iter])*100.0,
			(perf->p[perf_t_sys]/perf->p[perf_t_iter])*100.0,
			(perf->p[perf_t_solve]/perf->p[perf_t_iter])*100.0);
	  }
#endif
      
      }   
    //    return __eErr_no;
  
  
  };
  
};

