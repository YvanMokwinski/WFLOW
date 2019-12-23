#ifndef __header_LinearSystem_hpp__
#define __header_LinearSystem_hpp__

#include "ILinearSystem.hpp"
#include "eLinearSolver.h"
#include "eSparseFactorizationMethod.h"
#include "eSparseIterativeMethod.h"
#include "eSparsePreconditioningMethod.h"
#include "SparseFactorization.h"
#include "SparseIterative.h"


class LinearSystemIterative : public ILinearSystem
{
protected:

  pSparseIterative	m_iterative;


public:

  virtual ~LinearSystemIterative()
  {
    this->m_iterative 		= SparseIterative_kill(this->m_iterative);
  };

  
  LinearSystemIterative(cst_eSparseFactorizationMethod		factorization_method_,
			cst_eSparseIterativeMethod		iterative_method_,
			cst_eSparsePreconditioningMethod	preconditioning_method_)
  {

    this->m_iterative	= SparseIterative_new(iterative_method_,
					      preconditioning_method_);
    
  };
  
  void 		Precompute(pSparse S_)
  {
#if 0
    SparseIterative_precompute	(this->m_factorization,S_);
#endif
  };

  
  static void routine_SparseIterativeDriver(void * 		usrptr_,
					    const char * 	tr_,
					    pR 			y_,
					    cst_pR 		x_)
  {
    if (tr_[0] == 'N')
      {
	static const R r1=1.0;
	static const R r0=0.0;
	Sparse_gemv(&r1,
		    (cst_pSparse)usrptr_,
		    x_,
		    &r0,
		    y_);      
      }
    else
      {
	fprintf(stderr,"not implemented\n");
      }
  };
  
  void 		Compute(pSparse		B_,
			cst_pR		y_,
			pR 		x_)
  {
    
    const I n = Sparse_get_n(B_);
    Err err = __eErr_no;
    I rwork_n = 0;
    SparseIterative_get_memory	(this->m_iterative,
				 &n,
				 &rwork_n,
				 &err);
    pR rwork = (pR)malloc(sizeof(R)*rwork_n);
    SparseIterativeDriver A;
    SparseIterativeDriver_def(&A,
			      routine_SparseIterativeDriver,
			      B_);
	  
    const  R residual_norm = SparseIterative_compute	(this->m_iterative,
							 &n,
							 x_,
							 y_,
							 &A,
							 NULL,
							 &rwork_n,
							 rwork);
    free(rwork);
  };
  
};

#endif
