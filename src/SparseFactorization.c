#include "SparseFactorization.h"

#include <stdlib.h>


static void SparseFactorization_clr(pSparseFactorization self_)
{
  memset(self_,0,sizeof(SparseFactorization));
}


pSparseFactorization   SparseFactorization_kill(pSparseFactorization  self_)
{
  if (self_)
    {
#if 0
      self_->extern_umfpack 	= ExternUmfpack_kill(self_->extern_umfpack);
#endif
      self_->extern_sparsekit 	= ExternSparsekit_kill(self_->extern_sparsekit);
      self_->extern_pardiso 	= ExternPardiso_kill(self_->extern_pardiso);
      SparseFactorization_clr(self_);
      free(self_);
    }
  return NULL;
}


pSparseFactorization   SparseFactorization_new(const eSparseFactorizationMethod 	factorization_method_,
					       pErr					err_)
{

  err_[0] = __eErr_no;
  pSparseFactorization  self = (pSparseFactorization)calloc((size_t)1,sizeof(SparseFactorization));

  if (!self) { err_[0] = __eErr_memory;  return NULL; }

  self->factorization_method 	=  factorization_method_;

  switch(factorization_method_)
    {
    case __eSparseFactorizationMethod_UMFPACKLU:
      {	
#if 0
	self->extern_umfpack = ExternUmfpack_new();
#endif
	break;
      }

    case __eSparseFactorizationMethod_PARDISOLU:
      {
	self->extern_pardiso = ExternPardiso_new();	
	break;
      }
      
    case __eSparseFactorizationMethod_ILU0:
    case __eSparseFactorizationMethod_MILU0:
    case __eSparseFactorizationMethod_ILUT:
    case __eSparseFactorizationMethod_ILUTP:
    case __eSparseFactorizationMethod_ILUK:
    case __eSparseFactorizationMethod_ILUD:
      {	
	self->extern_sparsekit = ExternSparsekit_new();
	break;
      }

    case __eSparseFactorizationMethod_CHOL:
      {
	break;
      }
      
    case __eSparseFactorizationMethod_ERROR:
    case __eSparseFactorizationMethod_ALL:
      {
	err_[0] = __eErr_switch;
	break;
      }
    }

  if (err_[0])
    return SparseFactorization_kill(self);

  return self;   
}


void SparseFactorization_precompute(pSparseFactorization  	self_,
				    pSparse			S_)
{
  self_->matrix = S_;
  switch(self_->factorization_method)
    {            
    case __eSparseFactorizationMethod_UMFPACKLU:
      {
#if 0
	ExternUmfpack_precompute(self_->extern_umfpack,S_);
#endif
	break;
      }

    case __eSparseFactorizationMethod_PARDISOLU:
      { 
	ExternPardiso_precompute(self_->extern_pardiso,S_);
	break;
      }

    case __eSparseFactorizationMethod_ILU0:
    case __eSparseFactorizationMethod_ILUT:
    case __eSparseFactorizationMethod_ILUTP:
    case __eSparseFactorizationMethod_ILUK:
    case __eSparseFactorizationMethod_ILUD:
    case __eSparseFactorizationMethod_MILU0:      
      {
	ExternSparsekit_precompute(self_->extern_sparsekit,S_);
	break;
      }

    case __eSparseFactorizationMethod_CHOL:
      {
	break;
      }

    case __eSparseFactorizationMethod_ERROR:
    case __eSparseFactorizationMethod_ALL:
      { 
	break;
      }
    }
  return;
}



void  SparseFactorization_compute(void   * 	self_,
				  const char * 	tr_,
				  pR 		sol_,
				  cst_pR 	rhs_)  
{
  pSparseFactorization self = (pSparseFactorization)self_;
  switch(self->factorization_method)
    {

    case __eSparseFactorizationMethod_UMFPACKLU:
      {	
#if 0
	ExternUmfpack_compute(self->extern_umfpack,tr_,sol_,rhs_);
#endif
	break;
      }
      
    case __eSparseFactorizationMethod_PARDISOLU:
      {
	ExternPardiso_compute(self->extern_pardiso,tr_,sol_,rhs_);
	break;
      }

    case __eSparseFactorizationMethod_ILU0:
    case __eSparseFactorizationMethod_MILU0:
    case __eSparseFactorizationMethod_ILUT:
    case __eSparseFactorizationMethod_ILUTP:
    case __eSparseFactorizationMethod_ILUK:
    case __eSparseFactorizationMethod_ILUD:
      {
	ExternSparsekit_compute(self->extern_sparsekit,tr_,sol_,rhs_);
	break;
      }

    case __eSparseFactorizationMethod_CHOL:
      {
	break;
      }

    case __eSparseFactorizationMethod_ERROR:
    case __eSparseFactorizationMethod_ALL:
      {
	break;
      }      
    }
  return;
}

