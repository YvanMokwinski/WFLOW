#include "ExternPardiso.h"

#ifndef __MNS_WITH_PARDISO__  
#error __MNS_WITH_PARDISO__  
#endif

#if ( (__MNS_WITH_PARDISO__!=0) && (__MNS_WITH_PARDISO__!=1) )
#error wrong value for __MNS_WITH_PARDISO__  
#endif

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"

void 		ExternPardiso_precompute	(pExternPardiso 	const self_,
						 pSparse		S_)
{

  
#if __MNS_WITH_PARDISO__  

  self_->matrix = S_;
  self_->N = Sparse_get_n(S_);
  self_->pardiso_phase = 11;
  pardiso (&self_->pardiso_pt[0],
	   &self_->pardiso_maxfct,
	   &self_->pardiso_mnum,
	   &self_->pardiso_mtype,
	   &self_->pardiso_phase,
	   &self_->N,
	   Sparse_get_ownx(S_),
	   Sparse_get_ownb(S_),
	   Sparse_get_owni(S_),
#if 0
	   (pI)mat_->b,
	   (pI)mat_->i,
#endif
	   self_->pardiso_idum,
	   &self_->pardiso_nrhs,
	   self_->pardiso_iparm,
	   &self_->pardiso_msglvl,
	   &self_->pardiso_ddum,
	   &self_->pardiso_ddum,
	   &self_->pardiso_error);
  if (self_->pardiso_error != 0) 
    {
      fprintf(stderr,"nsSELFprecompute:pardiso failed during  symbolic factorization\n");
      return;
    }		
  self_->pardiso_phase = 22;
  pardiso (&self_->pardiso_pt[0],
	   &self_->pardiso_maxfct,
	   &self_->pardiso_mnum,
	   &self_->pardiso_mtype,
	   &self_->pardiso_phase,
	   &self_->N,
	   Sparse_get_ownx(S_),
	   Sparse_get_ownb(S_),
	   Sparse_get_owni(S_),
#if 0
	   mat_->own_x,
	   (pI)mat_->b,
	   (pI)mat_->i,
#endif
	   self_->pardiso_idum,
	   &self_->pardiso_nrhs,
	   self_->pardiso_iparm,
	   &self_->pardiso_msglvl,
	   &self_->pardiso_ddum,
	   &self_->pardiso_ddum,
	   &self_->pardiso_error);
  if (self_->pardiso_error != 0) 
    {
      fprintf(stderr,"nsSELFprecompute:pardiso failed during  numerical factorization\n");
      return;
    }		
#else

#endif  
}



void  		ExternPardiso_compute		(void   * 		self__,
						 const char * 		tr_,
						 pR 			sol_,
						 cst_pR 		rhs_)
{
  pExternPardiso self_ = (pExternPardiso)self__;  
#if __MNS_WITH_PARDISO__  

  self_->pardiso_phase = 33;
  self_->pardiso_iparm[7] = 2;  /* Max numbers of iterative refinement steps. */ 
  pardiso (&self_->pardiso_pt[0],	
	   &self_->pardiso_maxfct,
	   &self_->pardiso_mnum,
	   &self_->pardiso_mtype,
	   &self_->pardiso_phase,
	   &self_->N,
	   Sparse_get_ownx(self_->matrix),
	   Sparse_get_ownb(self_->matrix),
	   Sparse_get_owni(self_->matrix),

#if 0
	   &self_->matrix->n,	
	   self_->matrix->own_x,	
	   (pI)self_->matrix->b,	
	   (pI)self_->matrix->i,	
#endif

	   self_->pardiso_idum,
	   &self_->pardiso_nrhs,
	   self_->pardiso_iparm,
	   &self_->pardiso_msglvl,
	   (void*)rhs_,	
	   sol_,
	   &self_->pardiso_error);
  if (self_->pardiso_error != 0) 
    {
      fprintf(stderr,"ns_linsys_compute:pardiso failed during  numerical factorization\n");
      return;
    }

#else 

		
#endif
}



pExternPardiso 	ExternPardiso_new(pErr err_)
{
  
#if __MNS_WITH_PARDISO__  
  pExternPardiso self = (pExternPardiso)calloc(1,sizeof(ExternPardiso));  
  self->pardiso_mtype 		= 11 ;        /* Real unsymmetric matrix */
  self->pardiso_nrhs		= 1;
  
  { I i;
    for (i=0;i<64;++i)
      {
	self->pardiso_iparm[i]=(I)0;
      } }
  
  { I i;
    for (i=0;i<64;++i)
      {
	self->pardiso_pt[i]=NULL;
      } }
  
  { I i;
    for (i=0;i<64;++i)
      {
	self->pardiso_idum[i]=(I)0;
      } }
  
  self->pardiso_iparm[0] = 0; /* No solver default */
  self->pardiso_iparm[1] = 3; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  self->pardiso_iparm[2] = 16;
  self->pardiso_iparm[3] = 61; /* No iterative-direct algorithm                 ACCELERATION PARDISO KRYLOV */
  self->pardiso_iparm[4] = 0; /* No user fill-in reducing permutation */
  self->pardiso_iparm[5] = 0; /* Write solution into x */
  self->pardiso_iparm[6] = 0; /* Not in use */
  self->pardiso_iparm[7] = 4; /* Max numbers of iterative refinement steps */
  self->pardiso_iparm[8] = 0; /* Not in use */
  self->pardiso_iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  self->pardiso_iparm[10] = 1; /* Use nonsymmetric permutation and scaling MK */
  self->pardiso_iparm[11] = 0; /* Not in use */
  self->pardiso_iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  self->pardiso_iparm[13] = 0; /* Output: Number of perturbed pivots */
  self->pardiso_iparm[14] = 0; /* Not in use */
  self->pardiso_iparm[15] = 0; /* Not in use */
  self->pardiso_iparm[16] = 0; /* Not in use */
  self->pardiso_iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  self->pardiso_iparm[18] = -1; /* Output: Mflops for LU factorization */
  self->pardiso_iparm[19] = 0; /* Output: Numbers of CG Iterations */
  self->pardiso_maxfct 	= 1; /* Maximum number of numerical factorizations. */
  self->pardiso_mnum = 1; /* Which factorization to use. */
  self->pardiso_error = 0;/* Initialize error flag */
  self->pardiso_mnum   = 1; /* Which factorization to use. */
  self->pardiso_msglvl = 0; /* Print statistical information  */
  return self;
#else
  err_[0] = __eErr_user;
  return NULL;
#endif
}


pExternPardiso ExternPardiso_kill(pExternPardiso self_)
{
#if __MNS_WITH_PARDISO__
  if (self_)
    {
      /* Release internal memory. */
      self_->pardiso_phase 	= -1;                 
      pardiso (self_->pardiso_pt, 
	       &self_->pardiso_maxfct,
	       &self_->pardiso_mnum,
	       &self_->pardiso_mtype,
	       &self_->pardiso_phase,
	       &self_->N, 
	       &self_->pardiso_ddum[0], 
	       Sparse_get_ownb(self_->matrix), 
	       Sparse_get_owni(self_->matrix), 
	       self_->pardiso_idum, 
	       &self_->pardiso_nrhs,
	       self_->pardiso_iparm, 
	       &self_->pardiso_msglvl,
	       &self_->pardiso_ddum[0],
	       &self_->pardiso_ddum[0],
	       &self_->pardiso_error);
    }
#endif
  return NULL;
}

