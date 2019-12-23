#include "ExternSparsekit.h"

#ifndef __MNS_WITH_SPARSEKIT__  
#error __MNS_WITH_SPARSEKIT__  
#endif

#if ( (__MNS_WITH_SPARSEKIT__!=0) && (__MNS_WITH_SPARSEKIT__!=1) )
#error wrong value for __MNS_WITH_SPARSEKIT__  
#endif

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"




#ifdef FIELD
#error FIELD already defined
#else
#define FIELD 		((pExternSparsekitImpl)self_)->
#endif

#ifdef CONST_FIELD
#error CONST_FIELD already defined
#else
#define CONST_FIELD 	((cst_pExternSparsekitImpl)self_)->
#endif


void 		ExternSparsekit_precompute	(pExternSparsekit 	const self_,
						 pSparse		S_)
{
#ifdef TOTOTOTOTOTOTTO

#if __MNS_WITH_SPARSEKIT__  
    case __eSparseFactorizationMethod_ILU0:
      {
	
#if __MK_WITH_SPARSEKIT__
	{
	  const I len = ( (N+3)*(10+2) + (10+1)*10/2) * 10;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 	= calloc(len,sizeof(R));
#if 0
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_extern_lib;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
#endif
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
	      }
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 	= calloc(2*N,sizeof(I));
#if 0
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
#endif
	    }	
	  if (! self_->preconditioner)
	    {
	      self_->preconditioner =  Sparse_get_new(__ensSP_intern_format_row_fortran,
							 N,
							 N,
							 N+Sparse_get_nc(mat_),
							 NULL,NULL,NULL,
							 err_);
#if 0
	      if (err_[0])

		{
		  ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		  break;
		}
#endif
	    }
	}
	mok_use_sparsekit_ilu0(&N,
			       Sparse_get_x(mat_),
			       Sparse_get_i(mat_),
			       Sparse_get_b(mat_),
			       Sparse_get_x(self_->preconditioner),
			       Sparse_get_i(self_->preconditioner),
			       Sparse_get_b(self_->preconditioner),
			       self_->sparsekit_iwork,
			       &self_->sparsekit_ierr);

	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif
      }
    case __eSparseFactorizationMethod_milu0:      
      {
#if __MK_WITH_SPARSEKIT__
	{
	  const I len = ( (N+3)*(10+2) + (10+1)*10/2) * 10;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 	= calloc(len,sizeof(R));
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
	      }
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 	= calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }	
	  if (! self_->preconditioner)
	    {
	      self_->preconditioner =  Sparse_get_new(__ensSP_intern_format_row_fortran,
					       N,
					       N,
					       N+Sparse_get_nc(mat_),
					       NULL,NULL,NULL,
					       err_);
	      if (err_[0])
		{
		  ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		  break;
		}
	    }
	}
	mok_use_sparsekit_milu0(&N,
				Sparse_get_x(mat_),
				Sparse_get_i(mat_),
				Sparse_get_b(mat_),
				Sparse_get_x(self_->preconditioner),
				Sparse_get_i(self_->preconditioner),
				Sparse_get_b(self_->preconditioner),
				self_->sparsekit_iwork,
				&self_->sparsekit_ierr);
	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif
      }

    case __eSparseFactorizationMethod_ILUT:
      {

#if __MK_WITH_SPARSEKIT__
	{
	  const I len = N;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 	= calloc(len,sizeof(R));
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
		}
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 	= calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }	
	if (! self_->preconditioner)
	  {
	    self_->preconditioner =  Sparse_get_new(__ensSP_intern_format_row_fortran,
					N,
					N,
					N+2*N*self_->sparsekit_lfil,
					NULL,
					NULL,
					NULL,
					err_);
	    if (err_[0])
	      {
		ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		break;
	      }
	  }	
	}
	I iwk = 2*N*self_->sparsekit_lfil;
	mok_use_sparsekit_ilut(&N,
			       (mat_->intern_format!=__ensSP_intern_format_row)?cst_mkSP_trvalue(mat_):cst_Sparse_get_x(mat_),
			       cst_Sparse_get_i(mat_),
			       cst_Sparse_get_b(mat_),
			       &self_->sparsekit_lfil,
			       &self_->sparsekit_tol,
			       Sparse_get_x(self_->preconditioner),
			       Sparse_get_i(self_->preconditioner),
			       Sparse_get_b(self_->preconditioner),
			       &iwk,
			       self_->sparsekit_rwork,
			       self_->sparsekit_iwork,
			       &self_->sparsekit_ierr);
	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	if (err_[0])
	  {
	    ns_errmsg("mok_use_sparsekit_ilut failed");
	    if (self_->sparsekit_ierr>0)
	      fprintf(stderr,"\tzero pivot encountered at step number ierr="ifmt"\n",self_->sparsekit_ierr);
	    else if (self_->sparsekit_ierr==-1)
	      fprintf(stderr,"\tinput matrix may be wrong\n");
	    else if (self_->sparsekit_ierr==-2)
	      fprintf(stderr,"\tthe matrix L overflows\n");
	    else if (self_->sparsekit_ierr==-3)
	      fprintf(stderr,"\tthe matrix U overflows\n");
	    else if (self_->sparsekit_ierr==-4)
	      fprintf(stderr,"\tIllegal value for lfil\n");
	    else if (self_->sparsekit_ierr==-5)
	      fprintf(stderr,"\tzero row encountered\n");
	  }
	mkSP_rowfortran2row(mat_);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif

      }
    case  __eSparseFactorizationMethod_ILUTP:
      {
#if __MK_WITH_SPARSEKIT__
	{
	  const I len = ( (N+3)*(10+2) + (10+1)*10/2) * 10;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 	= calloc(len,sizeof(R));
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
		}
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 	= calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }	
	  if (! self_->preconditioner)
	    {
	      self_->preconditioner =  Sparse_get_new(__ensSP_intern_format_row_fortran,
					       N,
					       N,
					       N+2*N*self_->sparsekit_lfil,
					       NULL,NULL,NULL,err_);
	      if (err_[0])
		{
		  ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		  break;
		}
	    }	

	  if  ( (2*N>self_->sparsekit_iperm_n) OR (! self_->sparsekit_iperm) )
	    {
	      if (self_->sparsekit_iperm)
		{
		  free(self_->sparsekit_iperm);
		}
	      self_->sparsekit_iperm_n = 2*N;
	      self_->sparsekit_iperm 	 = calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iperm)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }	

	  self_->sparsekit_permtol = (R)0.1;
	  
	}

	mok_use_sparsekit_ilutp(&N,
				Sparse_get_x(mat_),
				Sparse_get_i(mat_),
				Sparse_get_b(mat_),
				&self_->sparsekit_lfil,
				&self_->sparsekit_tol,
				&self_->sparsekit_permtol,
				&N,
				Sparse_get_x(self_->preconditioner),
				Sparse_get_i(self_->preconditioner),
				Sparse_get_b(self_->preconditioner),
				&self_->sparsekit_rwork_n,
				self_->sparsekit_rwork,
				self_->sparsekit_iwork,
				self_->sparsekit_iperm,
				&self_->sparsekit_ierr);
	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif
      }

    case __eSparseFactorizationMethod_ILUK:
      {

#if __MK_WITH_SPARSEKIT__
	{
	  const I len = ( (N+3)*(10+2) + (10+1)*10/2) * 10;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 		= calloc(len,sizeof(R));
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
		}
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 		= calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }	
	  if (! self_->preconditioner)
	    {
	      self_->preconditioner = Sparse_get_new(__ensSP_intern_format_row_fortran,
					      N,
					      N,
					      N+2*N*self_->sparsekit_lfil,
					      NULL,NULL,NULL,err_);
	      if (err_[0])
		{
		  ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		  break;
		}
	    }		  
	}

	mok_use_sparsekit_iluk(&N,
			       Sparse_get_x(mat_),
			       Sparse_get_i(mat_),
			       Sparse_get_b(mat_),
			       &self_->sparsekit_lfil,
			       Sparse_get_x(self_->preconditioner),
			       Sparse_get_i(self_->preconditioner),
			       Sparse_get_b(self_->preconditioner),
			       self_->sparsekit_levs,
			       &self_->sparsekit_rwork_n,
			       self_->sparsekit_rwork,
			       self_->sparsekit_iwork,
			       &self_->sparsekit_ierr);
	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif
      }
    case __eSparseFactorizationMethod_ilud:
      {
#if __MK_WITH_SPARSEKIT__
	{
	  const I len = ( (N+3)*(10+2) + (10+1)*10/2) * 10;
	  if  ( (len>self_->sparsekit_rwork_n) OR (! self_->sparsekit_rwork) )
	    {
	      if (self_->sparsekit_rwork)
		{
		  free(self_->sparsekit_rwork);
		}
	      self_->sparsekit_rwork_n 	= len;
	      self_->sparsekit_rwork 		= calloc(len,sizeof(R));
	      if (! self_->sparsekit_rwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}
	    }
	  if  ( (2*N>self_->sparsekit_iwork_n) OR (! self_->sparsekit_iwork) )
	    {
	      if (self_->sparsekit_iwork)
		{
		  free(self_->sparsekit_iwork);
		}
	      self_->sparsekit_iwork_n 	= 2*N;
	      self_->sparsekit_iwork 		= calloc(2*N,sizeof(I));
	      if (! self_->sparsekit_iwork)
		{
		  err_[0] = __ens_err_memory;
		  ns_errmsg("SparseFactorization_compute:calloc failed");
		  break;
		}		  
	    }
	  if (! self_->preconditioner)
	    {
	      self_->preconditioner =  Sparse_get_new(__ensSP_intern_format_row_fortran,
					       N,
					       N,
					       N+2*N*self_->sparsekit_lfil,
					       NULL,NULL,NULL,err_);
	      if (err_[0])
		{
		  ns_errmsg("SparseFactorization_compute:sparsekit:Sparse_get_new failed");
		  break;
		}
	    }	
	}
	mok_use_sparsekit_ilud(&N,
			       Sparse_get_x(mat_),
			       Sparse_get_i(mat_),
			       Sparse_get_b(mat_),
			       &self_->sparsekit_alph,
			       &self_->sparsekit_tol,
			       Sparse_get_x(self_->preconditioner),
			       Sparse_get_i(self_->preconditioner),
			       Sparse_get_b(self_->preconditioner),
			       &self_->sparsekit_rwork_n,
			       self_->sparsekit_rwork,
			       self_->sparsekit_iwork,
			       &self_->sparsekit_ierr);
	err_[0] = (self_->sparsekit_ierr) ? (__ens_err_extern_lib) : (__ens_err_no);
	break;
#else
	ns_err("SparseFactorization_compute:sparsekit not available");
	break;
#endif
      }
#else

#endif  

#endif
}


void  		ExternSparsekit_compute		(void   * 		self__,
						 const char * 		tr_,
						 pR 			sol_,
						 cst_pR 		rhs_)
{

#if __MNS_WITH_SPARSEKIT__  

  
#if 0
  pExternSparsekit self_ = (pExternSparsekit)self__;
  I N = Sparse_get_n(self_->preconditioner);
  mok_use_sparsekit_lusol(&N, 
			  sol_,
			  rhs_,
			  Sparse_get_x(self_->preconditioner), 
			  Sparse_get_i(self_->preconditioner), 
			  Sparse_get_b(self_->preconditioner));
    
#endif
  
#else 

		
#endif
}



pExternSparsekit 	ExternSparsekit_new(pErr err_)
{ 
#if __MNS_WITH_SPARSEKIT__  

  pExternSparsekit self = (pExternSparsekit)calloc(1,sizeof(ExternSparsekit));  

  self->sparsekit_lfil 		= (I)7;
  self->sparsekit_alph 		= (R)0.5;
  self->sparsekit_tol 		= (R)1.0e-5;
  
  return self;
  
#else
  
  err_[0] = __eErr_user;
  return NULL;
  
#endif
}

pExternSparsekit ExternSparsekit_kill(pExternSparsekit self_)
{
#if __MNS_WITH_SPARSEKIT__
  if (self_)
    {

  self_->preconditioner = Sparse_kill(self_->preconditioner);
  if (self_->sparsekit_iperm) 	free(self_->sparsekit_iperm);
  if (self_->sparsekit_iwork) 	free(self_->sparsekit_iwork);
  if (self_->sparsekit_rwork) 	free(self_->sparsekit_rwork);
  if (self_->sparsekit_levs) 	free(self_->sparsekit_levs);
  if (self_->sparsekit_fpar) 	free(self_->sparsekit_fpar);    

    }
#endif
  return NULL;
}

