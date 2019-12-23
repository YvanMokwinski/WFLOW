#include "SparseIterative.h"
#include "Blas.h"
void 	SparseIterativeDriver_def		(pSparseIterativeDriver 		self_,
						 __routine_SparseIterativeDriver 	routine_,
						 void * 				usrptr_)
{
  self_->routine = routine_;
  self_->usrptr = usrptr_;
}

extern void gmres_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void fgmres_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void tfqmr_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void bcgstab_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void dbcg_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void cg_		(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void bcg_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void cgnr_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void fom_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);
extern void dqgmres_	(cst_pI n,cst_pR  rhs_,pR sol,pI ipar_,pR fpar_,pR rwork_);



static void 	SparseIterative_clr(pSparseIterative self_)
{ 
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  self_->m_iterative_method = __eSparseIterativeMethod_ERROR;
  self_->m_preconditioning_method = __eSparsePreconditioningMethod_ERROR;
  { I i;
    for (i=0;i<64;++i)
      {
	self_->m_sparsekit_ipar[i] = ((I)0);
      } }
  { I i;
    for (i=0;i<64;++i)
      {
	self_->m_sparsekit_fpar[i] = ((R)0);
      } }

}

#if __MK_WITH_SPARSEKIT2__
#endif

R distdot_(cst_pI n_,cst_pR x_,cst_pI xoff_,cst_pR y_,cst_pI yoff_)
{

#ifndef NDEBUG
  DebugVerif(n_);
  DebugVerif(n_[0]>0);
  DebugVerif(x_);
  DebugVerif(xoff_);
  DebugVerif(xoff_[0]>0);
  DebugVerif(y_);
  DebugVerif(yoff_);
  DebugVerif(yoff_[0]>0);
#endif

  return Blas_ddot(n_,x_,xoff_,y_,yoff_);
}

void SparseIterative_checkError	(cst_pSparseIterative 	self_,
				 STR 			msg_)
{
  const I err = self_->m_sparsekit_ipar[0];
  if (err<0)
    {
      switch(err)
	{
	case -1:
	  {
	    sprintf(msg_,"SparseIterative solver has iterated too many times");
	    break;
	  }
	case -2:
	  {
	    sprintf(msg_,"SparseIterative solver was not given enough work space");
	    break;
	  }
	case -3:
	  {
	    sprintf(msg_,"SparseIterative solver is facing a break-down");
	    break;
	  }
	case -4:
	  {
	    sprintf(msg_,"SparseIterative solver failed on switch");
	    break;
	  }
	default:	  
	  {
	    sprintf(msg_,"SparseIterative solver failed");
	    break;
	  }
	}
    }
}

pSparseIterative SparseIterative_kill(pSparseIterative 	self_)
{
  if (self_)
    {
      SparseIterative_clr(self_);
      free(self_);
    }
  return NULL;
}


pSparseIterative SparseIterative_new(const eSparseIterativeMethod 		iterative_method_,
				     const eSparsePreconditioningMethod 	preconditioning_method_)
{
  pSparseIterative self = (pSparseIterative)calloc(1,sizeof(SparseIterative));

  SparseIterative_clr(self);

  self->m_iterative_method 		= iterative_method_;
  self->m_preconditioning_method 		= preconditioning_method_;

  self->m_sparsekit_ipar[3-1] 		= (I)1; /* stopping criteria */
  self->m_sparsekit_ipar[5-1] 		= (I)20;
  self->m_sparsekit_ipar[6-1] 		= (I)160;
  self->m_sparsekit_fpar[0]		= (R)1.0e-6;
  self->m_sparsekit_fpar[1]		= (R)1.0e-10;

  switch(preconditioning_method_)
    {
    case __eSparsePreconditioningMethod_NO:
      {
	self->m_sparsekit_ipar[2-1] 		= (I)0; 
	break;
      }
    case __eSparsePreconditioningMethod_LEFT:
      {
	self->m_sparsekit_ipar[2-1] 		= (I)1; 
	break;
      }
    case __eSparsePreconditioningMethod_RIGHT:
      {
	self->m_sparsekit_ipar[2-1] 		= (I)2; 
	break;
      }
    case __eSparsePreconditioningMethod_ERROR:
    case __eSparsePreconditioningMethod_ALL:
      {
	self->m_sparsekit_ipar[0] = -4;
	break;
      }
    }

  switch(iterative_method_)
    {
    case __eSparseIterativeMethod_GMRES:
    case __eSparseIterativeMethod_FGMRES:
    case __eSparseIterativeMethod_TFQMR:
    case __eSparseIterativeMethod_BCGSTAB:
    case __eSparseIterativeMethod_DBCG:
    case __eSparseIterativeMethod_CG:
    case __eSparseIterativeMethod_CGNR:
    case __eSparseIterativeMethod_BCG:
    case __eSparseIterativeMethod_FOM:
    case __eSparseIterativeMethod_DQGMRES:
      {
	break;
      }
    case __eSparseIterativeMethod_ERROR:
    case __eSparseIterativeMethod_ALL:
      {
	self->m_sparsekit_ipar[0] = ((I)-4);
	break;
      }
    }

  return self;
}


void SparseIterative_get_memory(cst_pSparseIterative			self_,
				cst_pI					n_,
				pI					rwork_n_,
				pErr 					err_)
{
  const I m 	= self_->m_sparsekit_ipar[5-1];
  const I n 	= n_[0];
  I N 		= 0;
  err_[0] 	= __eErr_no;
  /*
    GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
    FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
    TFQMR   == 11 * n
    BCGSTAB == 8 * n
    DBCG    == 11 * n
    CG      == 5 * n
    CGNR    == 5 * n
    BCG     == 7 * n
    FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
    DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
  */
  switch(self_->m_iterative_method)
    {
    case __eSparseIterativeMethod_GMRES:
      {
	N = (n+3)*(m+2) + ((m+1)*m)/2;
	break;
      }
    case __eSparseIterativeMethod_FGMRES:
      {
	N = 2*(n*(2*m+1) + ((m+1)*m)/2 + 3*m + 2);
	break;
      }
    case __eSparseIterativeMethod_TFQMR:
      {
	N = 11*n;
	break;
      }
    case __eSparseIterativeMethod_BCGSTAB:
      {
	N = 8*n;
	break;
      }
    case __eSparseIterativeMethod_DBCG:
      {
	N = 11*n;
	break;
      }
    case __eSparseIterativeMethod_CG:
      {
	N = 5*n;
	break;
      }
    case __eSparseIterativeMethod_CGNR:
      {
	N = 5*n;
	break;
      }
    case __eSparseIterativeMethod_BCG:
      {
	N = 7*n;
	break;
      }
    case __eSparseIterativeMethod_FOM:
      {
	N = (n+3)*(m+2) + ((m+1)*m)/2;
	break;
      }
    case __eSparseIterativeMethod_DQGMRES:
      {
	N = n + (m+1)*(2*n+4);
	break;
      }
    case __eSparseIterativeMethod_ERROR:
    case __eSparseIterativeMethod_ALL:
      {
	err_[0] 	= __eErr_switch;
	break;
      }
    }
  rwork_n_[0] = N;
}


void SparseIterative_kernel(pSparseIterative 	self_,
			    cst_pI		n_,
			    pR 			y_,
			    cst_pR 		x_,
			    pR 			rwork_)
{

  switch(self_->m_iterative_method)
    {
    case __eSparseIterativeMethod_GMRES:
      {
	gmres_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_FGMRES:
      {
	fgmres_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_TFQMR:
      {
	tfqmr_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_BCGSTAB:
      {
	bcgstab_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_DBCG:
      {
	dbcg_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_CG:
      {
	cg_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_CGNR:
      {
	cgnr_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_BCG:
      {
	bcg_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_FOM:
      {
	fom_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_DQGMRES:
      {
	dqgmres_(n_,x_,y_,self_->m_sparsekit_ipar,self_->m_sparsekit_fpar,rwork_);
	break;
      }
    case __eSparseIterativeMethod_ALL:
    case __eSparseIterativeMethod_ERROR:
      {
	self_->m_sparsekit_ipar[0] = -4;
	break;
      }
    }  

}



R 			SparseIterative_compute		(pSparseIterative			self_,
							 cst_pI					n_,
							 pR 					y_,
							 cst_pR 				x_,
							 pSparseIterativeDriver			A_,
							 pSparseIterativeDriver			P_,
							 cst_pI					rwork_n_,
							 pR					rwork_)
{
  self_->m_sparsekit_ipar[4-1] 	= rwork_n_[0];
  self_->m_sparsekit_ipar[1-1] 	= (I)0;  
  I we_continue			= 1;

  while (we_continue)
    {	
      SparseIterative_kernel(self_,n_,y_,x_,rwork_);
      switch (self_->m_sparsekit_ipar[0])
	{
	case 1:
	  A_->routine(A_->usrptr,"N",&rwork_[self_->m_sparsekit_ipar[9-1]-1],&rwork_[self_->m_sparsekit_ipar[8-1]-1] );
	  break;
	case 2:
	  A_->routine(A_->usrptr,"T",&rwork_[self_->m_sparsekit_ipar[9-1]-1],&rwork_[self_->m_sparsekit_ipar[8-1]-1] );
	  break;
	case 3:
	case 5:	  
	  P_->routine(P_->usrptr,"N",&rwork_[self_->m_sparsekit_ipar[9-1]-1],&rwork_[self_->m_sparsekit_ipar[8-1]-1] );
	  break;
	case 4:
	case 6:	  
	  P_->routine(P_->usrptr,"T",&rwork_[self_->m_sparsekit_ipar[9-1]-1],&rwork_[self_->m_sparsekit_ipar[8-1]-1] );
	  break;
	default:
	  we_continue = 0;
	  break;
	}        
    }

  A_->routine(A_->usrptr,"N",rwork_,y_ );
    
  { static const I negal1 	= (I)1;
    static const R mregal1 	= (R)-1.0;
    Blas_daxpy(n_,&mregal1,x_,&negal1,rwork_,&negal1);  
    self_->m_sparsekit_fpar[5-1] = Blas_dnrm2(n_,rwork_,&negal1); }

  printf("residual %e niter "ifmt"\n",self_->m_sparsekit_fpar[5-1],self_->m_sparsekit_ipar[7-1]);
  return self_->m_sparsekit_fpar[5-1];
  
}


#undef FIELD

