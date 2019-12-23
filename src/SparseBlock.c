#include "SparseBlock.h"
#include <stdlib.h>

#ifndef NDEBUG
#define wrong_param_nsSPARSE_BLOCK(_cond,_i,_s) if ( (_cond) ) { fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); exit(1);} ((void)0)
#endif

#ifndef NDEBUG

cst_pSparse 		cst_SparseBlock_get_parent	(cst_pSparseBlock const self_) 	
{
  return self_->self_parent;
}

pSparse 		SparseBlock_get_parent		(pSparseBlock const self_) 	
{
  return self_->self_parent;
}

I 			SparseBlock_get_ni		(cst_pSparseBlock const self_) 	
{
  return self_->n_i;
}

I 			SparseBlock_get_nj		(cst_pSparseBlock const self_) 	
{
  return self_->n_j;
}

I 			SparseBlock_get_starti		(cst_pSparseBlock const self_) 	
{
  return self_->start_i;
}

I 			SparseBlock_get_startj		(cst_pSparseBlock const self_) 	
{
  return self_->start_i;
}

cst_pI 			cst_SparseBlock_get_bn		(cst_pSparseBlock const self_) 	
{
  return self_->bn;
}

cst_pI 			cst_SparseBlock_get_b		(cst_pSparseBlock const self_) 	
{
  return self_->b;
}

pI 			SparseBlock_get_bn		(pSparseBlock const self_) 	
{
  return self_->bn;
}

pI 			SparseBlock_get_b		(pSparseBlock const self_) 	
{
  return self_->b;
}

#endif


void SparseBlock_clr(pSparseBlock const self_)
{
  memset(self_,0,sizeof(SparseBlock));
}


void SparseBlock_free(pSparseBlock const self_)
{
  if (self_)
    {
      if (self_->b)
	{
	  free(self_->b);
	}
      if (self_->bn)
	{
	  free(self_->bn);
	}
      SparseBlock_clr(self_);
    }
}

void SparseBlock_clear(pSparseBlock const self_)
{
#if __ns_debug__
  wrong_param_nsSPARSE_BLOCK(! self_,1,SparseBlock_clear);
#endif
  pSparse A_ = SparseBlock_get_parent(self_);
  pR 	  Ax = Sparse_get_ownx(A_);
  { I i;
    for (i=0;i<self_->n_i;++i)
      {
	const I k = self_->b[i];
	{ I j;
	  for (j=0;j<self_->bn[i];++j)
	    {
	      Ax[k+j]=((R)0.0);
	    } }
      } }  
}


void SparseBlock_def_thread(pSparseBlock const	self_,
			    pSparse const	A_,
			    cst_pI		start_i_,
			    cst_pI		ni_,
			    cst_pI		start_j_,
			    cst_pI		nj_,
			    const I 		iproc_)
{
#ifndef NDEBUG
#if 0
  wrong_param_nsSPARSE_BLOCK(! self_,1,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! A_,2,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! start_i_,3,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_i_[0]+1<1,3,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_i_[0]+1>A_->n,3,_SparseBlock_def);

  wrong_param_nsSPARSE_BLOCK(! ni_,4,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! start_j_,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_j_[0]+1<1,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_j_[0]+1>A_->m,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! nj_,6,_SparseBlock_def);
#endif
#endif
  self_->self_parent 		= A_;
  self_->start_i 		= start_i_[0];
  self_->n_i 			= ni_[0];
  self_->start_j 		= start_j_[0];
  self_->n_j 			= nj_[0];
  self_->b 			= (pI)calloc((ni_[0]+1),sizeof(I));
  self_->bn 			= (pI)calloc((ni_[0]+1),sizeof(I));
  
  const L A_fortran_index 	= Sparse_is_fortran_index(A_);
  cst_pI Ai 			= Sparse_get_i(A_);
  cst_pI Ab 			= Sparse_get_b(A_);

  if (A_fortran_index)
    {
      { I i;
	for (i=start_i_[0];i<start_i_[0]+ni_[0];++i)
	  {
	    { I j;
	      for (j=Ab[i];j<Ab[i+1];++j)
		{
		  if (Ai[j-1]-1>=start_j_[0])
		    break;
		}
#if 0
	      /* 
		 mauvais test 
	      */
#if __ns_debug__
	      if (j>=Ab[i+1])
		{
		  fprintf(stderr,"SparseBlock_def:wrong algorithm line "ifmt" no element >="ifmt"\n",i,start_j_[0]);
		  exit(1);
		}
#endif
#endif
	      self_->b[i-start_i_[0]]=j-1;
	      for (;j<Ab[i+1];++j)
		{
		  if (Ai[j-1]-1>=nj_[0]+start_j_[0])
		    break;
		}
	      self_->bn[i-start_i_[0]]=(j-1)-self_->b[i-start_i_[0]]; }
	  } }
    }
  else
    {
      { I i;
	for (i=start_i_[0];i<start_i_[0]+ni_[0];++i)
	  {
	    { I j;
	      for (j=Ab[i];j<Ab[i+1];++j)
		{
		  if (Ai[j]>=start_j_[0])
		    break;
		}
#ifndef NDEBUG
	      if (j>=Ab[i+1])
		{
		  fprintf(stderr,"SparseBlock_def:wrong algorithm\n");
		  exit(1);
		}
#endif
	      self_->b[i-start_i_[0]]=j;
	      for (;j<Ab[i+1];++j)
		{
		  if (Ai[j]>=nj_[0]+start_j_[0])
		    break;
		}
	      self_->bn[i-start_i_[0]]=j-self_->b[i-start_i_[0]]; }
	  } }
    }
}






void SparseBlock_def(pSparseBlock const	self_,
		     pSparse const	A_,
		     cst_pI 		start_i_,
		     cst_pI 		ni_,
		     cst_pI 		start_j_,
		     cst_pI 		nj_)
{
#ifndef NDEBUG
#if 0
  wrong_param_nsSPARSE_BLOCK(! self_,1,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! A_,2,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! start_i_,3,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_i_[0]+1<1,3,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_i_[0]+1>A_->n,3,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! ni_,4,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! start_j_,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_j_[0]+1<1,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(start_j_[0]+1>A_->m,5,_SparseBlock_def);
  wrong_param_nsSPARSE_BLOCK(! nj_,6,_SparseBlock_def);
#endif
#endif
  self_->self_parent 	= A_;
  self_->start_i 	= start_i_[0];
  self_->n_i 		= ni_[0];
  self_->start_j 	= start_j_[0];
  self_->n_j 		= nj_[0];
  self_->b 		= (pI)calloc((ni_[0]+1),sizeof(I));
  self_->bn 		= (pI)calloc((ni_[0]+1),sizeof(I));

  const L A_fortran_index 	= Sparse_is_fortran_index(A_);
  cst_pI Ai 			= Sparse_get_i(A_);
  cst_pI Ab 			= Sparse_get_b(A_);


  if (A_fortran_index)
    {
      { I i;
	for (i=start_i_[0];i<start_i_[0]+ni_[0];++i)
	  {
	    { I j;
	      for (j=Ab[i];j<Ab[i+1];++j)
		{
		  if (Ai[j-1]-1>=start_j_[0])
		    break;
		}
#if 0
	      /* 
		 mauvais test 
	      */
#if __ns_debug__
	      if (j>=Ab[i+1])
		{
		  fprintf(stderr,"SparseBlock_def:wrong algorithm line "ifmt" no element >="ifmt"\n",i,start_j_[0]);
		  exit(1);
		}
#endif
#endif
	      self_->b[i-start_i_[0]]=j-1;
	      for (;j<Ab[i+1];++j)
		{
		  if (Ai[j-1]-1>=nj_[0]+start_j_[0])
		    break;
		}
	      self_->bn[i-start_i_[0]]=(j-1)-self_->b[i-start_i_[0]];
#if 0
	      if (self_->bn[i-start_i_[0]]==0)
		{
		  printf("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee starti = %lld, startj %lld   %lld/%lld  %lld\n",start_i_[0],start_j_[0],i-start_i_[0],self_->n_i,self_->bn[i-start_i_[0]]);
		  { I j;
		    for (j=Ab[i];j<Ab[i+1];++j)
		      {
			printf("%lld\n",Ai[j-1]-1);
		      } }
		  exit(1);
		}
#endif
	    }
	  } }
    }
  else
    {
      { I i;
	for (i=start_i_[0];i<start_i_[0]+ni_[0];++i)
	  {
	    { I j;
	      for (j=Ab[i];j<Ab[i+1];++j)
		{
		  if (Ai[j]>=start_j_[0])
		    break;
		}
#ifndef NDEBUG
	      if (j>=Ab[i+1])
		{
		  fprintf(stderr,"SparseBlock_def:wrong algorithm\n");
		  exit(1);
		}
#endif
	      self_->b[i-start_i_[0]]=j;
	      for (;j<Ab[i+1];++j)
		{
		  if (Ai[j]>=nj_[0]+start_j_[0])
		    break;
		}
	      self_->bn[i-start_i_[0]]=j-self_->b[i-start_i_[0]]; }
	  } }
    }
}


void  SparseBlock_ass_blank(pSparseBlock const 	block_,
			    cst_pI		ni_,
			    cst_pI		ddli_,
			    cst_pI		nj_,
			    cst_pI		ddlj_,
			    cst_pR		locmat_,
			    cst_pI		locmatoff_,
			    pI			blank_)
{
  pSparse A_ = SparseBlock_get_parent(block_);
  const L A_fortran_index 	= Sparse_is_fortran_index(A_);
  cst_pI 	Ai		= Sparse_get_i(A_);
  pR 		Ax		= Sparse_get_ownx(A_);

#if 0
  cst_pI Ab 			= Sparse_get_b(A_);
#if __ns_debug__
  ns_sparse_verif_lexicorder(A_);
#endif
#endif
#if 0
  {  I i;
  for (i=0;i<4;++i)
    {
      
      printf("ffffff %lld\n",block_->bn[i]);
    }}

#endif
  { I i,j;
    for (i=0;i<ni_[0];++i)
      {
	const I q = ddli_[i]-block_->start_i;	
	{ I k = block_->b[q];
	  for (j=0;j<block_->bn[q];++j)
	    {
	      if (A_fortran_index)
		{
		  blank_[Ai[k+j]-1]=k+j+1;
		}
	      else
		{
		  blank_[Ai[k+j]]=k+j+1;
		}
	    } }	
	for (j=0;j<nj_[0];++j)
	  {
	    const I jddl 	= ddlj_[j];	    
	    const I jddl_at = blank_[jddl]-1;

#if __ns_debug__
	    if (jddl_at<0) {printf("SparseBlock_ass_blank erreur jddl = "ifmt" blank = "ifmt" "ifmt" "ifmt"\n",jddl,blank_[jddl],block_->bn[q],q);
	      I k;
	      for (k=0;k<nj_[0];++k)
		{
		  printf("ddlj "ifmt"\n",ddlj_[k]);
		}
	      for (k=0;k<ni_[0];++k)
		{
		  printf("ddli "ifmt"\n",ddli_[k]);
		}
	      exit(1);}
#endif

	    Ax[jddl_at] += locmat_[j*locmatoff_[0]+i];
	  }	
	{ I k = block_->b[q];
	  for (j=0;j<block_->bn[q];++j)
	    {
	      if (A_fortran_index)
		blank_[Ai[k+j]-1]=(I)0;
	      else
		blank_[Ai[k+j]]=(I)0;
	    } }
      } }
}




void SparseBlock_gemv(cst_pR 			a_,
		      cst_pSparseBlock const 	block_,
		      cst_pR 			global_rhs_,
		      cst_pR 			b_,
		      pR 			global_sol_)
{
#ifndef NDEBUG
  wrong_param_nsSPARSE_BLOCK(! a_,1,SparseBlock_gemv);
  wrong_param_nsSPARSE_BLOCK(! block_,2,SparseBlock_gemv);
  wrong_param_nsSPARSE_BLOCK(! block_->self_parent,2,SparseBlock_gemv);
  wrong_param_nsSPARSE_BLOCK(! global_rhs_,3,SparseBlock_gemv);
  wrong_param_nsSPARSE_BLOCK(! b_,4,SparseBlock_gemv);
  wrong_param_nsSPARSE_BLOCK(! global_sol_,5,SparseBlock_gemv);
#endif
  cst_pSparse A_ = cst_SparseBlock_get_parent(block_);
  const L A_fortran_index 	= Sparse_is_fortran_index(A_);
  cst_pR	Ax		= Sparse_get_x(A_);
  cst_pI	Ai		= Sparse_get_i(A_);

  if (A_fortran_index)
    {
#if 0
      SparseBlock_gemv_threading(a_,
				     block_,
				     global_rhs_,
				     b_,
				     global_sol_);
      return;
#endif
      if (b_[0]==((R)0.0))
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*global_rhs_[Ai[k+j]-1];
		  }
		global_sol_[block_->start_i+i] = a_[0] * s;
	      } }
	}
      else
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*global_rhs_[Ai[k+j]-1];
		  }
		global_sol_[block_->start_i+i] = global_sol_[block_->start_i+i] * b_[0] + a_[0] * s;
	      } }
	}
    }
  else
    {
      if (b_[0]==((R)0.0))
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*global_rhs_[Ai[k+j]];
		  }
		global_sol_[block_->start_i+i] = a_[0] * s;
	      } }
	}
      else
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*global_rhs_[Ai[k+j]];
		  }
		global_sol_[block_->start_i+i] = global_sol_[block_->start_i+i] * b_[0] + a_[0] * s;
	      } }
	}
    }
}














void SparseBlock_gemv_sub(cst_pR 			a_,
			  cst_pSparseBlock const 	block_,
			  cst_pR 			rhs_,
			  cst_pR 			b_,
			  pR 				sol_)
{
#ifndef NDEBUG
  wrong_param_nsSPARSE_BLOCK(! a_,1,SparseBlock_gemv_sub);
  wrong_param_nsSPARSE_BLOCK(! block_,2,SparseBlock_gemv_sub);
  wrong_param_nsSPARSE_BLOCK(! block_->self_parent,2,SparseBlock_gemv_sub);
  wrong_param_nsSPARSE_BLOCK(! rhs_,3,SparseBlock_gemv_sub);
  wrong_param_nsSPARSE_BLOCK(! b_,4,SparseBlock_gemv_sub);
  wrong_param_nsSPARSE_BLOCK(! sol_,5,SparseBlock_gemv_sub);
#endif
  cst_pSparse 	A_ 		= cst_SparseBlock_get_parent(block_);
  const L 	A_fortran_index = Sparse_is_fortran_index(A_);
  cst_pR	Ax		= Sparse_get_x(A_);
  cst_pI	Ai		= Sparse_get_i(A_);

  if (A_fortran_index)
    {
#if 0
      SparseBlock_gemv_threading(a_,
				     block_,
				     global_rhs_,
				     b_,
				     global_sol_);
      return;
#endif
      if (b_[0]==((R)0.0))
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*rhs_[Ai[k+j]-1-block_->start_j];
		  }
		sol_[i] = a_[0] * s;
	      } }
	}
      else
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*rhs_[Ai[k+j]-1-block_->start_j];
		  }
		sol_[i] = sol_[i] * b_[0] + a_[0] * s;
	      } }
	}
    }
  else 
    {
      if (b_[0]==((R)0.0))
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*rhs_[Ai[k+j]-block_->start_j];
		  }
		sol_[i] = a_[0] * s;
	      } }
	}
      else
	{
	  { I i,j;
	    R s;
	    for (i=0;i<block_->n_i;++i)
	      {
		I k = block_->b[i];
		s = ((R)0.0);
		for (j=0;j<block_->bn[i];++j)
		  {
		    s+=Ax[k+j]*rhs_[Ai[k+j]-block_->start_j];
		  }
		sol_[i] = sol_[i] * b_[0] + a_[0] * s;
	      } }
	}
    }
}
