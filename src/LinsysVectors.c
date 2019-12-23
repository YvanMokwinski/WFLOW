#include "LinsysVectors.h"
#include "Blas.h"

#include <stdio.h>
#include <string.h>
#ifndef NDEBUG
#define wrong_param_Linsys(_cond,_i,_s) if ( (_cond) ) { fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); exit(1);} ((void)0)
#endif

I	LinsysVectors_get_n		(cst_pLinsysVectors const 	self_)
{
  return self_->m_N;
}


void * RESTRICT		LinsysVectors_get_friend_corr	(pLinsysVectors const 		self_) 
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_friend_corr);
#endif
  return self_->m_c;
}

const void * RESTRICT	LinsysVectors_get_corr		(cst_pLinsysVectors const 	self_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_corr);
#endif
  return self_->m_c;
}



void * RESTRICT		LinsysVectors_get_friend_grad	(pLinsysVectors const 		self_) 
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_friend_grad);
#endif
  return self_->m_grad;
}

const void * RESTRICT	LinsysVectors_get_grad		(cst_pLinsysVectors const 	self_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_grad);
#endif
  return self_->m_grad;
}


void * RESTRICT		LinsysVectors_get_friend_rhs	(pLinsysVectors const 		self_) 
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_friend_rhs);
#endif
  return self_->m_rhs;
}

const void * RESTRICT	LinsysVectors_get_rhs		(cst_pLinsysVectors const 	self_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_rhs);
#endif
  return self_->m_rhs;
}


void * RESTRICT		LinsysVectors_get_friend_x	(pLinsysVectors const 		self_) 
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_friend_x);
#endif
  return self_->m_x;
}

const void * RESTRICT	LinsysVectors_get_x		(cst_pLinsysVectors const 	self_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_x);
#endif
  return self_->m_x;
}

void * RESTRICT		LinsysVectors_get_friend_xi	(pLinsysVectors const 		self_,
							 cst_pI				istep_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_friend_xi);
  wrong_param_Linsys(!istep_,2,LinsysVectors_get_friend_xi);
  wrong_param_Linsys(istep_[0]+1<1,2,LinsysVectors_get_friend_xi);
  wrong_param_Linsys(istep_[0]>=self_->m_nstep,2,LinsysVectors_get_friend_xi);
#endif
  return &self_->m_x[ self_->m_N * istep_[0] ];
}

const void * RESTRICT	LinsysVectors_get_xi		(cst_pLinsysVectors const 	self_,
							 cst_pI 			istep_)
{ 
#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_get_xi);
  wrong_param_Linsys(!istep_,2,LinsysVectors_get_xi);
  wrong_param_Linsys(istep_[0]>0,2,LinsysVectors_get_xi);
  wrong_param_Linsys(istep_[0]<self_->m_nstep,2,LinsysVectors_get_xi);
#endif
  return &self_->m_x[ self_->m_N * istep_[0] ];
}



pLinsysVectors 	LinsysVectors_kill		(pLinsysVectors const 	self_)
{
  if (self_)
    {
      if (self_->m_x) 	free(self_->m_x);
      if (self_->m_rhs) 	free(self_->m_rhs);
      memset(self_,0,sizeof(LinsysVectors));
      free(self_);
    }
  return NULL;
}

pLinsysVectors 	LinsysVectors_new		(cst_pI N_,
						 cst_pI nstep_)
{

#ifndef NDEBUG
  wrong_param_Linsys(N_==NULL,1,LinsysVectors_new);
#endif

  pLinsysVectors self = (pLinsysVectors)calloc(1,sizeof(LinsysVectors));

  self->m_N 	= N_[0];
  self->m_nstep 	= nstep_[0];

  self->m_x 	= (pR)malloc(self->m_nstep*self->m_N*sizeof(R));
  self->m_rhs 	= (pR)malloc(3*self->m_N*sizeof(R));
  self->m_c 	= &self->m_rhs[self->m_N];
  self->m_grad 	= &self->m_rhs[self->m_N*2];

  return self;
}


void LinsysVectors_update(pLinsysVectors const 		self_,
			  cst_pI			nstep_,
			  cst_pI			start_,
			  cst_pI			start_N_)
{

#ifndef NDEBUG
  wrong_param_Linsys(!self_,1,LinsysVectors_update);
  wrong_param_Linsys(!nstep_,2,LinsysVectors_update);
  wrong_param_Linsys(!start_,3,LinsysVectors_update);
  wrong_param_Linsys(!start_N_,4,LinsysVectors_update);
#endif

  static const I nequal1	= ((I)1);	
  const I N 			= self_->m_N;
  const I nstep 		= nstep_[0];

  { I istep;
    for (istep=1;istep<nstep;++istep)
      {
	Blas_dcopy(start_N_,
		   &self_->m_x[start_[0] +  N * (nstep - 1 - istep) ],
		   &nequal1,
		   &self_->m_x[start_[0] + N * (nstep - istep) ],
		   &nequal1);
      } }

}

