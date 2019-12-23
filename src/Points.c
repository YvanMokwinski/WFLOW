#include "Points.h"
#include <stdio.h>
#include <stdlib.h>


pPoints Points_kill(pPoints const self_)
{
  if (self_)
    {
#if 0
      self_->topologyIds 	= IMatrix_kill(self_->topologyIds);
#endif
      self_->xyz 		= Matrix_kill(self_->xyz);
      self_->dim 		= __eDim_ERROR;

      free(self_);
    }
  return NULL;
}

pPoints Points_new(cst_eDim 	dim_,
		   const I	nb_)
{
#ifndef NDEBUG
  DebugVerif(nb_>0);
  DebugVerif( (dim_>__eDim_ERROR) AND (dim_<__eDim_ALL) );
#endif
  pPoints self 		= (pPoints)calloc(1,sizeof(Points));
  self->xyz 		= Matrix_malloc(dim_,nb_);
#if 0
  self->topologyIds 	= IMatrix_malloc(1,nb_);
#endif
  self->dim 		= dim_;
  return self;
}

I 	Points_get_n		(cst_pPoints	const 	self_) 
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return Matrix_get_ncol(self_->xyz);
}


eDim 	Points_get_dim		(cst_pPoints	const 	self_) 
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->dim;
}


void 	Points_get		(cst_pPoints	const 	self_,
				 const I 		ith_,
				 pR		const 	x_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  Matrix_copy_col(self_->xyz,
		  ith_,
		  x_,
		  1);
}


void 	Points_set		(pPoints	const 	self_,
				 const I 		ith_,
				 cst_pR		const 	x_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  Matrix_set_col(self_->xyz,
		  ith_,
		  x_);
}

void 	Points_fprintf_ascii	(cst_pPoints const 	self_,
				 FILE * 		out_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(out_);
#endif
  Matrix_fprintf_ascii(self_->xyz,
		       Transpose,
		       out_);
}
