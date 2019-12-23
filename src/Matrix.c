#include "Matrix.h"
#include "Blas.h"
#include <stdio.h>
#include <stdlib.h>

void Matrix_copy_row(cst_pMatrix 	const 	self_,
		     const I 			irow_,
		     pR 		const	row_,
		     const I 			rowoff_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(irow_<self_->nrow);
  DebugVerif(irow_+1>0);
  DebugVerif(rowoff_>=1);
#endif  
  const I ncol = self_->ncol;  
  { I i;
    for (i=0;i<ncol;++i)
      {
	row_[rowoff_*i] = self_->x[self_->off*i+irow_];
      } }
}

void Matrix_copy_col(cst_pMatrix 	const 	self_,
		     const I 			icol_,
		     pR 		const	col_,
		     const I 			coloff_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(icol_<self_->ncol);
  DebugVerif(icol_+1>0);
  DebugVerif(coloff_>=1);
#endif  
  const I nrow = self_->nrow;  
  { I i;
    for (i=0;i<nrow;++i)
      {
	col_[coloff_*i] = self_->x[self_->off*icol_+i];
      } }
}

void Matrix_set_col(pMatrix 	const 	self_,
		    const I 		icol_,
		    cst_pR 	const	col_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(icol_<self_->ncol);
  DebugVerif(icol_+1>0);
  DebugVerif(col_);
#endif  
  { I i;
    for (i=0;i<self_->nrow;++i)
      {
	self_->own_x[self_->off*icol_+i] = col_[i];
      } }
}


void 	Matrix_fprintf_ascii	(cst_pMatrix	const 	self_,
				 cst_eTranspose 	transpose_,
				 FILE * 		out_)
{
  const I ncol 	= self_->ncol;
  const I nrow 	= self_->nrow;
  const I off 	= self_->off;
  cst_pR  x 	= self_->x;
  switch(transpose_)
    {
    case NoTranspose:
      {
	{ I i;
	  for (i=0;i<nrow;++i)
	    {
	      { I j;
		for (j=0;j<ncol;++j)
		  {
		    fprintf(out_," "rfmt"",x[off*j+i]);
		  } }
	      fprintf(out_,"\n");
	    } }	
	break;
      }
    case Transpose:
      {
	{ I j;
	  for (j=0;j<ncol;++j)
	    {
	      { I i;
		for (i=0;i<nrow;++i)
		  {
		    fprintf(out_," "rfmt"",x[off*j+i]);
		  } }
	      fprintf(out_,"\n");
	    } }
	break;
      }
    case __eTranspose_ERROR:
    case __eTranspose_ALL:
      {
	break;
      }

    }

}


I 	Matrix_get_off	(cst_pMatrix	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->off;
}



I 	Matrix_get_nrow(cst_pMatrix	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->nrow;
}


I 	Matrix_get_ncol(cst_pMatrix	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->ncol;
}

pR   	Matrix_at		(pMatrix	const 	self_,
				 const I		irow_,
				 const I		icol_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->own_x[icol_*self_->off + irow_];
}


pR   Matrix_get_col(pMatrix		const 	self_,
		    const I 		icol_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->own_x[icol_*self_->off];
}

void Matrix_MM(cst_pMatrix 	const self_,
	       cst_pR		const xother_,
	       cst_pMatrix 	const other_,
	       cst_pR		const xresult_,
	       pMatrix 		const result_)
{  
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(xother_);
  DebugVerif(other_);
  DebugVerif(xresult_);
  DebugVerif(result_);
  DebugVerif(self_->ncol == other_->nrow);
  DebugVerif(self_->nrow == result_->nrow);
  DebugVerif(other_->ncol == result_->ncol);
#endif
  Blas_dgemm("N",
	     "N",
	     &result_->nrow,
	     &result_->ncol,
	     &self_->ncol,
	     xother_,
	     self_->x,
	     &self_->off,
	     other_->x,
	     &other_->off,
	     xresult_,
	     result_->x,
	     &result_->off);
}

pMatrix Matrix_kill(pMatrix const self_)
{
  if (self_)
    {
      if (self_->own_x)
	{
	  free(self_->own_x);
	}
      free(self_);
    }
  return NULL;
}

pMatrix Matrix_new(const I	nrow_,
		   const I 	ncol_,
		   cst_pR const	x_,
		   const I 	off_)
{
  pMatrix self 	= (pMatrix)calloc(1,sizeof(Matrix));
  self->x 	= x_;
  self->nrow 	= nrow_;
  self->ncol 	= ncol_;
  self->off 	= off_;
  return self;
}


pMatrix Matrix_malloc(const I	nrow_,
		      const I 	ncol_)
{
  pMatrix self 	= (pMatrix)calloc(1,sizeof(Matrix));
  self->nrow 	= nrow_;
  self->ncol 	= ncol_;
  self->own_x 	= (pR)malloc(sizeof(R)*nrow_*ncol_);
  self->x 	= self->own_x;
  self->off 	= nrow_;
  return self;
}





#if 0
void Matrix_MtMt(cst_pMatrix 	const self_,
		 cst_pR		const xother_,
		 cst_pMatrix 	const other_,
		 cst_pR		const xresult_,
		 pMatrix 		const result_)
{  
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(xother_);
  DebugVerif(other_);
  DebugVerif(xresult_);
  DebugVerif(result_);
  DebugVerif(self_->nrow == other_->ncol);
  DebugVerif(self_->ncol == result_->nrow);
  DebugVerif(other_->nrow == result_->ncol);
#endif
  Blas_dgemm("T",
	     "T",
	     &result_->nrow,
	     &result_->ncol,
	     &self_->nrow,
	     xother_,
	     self_->x,
	     &self_->off,
	     other_->x,
	     &other_->off,
	     xresult_,
	     result_->x,
	     &result_->off);
}

void Matrix_MtM(cst_pMatrix 	const self_,
	       cst_pR		const xother_,
	       cst_pMatrix 	const other_,
	       cst_pR		const xresult_,
	       pMatrix 		const result_)
{  
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(xother_);
  DebugVerif(other_);
  DebugVerif(xresult_);
  DebugVerif(result_);
  DebugVerif(self_->nrow == other_->nrow);
  DebugVerif(self_->ncol == result_->nrow);
  DebugVerif(other_->ncol == result_->ncol);
#endif
  Blas_dgemm("T",
	     "N",
	     &result_->nrow,
	     &result_->ncol,
	     &self_->nrow,
	     xother_,
	     self_->x,
	     &self_->off,
	     other_->x,
	     &other_->off,
	     xresult_,
	     result_->x,
	     &result_->off);
}


void Matrix_MMt(cst_pMatrix 	const self_,
	       cst_pR		const xother_,
	       cst_pMatrix 	const other_,
	       cst_pR		const xresult_,
	       pMatrix 		const result_)
{  
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(xother_);
  DebugVerif(other_);
  DebugVerif(xresult_);
  DebugVerif(result_);
  DebugVerif(self_->ncol == other_->ncol);
  DebugVerif(self_->nrow == result_->nrow);
  DebugVerif(other_->nrow == result_->ncol);
#endif
  Blas_dgemm("N",
	     "T",
	     &result_->nrow,
	     &result_->ncol,
	     &self_->ncol,
	     xother_,
	     self_->x,
	     &self_->off,
	     other_->x,
	     &other_->off,
	     xresult_,
	     result_->x,
	     &result_->off);
}

#endif
