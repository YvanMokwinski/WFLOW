#ifndef __header_Matrix_h__
#define __header_Matrix_h__

#include "eTranspose.h"
typedef struct
{
  I 		nrow;
  I 		ncol;
  I 		off;
  cst_pR 	x;
  pR 		own_x;
} Matrix,*RESTRICT pMatrix;

typedef const Matrix * RESTRICT cst_pMatrix;

void 	Matrix_fprintf_ascii	(cst_pMatrix	const 	self_,
				 cst_eTranspose 	tranpose_,
				 FILE * 		out_);

I 	Matrix_get_nrow		(cst_pMatrix	const 	self_);
I 	Matrix_get_ncol		(cst_pMatrix	const 	self_);
I 	Matrix_get_off		(cst_pMatrix	const 	self_);

pR   	Matrix_get_col		(pMatrix	const 	self_,
				 const I		icol_);

pR   	Matrix_at		(pMatrix	const 	self_,
				 const I		irow_,
				 const I		icol_);

void 	Matrix_copy_row		(cst_pMatrix 	const 	self_,
				 const I 		irow_,
				 pR 		const	row_,
				 const I 		rowoff_);

void 	Matrix_copy_col		(cst_pMatrix const 	self_,
				 const I 		icol_,
				 pR 			col_,
				 const I 		coloff_);

void 	Matrix_set_col		(pMatrix 	const 	self_,
				 const I 		icol_,
				 cst_pR 	const	col_);


pMatrix Matrix_kill		(pMatrix 	const self_);


pMatrix Matrix_new		(const I		nrow_,
				 const I 		ncol_,
				 cst_pR 	const	x_,
				 const I 		off_);

pMatrix Matrix_malloc		(const I		nrow_,
				 const I 		ncol_);

void 	Matrix_MtMt		(cst_pMatrix 	const self_,
				 cst_pR		const xother_,
				 cst_pMatrix 	const other_,
				 cst_pR		const xresult_,
				 pMatrix 	const result_);


void 	Matrix_MMt		(cst_pMatrix 	const self_,
				 cst_pR		const xother_,
				 cst_pMatrix 	const other_,
				 cst_pR		const xresult_,
				 pMatrix 	const result_);


void 	Matrix_MtM		(cst_pMatrix 	const self_,
				 cst_pR		const xother_,
				 cst_pMatrix 	const other_,
				 cst_pR		const xresult_,
				 pMatrix 	const result_);


void 	Matrix_MM		(cst_pMatrix 	const self_,
				 cst_pR		const xother_,
				 cst_pMatrix 	const other_,
				 cst_pR		const xresult_,
				 pMatrix 	const result_);


#endif
