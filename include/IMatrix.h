#ifndef __header_IMatrix_h__
#define __header_IMatrix_h__

#include "eTranspose.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    I 		m_numRows;
    I 		m_numCols;
    I 		m_offset;
    cst_pI 	m_ix;
    pI 		m_own_ix;
  } IMatrix,*RESTRICT pIMatrix;

  typedef const IMatrix * RESTRICT pReadOnlyIMatrix;

  I 		ReadOnlyIMatrix_get_numRows	(pReadOnlyIMatrix	const 	self_);
  I 		ReadOnlyIMatrix_get_numCols	(pReadOnlyIMatrix	const 	self_);
  I 		ReadOnlyIMatrix_get_offset		(pReadOnlyIMatrix	const 	self_);

  I		ReadOnlyIMatrix_get		(pReadOnlyIMatrix	const 	self_,
						 const I			rowIndex_,
						 const I			colIndex_);

  cst_pI	ReadOnlyIMatrix_get_col		(pReadOnlyIMatrix	const 	self_,
						 const I			colIndex_);


  void 		ReadOnlyIMatrix_copy_col	(pReadOnlyIMatrix 	const 	self_,
						 const I 			colIndex_,
						 pI 			const	col_,
						 const I 			colOffset_);

  void 		ReadOnlyIMatrix_copy_row	(pReadOnlyIMatrix 	const 	self_,
						 const I 			rowIndex_,
						 pI 			const	row_,
						 const I 			rowOffset_);

  pI   		IMatrix_get_col			(pIMatrix		const 	self_,
						 const I			colIndex_);

  void 		IMatrix_clear			(pIMatrix 		const 	self_);
  void 		IMatrix_free			(pIMatrix 		const 	self_);

  void 		IMatrix_def			(pIMatrix 		const 	self_,
						 const I 			numRows_,
						 const I 			numCols_,
						 cst_pI 		const	ix_,
						 const I 			ixOffset_);
  
  pIMatrix 	IMatrix_kill			(pIMatrix 		const 	self_);
  
  pIMatrix 	IMatrix_new			(const I 			numRows_,
						 const I 			numCols_,
						 cst_pI 		const	ix_,
						 const I 			off_);

  pIMatrix 	IMatrix_malloc			(const I 			numRows_,
						 const I 			numCols_);


#ifdef __cplusplus
}
#endif

#endif
