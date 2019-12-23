#ifndef __header_Sparse_h__
#define __header_Sparse_h__

#include "SparseIterative.h"
#include "eSparseFormat.h"
#ifdef __cplusplus
extern "C"
{
#endif

/**\brief struct to define a sparse matrix */
typedef struct
{
  /**\brief number of rows */
  I 		n;
  /**\brief number of columns */
  I 		m;
  /**\brief number of coefficient */
  I 		nc;
  /**\brief row bandwith */
  I 		band_row;
  /**\brief col bandwith */
  I 		band_col;
  /**\brief format */
  I 		format;
  /**\brief begin index */
  cst_pI   	b;
  /**\brief indices */
  cst_pI   	i;
  /**\brief values */
  cst_pR   	x;
  /**\brief (modifiable) values */
  pR 		own_x;
  /**\brief (modifiable) indices */
  pI		own_i;
  /**\brief (modifiable) begin */
  pI  		own_b;
} Sparse, *RESTRICT pSparse;


typedef const Sparse * 	RESTRICT cst_pSparse;


void 		Sparse_gemv			(cst_pR 	const	a_,
						 cst_pSparse 	const	self_,
						 cst_pR 	const	rhs_,
						 cst_pR 	const	b_,
						 pR 		const	sol_);

pSparse 	Sparse_build			(const I 		n_,
						 const I 		m,
						 const I 		nc_,
						 pI 		const	own_b_,
						 pI 		const	own_i_,
						 pR 		const	own_x_);


pSparse 	Sparse_buildReference		(const I 		n_,
						 const I 		m,
						 const I 		nc_,
						 cst_pI		const	b_,
						 cst_pI		const	i_,
						 cst_pR 	const	x_);


pSparse 	Sparse_clone			(cst_pSparse 	const 	self_);

/**
   \brief Create a new Sparse matrix
   @param nrow_ number of rows
   @param ncol_ number of columns
   @param ncoeffs_ number of coefficients
*/
pSparse		Sparse_new			(const I 		nrow_,
						 const I 		ncol_,
						 const I 		ncoeffs_);


/**
   \brief Kill a Sparse matrix
   @param self_ self matrix
*/
pSparse 	Sparse_kill			(pSparse 	const self_);

/**
   \brief Is a fortran indexation
   @param self_ self matrix
   @return logical value
*/
L 		Sparse_is_fortran_index		(cst_pSparse const self_); 	

/**
   \brief Get the number of rows
   @param self_ self matrix
   @return number of rows
*/
I 		Sparse_get_n			(cst_pSparse 	const self_);

/**
   \brief Get the number of columns
   @param self_ self matrix
   @return number of columns
*/
I 		Sparse_get_m			(cst_pSparse 	const self_);

/**
   \brief Get the number of coefficients
   @param self_ self matrix
   @return the number of coefficients
*/
I 		Sparse_get_nc			(cst_pSparse 	const self_);

/**
   \brief Get pointer to constant index
   @param self_ self matrix
   @return pointer to index
*/
cst_pI 		Sparse_get_i			(cst_pSparse 	const self_);


/**
   \brief Get pointer to constant begin
   @param self_ self matrix
   @return pointer to begin
*/
cst_pI 		Sparse_get_b			(cst_pSparse 	const self_);

/**
   \brief Get address of constant values
   @param self_ self matrix
   @return pointer to constant values
*/
cst_pR 		Sparse_get_x			(cst_pSparse 	const self_);

/**
   \brief Get begin index
   @param self_ self matrix
   @return pointer to begin index
*/
pI 		Sparse_get_owni			(pSparse 	const self_);


/**
   \brief Get begin array
   @param self_ self matrix
   @return pointer to begin array
*/
pI 		Sparse_get_ownb			(pSparse 	const self_);

/**
   \brief Get address of values
   @param self_ self matrix
   @return pointer to values
*/
pR 		Sparse_get_ownx			(pSparse 	const self_);

/**
   \brief Clear the Sparse matrix
   @param self_ self matrix
*/
void 		Sparse_clr			(pSparse 	const self_);

/**
   \brief Clear values of the matrix
   @param self_ self matrix
*/
void 		Sparse_clear			(pSparse 	const self_);

/**
   \brief Convert coefficient to a Fortran indexation
   @param self_ self matrix
*/
void 		Sparse_fortran_indexation	(pSparse 	const self_);

/**
   \brief Convert coefficient to a C indexation
   @param self_ self matrix
*/
void 		Sparse_c_indexation		(pSparse 	const self_);

/**
   \brief Sort coefficient of a Sparse matrix
   @param self_ self matrix
*/
void 		Sparse_sort			(pSparse 	const self_);

/**
   \brief Write the structure associated to a Sparse matrix
   @param self_ self matrix
   @param filename_ name of the file
*/
void 		Sparse_spy			(cst_pSparse 	const 	self_,
						 cst_pS 		filename_,
						 ...);
/**
   \brief Write a Sparse matrix to file
   @param self_ self matrix
   @param filename_ name of the file
*/
void 		Sparse_write			(cst_pSparse 	const 	self_,
						 cst_pS			filename_,
						 ...);

void 		Sparse_dirichlet		(pSparse 	const 	self_,
						 const I 		n_,
						 cst_pI 		dir_,
						 const I 		dir_dec_);

  void 		Sparse_ass			(pSparse const	self_,
						 const I	i_,
						 const I	j_,
						 const R	a);

  
void       	Sparse_extractDiagonal		(cst_pSparse const 	self_,
						 pR 			values_);

  
#ifdef __cplusplus
}
#endif

#endif
