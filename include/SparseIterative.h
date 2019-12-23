#ifndef __header_SparseIterative_h__
#define __header_SparseIterative_h__

#include "eSparseIterativeMethod.h"
#include "eSparsePreconditioningMethod.h"
#include "Err.h"
#ifdef __cplusplus
extern "C"
{
#endif

/**
   \brief Type header of SparseIterative driver
*/
typedef void (*__routine_SparseIterativeDriver)(void * 		usrptr_,
						const char * 	tr_,
						pR 		y_,
						cst_pR 		x_);

/**
   \brief Structure to define a SparseIterativeDriver
*/
typedef struct
{
  __routine_SparseIterativeDriver 	routine;
  void * 				usrptr;
} SparseIterativeDriver,* RESTRICT pSparseIterativeDriver;


/**
   \brief Type a pointer to a constant SparseIterativeDriver
*/
typedef const SparseIterativeDriver * RESTRICT 	cst_pSparseIterativeDriver;



/**
   \brief Define a SparseIterativeDriver
   @param self_ self SparseIterativeDriver
   @param routine_ driver
   @param usrptr_ external user pointer
*/
void 	SparseIterativeDriver_def		(pSparseIterativeDriver 		self_,
						 __routine_SparseIterativeDriver 	routine_,
						 void * 				usrptr_);

/**
   \brief Free a SparseIterativeDriver
   @param self_ self SparseIterativeDriver
*/
void 	SparseIterativeDriver_free		(pSparseIterativeDriver 		self_);



/**
   \brief Structure to define a SparseIterative object
*/
typedef struct
{
  eSparseIterativeMethod 	m_iterative_method;
  eSparsePreconditioningMethod 	m_preconditioning_method;
  I				m_sparsekit_ipar[64];
  R 				m_sparsekit_fpar[64];
} SparseIterative,*RESTRICT pSparseIterative;


/**
   \brief Type a pointer to a constant SparseIterative object
*/
typedef const SparseIterative * RESTRICT 	cst_pSparseIterative;


/**
   \brief Create a new SparseIterative object
   @param self_ self object
   @param msg_  error msg
*/
void 			SparseIterative_checkError	(cst_pSparseIterative 	self_,
							 STR 			msg_);



/**
   \brief Get the required memory
   @param self_ 	self object
   @param n_  		dimension
   @param rwork_n_  	required memory
   @param err_  	error
*/
void 			SparseIterative_get_memory	(cst_pSparseIterative	self_,
							 cst_pI			n_,
							 pI			rwork_n_,
							 pErr 			err_);


/**
   \brief Kill a SparseIterative object
   @param self_ self SparseIterativeDriver
*/
pSparseIterative 	SparseIterative_kill		(pSparseIterative			self_);


/**
   \brief Create a new SparseIterative object
   @param iterative_method_ iterative method to use
   @param preconditioning_method_ preconditioning method to use with the iterative method 
*/
pSparseIterative 	SparseIterative_new		(const eSparseIterativeMethod		iterative_method_,
							 const eSparsePreconditioningMethod	preconditioning_method_);



/**
   \brief Solve linear system
   @param self_ self
   @param n_ dimension
   @param y_ solution
   @param x_ right hand side
   @param A_ iterative driver for the matrix of the linear system
   @param P_ iterative driver for the preconditioner
   @param rwork_n_ size of rwork_ array
   @param rwork_ real working array
   @return ?
*/
R 			SparseIterative_compute		(pSparseIterative			self_,
							 cst_pI					n_,
							 pR 					y_,
							 cst_pR 				x_,
							 pSparseIterativeDriver			A_,
							 pSparseIterativeDriver			P_,
							 cst_pI					rwork_n_,
							 pR					rwork_);

#ifdef __cplusplus
}
#endif

#endif
