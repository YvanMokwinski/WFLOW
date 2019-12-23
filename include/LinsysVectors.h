#ifndef __header_LinsysVectors_h__
#define __header_LinsysVectors_h__

#include "Config.h"
#include "Type.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
     \brief Empty structure for LinsysVectors object
  */
  typedef struct
  {
    I 		m_N;
    I 		m_nstep;
    pR 		m_x;
    pR 		m_rhs;
    pR 		m_c;
    pR 		m_grad;
  } LinsysVectors,*RESTRICT pLinsysVectors;

  /**
     \brief Type a pointer to a constant LinsysVectors object
  */
  typedef const LinsysVectors * RESTRICT cst_pLinsysVectors;

  const void * RESTRICT		LinsysVectors_get_rhs		(cst_pLinsysVectors const 	self_);
  void * RESTRICT		LinsysVectors_get_friend_rhs	(pLinsysVectors const 		self_);

  const void * RESTRICT		LinsysVectors_get_corr		(cst_pLinsysVectors const 	self_);
  void * RESTRICT		LinsysVectors_get_friend_corr	(pLinsysVectors const 		self_);

  const void * RESTRICT		LinsysVectors_get_grad		(cst_pLinsysVectors const 	self_);
  void * RESTRICT		LinsysVectors_get_friend_grad	(pLinsysVectors const 		self_);

  /**
     \brief Get the dimension of the LinsysVectors object
     @param self_ self object
     @return dimension
  */
  I				LinsysVectors_get_n		(cst_pLinsysVectors const 	self_);

  /**
     \brief Kill a LinsysVectors object
     @param self_ self object
     @return NULL
  */
  pLinsysVectors 		LinsysVectors_kill		(pLinsysVectors const 		self_);

  /**
     \brief Create a new LinsysVectors object
     @param N_ length of the vectors
     @param nstep_ number of vectors to allocate
     @return the new object
  */
  pLinsysVectors 		LinsysVectors_new		(cst_pI				N_,
								 cst_pI				nstep_);


  /**
     \brief Update a LinsysVectors object
     @param self_ self object
     @param nstep_ how many Xi vectors to update
     @param start_ starting index for updating
     @param start_N_ length after the starting index for updating
  */
  void 				LinsysVectors_update		(pLinsysVectors const 		self_,
								 cst_pI				nstep_,
								 cst_pI				start_,
								 cst_pI				start_N_);
  
  /**
     \brief Get address of the values of the X vector
     @param self_ self object
     @return address of the values of the X vector
  */
  void * RESTRICT		LinsysVectors_get_friend_x	(pLinsysVectors const 		self_);

  /**
     \brief Get address of the constant values of the X vector
     @param self_ self object
     @return address of the constant values of the X vector
  */
  const void * RESTRICT	LinsysVectors_get_x		(cst_pLinsysVectors const 	self_);

  /**
     \brief Get address of the values of one of the Xi vectors
     @param self_ self object
     @param istep_ which Xi vector to select
     @return address of the constant values of the X vector
  */
  void * RESTRICT		LinsysVectors_get_friend_xi	(pLinsysVectors const 		self_,
								 cst_pI				istep_);

  /**
     \brief Get address of the constant values of one of the Xi vectors
     @param self_ self object
     @param istep_ which Xi vector to select
     @return address of the constant values of the X vector
  */
  const void * RESTRICT	LinsysVectors_get_xi		(cst_pLinsysVectors const 	self_,
							 cst_pI 			istep_);


#ifdef __cplusplus
}
#endif

#endif
