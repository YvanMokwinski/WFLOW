#ifndef __header_ILinearSystem_hpp__
#define __header_ILinearSystem_hpp__

#include "Sparse.h"

class ILinearSystem
{
  
public:

  virtual void	Precompute	(pSparse 	S_) = 0;
  virtual void	Compute		(pSparse 	S_,
				 cst_pR		y_,
				 pR 		x_) = 0;

};

#endif
