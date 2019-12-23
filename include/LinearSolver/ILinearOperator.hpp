#ifndef __header_ILinearOperator_hpp__
#define __header_ILinearOperator_hpp__

namespace LinearSolver
{
  class ILinearOperator
  {
  public:
    virtual void Apply(const char * transpose_,
		       pR y_,
		       cst_pR x_) = 0;
  };
};

#endif
