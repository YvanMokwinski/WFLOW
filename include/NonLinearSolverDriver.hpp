#ifndef __header_NonLinearSolverDriver_hpp__
#define __header_NonLinearSolverDriver_hpp__
#error "should not be included"
#include "nsGLOBAL.h"
#include "Monitor.h"
#include "NonLinearSolverNewton.hpp"
#include "INonLinearSolverDriver.hpp"

class NonLinearSolverDriver : public INonLinearSolverDriver
{
 protected:
  void* m_usrptr;
public:
  
  NonLinearSolverDriver(void* usrptr_)
  {
    this->m_usrptr = usrptr_;
  };
  
  virtual void * GetUsrPtr()
  {
    return this->m_usrptr;
  };
  
  virtual void BuildSystem(cst_pI 		nx_,
			   cst_pR 		x_,
			   cst_pI 		xoff_)
  {
    ns_build_system((pGlobal)this->m_usrptr,
		    nx_,
		    x_,
		    xoff_);
  };
  
  virtual void BuildResidu(pR 			global_rhs_,
			   cst_pR 		x_,
			   cst_pR 		xi_,
			   cst_pR 		xii_)
  {
    ns_build_residu((pGlobal)this->m_usrptr,
		    global_rhs_,
		    x_,
		    xi_,
		    xii_);
  };
  
};

#endif
