#pragma once

class INonLinearSolverDriver
{
 protected:
  INonLinearSolverDriver() { };
public:
  
  virtual void * 	GetUsrPtr() = 0;
  
  virtual void 		BuildSystem(cst_pI 		nx_,
				    cst_pR 		x_,
				    cst_pI 		xoff_) = 0;
  
  virtual void 		BuildResidu(pR 			global_rhs_,
				    cst_pR 		x_,
				    cst_pR 		xi_,
				    cst_pR 		xii_) = 0;
  
};
