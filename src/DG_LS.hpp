#pragma once
#include "Direct/MKL/Pardiso.hpp"


struct DG_LS
{

  WLS::Direct::MKL::Pardiso * m_pardiso {};

  ~DG_LS()
  {
    if (this->m_pardiso)
      {
	delete this->m_pardiso;
	this->m_pardiso = nullptr;
      }
  };
  
  DG_LS(DG_JACOBIAN&J)
  {
  };

  void presolve(DG_JACOBIAN&J)
  {
    //
    // Call pardiso.
    //
    if (!this->m_pardiso)
      {
	std::cout << "precompute: symbolic phase ... " << std::endl;
	std::cout << "            n  = " << J.m_n << std::endl;
	std::cout << "            nc = " << J.m_nc << std::endl;
	this->m_pardiso = new WLS::Direct::MKL::Pardiso (J.m_n,
							 J.m_begin,
							 J.m_index,
							 J.m_values,
							 false);
      }
  };
  
  void solve(DG_JACOBIAN&J,
	     DG_VAR&F,
	     DG_VAR&Residual,
	     DG_VAR&E,
	     DG_VAR&LR)
  {
    bool hasFailed;
    m_pardiso->Compute(&hasFailed);
    if (hasFailed)
      {
	std::cerr << "compute failed: " << m_pardiso->GetErrorMessage()  << std::endl;
	exit(1);
      }
    
    {
      double* temporaryVector = nullptr;
      WLS::integer_t sizeOfTemporaryVector = m_pardiso->GetSizeOfTemporaryVector();
      if (sizeOfTemporaryVector>0)
	{
	  temporaryVector = new double[sizeOfTemporaryVector];
	}
    
      m_pardiso->Apply("No transpose",
		       E.matrix().x,
		       Residual.matrix().x,
		       sizeOfTemporaryVector,
		       temporaryVector,
		       &hasFailed);
    if (hasFailed)
      {
	std::cerr << "apply failed: " << m_pardiso->GetErrorMessage()  << std::endl;
	exit(1);
      }
      
      delete [] temporaryVector;
    }
    
  }
  
};
