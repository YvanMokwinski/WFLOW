#ifndef __header_LinsysVectors_hpp__
#define __header_LinsysVectors_hpp__

#include "Config.h"
#include "Type.h"

#include "Blas.h"
#include <stdio.h>
#include <string.h>

class LinsysVectors
{

private:  

  void ClearMembers()
  {
    this->m_N 		= 0;
    this->m_nstep 	= 0;
    this->m_x 		= NULL;
    this->m_rhs 	= NULL;
    this->m_c 		= NULL;
    this->m_grad 	= NULL;
  };

protected:

  I 		m_N;
  I 		m_nstep;
  pR 		m_x;
  pR 		m_rhs;
  pR 		m_c;
  pR 		m_grad;
  
public:

  I GetN() const
  {
    return this->m_N;
  };
  
  pR GetCorr() 
  { 
    return this->m_c;
  };

  cst_pR GetCorr() const
  { 
    return this->m_c;
  };

  pR GetGrad() 
  { 
    return this->m_grad;
  };

  cst_pR GetGrad() const
  { 
    return this->m_grad;
  };
  
  pR GetRhs() 
  { 
    return this->m_rhs;
  };

  cst_pR GetRhs() const
  { 
    return this->m_rhs;
  };

  pR GetX() 
  { 
    return this->m_x;
  };

  cst_pR GetX() const
  { 
    return this->m_x;
  };

  cst_pR GetXi(cst_pI istep_) const
  {
    return &this->m_x[ this->m_N * istep_[0] ];
  };
  
  pR GetXi(cst_pI istep_) 
  {
    return &this->m_x[ this->m_N * istep_[0] ];
  };

  virtual ~LinsysVectors()
  {
    
    if (NULL != this->m_x)
      {
	free(this->m_x);
      }
    
    if (NULL != this->m_rhs)
      {
	free(this->m_rhs);
      }
    
    ClearMembers();
  };

  LinsysVectors(cst_pI N_,
		cst_pI nstep_)
  {
    ClearMembers();
    this->m_N 		= N_[0];
    this->m_nstep 	= nstep_[0];
    
    this->m_x 		= (pR)malloc(this->m_nstep*this->m_N*sizeof(R));
    this->m_rhs 	= (pR)malloc(3*this->m_N*sizeof(R));
    this->m_c 		= &this->m_rhs[this->m_N];
    this->m_grad 	= &this->m_rhs[this->m_N*2];
    
    return;
  };
  
  void Update(cst_pI			nstep_,
	      cst_pI			start_,
	      cst_pI			start_N_)
  {
    
    static const I nequal1	= ((I)1);	
    const I N 			= this->m_N;
    const I nstep 		= nstep_[0];
    
    { I istep;
      for (istep=1;istep<nstep;++istep)
	{
	  Blas_dcopy(start_N_,
		     &this->m_x[start_[0] +  N * (nstep - 1 - istep) ],
		     &nequal1,
		     &this->m_x[start_[0] + N * (nstep - istep) ],
		     &nequal1);
	} }
    
  };
  
};



#endif
