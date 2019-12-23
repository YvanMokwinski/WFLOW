#ifndef __header_Variables_hpp__
#define __header_Variables_hpp__
#include "eVariable.h"


class Variables
{
protected:

  I m_numTimeSteps;
  pR * m_x[__eVariable_ALL];

public:

  inline pR 	Get(const eVariable 	variable_,
		    const I 		timeStep_)
  {
    return m_x[variable_][timeStep_];    
  };
  

  inline cst_pR 	Get(const eVariable 	variable_,
			    const I 		timeStep_) const
  {
    return m_x[variable_][timeStep_];    
  };
  
  inline void Clear()
  {
    { eVariable variable = __eVariable_ERROR;  
      for (++variable;variable<__eVariable_ALL;++variable)
	{
	  for (I timeStep=0;timeStep < m_numTimeSteps;++timeStep)
	    {
	      this->m_x[variable][timeStep] = NULL;
	    }
	} }
  };

  inline void 	Set(const eVariable 	variable_,
		    const I 		timeStep_,
		    pR 			x_)
  {
    this->m_x[variable_][timeStep_] = x_;
  };

  virtual ~Variables()
  {
    { eVariable variable = __eVariable_ERROR;  
      for (++variable;variable<__eVariable_ALL;++variable)
	{
	  delete [] this->m_x[variable];
	  this->m_x[variable] = NULL;
	} }
  };
  
  Variables(const I 	bound_,
	    const I 	numTimeSteps_,
	    pR		x_[],
	    const I	nddls_	[__eVariable_ALL],
	    const I	shifts_	[__eVariable_ALL])
    : m_numTimeSteps(numTimeSteps_)
  {
    
    
    { eVariable variable = __eVariable_ERROR;  
      for (++variable;variable<__eVariable_ALL;++variable)
	{
	  if ( (shifts_[variable]<0) || (shifts_[variable]>=bound_) )
	    {
	      fprintf(stderr,"Variables_def invalid shifts_[%s] = " ifmt "",
		      eVariable_get_string(variable),
		      shifts_[variable]);
	      exit(1);
	    }
	} }
    
    { eVariable variable = __eVariable_ERROR;  
      for (++variable;variable<__eVariable_ALL;++variable)
	{
	  this->m_x[variable] = new pR[numTimeSteps_];
	} }
    
    { eVariable variable = __eVariable_ERROR;  
      for (++variable;variable<__eVariable_ALL;++variable)
	{
	  
	  for (I timeStepIndex=0;timeStepIndex<numTimeSteps_;++timeStepIndex)
	    {
	      this->Set(variable,
			timeStepIndex,
			x_[timeStepIndex] + shifts_[variable]);
	    }
	} }
  };
  
};

#endif

