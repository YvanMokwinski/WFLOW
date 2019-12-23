#ifndef __header_LinearSolver_Iterative_MKL_Parameters_hpp__
#define __header_LinearSolver_Iterative_MKL_Parameters_hpp__

#include "LinearSolver/ILinearOperator.hpp"

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_rci.h"

namespace LinearSolver
{
  namespace Iterative
  {
    namespace MKL
    {

      /// <summary>
      /// Base class for Parameters of  MKL sparse iterative solvers.
      /// </summary>
      class Parameters
      {
	
      protected:
	
	/// <summary>
	/// The array of integer parameters.
	/// </summary>
	I m_ipar[128];
	
	/// <summary>
	/// The array of double parameters.
	/// </summary>
	R m_dpar[128];
        
      public:
	
	Parameters()
	{
	  for (int i=0;i<128;++i)
	    {
	      m_dpar[i] = ((R)0.0);
	    }
	  for (int i=0;i<128;++i)
	    {
	      m_ipar[i] = ((I)0);
	    }
	};

      public:
	
	virtual ~Parameters()
	{
	  for (int i=0;i<128;++i)
	    {
	      m_dpar[i] = ((R)0.0);
	    }
	  for (int i=0;i<128;++i)
	    {
	      m_ipar[i] = ((I)0);
	    }
	};

	
      public:
	
	/// <summary>
	/// The array of the integer parameters.
	/// </summary>
	pI GetParamIntegers()
	{
	  return this->m_ipar;
	};
	
	/// <summary>
	/// The array of the real parameters.
	/// </summary>
	pR GetParamReals()
	{
	  return this->m_dpar;
	};
       
      };
      
    };

  };
  
};
#endif
