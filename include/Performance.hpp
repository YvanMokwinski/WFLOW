#ifndef __header_Performance_hpp__
#define __header_Performance_hpp__

#include "ePerformance.h"



struct PerformanceOld
{
public:
  
  typedef enum Value
    {
      t0 = 0,
      t_residu,
      t_sys,
      t_solve,
      tot_residu,
      tot_sys,
      tot_solve,
      t_work,
      t_iter,
      tot_iter,
      ALL_VALUES
    } Type;
  
protected:
  
  R m_performances[__ePerformance_ALL];
  
public:

  inline PerformanceOld()
  {
    this->Reset();
  };

  inline void Reset()
  {
    { ePerformance performance = __ePerformance_ERROR;
      for (++performance;performance<__ePerformance_ALL;++performance)
	{
	  this->m_performances[performance] = ((R)0.0);
	} }     
  };

  
  inline R 	Get	(const ePerformance performance_) const
  {
    return this->m_performances[performance_];
  };
  
  inline void 	Set	(const ePerformance performance_,const R value_) 
  {
    this->m_performances[performance_] = value_;    
  };

  inline void 	Add	(const ePerformance performance_,const R value_) 
  {
    this->m_performances[performance_] += value_;    
  };

  
};



struct Performance
{
public:
  
  typedef enum Value
    {
      Initial = 0,
      Residu,
      System,
      Solve,
      Work,
      Iter,
      TotalResidu,
      TotalSystem,
      TotalSolve,
      TotalIter,
      ALL_VALUES
    } Type;
  
protected:
  
  R m_performances[ALL_VALUES];
  
public:

  inline Performance()
  {
    this->Reset();
  };

  inline void Reset()
  {
    for (int i=0;i<ALL_VALUES;++i)
      {
	this->m_performances[i] = ((R)0.0);
      }
  };

  inline R operator[](const Type field_) const
  {
    return this->m_performances[field_];
  };

  inline R& operator[](const Type field_)
  {
    return this->m_performances[field_];
  };

  
};

#endif

