#ifndef __HEADER_Shape_HPP__
#define __HEADER_Shape_HPP__

#include "eTopology.h"

struct FiniteElementFamily
{
private:
  FiniteElementFamily(){
  };
public:

  typedef enum   { Lagrange = 0,
		   Canonic = 1,
		   L2Orthonormal = 2,
		   LagrangeBubble = 3,
		   begin = Lagrange,
		   last = LagrangeBubble} ValueType;  
  
  static const int Size = 4;
  
  static inline ValueType GetValue(const int valueIndex) 
  {
    static FiniteElementFamily::ValueType s_values[] = {FiniteElementFamily::Lagrange,
							FiniteElementFamily::Canonic,
							FiniteElementFamily::L2Orthonormal,
							FiniteElementFamily::LagrangeBubble};
    return s_values[valueIndex]; 
  };

  
#if 0
  inline ValueType operator[](const int valueIndex) 
  {
    static ValueType s_values[] = {Lagrange,
				   Canonic,
				   L2Orthonormal,
				   LagrangeBubble};
    return s_values[valueIndex]; 
  };
#endif
  
};

inline FiniteElementFamily::ValueType & operator ++(FiniteElementFamily::ValueType & self_)
{
  self_=static_cast<FiniteElementFamily::ValueType>(self_+1); return self_;
};

#if 1

std::ostream& operator<<(std::ostream&s_,const FiniteElementFamily::ValueType& d_)
{
  switch(d_)
    {
    case FiniteElementFamily::Lagrange:
      {
	s_ << "Lagrange";
	break;
      }
      
    case FiniteElementFamily::Canonic:
      {
	s_ << "Canonic";
	break;
      }
      
    case FiniteElementFamily::L2Orthonormal:
      {
	s_ << "L2Orthonormal";
	break;
      }
      
    case FiniteElementFamily::LagrangeBubble:
      {
	s_ << "LagrangeBubble";
	break;
      }
      
    }
  return s_;
}
#endif


template <eTopology _topology,FiniteElementFamily::ValueType _family,int _degree> class Shape
{

public:

  static void Eval	(const char* 	transp,
			 const char* 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info);
  
  static void DEval	(cst_pI 	componentIndex,
			 const char* 	transp,
			 const char* 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info);
};



template <> class Shape<__eTopology_TRIANGLE,FiniteElementFamily::Lagrange,1>
{

public:

  static void Eval	(cst_pS 	transp,
			 cst_pS 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info)
  {
    std::cout << "ccccc" << std::endl; 
  };
  
  static void DEval	(cst_pI 	componentIndex,
			 cst_pS 	transp,
			 cst_pS 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info)
  {

  };

};


template <> class Shape<__eTopology_TRIANGLE,FiniteElementFamily::Lagrange,0>
{

public:

  static void Eval	(cst_pS 	transp,
			 cst_pS 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info)
  {

  };
  
  static void DEval	(cst_pI 	componentIndex,
			 cst_pS 	transp,
			 cst_pS 	transr,
			 cst_pI		n,
			 cst_pR 	pos,
			 cst_pI 	posoff,
			 pR 		res,
			 cst_pI 	resoff,
			 cst_pI 	info)
  {

  };
};


#endif
