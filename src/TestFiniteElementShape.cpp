
#include "Shape.hpp"

int main(int 		argc,
	 const char**	argv)
{

  for (FiniteElementFamily::ValueType finiteElementFamily  = FiniteElementFamily::begin;finiteElementFamily <= FiniteElementFamily::last;++finiteElementFamily)
    {
      std::cout << finiteElementFamily << std::endl;
    }
  
  for (int i=0;i<FiniteElementFamily::Size;++i)
    {
      std::cout << FiniteElementFamily::GetValue(i) << std::endl;
    }
  
  
  {
    static const char transp = 'N';
    static const char transr = 'N';
    static const I n = 2;
    R pos[] = {0.0,0.0,
	       1.0,0.0};
    static const I posoff = 2;
    R res[] = {0.0,0.0,0.0,
	       0.0,0.0,0.0};
    static const I resoff = 3;
    
    I info;

    
    Shape<__eTopology_TRIANGLE,FiniteElementFamily::Lagrange,1>::Eval(&transp,
							      &transr,
							      &n,
							      pos,
							      &posoff,
							      res,
							      &resoff,
							      &info);    


    


    
    return 0;
  }
}
