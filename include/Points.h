#ifndef __header_Points_h__
#define __header_Points_h__
#include "Matrix.h"
#include "IMatrix.h"
#include "eDim.h"
#if 0
#ifdef __cplusplus

class Points
{
 public:
  
  ~Points()
    {
      this->m_xyz = Matrix_kill(this->m_xyz);
    };
  
  Points(const unsigned int& 	dim_,
	 const size_t&		numPoints_)
    {
#ifndef NDEBUG
      DebugVerif(numPoints_>0);
      DebugVerif( dim_>0 );
#endif
      this->m_xyz = Matrix_malloc(dim_,numPoints_);
    };
  
  
  size_t 		GetNumPoints()  const
  {
    return Matrix_get_ncol(this->m_xyz);
  };
  
  unsigned int		GetDimension() const
  {
    return Matrix_get_nrow(this->m_xyz);
  };
  
  cst_pR 		operator[](const size_t&pointIndex_) const
  {
    return &this->m_xyz->x[this->m_xyz->off*pointIndex_];
  };

  pR 			operator[](const size_t&pointIndex_) 
  {
    return &this->m_xyz->own_x[this->m_xyz->off*pointIndex_];
  };
  
  void 			Get	(const size_t& 	pointIndex_,
				 pR		const 	x_)
  {
    Matrix_copy_col(this->m_xyz,index_,x_);
  };
  
  void 			Set	(const size_t&	pointIndex_,
				 cst_pR		const 	x_)
  {
    Matrix_set_col(this->m_xyz,index_,x_);
  };
  
 protected:
  pMatrix m_xyz;
};


namespace std
{
  template <class T> std::ostream& operator<<(std::ostream&s, const Points<T> & b)
    {
      const eDim dim = b.GetDimension();
      const size_t N = b.GetNumPoints();
      s << dim << " " << N << endl;
      for (size_t i=0;i<N;++i,s << endl)
	{	  
	  for (unsigned int j=0;j<dim;++j)
	    s << " " << b[i][j];	  
	}
      return s;
    };
};


#endif
#endif


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
  pMatrix	xyz;
  pIMatrix	topologyIds;
  eDim		dim;
} Points,*RESTRICT pPoints;

typedef const Points * RESTRICT cst_pPoints;



pPoints Points_kill		(pPoints 	const self_);
pPoints Points_new		(cst_eDim		dim_,
				 const I 		nb_);
eDim 	Points_get_dim		(cst_pPoints	const 	self_);
I 	Points_get_n		(cst_pPoints	const 	self_);
void 	Points_get		(cst_pPoints	const 	self_,
				 const I 		ith_,
				 pR		const 	x_);
void 	Points_set		(pPoints	const 	self_,
				 const I 		ith_,
				 cst_pR		const 	x_);


void 	Points_fprintf_ascii	(cst_pPoints const 	self_,
				 FILE * 		out_);

#ifdef __cplusplus
}
#endif

#endif
