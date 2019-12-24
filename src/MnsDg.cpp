

#include <signal.h>
#include <pthread.h>
#include <stdio.h>
#include <math.h>


#include "Monitor.h"
#include "Blas.h"
#include "Type.h"
#include "ns_sys.h"
#include "ns_mesh.h"
#include "mkS.h"
#include "ns_config_lapack.h"
#include "ns_constantes.h"
#include "ExternPardiso.h"
#include "ensBASIS.h"
#include "SmoothedHeaviside.h"
#include "WLA/include/matrix.hpp"



struct DG_JACOBIAN;
struct DG_VAR;

struct temp_jacvar
{
  const DG_JACOBIAN&j;
  const DG_VAR&x;
};

struct vector_handle;
struct matrix_handle;

struct temp_gemm
{
  const R a;
  const char transA;
  const char transB;
  const matrix_handle&A;
  const matrix_handle&B;
};




struct temp_gemv
{
  const R a;
  const char trans;
  const matrix_handle&A;
  const vector_handle&b;
};
struct temp_scal;
struct temp_transpose
{
  const char trans;
  const matrix_handle&A;
  inline temp_scal operator * (const R&v_) const;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
};

struct temp_scal
{
  const R a;
  const char trans;
  const matrix_handle&A;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
  inline temp_gemm operator * (const temp_transpose&v_) const;
};

inline temp_gemm temp_scal::operator * (const temp_transpose&v_) const
{
  return {a,trans,v_.trans,A,v_.A};
}

inline temp_gemv temp_scal::operator * (const vector_handle&v_) const
{
  return {a,trans,A, v_};
};
inline temp_gemm temp_scal::operator * (const matrix_handle&v_) const
{
  return {a, trans,'N',A, v_};
};

inline temp_scal temp_transpose::operator * (const R&v_) const
{
  return {v_,trans,A};
};
inline temp_gemv temp_transpose::operator * (const vector_handle&v_) const
{
  return {1.0,trans,A, v_};
};
  
inline temp_gemm temp_transpose::operator * (const matrix_handle&v_) const
{
  return {1.0, trans,'N', A, v_};
};

struct matrix_handle
{
  pR  x{};
  I   n{};
  I   m{};
  I   ld{};
  
  matrix_handle();
  matrix_handle(I n_, I m_,pR  x_,I ld_);


  inline matrix_handle& operator += (const matrix_handle& that)
  {
    I n1 = 1;    
    if (that.m != m)
      {
	std::cerr << "error dimension"  << std::endl;
      }
    R r1=1.0;
    I nn = (that.n < n) ? that.n : n;
    for (I j=0;j<m;++j)
      {
	daxpy(&nn,&r1, that.x + that.ld * j, &n1,x + ld * j,&n1);
      }
    return *this;
  }
  inline matrix_handle& operator -= (const matrix_handle& that)
  {
    I n1 = 1;    
    if (that.m != m)
      {
	std::cerr << "error dimension"  << std::endl;
      }
    R mr1=-1.0;
    I nn = (that.n < n) ? that.n : n;
    for (I j=0;j<m;++j)
      {
	daxpy(&nn,&mr1, that.x + that.ld * j, &n1,x + ld * j,&n1);
      }
    return *this;
  }

  inline matrix_handle& operator *= (double alpha_)
  {
    I n1 = 1;
    if (n == ld)
      {
	I N = n*m;
	dscal(&N,&alpha_, x, &n1);
      }
    else
      {
	for (I j=0;j<m;++j)
	  {
	    dscal(&n,&alpha_, x + ld * j, &n1);
	  }
      }
    return *this;
  }

  inline temp_gemm operator * (const matrix_handle&v_)const;  
  inline temp_gemv operator * (const vector_handle&v_)const;  
  inline temp_scal operator * (const R&v_)const;  
  inline temp_transpose transpose() const
  {
    return {'T',*this};
  };

  void clear()
  {
    if (ld == n)
      {
	for (I i=0;i<n*m;++i)
	  {
	    x[i] = 0.0;
	  }
      }
    else
      {
	for (I j=0;j<m;++j)
	  {
	    for (I i=0;i<n;++i)
	      {	    
		x[j*ld+i] = 0.0;
	      }
	  }
      }
  }

  inline matrix_handle& operator=(const temp_gemm&temp)
  {
    const R r0 = 0.0;
    I k = (temp.transA == 'N') ? temp.A.m : temp.A.n;
    Blas_dgemm(&temp.transA,
	       &temp.transB,
	       &n,
	       &m,
	       &k,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.B.x,
	       &temp.B.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };

};


temp_scal operator * (const R&v_,const matrix_handle &m) 
{
  return {v_,'N',m};
};

matrix_handle::matrix_handle(){};
matrix_handle::matrix_handle(I n_, I m_,pR  x_,I ld_) : x(x_), n(n_), m(m_), ld(ld_)
  {
  };

  inline temp_gemv matrix_handle::operator * (const vector_handle&v_) const
  {
    return {1.0,'N',*this, v_};
  };

inline temp_gemm matrix_handle::operator * (const matrix_handle&v_) const
  {
    return {1.0,'N','N',*this, v_};
  };

inline temp_scal matrix_handle::operator * (const R&v_) const
  {
    return {v_,'N',*this};
  };
  

struct vector_handle
{
  pR  x{};
  I   n{};
  I   ld{};
  vector_handle(){};
  vector_handle(I n_, pR  x_,I ld_)
    : x(x_),n(n_),ld(ld_)
  {};


  
  
  template <typename F>
  inline void apply(F f)
  {
    for (I j=0;j<n;++j)
      {
	pR e = x + j*ld;
	*e = f(*e);
      }    
  };
  
  inline vector_handle& operator=(const temp_gemv&temp)
  {
    const R r0 = 0.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };
  
  inline vector_handle& operator+=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };

  inline vector_handle& operator-=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    R a = -temp.a;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };


};

void vector_handle_def(struct vector_handle * h,I n_, pR  x_,I ld_)
{
  h->n  = n_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_def(struct matrix_handle * h,I n_, I m_,pR  x_,I ld_)
{
  h->n  = n_;
  h->m  = m_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_gemv_low(const struct matrix_handle * h,const char * trans,cst_pR a,cst_pR x,cst_pI xoff,cst_pR ry,pR y,cst_pI yoff)
{

  Blas_dgemv(trans,
	     &h->n,
	     &h->m,
	     a,
	     h->x,
	     &h->ld,
	     x,
	     xoff,
	     ry,
	     y,
	     yoff);

}

void matrix_handle_gemv(const struct matrix_handle * h,const char * trans,cst_pR a,const struct vector_handle * x,cst_pR ry,struct vector_handle * y)
{
  matrix_handle_gemv_low(h,trans,a,x->x,&x->ld,ry,y->x,&y->ld);
}


void matrix_handle_gemm(const struct matrix_handle * h,const char * transA,const char * transB,cst_pR a,const struct matrix_handle * x,cst_pR ry,struct matrix_handle * y)
{
  const I s = transA[0]=='N' ? h->m : h->n;
  Blas_dgemm(transA,transB,&y->n,&y->m,&s,a,h->x,&h->ld,x->x,&x->ld,ry,y->x,&y->ld);

}

void vector_handle_print(const struct vector_handle * h,FILE * f)
{
  for (I i=0;i<h->n;++i)
    {
      fprintf(f," " rfmt,h->x[h->ld*i]);	  
      fprintf(f,"\n");
    }

}

void matrix_handle_print(const struct matrix_handle * h,FILE * f)
{
  for (I i=0;i<h->n;++i)
    {
      for (I j=0;j<h->m;++j)
	{
	  fprintf(f," " rfmt,h->x[h->ld*j+i]);	  
	}
      fprintf(f,"\n");
    }

}

void matrix_handle_gesv(const struct matrix_handle * h,struct vector_handle * rhs,pI lcperm)
{
  I n1=1;
  I info_lapack;
  dgesv(&h->n,
	&n1,
	h->x,
	&h->ld,
	lcperm,
	rhs->x,
	&h->n,
	&info_lapack);
}




#define DGERR_MEMORY  2
#define DGERR_USER    3

#define DG_r_lc        			0
#define DG_r_n         			1

#define DG_ires_err 			0 
#define DG_ires_convergence 		1
#define DG_ires_iter_gauss_seidel 	2
#define DG_ires_required_iw_n 		3
#define DG_ires_required_rw_n 		4
#define DG_ires_n 			5

#define DG_rres_max 			0 
#define DG_rres_nrmL2 			1 
#define DG_rres_nrmLInf			2 
#define DG_rres_areaL1 			3 
#define DG_rres_jumpL2 			4 
#define DG_rres_johnson			5 
#define DG_rres_n 			6


  template <typename impl_t,
	    eTopology 	_faceShape,
	    unsigned int 	_degree>
  class treilli2d_t;
  
  template <eTopology _faceShape,unsigned int _degree> struct treilli2d_traits_t
  {
  };
  template <unsigned int _degree> struct treilli2d_traits_t<__eTopology_TRIANGLE,_degree>
  {
    static constexpr const unsigned int nnodes 		= ( (_degree+1)*(_degree+2))/2;
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim	= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 3;
  };

    
  template <eTopology _faceShape,unsigned int _degree> struct treilli2d_utils
  {	
  };
  
  template <unsigned int _degree> struct treilli2d_utils<__eTopology_TRIANGLE,_degree>
  {
    using traits_t = treilli2d_traits_t<__eTopology_TRIANGLE,_degree>;
    
    using faces_to_nodes_t 		= std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
    using nodes_int_t 			= std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
    using nodes_real_t 	= std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
    
    template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& 	w_,														nodes_real_t& 		p_) noexcept
    {
      
      static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t zero(0.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  unsigned int startedge0 = 3;
	  unsigned int startedge1 = 3 + numPointsOnEdge + numPointsOnEdge-1;  
	  unsigned int startedge2 = 3 + numPointsOnEdge*2 + numPointsOnEdge-1;
	  unsigned int startInterior = 3*numPointsOnEdge+3;
	  
	  p_[0][0] 	= zero;
	  p_[0][1] 	= zero;
	  p_[1][0] 	= one;
	  p_[1][1] 	= zero;
	  p_[2][0] 	= zero;
	  p_[2][1] 	= one;

	  //
	  // Third edge
	  //
	  for (unsigned int j=1;j<_degree;++j)
	    {
	      const auto wj 	= w_[j];
	      const auto wk 	= w_[_degree-j];
	      p_[startedge2][0] 	= zero;
	      p_[startedge2][1] 	= (three + two * wj - wk) / six;
	      //	  std::cout << p_[startedge2][1] << std::endl;
	      //	  std::cout << (three + two * wj - wk) / six << std::endl;
	      --startedge2;
	    }
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi 	= w_[i];
	      const auto wk 	= w_[_degree - i];
	      p_[startedge0][0] 	= (three + two * wi - wk) / six;
	      p_[startedge0][1] 	= zero;
	      ++startedge0;
	    }
	  
	  startedge0=3;
	  // 
	  // Second edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[startedge1][0] 	= p_[startedge0++][0]; 
	      p_[startedge1][1] 	= p_[++startedge2][1];
	      --startedge1;
	    }
	  //
	  // Interior
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi = w_[i];
	      for (unsigned int j=1;j<_degree-i;++j)
		{
		  const auto wj 	= w_[j];
		  const auto wk 	= w_[_degree-i-j];
		  p_[startInterior][0] 	= ( two* (one + wi) - (wj + wk) ) / six;
		  p_[startInterior][1] 	= ( two* (one + wj) - (wi + wk) ) / six;
		  startInterior++;
		}
	    }
	  
	};

	
	static inline void ComputeSubcnc(faces_to_nodes_t & subcnc_)
	{
	  unsigned int subCellIndex = 0;
#define _dec(_i,_j) (( (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	  for (unsigned int j=0;j<_degree;j++)
	    {
	    
	      for (unsigned int i=0;i<j;i++)
		{
		 
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-i);
		  subcnc_[subCellIndex][2] = _dec(i,j+1-i); 
		  ++subCellIndex;
		
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-1-i);
		  subcnc_[subCellIndex][2] = _dec(i+1,j-i);
		  ++subCellIndex;

		}
	    
	      subcnc_[subCellIndex][0] = _dec(j,0);
	      subcnc_[subCellIndex][1] = _dec(j+1,0);
	      subcnc_[subCellIndex][2] = _dec(j,1);
	      ++subCellIndex;
	    }
#undef _dec  
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  // COMPUTE GRID     
	  //	6 
	  //	3 7
	  //	1 4 8 
	  //	0 2 5 9
	  unsigned int nodeIndex = 0;
	  for (unsigned int i=0;i<_degree+1;i++)
	    {
	      for (unsigned int j=0;j<=i;j++)
		{
		  icoo_[nodeIndex][0] = j;
		  icoo_[nodeIndex][1] = i-j;
		  ++nodeIndex;
		}
	    }
	};

	
	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	  // VERTEX 0 
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;
	
	  // VERTEX 1 
	  icoo_[nodeIndex][0] = _degree;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;

	  // VERTEX 2
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = _degree;
	  ++nodeIndex;
	    
	  // EDGE 0 
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      icoo_[nodeIndex][0] = i+1;
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	  // EDGE 1 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree-(i+1); 
	      icoo_[nodeIndex][1] = i+1;
	      ++nodeIndex;
	    }
	  // EDGE 2 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree-(i+1);
	      ++nodeIndex;
	    }
	  // INTERIOR
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      for (unsigned int j=0;j<_degree-i-2;++j)
		{
		  icoo_[nodeIndex][0] = i + 1;
		  icoo_[nodeIndex][1] = j + 1;
		  ++nodeIndex;
		}
	    }
	};
	  
      };


#if 0
  template <unsigned int _degree> struct treilli2d_utils<FaceType::Quadrilateral,_degree>
      {
	using traits_t = treilli2d_traits_t<FaceType::Quadrilateral,_degree>;		
	using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

	template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& w_,													    nodes_real_t& 	p_) noexcept
	{
	  static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t mone(-1.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  p_[0][0] 	= mone;
	  p_[0][1] 	= mone;
	  p_[1][0] 	= one;
	  p_[1][1] 	= mone;
	  p_[2][0] 	= one;
	  p_[2][1] 	= one;
	  p_[3][0] 	= mone;
	  p_[3][1] 	= one;

	  unsigned int pointIndex = 4;
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[i];
	      p_[pointIndex][1] 	= mone;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= one;
	      p_[pointIndex][1] 	= w_[i];
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[_degree - i];
	      p_[pointIndex][1] 	= one;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= mone;
	      p_[pointIndex][1] 	= w_[_degree - i];
	      ++pointIndex;
	    }
	  
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  p_[pointIndex][0] = w_[i];
		  p_[pointIndex][1] = w_[j];	
		  ++pointIndex;
		} 
	    }
	};

	static inline void ComputeSubcnc(faces_to_nodes_t& subcnc_)
	{
	  unsigned int subCellIndex = 0;
	  // COMPUTE CNC
#define _dec(_i,_j)   (_degree+1) * (_i) + (_j) 
	  { 
	    for (unsigned int i=0;i<_degree;i++)
	      {
		for (unsigned int j=0;j<_degree;j++)
		  {		    
		    subcnc_[subCellIndex][0] = _dec( (i+1), (j+1) );
		    subcnc_[subCellIndex][1] = _dec( (i), (j+1) );
		    subcnc_[subCellIndex][2] = _dec( (i), (j) );
		    subcnc_[subCellIndex][3] = _dec( (i+1), (j) );
		    ++subCellIndex;
		  } 
	      }
	  }
#undef _dec
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  for (unsigned int i=0;i<_degree+1;++i)
	    {
	      for (unsigned int j=0;j<_degree+1;j++)
		{
		  icoo_[i*(_degree+1)+j][0] = i;	  
		  icoo_[i*(_degree+1)+j][1] = j;	
		}
	    }
	};


	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
	
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  /* EDGE 0 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = i;	  
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	
	  /* EDGE 1 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = _degree;	  
	      icoo_[nodeIndex][1] = i;
	      ++nodeIndex;
	    } 
	
	  /* EDGE 2 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree - i - 1;	  
	      icoo_[nodeIndex][1] = _degree;	
	      ++nodeIndex;
	    }

	  /* EDGE 3 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree - i - 1;	
	      ++nodeIndex;
	    }
	
	  /* INSIDE */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  icoo_[nodeIndex][0] = i;
		  icoo_[nodeIndex][1] = j;	
		  ++nodeIndex;
		} 
	    }
	  
	};

  };
#endif
      template <typename impl_t, eTopology _faceShape,unsigned int _degree>
      class treilli2d_t // : public CRTP< treilli2d_t<impl_t,_faceShape,_degree> > 
      {
      private: using traits_t 	= treilli2d_traits_t<_faceShape,_degree>;
	
      private: using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	
      protected: static constexpr const unsigned int s_numNodes 	= traits_t::nnodes;
      protected: static constexpr const unsigned int s_numElements 	= traits_t::nsubcells;
      protected: static constexpr const unsigned int s_dimension 	= traits_t::dim;
      protected: static constexpr const unsigned int s_degree 		= traits_t::degree;
      protected: static constexpr const unsigned int s_numNodesInFace	= traits_t::nnodesInFace;
	
      public: inline unsigned int 	degree() 		const noexcept { return s_degree; };
      public: inline unsigned int 	dimension() 		const noexcept { return s_dimension; };
      public: inline unsigned int 	nnodes() 		const noexcept { return s_numNodes; };    
      public: inline unsigned int 	nsubcells() 		const noexcept { return s_numElements; };
      public: inline eTopology	celltype()		const noexcept { return _faceShape; };
      public: inline unsigned int 	nnodesincell() 		const noexcept { return s_numNodesInFace; };    
	
      protected: faces_to_nodes_t m_subcnc;

      public: template <typename _float_type>
      inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				       const unsigned int dimensionIndex_) const noexcept
	{
	  return static_cast<const impl_t&>(*this).GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
	};
	
      public: inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
					       const unsigned int&localNodeIndex_) const noexcept
	{
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,subElementIndex_,(unsigned int)0,this->s_numElements-1);
//	  Debug::IsInRange(__TRACE__,localNodeIndex_,(unsigned int)0,this->s_numNodesInFace-1);
//#endif
	  return this->m_subcnc[subElementIndex_][localNodeIndex_];
	};
	
      protected: inline treilli2d_t(nodes_int_t& icoo) noexcept
	{
	  //
	  // Compute integer coordinates.
	  //
	  treilli2d_utils<_faceShape,_degree>::ComputeCoordinates(icoo);
	  std::array<unsigned int, (_degree+1)*(_degree+1)> perm;	  
	  nodes_int_t  ilagr;
	  {
	  //
	  // Compute integer coordinates.
	  //
	    treilli2d_utils<_faceShape,_degree>::ComputeCoordinates2(ilagr);	  
	    for (unsigned int i=0;i<s_numNodes;++i)
	      {	
		perm[ icoo[i][0] * (_degree+1) + icoo[i][1] ] = i+1;
	      } 
	  
	    treilli2d_utils<_faceShape,_degree>::ComputeSubcnc(this->m_subcnc);
	  
	    for (unsigned int subElementIndex=0;subElementIndex<s_numElements;++subElementIndex)
	      {
		for (unsigned int iv=0;iv<s_numNodesInFace;++iv)
		  {	  
		    const unsigned int l 				= m_subcnc[subElementIndex][iv];
		    const unsigned int i 				= ilagr[l][0];
		    const unsigned int j 				= ilagr[l][1];
		    m_subcnc[subElementIndex][iv] = perm[(_degree+1)*i+j]-1;
		  }
	      }
	  }
	
	};
    
	inline ~treilli2d_t() noexcept
	{
	};
    
      };
    

      template <eTopology _faceShape,unsigned int _degree>
      class Uniform : public treilli2d_t<Uniform<_faceShape,_degree> ,_faceShape,_degree>
      {
      
      private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
      private: nodes_int_t m_icoo;

      public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									       const unsigned int dimensionIndex_) const noexcept
	{
	  static constexpr const _float_type idegree = _float_type(1.0) / _float_type(_degree);
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
//	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
//#endif      
	  return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
	};
            
      public: inline Uniform() noexcept : treilli2d_t<Uniform<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
	{  	
	};
	
	inline ~Uniform() noexcept
	{
	};
    
      };


  template <eTopology _faceShape,unsigned int _degree>
  class Generator : public treilli2d_t<Generator<_faceShape,_degree> ,_faceShape,_degree>
  {
    
  private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
  protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
  protected: using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

  private: nodes_int_t m_icoo;
  private: nodes_real_t m_rcoo;
    
    //!
    //! @brief Get the coordinates.
    //!
  public: template <typename _float_type>
  inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				   const unsigned int dimensionIndex_) const noexcept
    {
      //#ifndef NDEBUG
      //	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
      //	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
      //#endif      
      return this->m_rcoo[nodeIndex_][dimensionIndex_];
    };
    
    
  public: inline Generator(const double * p) noexcept : treilli2d_t<Generator<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
    {
      //      Quadrature::Edge::Legendre<double,_degree-1> l;
      std::array<double,_degree+1> w;
      w[0]=-1.0;
      for (int i=0;i<_degree-1;++i)
	{
	  w[1+i] = p[i]; // l.GetPosition(i,0);
	}	  
      w[_degree]=1.0;
      treilli2d_utils<_faceShape,_degree>::LobattoPoints(w, 
						       this->m_rcoo);
      
    };
    
    inline ~Generator() noexcept
    {
    };
  };

  





cst_mkS 		mkS_derivative	(mkS shape_,const I idim)
{
  if (idim==0)
    {
      return mkS_dx(shape_);
    }
  else if (idim==1)
    {
      return mkS_dy(shape_);
    }
  else if (idim==2)
    {
      return mkS_dz(shape_);
    }
  else
    {
      fprintf(stderr,"mkS_derivative error idim = " ifmt "\n",idim);
      exit(1);
    }
  return NULL;
}


void dg_print_sol(mkS shape,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }

#if 1
    
    mkS_st shapeP;
    I degree = mkS_k(shape);
    Err err;
    mkS_definit	(&shapeP,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 degree,
		 __emk_discontinuous,
		 &err);

    double rst[64];
    double beval[64];
    
    R rwork[1024*2];
    I rwork_n = 1024*2;

    I np = mkS_n(&shapeP);
    I ee;
    mkS_lagrange_localspl_tria(&degree,
			       rst,
			       &n);
    
    mkS_basis(mkS_b(shape),
	      &np,	      
	      beval,
	      &n,
	      rst,
	      &n,
	      rwork,
	      &rwork_n,
	      &ee);
    
#if 0
    for (I i=0;i<n;++i)
      {
	for (I j=0;j<n;++j)
	  {
	    std::cout << " " <<beval[j*n+i] ;
	  }
	std::cout    << std::endl;

      }
#endif
#endif
    
    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*n);
    { I j;
      for (j=0;j<N;++j)
	{
	      double r1=1.0;
	      double r0=0.0;
	      I n1=1;
	      Blas_dgemv("T",
			 &n,
			 &n,
			 &r1,
			 beval,
			 &n,
			 &sol[n*j],
			 &n1,
			 &r0,
			 tmp,
			 &n1);
	  for (I i=0;i<n;++i)
	    {
	      fprintf(fil,"" rfmt "\n",tmp[i]);
	    }
	} } 
    fclose(fil); }
}

void dg_print_mod(mkS shape,I kmod,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }

#if 1
    mkS_st shapeP;
    I degree = mkS_k(shape);
    Err err;
    mkS_definit	(&shapeP,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 degree,
		 __emk_discontinuous,
		 &err);

    double rst[64];
    double beval[64];
    
    R rwork[1024*2];
    I rwork_n = 1024*2;

    I np = mkS_n(&shapeP);
    I ee;
    mkS_lagrange_localspl_tria(&degree,
			       rst,
			       &n);
    double tt[64];
    mkS_basis(mkS_b(shape),
	      &np,	      
	      beval,
	      &n,
	      rst,
	      &n,
	      rwork,
	      &rwork_n,
	      &ee);
    
#if 0
    for (I i=0;i<n;++i)
      {
	for (I j=0;j<n;++j)
	  {
	    std::cout << " " <<beval[j*n+i] ;
	  }
	std::cout    << std::endl;

      }
#endif
#endif
    
    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*n);
    { I j;
      for (j=0;j<N;++j)
	{
	      double r1=1.0;
	      double r0=0.0;
	      I n1=1;
	      for (I i=0;i<n;++i)
		{
		  tt[i] = sol[n*j+i];
		}

	      // 0 [0 0[ [1,n[
	      // 1 [0 1[ [3,n[
	      // 2 [0 3[ [6,n[
	      // 3 [0 6[ [10,n[
	      I low = kmod==0 ? 0 : kmod-1;	      
	      I start = 0;
	      I bound = (kmod==0) ? 0 : ((low+1)*(low+2))/2;
	      for (I i=start;i<bound;++i)
		{
		  tt[i] = 0.0;
		}
	      start = ((kmod+1)*(kmod+2))/2;

	      for (I i=start;i<n;++i)
		{
		  tt[i] = 0.0;
		}


	      
	      Blas_dgemv("T",
			 &n,
			 &n,
			 &r1,
			 beval,
			 &n,
			 tt,
			 &n1,
			 &r0,
			 tmp,
			 &n1);
		   
	  for (I i=0;i<n;++i)
	    {
	      fprintf(fil,"" rfmt "\n",tmp[i]);
	    }
	} } 
    fclose(fil); }
}

void dg_print_mesh(const ns_mesh*s_,const char * name_,...)
{
  //  const I numNodes = ns_mesh_get_numNodes(s_);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	} } 

    fprintf(fil,"Triangles\n" ifmt "\n",s_->nelm); 
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,((I)0));
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void dg_print_mesh(const ns_mesh*s_,cst_pI codelm,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	} } 

    fprintf(fil,"Triangles\n" ifmt "\n",s_->nelm); 
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void dg_print_mesh(I numNodes,I nelm,pI cnc,I cncoff,pR coo,pI codnodes,I codnodesoff,pI codelm,I codelmoff,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",nelm*3);
    { I i;
      for (i=0;i<nelm;++i)
	{
	  for (I k=0;k<3;++k)
	    {
	      cncelm[k] = cnc[cncoff*i+k];
	    }

	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",coo[cncelm[j]*2+0],coo[cncelm[j]*2+1],codnodes[codnodesoff * cncelm[j]]); } }
	} } 
    
    fprintf(fil,"Triangles\n" ifmt "\n",nelm); 
    { I i;
      for (i=0;i<nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i*codelmoff]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}

#include "CG_VAR.hpp"
#include "DG_VAR.hpp"
#include "DG_VIEW.hpp"


void dg_print_sol(DG_VIEW& 		dg_view_,
		  mkS shape_x,mkS shape,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  I nx = mkS_n(shape_x);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }
    //    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*dg_view_.m_nsubcells*3);
    WLA::vector_h f;
    R fx[128];
    WLA::vector_h::define(f,n,fx,1);
    WLA::vector_h f2;
    R fx2[128];
    WLA::vector_h::define(f2,nx,fx2,1);
    { I j;
      for (j=0;j<N;++j)
	{
	  for (I i=0;i<n;++i)
	    {
	      fx[i] = sol[n*j+i];
	    }
	  
	  dg_view_.f(f);
	  for (I k=0;k<dg_view_.m_nsubcells;++k)
	    {
	      dg_view_.f(k,f2);
	      for (I i=0;i<nx;++i)
		{
		  fprintf(fil,"" rfmt "\n",f2.x[i]);
		}
	    }
	} } 
    fclose(fil); }
}


void dg_print_mesh(const ns_mesh*	s_,
		   DG_VIEW& 		dg_view_,
		   const char * 	name_,...)
{
  
  //  const I numNodes = ns_mesh_get_numNodes(s_);

  
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }

    I nsubcells = dg_view_.m_nsubcells;
    I nelm = s_->nelm*nsubcells;
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*nsubcells*3);
    
    R cooelm[6];
    WLA::matrix_h xyz;
    WLA::matrix_h::define(xyz,3,2,cooelm,3);
    R cooelm2[6];
    WLA::matrix_h xyz2;
    WLA::matrix_h::define(xyz2,3,2,cooelm2,3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);	  
	  ns_mesh_cooelm(s_,
			 &i,
			 cooelm);	  
	  dg_view_.xyz(xyz);
	  
	  for (I k=0;k<nsubcells;++k)
	    {	      
	      dg_view_.xyz(k,xyz2);
	      {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",xyz2.x[xyz2.ld*0+j],xyz2.x[xyz2.ld*1+j],((I)0)); }}
	    }	  
	} } 

    
    fprintf(fil,"Triangles\n" ifmt "\n",nelm); 
    { I i;
      for (i=0;i<nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,((I)0));
	} }

    
    fprintf(fil,"End\n");						
    fclose(fil); }  
}

#define dim 		((I)2)
#define nfaceinelm 	((I)3)

int comp(const void * a ,const void * b)
{
  cst_pI a_ = (cst_pI)a;
  cst_pI b_ = (cst_pI)b;
  if (a_[0] < b_[0] )
    {
      return -1;
    }
  else if (a_[0] > b_[0] )
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

#include "DG_JACOBIAN.hpp"

DG_VAR& DG_VAR::operator = (const temp_jacvar&that)
  {
    that.j.gemv(that.x.m_values3.x,
		 that.x.m_values3.ld,
		 m_values3.x,
		 m_values3.ld);
    return *this;
  };


#include "DG_DATA.hpp"
#include "DG_HANDLE.hpp"


struct DG
{  

  
public:  
  typedef enum enum_DG
    {
    
      IA_lc=0,
      IA_lc_elm=1,
      IA_lc_face=2,
      IA_lc_elm_a=3,
      IA_lc_u=4,
      IA_lc_elm_u0=5,
      IA_lc_elm_u1=6,
      I_npoints=11,
      I_npoints_element=12,
      I_npoints_boundary=13,
      RA_lc=14,
      I_lc_len=16,
      RA_bmat=17,
      RA_bmatx=18,
      I_bmat_n=19,
      I_bmat_m=20,
      I_bmat_len=21,
      ERR=22,
      RA_EVAL_TETA_U=23,

      I_QFACE_N=28,
      I_QELM_N=29,

      I_TRIAL_FAMILY=30,
      I_TRIAL_DEGREE=31,
      I_TRIAL_NBASIS=32,

    
      I_TEST_FAMILY=33,
      I_TEST_DEGREE=34,
      I_TEST_NBASIS=35,
    
      I_TETA_FAMILY=36,
      I_TETA_DEGREE=37,
      I_TETA_NBASIS=38,
    
      I_TETA_U_FAMILY=39,
      I_TETA_U_DEGREE=40,
      I_TETA_U_NBASIS=41,
    
      I_TETA_A_FAMILY=42,
      I_TETA_A_DEGREE=43,
      I_TETA_A_NBASIS=44,

      I_n=64

    } info_t;

struct INFO
{
  I 		dg_iinfo[DG::I_n];
  R		dg_rres	[DG_rres_n];
  I		dg_ires	[DG_ires_n];
  I 		dg_rinfo_n;
  I 		dg_iinfo_n;
  pR		dg_rinfo;
  I 		dg_rwork_n;
  pR		dg_rwork;
  I 		dg_iwork_n;
  pI		dg_iwork;
  INFO()
  {
    dg_iinfo[DG::I_QELM_N]	= 10;
    dg_iinfo[DG::I_QFACE_N]  	= 10;  
    dg_rinfo_n			= (I)640000;
    dg_rinfo 			= (pR)malloc(sizeof(R)*dg_rinfo_n);
    dg_rwork_n			= (I)1280000;
    dg_rwork 			= (pR)malloc(sizeof(R)*dg_rwork_n);
    dg_iwork_n			= (I)1280000;
    dg_iwork 			= (pI)malloc(sizeof(I)*dg_iwork_n);
  };
  
  static void define(INFO*self)
  {
    //    printf("I_QFACE_N %d\n",DG::I_QFACE_N);
    self->dg_iinfo_n = DG::I_n;
    self->dg_iinfo[DG::I_QELM_N]	= 10;
    self->dg_iinfo[DG::I_QFACE_N]  	= 10;  
    self->dg_rinfo_n			= (I)640000;
    self->dg_rinfo 			= (pR)malloc(sizeof(R)*self->dg_rinfo_n);
    self->dg_rwork_n			= (I)1280000;
    self->dg_rwork 			= (pR)malloc(sizeof(R)*self->dg_rwork_n);
    self->dg_iwork_n			= (I)1280000;
    self->dg_iwork 			= (pI)malloc(sizeof(I)*self->dg_iwork_n);
  };
  
};
  
public: static DG_DATA * create_data(cst_mkS 		s_teta_a_,
				     cst_mkS 		s_teta_u_,
				     cst_mkS 		s_teta_,
				     cst_mkS 		s_test_,
				     cst_mkS 		s_trial_)
  {
    return new DG_DATA(mkS_n(s_teta_),
		       mkS_n(s_trial_),
		       mkS_n(s_test_),
		       mkS_n(s_teta_u_));
  };
  
public: static void dgadvection_init_with_quadrature(cst_mkS 		s_teta_a_,
						     cst_mkS 		s_teta_u_,
						     cst_mkS 		s_teta_,
						     cst_mkS 		s_test_,
						     cst_mkS 		s_trial_,
						     
						     cst_pI 		qelm_n_,
						     cst_pR 		qelm_p_,
						     cst_pR 		qelm_w_,
						     
						     cst_pI 		qface_n_,
						     cst_pR 		qface_p_,
						     cst_pR 		qface_w_,
						     
						     cst_pI           	iinfo_n_,
						     pI			iinfo_,
						     cst_pI		rinfo_n_,
						     pR 		rinfo_,
						     cst_pI 		rwork_n_,
						     pR			rwork_)
  {
  
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_n_[0] < DG::I_n)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_N]<1\n");
	exit(1);
      }
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }
  
    iinfo_[DG::I_TETA_A_NBASIS]  = (s_teta_a_) ? mkS_n(s_teta_a_):((I)0);
    iinfo_[DG::I_TETA_U_NBASIS]	 = mkS_n(s_teta_u_);
    iinfo_[DG::I_TETA_NBASIS]	 = mkS_n(s_teta_);
    iinfo_[DG::I_TRIAL_NBASIS]	 = mkS_n(s_trial_);
    iinfo_[DG::I_TEST_NBASIS]	 = mkS_n(s_test_);
    iinfo_[DG::I_QELM_N] 	 = qelm_n_[0];  

    if (ferr_[0]) 
      {
	return;
      }
  
    if (iinfo_[DG::I_TRIAL_NBASIS]!=iinfo_[DG::I_TEST_NBASIS])   
      { 
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_TEST_NBASIS]!=IINFO[DG::I_TRIAL_NBASIS]\n");    
      }

    mkS_st s_teta_a;
    mkS_st s_teta_u;
    mkS_st s_teta;
    mkS_st s_test;
    mkS_st s_trial;

    mkS_copy(&s_teta_a,s_teta_a_);
    mkS_copy(&s_teta_u,s_teta_u_);
    mkS_copy(&s_teta,s_teta_);
    mkS_copy(&s_test,s_test_);
    mkS_copy(&s_trial,s_trial_);

    const I s_teta_a_degree = mkS_k(s_teta_a_);
    const I s_teta_u_degree = mkS_k(s_teta_u_);

    const I trial_n	= iinfo_[DG::I_TRIAL_NBASIS];
    const I test_n	= iinfo_[DG::I_TEST_NBASIS];
    const I teta_n	= iinfo_[DG::I_TETA_NBASIS];
    const I teta_a_n 	= iinfo_[DG::I_TETA_A_NBASIS];
    const I teta_u_n 	= iinfo_[DG::I_TETA_U_NBASIS];

#if 0
    fprintf(stdout,"  Discontinuous Galerkin Definition\n");
    fprintf(stdout,"  trial_n   " ifmt "\n",trial_n);
    fprintf(stdout,"  test_n    " ifmt "\n",test_n);
    fprintf(stdout,"  teta_n    " ifmt "\n",teta_n);
    fprintf(stdout,"  teta_a_n  " ifmt "\n",teta_a_n);
    fprintf(stdout,"  teta_u_n  " ifmt "\n",teta_u_n);
#endif
    
    const I trial_nXtrial_n 	= trial_n * trial_n;

    
    const I npoints_element	= teta_a_n + dim * teta_u_n;
    const I npoints_boundary	= nfaceinelm * qface_n_[0];
    const I npoints 	 	= npoints_element + npoints_boundary;  
    const I bmat_n 		= trial_n * test_n + teta_n * test_n;
    //    const I bmat_m 		= npoints;
    
    const I bmatx_size 		= bmat_n*(nfaceinelm * qface_n_[0] * nfaceinelm);
    const I evalu_size 		= teta_u_n * (dim*teta_u_n + nfaceinelm * qface_n_[0]);  
    if (rinfo_n_[0] < 2 * npoints + bmat_n*npoints + bmatx_size + evalu_size)
      {
	fprintf(stderr,"*** DGI too small rinfo_ array (" ifmt "<" ifmt ")\n",
		rinfo_n_[0],
		2*npoints+bmat_n*npoints+bmatx_size+evalu_size);
	ferr_[0] = DGERR_MEMORY;
	return;
      }
  
    /* on recopie la quadrature pour les faces */  
    /* on recopie les interpolants */
    
    iinfo_[DG::I_npoints]     		= npoints;
    iinfo_[DG::I_npoints_element]  	= npoints_element;
    iinfo_[DG::I_npoints_boundary] 	= npoints_boundary;
    
    iinfo_[DG::I_bmat_n] 		= bmat_n;
    iinfo_[DG::I_bmat_m] 		= iinfo_[DG::I_npoints];
    iinfo_[DG::I_bmat_len] 		= iinfo_[DG::I_bmat_n]*iinfo_[DG::I_npoints];

    iinfo_[DG::I_lc_len] 		= 2*iinfo_[DG::I_npoints];
  
    /* matrice de flux */  
    /* les points d evaluation en coordonnees locales */
    iinfo_[DG::IA_lc]        	= 0;
    iinfo_[DG::IA_lc_elm]    	= 0;
    iinfo_[DG::IA_lc_elm_a]  	= 0;
    iinfo_[DG::IA_lc_elm_u0] 	= iinfo_[DG::I_TETA_A_NBASIS];
    iinfo_[DG::IA_lc_elm_u1] 	= iinfo_[DG::I_TETA_A_NBASIS]+iinfo_[DG::I_TETA_U_NBASIS];
    iinfo_[DG::IA_lc_face]   	= npoints_element;
    iinfo_[DG::IA_lc_u]      	= iinfo_[DG::I_TETA_A_NBASIS];
  
    iinfo_[DG::RA_bmat] 		= 2*npoints;
    iinfo_[DG::RA_bmatx] 		= iinfo_[DG::RA_bmat]  + iinfo_[DG::I_bmat_len];
    iinfo_[DG::RA_EVAL_TETA_U]    = iinfo_[DG::RA_bmatx] + bmatx_size;

    

    const I npoints_velocity_involved 	= iinfo_[DG::I_npoints]-iinfo_[DG::I_TETA_A_NBASIS];  
    
    pR bmat_x = &rinfo_[iinfo_[DG::RA_bmat]];  
    const I bmat_ld = bmat_n;
    pR teta_a_trial_test 	= &bmat_x[0];

    pR teta_u_nabla_trial_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_trial_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld];
      }
    
    //    pR teta_u_dxtrial_test 	= &bmat_x[teta_a_n*bmat_ld];
    //    pR teta_u_dytrial_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_ld];
    pR teta_u_nabla_teta_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_teta_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld+trial_nXtrial_n];
      }

    pR teta_a_teta_test  	= &bmat_x[trial_nXtrial_n];
    //    pR teta_u_dxteta_test 	= &bmat_x[teta_a_n*bmat_n+trial_nXtrial_n];
    //    pR teta_u_dyteta_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_n+trial_nXtrial_n];
    
    if (teta_a_n>0)
      {
	
	//
	// Form the matrix.
	//
	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_trial),
		mkS_b(&s_test),
		teta_a_trial_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);

#if 0
	I i;for (i=0;i<36;++i) { printf("hh %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }      
#endif

	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_teta),
		mkS_b(&s_test),
		teta_a_teta_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
#if 0
	for (i=0;i<36;++i) { printf("gg %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }
#endif
      
      }

  
    if (teta_u_n>0)
      {
	
	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_trial,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_trial_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }

	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_teta,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_teta_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }
      }  

    /* -------------------------------------------------------------------------------------- */
    /* LOCAL COORDINATES -------------------------------------------------------------------  */
    /* -------------------------------------------------------------------------------------- */            
    { pR lc = &rinfo_[iinfo_[DG::IA_lc]];

      if (iinfo_[DG::I_TETA_A_NBASIS]>0)
	{
	  mkS_lagrange_localspl_tria(&s_teta_a_degree,
				     &lc[iinfo_[DG::IA_lc_elm_a]],
				     &npoints);
	}
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u0]],
					 &npoints);
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u1]],
					 &npoints);
    
      mkS_bmapping		(qface_n_[0],
				 &lc[iinfo_[DG::IA_lc_face]],
				 &npoints,
				 qface_p_);
    
      /* -------------------------------------------------------------------------------------- */
      /* CALCUL DE L EVALUATION DE TETA_U POUR TOUS LES POINTS D INTEGRATION ELM + 3*NFACES */
      /* -------------------------------------------------------------------------------------- */            
      mkS_basis(mkS_b(&s_teta_u),
		&npoints_velocity_involved,
		&rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		&iinfo_[DG::I_TETA_U_NBASIS],
		&rinfo_[iinfo_[DG::IA_lc_u]],
		&iinfo_[DG::I_npoints],// DG::I_npoints			   			   
		rwork_,
		rwork_n_,
		ferr_);  }
  
#if 0
    {
      I i,j;pR lc = &rinfo_[iinfo_[DG::IA_lc]];
      for (i=0;i<npoints_velocity_involved;++i)
	{
	  printf("allo %e %e " ifmt "\n",lc[iinfo_[DG::IA_lc_face]+i],lc[iinfo_[DG::IA_lc_face]+npoints+i],iinfo_[DG::I_TETA_U_NBASIS]);
	  for (j=0;j<iinfo_[DG::I_TETA_U_NBASIS];++j)
	    {
	      printf("%e\n",rinfo_[iinfo_[DG::RA_EVAL_TETA_U] + iinfo_[DG::I_TETA_U_NBASIS]*i+j]);
	    }
	}
      exit(1);}
#endif


    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD */
    /* -------------------------------------------------------------------------------------- */        

    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_trial),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n],
	     bmat_n,
	     ferr_);
  
    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_teta),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n+trial_nXtrial_n],
	     bmat_n,
	     ferr_);

  
    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD CONTRAPOSEE AVEC LE BORD */
    /* -------------------------------------------------------------------------------------- */            

    { pR bmatflux = &rinfo_[iinfo_[DG::RA_bmatx]];

      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_trial),
		   mkS_b(&s_test),
		   bmatflux,
		   bmat_n,
		   ferr_);    
    
      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_teta),
		   mkS_b(&s_test),
		   &bmatflux[trial_nXtrial_n],
		   bmat_n,
		   ferr_); }
    return;
  };

  static void define(DG_HANDLE* 	handle_,
		     cst_mkS 		s_teta_a_,
		     cst_mkS 		s_teta_u_,
		     cst_mkS 		s_teta_,
		     cst_mkS 		s_test_,
		     cst_mkS 		s_trial_,
		     cst_pI         	iinfo_n_,
		     pI			iinfo_,		  
		     cst_pI		rinfo_n_,
		     pR 		rinfo_,
		     cst_pI 		rwork_n_,
		     pR			rwork_)
  {
    
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }

    const I qelm_n_	[1] 	= {iinfo_[DG::I_QELM_N]*iinfo_[DG::I_QELM_N]};
    const I qface_n_	[1] 	= {iinfo_[DG::I_QFACE_N]};
    
    R qelm_p_	[64*64*2];
    R qelm_w_	[64*64];
    R qface_p_	[64];
    R qface_w_	[64];
#if 0  
    fprintf(stdout,"generate quadrature " ifmt " " ifmt "\n",
	    qelm_n_[0],
	    qface_n_[0]);
#endif  
    {
      const I rwork_n = (rwork_n_[0]>3*iinfo_[DG::I_QELM_N])?rwork_n_[0]-3*iinfo_[DG::I_QELM_N]:(I)0;
      if (rwork_n<1)
	{
	  fprintf(stderr,"*** MOK_DGI_INIT ERROR TOO SMALL RWORK_N\n");
	  ferr_[0] = DGERR_MEMORY;
	  return;
	}
      mkQ_legendre_interval_(qface_n_,
			     qface_p_,
			     &negal1,
			     qface_w_,
			     &negal1,
			     rwork_,
			     rwork_n_,
			     ferr_);
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:first mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      mkQ_legendre_interval_(&iinfo_[DG::I_QELM_N],
			     rwork_,
			     &negal1,
			     &rwork_[2*iinfo_[DG::I_QELM_N]],
			     &negal1,
			     &rwork_[3*iinfo_[DG::I_QELM_N]],
			     rwork_n_,
			     ferr_);    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:second mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      I qelm_n=0;
      mkQ_collapse_		(&iinfo_[DG::I_QELM_N],
				 rwork_,
				 &negal1,
				 &rwork_[2*iinfo_[DG::I_QELM_N]],
				 &negal1,
				 &qelm_n,
				 qelm_p_,
				 qelm_n_,
				 qelm_w_,
				 &negal1,
				 ferr_);
    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:mkQ_legendre_interval_ failed\n");
	  return;
	}
    }
#if 0
    fprintf(stdout,"generate quadrature done\n");
#endif


    dgadvection_init_with_quadrature(s_teta_a_,
				     s_teta_u_,
				     s_teta_,
				     s_test_,
				     s_trial_,
		  
				     qelm_n_,
				     qelm_p_,
				     qelm_w_,
		  
				     qface_n_,
				     qface_p_,
				     qface_w_,

				     iinfo_n_,
				     iinfo_,

				     rinfo_n_,
				     rinfo_,

				     rwork_n_,
				     rwork_);


    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	na      = &iinfo_[DG::I_TETA_A_NBASIS];
    const I bmat_n 		= iinfo_[DG::I_bmat_n];  
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  


    
    pR brhs      		= &rwork_[0];
    pR brhs_uelm 		= &rwork_[1];
    I  npts_involving_u 	= bmat_m-na[0];/*nTot-teta_a_n*/    
    pR tmpbrhs        		= &rwork_[bmat_m]; 
    pR tmpbrhs_uelm   		= &tmpbrhs[na[0]]; 
    pR tmpbrhs_u 		= tmpbrhs_uelm;
    const I tmpbrhs_ufaceoff 	= bmat_m;
    
    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    //    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    matrix_handle_def(&handle_->m_BMAT,
		      bmat_n,
		      bmat_m,(pR)bmat,
		      bmat_n);
    

    matrix_handle_def(&handle_->m_mat_tmpbrhs_uelm,
		      nu[0],
		      dim,
		      tmpbrhs_uelm,
		      tmpbrhs_ufaceoff);

    matrix_handle_def(&handle_->m_UVWDOFS,
		      nu[0],
		      dim,
		      brhs_uelm,
		      nu[0]);

    matrix_handle_def(&handle_->m_brhs_uvw,
		      npts_involving_u,
		      dim,
		      tmpbrhs_u,
		      tmpbrhs_ufaceoff);
    
    vector_handle_def(&handle_->m_BRHS,
		      bmat_m,
		      brhs,
		      1);

    matrix_handle_def(&handle_->m_EVALU,
		      nu[0],
		      npts_involving_u,
		      &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		      iinfo_[DG::I_TETA_U_NBASIS]);
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {    
	matrix_handle_def(&handle_->m_mat_tmpbrhs_uface[localFaceIndex],
			  qface_n_[0],
			  dim,
			  &tmpbrhs[na[0]+dim*nu[0]+localFaceIndex*qface_n_[0]],
			  tmpbrhs_ufaceoff);
      }
    
    
  };


  
  static void solve(DG_HANDLE*  handle_,
		    DG_DATA*    data_,
		    cst_pR 	xa_,
		    cst_pR 	xu_,
		    cst_pR	rhs_,
		    cst_pI 	rhsoff_,
		    
		    cst_pI  	cnc_u_,
		    cst_pI  	cncoff_u_,
		    cst_pR 	data_u_,		       
		    cst_pR	data_v_,
		    
		    cst_pR	sol_,
		    cst_pI 	soloff_,
		    
		    pR		corr_,
		    cst_pI 	corroff_,
		    
		    cst_pR 	t_,
		    cst_pI 	nelm_,
		    cst_pR 	coo_,
		    cst_pI 	cooff_,
		    cst_pI 	cnc_,
		    cst_pI 	cncoff_,
		    cst_pI 	adj_,
		    cst_pI 	adjoff_,
		    cst_pI 	vcod_,
		    cst_pI 	noboundary_cod_,
		    
		    cst_pI 	rwork_n_,
		    pR  	rwork_,
		    cst_pI 	iwork_n_,
		    pI	iwork_,
		    
		    cst_pR 	rinfo_,
		    const I  	iinfo_[DG::I_n],
		    R 		rres_[DG_rres_n],
		    I 		ires_[DG_ires_n])
  {
    
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    cst_pR bmat 	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1,
      n2 	= (I)2;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/

    I marker_size	= 0;
  
    pI
      lcperm 	= NULL,
      perm 	= NULL,
      graph 	= NULL,
      marker 	= NULL,
      marker_begin= NULL,
      marker_flag = NULL;

    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  
    R mx;  
    I
      id,
      ielm,
      iter_gauss_seidel,
      b1,
      b2,
      pp,
      n,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R
      jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1],
      lcsol[21];
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*nu[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*nu[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/



    //
    //
    //

#if 0
    R 	matflux[21*(21+21)];/*trial_n*(trial_n+teta_n)*/
    pR 		matflux_corr 	= &matflux[0];
    pR 		matflux_residu 	= &matflux[trial_n[0]*trial_n[0]];

    matrix_handle_def(&handle_->m_FC,trial_n[0],trial_n[0],matflux_corr,trial_n[0]);
    matrix_handle_def(&handle_->m_FR,trial_n[0],teta_n[0],matflux_residu,trial_n[0]);
#endif

    // R lcrhs[21*2];
    const I renum_iwork_n = 3*nelm_[0] + 4*(nelm_[0]+1) + trial_n[0];
    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
  
    const I tmpbrhs_ufaceoff 	= bmat_m;
    pR tmpbrhs        	= &rwork_[bmat_m]; 
    pR tmpbrhs_ufaces 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface0 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface1 	= &tmpbrhs[1+2*nu[0]+qface_n[0]];
    pR tmpbrhs_uface2 	= &tmpbrhs[1+2*nu[0]+2*qface_n[0]];
    
    //  pR tmpbrhs_uelm   	= &tmpbrhs[1]; 
    //  pR tmpbrhs_u	= tmpbrhs_uelm;
    /*_____________________________________________________*/
    if (renum_iwork_n>iwork_n_[0])
      {
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_iw_n] = renum_iwork_n;
	fprintf(stderr,"too small iwork_n_ " ifmt " " ifmt "\n",renum_iwork_n,iwork_n_[0]);
	return;
      }

    lcperm  	= &iwork_[0];
    perm  	= &iwork_[trial_n[0]];
    graph 	= &iwork_[(nelm_[0]+1) + trial_n[0] ];
    marker 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]];
    marker_begin 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ (nelm_[0]+1)];
    marker_flag 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ 2*(nelm_[0]+1)];
    /*_____________________________________________________*/

    goto __state_renum;

    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  
  __state_renum:
    b1=1;
    b2=0;    
    { I j;for (j=trial_n[0];j<renum_iwork_n;++j) iwork_[j]=(I)0;}  
    /* */
    /* on donne un sens a l adjacence */
    { I jelm;
      for (jelm=0;jelm<nelm_[0];++jelm)
	{        
	  /* compute jaface and nrmelm */
	  neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];
	  cnc[0]     = cnc_[jelm*cncoff_[0]+0];cnc[1]     = cnc_[jelm*cncoff_[0]+1];cnc[2]     = cnc_[jelm*cncoff_[0]+2];
	  cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	  cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];

	  jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	  jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	  jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  

	  longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];  
	  data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];				
	  x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;						  
	  data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	  x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;  
	  data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	  x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	  jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	  /*____ INFO ELEMENT DONE */ 	
	  /*___________________________________________________________________________________________________________________*/  
	  /* 
	     eval u faces
	  */      
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
  
	  /*___________________________________________________________________________________________________________________*/
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);


	  
	  Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	
	  { I k;
	    for (k=0;k<qface_n[0];++k) if (uface0[k]<((R)0.0)) break;
	    graph[jelm*3+0] = (neids[0]) ? ( (k<qface_n[0])?-neids[0]:neids[0] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	  { I k;
	    for (k=0;k<qface_n[0];++k)if (uface1[k]<((R)0.0))break;
	    graph[jelm*3+1] = (neids[1]) ? ( (k<qface_n[0])?-neids[1]:neids[1] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	  { I k;
	    for (k=0;k<qface_n[0];++k)if (uface2[k]<((R)0.0))break;
	    graph[jelm*3+2] = (neids[2]) ? ( (k<qface_n[0])?-neids[2]:neids[2] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }	
	} }
    /* on met tous ceux qui n'ont pas d antecedent en premier*/
    { I jelm;
      for (jelm=0;jelm<nelm_[0];++jelm)
	{ 
	  if (graph[jelm*3]<-nelm_[0]) {perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	  if (graph[jelm*3+1]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	  if (graph[jelm*3+2]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	} }
    /*________________________________________________________________________________________________________________*/  
    pp	= b2;
    n	= (I)1;
  __state_renum_do1:
    marker_size=0;
    /* INJECTION */
    n=0;
    { I i;	
      for (i=b1;i<=b2;++i)
	{
	  { I j;
	    for (j=0;j<3;++j)
	      {
		const I s = graph[(perm[i]-1)*3+j];
		if ( (s>0) AND ( (NOT marker[s]) AND (NOT marker_flag[s]) ))
		  {
		    marker_size=marker_size+1; marker_begin[marker_size] = s;
		    marker[s] = marker_size;n=n+1;
		  }
	      } }
	} }
    
    /* FIN INJECTION */
    { I i;	    
      for (i=1;i<=marker_size;++i)
	{
	  marker_flag[marker_begin[i]]=1;perm[b2+i]=marker_begin[i];
	  marker[marker_begin[i]]=0;marker_begin[i]=0;
	} }
  
    pp+=n;b1 = b2+1;b2 = b2+n;
    if (b1>b2)
      {
	if (pp==nelm_[0])
	  goto __state_renum_do1_done;
	else
	  {
	    /* on cherche un element*/
	    { I jelm;
	      for (jelm=1;jelm<=nelm_[0];++jelm)
		{
		  if (NOT marker_flag[jelm])
		    {
		      marker_flag[jelm]=(I)1;
		      perm[b1] = jelm;
		      ++n;b2 = b1;
		      break;
		    }
		} }
	  }
      }
    if (n>0)
      {
	/* on continue a injecter */
	goto __state_renum_do1;
      }
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_renum_do1_done:
    if (b2!=nelm_[0])
      {      
	fprintf(stderr,"*** DGERR RENUM FAILED " ifmt "!=" ifmt "",b2,nelm_[0]);
	goto __state_error;
      }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    iter_gauss_seidel = 0;
  __state_next_iter_gauss_seidel:

    fprintf(stdout,"next iter gauss seidel \n");
    mx = (R)0.0;
    /*___________________________________________________________________________________________________________________*/  
    ielm=0;      
  __state_next_elm:
    id 		= perm[ielm+1]-1;
    neids[0]	= adj_[id*adjoff_[0]];
    neids[1]	= adj_[id*adjoff_[0]+1];
    neids[2]	= adj_[id*adjoff_[0]+2];

    cnc[0]        = cnc_[id*cncoff_[0]+0];
    cnc[1]        = cnc_[id*cncoff_[0]+1];
    cnc[2]        = cnc_[id*cncoff_[0]+2];

#if 0    
    vcod[0] 	= vcod_[cnc[0]-1];
    vcod[1] 	= vcod_[cnc[1]-1];
    vcod[2] 	= vcod_[cnc[2]-1];
    codface[0] 	= (neids[0])? noboundary_cod_[0] :  MAX(vcod_[0],vcod_[1]);
    codface[1] 	= (neids[1])? noboundary_cod_[0] :  MAX(vcod_[1],vcod_[2]);
    codface[2] 	= (neids[2])? noboundary_cod_[0] :  MAX(vcod_[2],vcod_[0]);
#endif
    
    neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
    /*---*/
    if (neids[0])
      {	if (neids[0]<0) neids[0]=-neids[0];
	if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
    if (neids[1])
      {	if (neids[1]<0) neids[1]=-neids[1];
	if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
    if (neids[2])
      {	if (neids[2]<0) neids[2]=-neids[2];
	if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
    /*---*/
    cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
    cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
    cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
    cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
    cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
    cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
    jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
    jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
    jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
    longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
    jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
    data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
    data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
    data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
    x = data_->mat_belm.x[2];
    data_->mat_belm.x[2]=data_->mat_belm.x[1];
    data_->mat_belm.x[1]=x;
    
    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
    /*____ INFO ELEMENT DONE */ 
    /*____ COPY LOCAL RHS */ 
    {I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	


    if (data_a_)
      {
	for (I k=0;k<teta_a_n[0];++k)
	  {
	    brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	  }
      }
    else
      {
	brhs_a[0] = jacelm[0] * xa_[0];
      }
      
  
    /*___________________________________________________________________________________________________________________*/  
    /* 
       eval u : elm and faces
    */  
    /*___________________________________________________________________________________________________________________*/  

    //
    // Copy the velocity.
    //
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }

    //
    // Wrap the solution.
    //
    vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
  
    //
    // Evaluate u at all the required coordinates.
    //
    handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
    
    /*____ TRANSFORM ON UELM */
    handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

    //
    // On the boundary.
    //
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
      }
    
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex].apply([xu_](const R r)
					       {
						 return (r<0.0) ? -r*xu_[0] : ((R)0.0);
					       });
      }
      
    /*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	if (neids[localFaceIndex])
	  {
	    Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);

	    
	    vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
	    vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);

	    
	    data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
	    data_->m_local_rhs += data_->m_flux_part_sol  * data_->vec_neisol;	 
	  }
      }
  
  
    /*--- COMPUTE MATRICES  */
    data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
  
    /*--- COMPUTE LOCAL RESIDUAL  */
    data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
  
    /*--- COMPUTE LOCAL CORRECTION   */
#if 0
    std::cout << "ielm " << ielm << std::endl;
    matrix_handle_print(&data_->m_local_matrix_part_corr,stdout);
    std::cout << " rhs " << ielm << std::endl;
    vector_handle_print(&data_->m_local_rhs,stdout);
#endif    
    matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
  
    /*--- COMPUTE NRMS OF LOCAL CORRECTION */
    { I j; for (j=0;j<test_n[0];++j) {lcsol[j] = corr_[corroff_[0]*id+j];}}
    { I j;for (j=0;j<test_n[0];++j) corr_[corroff_[0]*id+j] = data_->m_local_rhs.x[j];}
    { I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j]-=lcsol[j];}
    { R xx=(R)0.0; {I j;for (j=0;j<test_n[0];++j) xx+=data_->m_local_rhs.x[j]*data_->m_local_rhs.x[j];}
      if (mx<xx) mx = xx; }		    
  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    ++ielm;
    if (ielm<nelm_[0])
      {
	goto  __state_next_elm;
      }  
    else
      {
	goto  __state_next_elm_done;
      }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_next_elm_done:
    {
      fprintf(stdout,"*** DGINFO time %.5f iter " ifmt "/" ifmt " tol = %8.15e\n",t_[0],iter_gauss_seidel,nelm_[0],mx);
      if ( (mx>((R)1.0e-24)) AND (iter_gauss_seidel<nelm_[0]) )
	{
	  ++iter_gauss_seidel;
	  goto __state_next_iter_gauss_seidel;
	}
      else
	{	
	  goto __state_computation_done;
	}
    }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_computation_done:
    {    
      ires_[DG_ires_iter_gauss_seidel] 	= iter_gauss_seidel;
      rres_[DG_rres_max] 			= mx;
      if (iter_gauss_seidel<nelm_[0])
	{
	  ires_[DG_ires_convergence] = (I)1;
	}
      else
	{
	  ires_[DG_ires_convergence] = (I)0;
	}

      rres_[DG_rres_areaL1]	= 0.0; 
      rres_[DG_rres_nrmL2]	= 0.0; 
      rres_[DG_rres_nrmLInf]	= 0.0; 
      rres_[DG_rres_jumpL2]	= 0.0; 
      rres_[DG_rres_johnson]	= 0.0;   
      { I jelm;
	for (jelm=0;jelm<nelm_[0];++jelm)
	  { 
	    /**/
	    neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];

	    cnc[0]        = cnc_[jelm*cncoff_[0]+0];
	    cnc[1]        = cnc_[jelm*cncoff_[0]+1];
	    cnc[2]        = cnc_[jelm*cncoff_[0]+2];	  
	    cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	    cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	    cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	    cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	    cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	    cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
	    jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	    jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	    jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	    longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	    jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  
	    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];data_->mat_belm.x[1] = cooelm[0]-cooelm[2];data_->mat_belm.x[2] = cooelm[3]-cooelm[4];data_->mat_belm.x[3] = cooelm[1]-cooelm[0];  	  
	    x = data_->mat_belm.x[2];data_->mat_belm.x[2]=data_->mat_belm.x[1];data_->mat_belm.x[1]=x;	  
	    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	    /**/	  

	    Blas_dgemv("N",trial_n,trial_n,jacelm,bmat,trial_n,&corr_[jelm*corroff_[0]],&n1,&r0,data_->m_local_rhs.x,&n1);

	    rres_[DG_rres_nrmL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);

	    { I j;
	      for (j=0;j<trial_n[0];++j)
		{
		  rres_[DG_rres_areaL1] += data_->m_local_rhs.x[j];
		  if (rres_[DG_rres_nrmLInf]<nsFABS(corr_[corroff_[0]*jelm+j]))
		    {
		      rres_[DG_rres_nrmLInf]=nsFABS(corr_[corroff_[0]*jelm+j]);
		    }
		} }

	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	  
	    /* u au points de gauss de l element */
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);
	  
	    /* u.n sur les faces */
	    Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	  
	    /* mise en a zero */
	    { I j;for (j=0;j<qface_n[0];++j) uface0[j] = (uface0[j]<0.0) ? uface0[j]*(-xu_[0]) : (R)0.0;}       
	    { I j;for (j=0;j<qface_n[0];++j) uface1[j] = (uface1[j]<0.0) ? uface1[j]*(-xu_[0]) : (R)0.0;}	
	    { I j;for (j=0;j<qface_n[0];++j) uface2[j] = (uface2[j]<0.0) ? uface2[j]*(-xu_[0]) : (R)0.0;}
	  	  
	    { I k;
	      for (k=0;k<3;++k)
		{
		  if ( neids[k] )
		    { 
		      /* 
			 on forme la matrice de flux pour l'arete
		      */
		      Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[k][neids_face[k]],&bmat_n,(k==2)?uface2:((k==1)?uface1:uface0),&n1,&r0,data_->m_flux_memory,&n1);
		      
		      /* 
			 on calcule le produit scalaire
		      */
		      Blas_dgemv("N",trial_n,trial_n,&jacface[k],data_->m_flux_part_corr.x,trial_n,&corr_[(neids[k]-1)*corroff_[0]],&n1,&r1,data_->m_local_rhs.x,&n1);	   
		      rres_[DG_rres_jumpL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);/* saut */	
		    } 
		} }
	  
	  } }

      printf("compute nrms\n"
	     "   areaL1  : %e\n"
	     "   nrmL2   : %e\n"
	     "   nrmLInf : %e\n"
	     "   jumpL2  : %e\n"
	     "   jhonson : %e\n",
	     rres_[DG_rres_areaL1],
	     rres_[DG_rres_nrmL2],
	     rres_[DG_rres_nrmLInf],
	     rres_[DG_rres_jumpL2],
	     rres_[DG_rres_johnson]);
    
      return;
    }  
    /*_________________________________________________________________________________________________________________*/    
  __state_error:
    {
      printf("out\n");
      ires_[DG_ires_err] = (I)1;
      return;
    }
  }

// x_{n+1} = (D-L)^{-1} * U * x_n + (D-L)^{-1} * b_;
// x_{n+1} - x_n = ( (D-L)^{-1} * U - I) x_n + (D-L)^{-1} * b_;







  
  static DG_JACOBIAN* create_jacobian(DG_HANDLE* 	handle_,				      
				      DG::INFO * 	dgi,
				      DG_DATA*   	data_,
				      cst_pR 		xa_,
				      cst_pR 		xu_,
				      const CG_VAR&	velocity_,
				      cst_pI 		nelm_,
				      cst_pR 		coo_,
				      cst_pI 		cooff_,
				      cst_pI 		cnc_,
				      cst_pI 		cncoff_,
				      cst_pI 		adj_,
				      cst_pI 		adjoff_,
				      cst_pI 		vcod_,
				      cst_pI 		noboundary_cod_)
  {
    cst_pI 	rwork_n_ 	= &dgi->dg_rwork_n;      
    pR  	rwork_ 		= &dgi->dg_rwork[0];
    //    cst_pI 	iwork_n_	= &dgi->dg_iwork_n;
    //    pI		iwork_		= &dgi->dg_iwork[0];			      
    cst_pR 	rinfo_ 		= &dgi->dg_rinfo[0];
    pI  	iinfo_ 		= &dgi->dg_iinfo[0];
    //    pR 		rres_ 		= &dgi->dg_rres[0];
    pI 		ires_		= &dgi->dg_ires[0];


    
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    
    DG_JACOBIAN* jacobian = new DG_JACOBIAN(nelm_[0],
					    nfaceinelm,
					    trial_n[0],
					    adj_,
					    adjoff_);
    
    //    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    //    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];

    //    cst_pI	teta_n      = &iinfo_[DG::I_TETA_NBASIS];
    //    cst_pI	test_n      = &iinfo_[DG::I_TEST_NBASIS];
    //    cst_pI	teta_u_n      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    //const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    //    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    //    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    //    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      //      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/


    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  

    I
      id,
      ielm,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1];
      
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //
#if 0
#if __mk_debug__
    if (bmat_m!=1 + 2*teta_u_n[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*teta_u_n[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/

    //
    //
    //
    // R lcrhs[21*2];

    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return NULL;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
      
    /*___________________________________________________________________________________________________________________*/  
    for (ielm = 0; ielm < nelm_[0];++ielm)
      {
	id = ielm;

	//	std::cout << "ielm " << ielm << std::endl;
	neids[0]	= adj_[id*adjoff_[0]];
	neids[1]	= adj_[id*adjoff_[0]+1];
	neids[2]	= adj_[id*adjoff_[0]+2];
	cnc[0]        = cnc_[id*cncoff_[0]+0];
	cnc[1]        = cnc_[id*cncoff_[0]+1];
	cnc[2]        = cnc_[id*cncoff_[0]+2];

	
	neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
	/*---*/
	if (neids[0])
	  {	if (neids[0]<0) neids[0]=-neids[0];
	    if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
	if (neids[1])
	  {	if (neids[1]<0) neids[1]=-neids[1];
	    if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
	if (neids[2])
	  {	if (neids[2]<0) neids[2]=-neids[2];
	    if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
	/*---*/
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
	jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

	data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
	data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
	data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
	data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
	x = data_->mat_belm.x[2];
	data_->mat_belm.x[2]=data_->mat_belm.x[1];
	data_->mat_belm.x[1]=x;
    
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  

	//
	//
	//
	if (data_a_)
	  {
	    for (I k=0;k<teta_a_n[0];++k)
	      {
		brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	      }
	  }
	else
	  {
	    brhs_a[0] = jacelm[0] * xa_[0];
	  }
  
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u : elm and faces
	*/  
	/*___________________________________________________________________________________________________________________*/  

	//
	// Copy the velocity.
	//
	velocity_.dofselm(id,0,data_->uvw_ldofs.x,1);
	velocity_.dofselm(id,1,data_->uvw_ldofs.x+data_->uvw_ldofs.ld,1);
#if 0
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }

	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
#endif
	
	//
	// Evaluate u at all the required coordinates.
	//
	handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
	
	/*____ TRANSFORM ON UELM */
	handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

	//
	// On the boundary.
	//
	for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
	  }
    
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex].apply([xu_](const R r)
						   {
						     return (r<0.0) ? -r*xu_[0] : ((R)0.0);
						   });
	  }
	
	/*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    if (neids[localFaceIndex])
	      {
		//		printf("  flux edge " ifmt " \n",localFaceIndex);
		double mr1= -1.0;
		Blas_dgemv("N",&bmat_n,qface_n,&mr1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);
		// vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
		// vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);		
		jacobian->addelm(id,neids[localFaceIndex]-1,0.0,data_->m_flux_part_corr.x,data_->m_flux_part_corr.ld);
		
		//
		//
		//
		//		data_->m_flux_part_corr;
		//
		//		data_->m_flux_part_sol;

		// data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
		//	vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);
		//		data_->m_local_rhs -= data_->m_flux_part_sol  * data_->vec_neisol;	 
	      }
	  }
	
	/*--- COMPUTE MATRICES  */
	{
	  // matrix_handle msub;
	  // vector_handle rsub;
	  
	  // matrix_handle_def(&msub,handle_->m_BMAT.n,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BMAT.x,handle_->m_BMAT.ld);
	  // vector_handle_def(&rsub,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BRHS.x,1);
	  
	  data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
	  
	  // data_->m_local_matrices = msub * rsub;
	  jacobian->addelm(id,id,0.0,data_->m_local_matrix_part_corr.x,data_->m_local_matrix_part_corr.ld);

	  
	  
	  //	  vector_handle_print(&data_->m_local_matrices,stdout);
	  //	  matrix_handle hh;
	  //	  matrix_handle_def(&hh,trial_n[0],trial_n[0],data_->m_local_matrix_part_corr.x,trial_n[0]);
	  //	  printf("$$$$ " ifmt "\n",data_->m_local_matrix_part_corr.ld);
	  //	  matrix_handle_print(&hh,stdout);
	  //	  exit(1);

	  
	  //	  vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
	  //	  data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	  //	  printf("addelmcorr\n");
	}
	

	//	{I j;for (j=0;j<test_n[0];++j) printf(" sisis" ifmt " %e\n",j, rhs_[rhsoff_[0]*id+j]);}
	//
	// aE + u.nabla(E) + h(u,E) = -aF - u.nabla(F) - h(u,F)
	// JACOBIAN(trial,test) E = -R(F,test)
	// 
	// corr P3 F P2  Operator P2-P3 Operator P3-P3
	//
	// P0
	//   P1\P0
	//        P2\P0
	//
	//
	//
	//
	//
	//
	//
	//
	//

	//	data_->m_local_matrix_part_corr;
	/*--- COMPUTE LOCAL RESIDUAL  */
	// data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	/*--- COMPUTE LOCAL CORRECTION   */
	//	matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
	// matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
      }
    return jacobian;
  }


  static void compute_residual(DG_HANDLE* 	handle_,				      
			       DG::INFO * 	dgi,
			       DG_DATA*   	data_,
			       cst_pR 		xa_,
			       cst_pR 		xu_,
			       const DG_VAR& 	sol_,
			       const CG_VAR& 	velocity_,
			       DG_VAR& 		residual_,
			       cst_pI 		rwork_n_,
			       pR  		rwork_)
  {
    ns_mesh* mesh = residual_.mesh();
    I nelm_[1];
    nelm_[0] = mesh->nelm;
    cst_pI adj_ = mesh->adj;
    I adjoff_[1];
    adjoff_[0] =3;
    cst_pI cnc_ = mesh->cnc;
    I cncoff_[1];
    cncoff_[0]=6;
    cst_pR coo_ = mesh->coo;
    I cooff_[1];
    cooff_[0]=2;
    cst_pI	teta_n      	= &dgi->dg_iinfo[DG::I_TETA_NBASIS];
    //  cst_pI	test_n  = &dgi->dg_iinfo[DG::I_TEST_NBASIS];
    //  cst_pI	teta_u_n      = &dgi->dg_iinfo[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &dgi->dg_iinfo[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &dgi->dg_iinfo[DG::I_QFACE_N];
    const I bmat_n 		= dgi->dg_iinfo[DG::I_bmat_n];
    const I bmat_m 		= dgi->dg_iinfo[DG::I_bmat_m];  
    cst_pR bmatflux 		= &dgi->dg_rinfo[dgi->dg_iinfo[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[dgi->dg_iinfo[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      //      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/


    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  

    I
      id,
      ielm,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1];
      
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

    dgi->dg_ires[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	dgi->dg_ires[DG_ires_err] = (I)1;
	dgi->dg_ires[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[dgi->dg_iinfo[DG::IA_lc_face]];
    uface1 	= &rwork_[dgi->dg_iinfo[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[dgi->dg_iinfo[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
      
    /*___________________________________________________________________________________________________________________*/  
    for (ielm = 0; ielm < nelm_[0];++ielm)
      {
	id = ielm;
	residual_.dofselm(id,0,data_->m_local_rhs.x,1);
	//	{I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	
	
	
	//	std::cout << "ielm " << ielm << std::endl;
	neids[0]	= adj_[id*adjoff_[0]];
	neids[1]	= adj_[id*adjoff_[0]+1];
	neids[2]	= adj_[id*adjoff_[0]+2];
	cnc[0]        = cnc_[id*cncoff_[0]+0];
	cnc[1]        = cnc_[id*cncoff_[0]+1];
	cnc[2]        = cnc_[id*cncoff_[0]+2];

	
	neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
	/*---*/
	if (neids[0])
	  {	if (neids[0]<0) neids[0]=-neids[0];
	    if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
	if (neids[1])
	  {	if (neids[1]<0) neids[1]=-neids[1];
	    if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
	if (neids[2])
	  {	if (neids[2]<0) neids[2]=-neids[2];
	    if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
	/*---*/
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
	jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

	data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
	data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
	data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
	data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
	x = data_->mat_belm.x[2];
	data_->mat_belm.x[2]=data_->mat_belm.x[1];
	data_->mat_belm.x[1]=x;
    
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  

	//
	//
	//
	if (data_a_)
	  {
	    for (I k=0;k<teta_a_n[0];++k)
	      {
		brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	      }
	  }
	else
	  {
	    brhs_a[0] = jacelm[0] * xa_[0];
	  }
  
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u : elm and faces
	*/  
	/*___________________________________________________________________________________________________________________*/  

	//
	// Copy the velocity.
	//
	velocity_.dofselm(id,0,data_->uvw_ldofs.x,1);
	velocity_.dofselm(id,1,data_->uvw_ldofs.x+data_->uvw_ldofs.ld,1);
	//
	// Evaluate u at all the required coordinates.
	//
	handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
	
	/*____ TRANSFORM ON UELM */
	handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

	//
	// On the boundary.
	//
	for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
	  }
    
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex].apply([xu_](const R r)
						   {
						     return (r<0.0) ? -r*xu_[0] : ((R)0.0);
						   });
	  }
	
	const WLA::matrix_h& solmat = sol_.matrix();
	/*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    if (neids[localFaceIndex])
	      {
		double mr1= -1.0;
		Blas_dgemv("N",&bmat_n,qface_n,&mr1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);

		vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&solmat.x[(neids[localFaceIndex]-1)*solmat.ld], 1);
		//		vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_.m_values3.x[(neids[localFaceIndex]-1)*sol_.m_values3.ld], 1);
		data_->m_local_rhs -= data_->m_flux_part_sol  * data_->vec_neisol;	 
	      }
	    
	  }
	
	/*--- COMPUTE MATRICES  */
	{
	  data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
	  vector_handle_def(&data_->hsol,teta_n[0],(pR)&solmat.x[id*solmat.ld],1);
	  data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	}

	residual_.setdofselm(id,0,data_->m_local_rhs.x,1);
	//	{I j;for (j=0;j<test_n[0];++j) rhs_[rhsoff_[0]*id+j] = data_->m_local_rhs.x[j];}	
      }

  }
  
};



#include "DG_BOUNDARY_CONDITION.hpp"

struct DG_OPERATOR
{
  DG_JACOBIAN* 			m_jacobian;
  DG_BOUNDARY_CONDITION* 	m_bc;
  DG_DATA* 			m_dgdata;
  DG_HANDLE* 			m_dghandle;
  DG::INFO* 			m_dginfo;
  
  DG_OPERATOR(mkS 		shape_X,
	      mkS 		shape_A,
	      mkS 		shape_F,
	      mkS 		shape_E,
	      mkS 		shape_Test,
	      cst_pR 		xa_,
	      cst_pR 		xu_,
	      CG_VAR&		velocity_,
	      cst_pI 		nelm_,
	      cst_pR 		coo_,
	      cst_pI 		cooff_,
	      cst_pI 		cnc_,
	      cst_pI 		cncoff_,
	      cst_pI 		adj_,
	      cst_pI 		adjoff_,
	      cst_pI 		vcod_,
	      cst_pI 		noboundary_cod_)
  {
    mkS_st shape_U_st;
    memcpy(&shape_U_st,velocity_.shape(),sizeof(mkS_st));
    mkS shape_U = &shape_U_st;
    
    double iw[6];
    double ip[6];
    iw[0] = 1.713244923791704e-1;
    iw[1] = 3.607615730481386e-1;
    iw[2] = 4.679139345726911e-1;
    iw[3] = 4.679139345726911e-1;
    iw[4] = 3.607615730481386e-1;
    iw[5] = 1.713244923791704e-1;
  
    ip[0] = -9.32469514203152e-1;
    ip[1] = -6.612093864662645e-1;
    ip[2] = -2.386191860831969e-1;
    ip[3] = 2.386191860831969e-1;
    ip[4] = 6.612093864662645e-1;
    ip[5] = 9.32469514203152e-1;
    //    I rhsoff 	= mkS_n(shape_E);
    //    I soloff 	= mkS_n(shape_F);
    //    I corroff 	= mkS_n(shape_E);
    //    I adjoff=3;
    //    I cncoff=6;

    this->m_bc 		= (DG_BOUNDARY_CONDITION *)calloc(1,sizeof(DG_BOUNDARY_CONDITION));
    this->m_dghandle 	= (DG_HANDLE *)calloc(1,sizeof(DG_HANDLE));
    this->m_dginfo 	= (DG::INFO *)calloc(1,sizeof(DG::INFO));
    //    printf("info define\n");
    DG::INFO::define(this->m_dginfo);

    //    printf("bc define\n");
    DG_BOUNDARY_CONDITION::define(this->m_bc,
				  shape_E,
				  shape_U,
				  shape_X,
				  6,iw,ip);
    //    printf("solver define\n");
    DG::define(this->m_dghandle,
	       shape_A,
	       shape_U,
	       shape_F,
	       shape_E,
	       shape_E,
	       &this->m_dginfo->dg_iinfo_n,
	       this->m_dginfo->dg_iinfo,
	       &this->m_dginfo->dg_rinfo_n,
	       this->m_dginfo->dg_rinfo,
	       &this->m_dginfo->dg_rwork_n,
	       this->m_dginfo->dg_rwork);
    
    this->m_dgdata 	= DG::create_data	(shape_A,
						 shape_U,
						 shape_F,
						 shape_E,
						 shape_Test);
    
    this->m_jacobian 	= DG::create_jacobian	(this->m_dghandle,
						 this->m_dginfo,
						 this->m_dgdata,
						 xa_,
						 xu_,
						 velocity_,
						 nelm_,
						 coo_,
						 cooff_,
						 cnc_,
						 cncoff_,
						 adj_,
						 adjoff_,
						 vcod_,
						 noboundary_cod_);
    
  };

  void compute_residual(const DG_VAR& sol_,
			const CG_VAR& velocity_,
			DG_VAR& residual_,
			cst_pI rwork_n_,
			pR rwork_)
  {
    R xa_[1] = {0.0};
    R xu_[1] = {1.0};
    DG::compute_residual(this->m_dghandle,
			 this->m_dginfo,
			 this->m_dgdata,
			 xa_,
			 xu_,					   
			 sol_,
			 velocity_,
			 residual_,
			 rwork_n_,
			 rwork_);
  };
  
};

I jj{};

void compute_nrm(DG_JACOBIAN*jacobian,cst_pR jacelms_,I degree,pR 		x_)
{
  {
    I n1 = 1;
    double h = dnrm2(&jacobian->m_n,x_,&n1);
    printf(" " ifmt " %e",jj++,h);
    {
      for (I ideg=0;ideg<=degree;++ideg)
	{
	  double nrm=0.0;
	  for (I ielm=0;ielm<jacobian->nelm();++ielm)
	    {
	      double h1 = 0.0;
	      I start = ((ideg) * (ideg+1) )/2;
	      I bound = ((ideg+1) * (ideg+2) )/2;
	      for (I j = start;j<bound;++j)
		{
		  double x = x_[ielm*jacobian->size_block()+j];
		  h1 += x*x;
		}
	      h1 *= jacelms_[ielm];
	      nrm+=h1;
	    }
	  printf(" %e",sqrt(nrm));
	}
      printf("\n");
    }	
  }  
}

void linear_residual(DG_OPERATOR*	dg,
		     pR 		x_,
		     pR 		b_,
		     pR 		r_)
{
  R mr1=-1.0;
  R r1=1.0;
  I n1=1;
  DG_JACOBIAN*jacobian 	= dg->m_jacobian;
  jacobian->gemv(x_,r_);
  dscal(&jacobian->m_n,&mr1, r_, &n1);
  daxpy(&jacobian->m_n,&r1,b_,&n1,r_,&n1);
}


void down(DG_OPERATOR * dg[],I degree,cst_pR y_,pR x_)
{
  I k = degree-1;
  {
    for (I ielm=0;ielm<dg[degree]->m_jacobian->nelm();++ielm)
      {
	// I start = ((k) * (k+1) )/2;
	I bound = ((k+1) * (k+2) )/2;
	//		  for (I j = start;j<bound;++j)
	for (I j = 0;j<bound;++j)
	  {
	    x_[ielm*dg[k]->m_jacobian->size_block()+j] = y_[ielm*dg[degree]->m_jacobian->size_block()+j];
	  }
      }    
  }  
}

//
// What is the shape of the residual? 
//
//
void down_residual(DG_OPERATOR * dg[],I degree,cst_pR y_,pR x_)
{
  I k = degree-1;
  {
    for (I ielm=0;ielm<dg[degree]->m_jacobian->nelm();++ielm)
      {
	// I start = ((k) * (k+1) )/2;
	I bound = ((k+1) * (k+2) )/2;
	//		  for (I j = start;j<bound;++j)
	for (I j = 0;j<bound;++j)
	  {
	    x_[ielm*dg[k]->m_jacobian->size_block()+j] = y_[ielm*dg[degree]->m_jacobian->size_block()+j] ;
	  }
      }    
  }
  
}

void up(DG_OPERATOR * dg[],I yk_,cst_pR y_,I xk_,pR x_)
{
  {
    for (I ielm=0;ielm<dg[yk_]->m_jacobian->nelm();++ielm)
      {
	double h1 = 0.0;
	// I start = ((k) * (k+1) )/2;
	I bound = ((yk_+1) * (yk_+2) )/2;
	//		  for (I j = start;j<bound;++j)
	for (I j = 0;j<bound;++j)
	  {
	    x_[ielm*dg[xk_]->m_jacobian->size_block()+j] += y_[ielm*dg[yk_]->m_jacobian->size_block()+j];
	  }
#if 0
	I bound2 = ((xk_+1) * (xk_+2) )/2;
	for (I j = bound;j<bound2;++j)
	  {
	    x_[ielm*dg[xk_]->m_jacobian->size_block()+j] =0.0;
	  }
#endif
	//
      }    
  }  
}


#include "DG_LINEAR_SOLVER.hpp"
#include "DG_LINEAR_SOLVER_JACOBI.hpp"
#include "DG_LS.hpp"

void eval_u(double x,
	    double y,
	    double u[])
{
  u[0] = 10.0*y*y-12.0*x+1.0;  
  u[1] = 1.0 + y;
};




int main(int 		argc,
	 const char**	argv)
{
  const char * ifilename = argv[1];
  I degree = atol(argv[2]);
  char bbname[256];
  sprintf(   bbname,"%s",argv[4]);

  Err err;
  ns_mesh * mesh = (ns_mesh*)calloc(1,sizeof(ns_mesh));
  {
    STR errmsg;
    ns_mesh_read(mesh,
		 errmsg,
		 &err,
		 ifilename);
  }
  
  std::cout << " degree "  << degree << std::endl;

  I degree_A = 0;
  I degree_X = 1;
  I degree_U = 1;
  I degree_F = degree;
  I degree_E = degree;
  
  mkS_st shape_A;
  mkS_st shape_F;
  mkS_st shape_U;
  mkS_st shape_E;
  mkS_st shape_X;
  mkS_st shape_EE[32];
  mkS_st shape_FF[32];
  
  mkS_definit	(&shape_X,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,degree_X,__emk_discontinuous,&err);
  mkS_definit	(&shape_A,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,degree_A,__emk_discontinuous,&err);
  mkS_definit	(&shape_U,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,degree_U,__emk_discontinuous,&err);
  mkS_definit	(&shape_F,__eTopology_TRIANGLE,__emkS_FAMILY_l2orthonormal,degree_F,__emk_discontinuous,&err);
  mkS_definit	(&shape_E,__eTopology_TRIANGLE,__emkS_FAMILY_l2orthonormal,degree_E,__emk_discontinuous,&err);
  
  for (I i=0;i<=degree;++i)
    {
      mkS_definit	(&shape_EE[i],
			 __eTopology_TRIANGLE,
			 __emkS_FAMILY_l2orthonormal,
			 i,
			 __emk_discontinuous,
			 &err);
      
      mkS_definit	(&shape_FF[i],
			 __eTopology_TRIANGLE,
			 __emkS_FAMILY_l2orthonormal,
			 i,
			 __emk_discontinuous,
			 &err);
    }

  //
  // Create variables.
  //

  //
  // Fluid
  //
  DG_VAR F(mesh,
	   &shape_F);

  //
  // Velocity
  //
  CG_VAR U(mesh,2,
	   &shape_U);

  //
  // Interpolate velocity.
  //
  U.ih(eval_u);

  DG_OPERATOR*dg[32];    
  const I noboundary_cod  	= 100;  
  const I cooff = 2;
  const I cncoff = 6;
  const I adjoff = 3;

  //
  // A0 corr0 = R(F0)
  //
  // A0 E  corr0   = R(F0)
  // F  Q1 corr1\0 = R(F1\F0)
  //
  //
  //
  //  cst_pI cnc_u = mesh->cnc;
  double xa = 0.0, xu = 1.0;
  dg[degree] = new DG_OPERATOR(&shape_X,
			       &shape_A,				
			       &shape_FF[degree],
			       &shape_EE[degree],
			       &shape_EE[degree],
			       &xa,
			       &xu,
			       U,
			       &mesh->nelm,
			       mesh->coo,
			       &cooff,
			       mesh->cnc,
			       &cncoff,
			       mesh->adj,
			       &adjoff,
			       mesh->cod,
			       &noboundary_cod);
#define SUB 0
  dg[SUB] = new DG_OPERATOR(&shape_X,
			  &shape_A,				
			  &shape_FF[SUB],
			  &shape_EE[SUB],
			  &shape_EE[SUB],
			  &xa,
			  &xu,
			  U,
			  &mesh->nelm,
			  mesh->coo,
			  &cooff,
			  mesh->cnc,
			  &cncoff,
			  mesh->adj,
			  &adjoff,
			  mesh->cod,
			  &noboundary_cod);

  DG_VAR F0(mesh,&shape_FF[SUB]);
  DG_VAR E0(mesh,&shape_EE[SUB]);  
  DG_VAR Residual0	(mesh, &shape_EE[SUB]);
  DG_VAR LinearResidual0(mesh, &shape_EE[SUB]);
  
  dg[SUB]->m_bc->boundary_condition(U,
				  Residual0);
  dg[SUB]->compute_residual(F0,
			  U,
			  Residual0,
			  &dg[SUB]->m_dginfo->dg_rwork_n,
			  dg[SUB]->m_dginfo->dg_rwork);

  //
  // P0  F       x0  r0
  // E  Q1       x1  r1
  //
  // P0  F       x0  r0
  // E  Q1       x1  r1
  //
  // P0  F       x0  r0()
  // E  Q1       x1  r1
  //
  //
  // x0 = inv(G0) (r0 - F x1)
  // (G1/0 - Einv(G0) F ) x1 = inv(G0) ( r0 - F x1) - E inv(G0) r0
  //

  
  // (G1/0 ) x1 = inv(G0) ( r0 - F x1) - E inv(G0) r0
  //
  // r0 - r0 - F x1 G0 = 0
  // r1 - Einv(G0)r0  - G1/0 x1 = 0
  //  
  //  DG_LINEAR_SOLVER linear_solver0(*dg[SUB]->m_jacobian);

  {
    DG_LS linear_solver0(*dg[SUB]->m_jacobian);
    linear_solver0.presolve(*dg[SUB]->m_jacobian);
    linear_solver0.solve(*dg[SUB]->m_jacobian,
			 F0,
			 Residual0,
			 E0,
			 LinearResidual0);
  }
  
  F0 += E0;
  F += E0;
  
  {
    bbname[strlen(bbname)-3]='\0';    
    DG_VIEW view(&shape_X,&shape_F);
    printf("# Export %s\n", bbname);  
    dg_print_mesh( mesh,
		   view,
		   bbname);
    dg_print_sol  (view,
		   &shape_X,
		   &shape_F,
		   mesh->nelm,
		   F.matrix().x,
		   bbname);
  }

  DG_VAR E(mesh,&shape_E);  
  DG_VAR Residual	(mesh, &shape_EE[degree]);
  DG_VAR LinearResidual	(mesh, &shape_EE[degree]);
  //  DG_VAR Correction	(mesh, &shape_EE[degree]);
  DG_VAR Correction0	(mesh, &shape_EE[degree-1]);
  
  dg[degree]->m_bc->boundary_condition(U,
				       Residual);
  dg[degree]->compute_residual(F,
			       U,
			       Residual,
			       &dg[degree]->m_dginfo->dg_rwork_n,
			       dg[degree]->m_dginfo->dg_rwork);

  DG_LINEAR_SOLVER linear_solver(*dg[degree]->m_jacobian);
  linear_solver.presolve(*dg[degree]->m_jacobian);
  linear_solver.solve(*dg[degree]->m_jacobian,
		      F,
		      Residual,
		      E,
		      LinearResidual);
  F += E;

  
#if 0
  Residual.clear();
  dg[degree]->m_bc->boundary_condition(U,Residual);
  dg[degree]->compute_residual(F,U,
			       Residual,
			       &dg[degree]->m_dginfo->dg_rwork_n,
			       dg[degree]->m_dginfo->dg_rwork);
  linear_solver.presolve(*dg[degree]->m_jacobian);
  E.clear();
  linear_solver.solve(*dg[degree]->m_jacobian,
		      F,
		      Residual,
		      E,
		      LinearResidual);
  F+=E;
#endif

  {
    char namerr[256];
    sprintf(namerr,"%s",bbname);
    
    DG_VIEW view(&shape_X,&shape_F);
    
    printf("# Export %s\n",namerr);
    dg_print_mesh( mesh,
		   view,
		   namerr);  
    
    dg_print_sol	(view,
			 &shape_X,
			 &shape_F,
			 mesh->nelm,
			 F.matrix().x,
			 namerr);

  }
  
  //
  // F = F0
  // E = F 
  //
  {
    WLA::matrix_h& fmat	 = F.matrix();
    WLA::matrix_h& fmat0 = F0.matrix();
    for (I j=0;j<mesh->nelm;++j)
      {
	
	for (I i=0;i<fmat0.m;++i)
	  {
	    fmat.x[fmat.ld*j+i] -= fmat0.x[fmat0.ld * j + i];
	  }
	
	for (I i=fmat0.m;i<fmat.m;++i)
	  {
	    fmat.x[fmat.ld*j+i] = 0.0;
	  }
	
      }
  }
    
  {
    char namerr[256];

    sprintf(namerr,"%s.error",bbname);

    DG_VIEW view(&shape_X,&shape_F);

    printf("# Export %s\n",namerr);

    dg_print_mesh(mesh,
		  view,
		  namerr);  
    
    dg_print_sol(view,
		 &shape_X,
		 &shape_F,
		 mesh->nelm,
		 F.matrix().x,
		 namerr);    
  }
  return 0;
}




#if 0
void correction_scheme(DG_OPERATOR*	dg[],
		       I 		degree,
		       I 		mxiter,
		       DG_VAR& 		x_,
		       DG_VAR&		b_,
		       CG_VAR& 		velocity_)
{
  //  DG_VAR Residual	(velocity_->mesh, x_.shape(&shape_EE[degree]);

  //
  // 1. Pre-smoothing and get x^{p*}
  // 2. Approximate x^{p*} in x^{q*} 
  // 3. Solve with sol = x^{q*} and rhs=rhs_q(x^{p*} 
  // 4. injection the correction.
  // 5. Post smoothing
  //
  
  pR b  = b_.m_memory;
  pR x  = x_.m_memory;
  pR r  = (pR)malloc(sizeof(R)*dg[degree]->m_jacobian->m_n);
  pR x0 = (pR)malloc(sizeof(R)*dg[degree-1]->m_jacobian->m_n);
  pR r0 = (pR)malloc(sizeof(R)*dg[degree-1]->m_jacobian->m_n);

  // A0
  //
  // Approx [P0 P1]
  //
  // 
  //
  for (I iter=0;iter<10;++iter)
    {
      // 1/ Smooth
      relax(dg[degree],degree,3,x,b,r);
      // 2/ r = b_ - A_x and project.
      linear_residual(dg[degree],x,b,r);  
      //   compute_nrm(dg[degree]->m_jacobian,degree,r);
      
      down(dg,degree,x,x0);      
      //
      // Solve at degree-1
      //


#if 0
      dg[degree-1]->m_bc->boundary_condition(velocity,
					     Residual);
      dg[degree-1]->compute_residual(F,
				     velocity,
				     Residual,
				     &dg[degree-1]->m_dginfo->dg_rwork_n,
				     dg[degree-1]->m_dginfo->dg_rwork);
#endif 
      // compute_nrm(dg[degree-1]->m_jacobian,degree-1,r0);
      
      // 3/ Solve A0 e0 = r0
      relax(dg[degree-1],degree-1,200,x0,r0,r,false);
      //      printf("SUBCORR\n");
      compute_nrm(dg[degree-1]->m_jacobian,degree-1,x0);
      // 4/ x_ += e0
      up(dg,degree-1,x0,degree,x);
      
      linear_residual(dg[degree],x,b,r);  
      //      printf("Step3\n");
      compute_nrm(dg[degree]->m_jacobian,degree,r);
      
      //      printf("Step33\n");
      relax22(dg[degree],degree,3,x,b,r);
    }
  
}

#endif

#if 0
void solve_gauss_seidel(DG_JACOBIAN&	jacobian,
			pR 		jacelms_,
			I 		mxiter,
			I 		degree,
			pR 		x_,
			pR 		b_)
{  
  jacobian.inverse_diagonal();
  I size_block = jacobian.size_block();
  jacobian.inverse_lower_gemv(b_,&size_block);
  
  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);

  
  for (I i=0;i<mxiter;++i)
    {
      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);
      
      //
      // Compute the correction.
      //
      jacobian.extra_gemv_upper(tmpx,tmp2);
      jacobian.inverse_lower_gemv(tmpx,&size_block);      

      daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);
      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<jacobian.nelm();++ielm)
	      {
		double h1 = 0.0;
		I start = ((ideg) * (ideg+1) )/2;
		I bound = ((ideg+1) * (ideg+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*size_block+j];
		    h1 += x*x;
		  }
		h1 *= jacelms_[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }

      if (h < 1.0e-12)
	{
	  break;
	}
    }
  //
  //
  // corr = extra * corr;
  // corr = diag^{-1} * corr
  // corr += rhs;
  //
  //
}

// jacobian.inverse_diagonal();
///  jacobian.inverse_lower_gemv(b_,&jacobian.size_block());


void solve_jacobi(DG_JACOBIAN&	jacobian,
		  cst_pR 	jacelms_,
		  I 		mxiter,
		  I 		degree,
		  pR 		x_,
		  pR 		b_)
{  
  jacobian.inverse_diagonal();
  I size_block = jacobian.size_block();
  jacobian.inverse_diagonal_gemv(b_,&size_block);
  
  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  
  for (I i=0;i<mxiter;++i)
    {
      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);
      
      //
      // Compute the correction.
      //
      jacobian.extra_gemv(tmpx,tmp2);
      jacobian.inverse_diagonal_gemv(tmpx,&size_block);      

      daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);
      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<jacobian.nelm();++ielm)
	      {
		double h1 = 0.0;
		I start = ((ideg) * (ideg+1) )/2;
		I bound = ((ideg+1) * (ideg+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*size_block+j];
		    h1 += x*x;
		  }
		h1 *= jacelms_[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }

      if (h < 1.0e-12)
	{
	  break;
	}
    }
  //
  //
  // corr = extra * corr;
  // corr = diag^{-1} * corr
  // corr += rhs;
  //
  //
}


void solve_jacobi_truncated(DG_JACOBIAN&	jacobian,
			    cst_pR      jacelms_,
			    I 		mxiter,
			    I 		degree,
			    pR 		x_,
			    pR 		b_)
{  
  jacobian.inverse_diagonal();
  //  jacobian.inverse_lower_gemv(b_,&jacobian.size_block());
  
  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  I size_block = jacobian.size_block();
  for (I i=0;i<mxiter;++i)
    {
      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);
      
      //
      // Compute the correction.
      //
      jacobian.extra_gemv(tmpx,tmp2);
      jacobian.inverse_diagonal_gemv(tmpx,&size_block);      

      daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);
      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<jacobian.nelm();++ielm)
	      {
		double h1 = 0.0;
		I start = ((ideg) * (ideg+1) )/2;
		I bound = ((ideg+1) * (ideg+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*size_block+j];
		    h1 += x*x;
		  }
		h1 *= jacelms_[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }

      if (h < 1.0e-12)
	{
	  break;
	}
    }
  //
  //
  // corr = extra * corr;
  // corr = diag^{-1} * corr
  // corr += rhs;
  //
  //
  free(tmp2);
  free(tmpx);
}





void solve_gauss_seidel_select(DG_JACOBIAN&jacobian,cst_pR jacelms_,I q,I degree,pR x_,pR b_)
{
  jacobian.inverse_diagonal();
  I size_block = jacobian.size_block();
  jacobian.inverse_lower_gemv(b_,&size_block);

  I n1=1;
  //  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)calloc(jacobian.m_n,sizeof(R));
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  
  for (I i=0;i<10;++i)
    {
      for (I j=0;j<jacobian.m_n;++j) tmpx[i]=0.0;
      
      //      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);
      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*size_block+j] = x_[ielm*size_block + j];
	    }
	}
      //
      // Compute the correction.
      //
      jacobian.extra_gemv_upper(tmpx,tmp2);
      
      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  
	  for (I j = 0;j<size_block;++j)
	    {
	      if (j<start || j>=bound)
		{
		  tmpx[ielm*size_block+j] = 0.0;
		}
	    }
	}

      jacobian.inverse_lower_gemv(tmpx,&size_block);      
      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  
	  for (I j = 0;j<size_block;++j)
	    {
	      if (j<start || j>=bound)
		{
		  tmpx[ielm*size_block+j] = 0.0;
		}
	    }
	}

      // daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*size_block+j] -= x_[ielm*size_block + j];
	    }
	}

      

      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*size_block+j] += b_[ielm*size_block + j];
	    }
	}
      //      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);

      //
      //
      //
      // daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);      
      for (I ielm=0;ielm<jacobian.nelm();++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      x_[ielm*size_block + j] += tmpx[ielm*size_block+j];
	    }
	}

      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	// for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<jacobian.nelm();++ielm)
	      {
		double h1 = 0.0;
		I start = ((q) * (q+1) )/2;
		I bound = ((q+1) * (q+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*size_block+j];
		    h1 += x*x;
		  }
		h1 *= jacelms_[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }
#if 0      
      double mr1=-1.0;
      daxpy(&jacobian.m_n,&mr1,corr_,&n1,tmp2,&n1);
      double h = dnrm2(&jacobian.m_n,tmp2,&n1);
      double h0 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*size_block];
	  h0 += mesh->jacelm[ielm]*x*x;
	}
      double h1 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*size_block+1];
	  double y = tmp2[ielm*size_block+2];
	  h1 += mesh->jacelm[ielm]*(x*x + y*y);
	}
      double h2 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*size_block+3];
	  double y = tmp2[ielm*size_block+4];
	  double z = tmp2[ielm*size_block+5];
	  h2 += mesh->jacelm[ielm]*(x*x + y*y + z*z);
	}
      double h3 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*size_block+6];
	  double y = tmp2[ielm*size_block+7];
	  double z = tmp2[ielm*size_block+8];
	  double t = tmp2[ielm*size_block+9];
	  h3 += mesh->jacelm[ielm]*(x*x + y*y + z*z + t*t);
	}
      printf(" " ifmt " %e %e %e %e %e\n",i,h,sqrt(h0),sqrt(h1),sqrt(h2),sqrt(h3));
#endif

      if (h < 1.0e-12)
	{
	  break;
	}
      
	}
      //
      //
      // corr = extra * corr;
      // corr = diag^{-1} * corr
      // corr += rhs;
      //
      //
    }








void solve_gauss_seidel_mg(DG_JACOBIAN&jacobian,cst_pR jacelms_,I degree,pR x_,pR b_)
{
  jacobian.inverse_diagonal();
  I size_block = jacobian.size_block();
  jacobian.inverse_lower_gemv(b_,&size_block);
  
  //  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpb = (pR)malloc(sizeof(R)*jacobian.m_n);
  //  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  //  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  I n1=1;
  dcopy(&jacobian.m_n,b_,&n1,tmpb,&n1);
  for (I ielm=0;ielm<jacobian.nelm();++ielm)
    {
      for (I j = 1;j<size_block;++j)
	{
	  tmpb[ielm*size_block+j]=0.0;
	}
    }
  
  solve_gauss_seidel_select(jacobian,
			    jacelms_,
			    0,
			    degree,
			    x_,
			    tmpb);
}

#endif
