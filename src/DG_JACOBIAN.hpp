#pragma once
#include "Sparse/Matrix.hpp"
#include "DiagonalBlockMatrix.hpp"

struct DG_JACOBIAN
{
private:
  
  I  m_nelm;
  I  m_nfaceinelm;
  I  m_size_block;
  I  m_size_blockXsize_block;

public:
  
  WLS::Sparse::Matrix<double>* 	m_matrix;
  WLS::Sparse::SparsityPattern* m_sparsityPattern;
  WLS::DiagonalBlockMatrix* 	m_blockJacobi;
  
#if 1
  
  I  	m_n;
  I  	m_nc;
  pI 	m_begin;
  pI 	m_index;
  
  pR 	m_values;
  pR 	m_original_values;
  
  pR   	m_idiagonal {};
  pR   	m_diagonal  {};
  
  bool 	has_inverse_diagonal{};

#endif
  
public:

  inline I  nelm() const { return this->m_nelm; };
  inline I  size_block() const { return m_size_block; };
  inline I  size_blockXsize_block()const{return m_size_blockXsize_block;};
  
  void spy(const char * filename)
  {
    FILE * f = fopen(filename,"w");
    for (I i = 0;i<m_n;++i)
      {
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    fprintf(f," " ifmt " " ifmt " " rfmt"\n",m_n - i, m_index[at] + 1,m_values[at]);
	    //  fprintf(f," " ifmt " " ifmt "\n",i, m_index[at]);
	  }
	fprintf(f,"\n");
      }
    fclose(f);
  };
  
  void gemv(cst_pR 	x,
	    pR 		y)
  {
    for (I i=0;i<m_n;++i)
      {
	R s = 0.0;
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    s+=x[m_index[at]] * m_original_values[at];
	  }
	y[i] = s;
      }
  };

  void gemv(cst_pR 	x,
	    I 		xoff,
	    pR 		y,
	    I 		yoff) const
  {
    for (I i=0;i<m_n;++i)
      {
	I row = i;
	I ielm = row / m_size_block;
	I iloc = row % m_size_block;
	R s = 0.0;
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    I col = m_index[at];
	    I jelm = col / m_size_block;
	    I jloc = col % m_size_block;
	    s+=x[jelm * xoff + jloc] * m_values[at];
	  }
	y[ielm*yoff + iloc] = s;
      }
  };

  

  void inverse_diagonal()
  {
    
    if (has_inverse_diagonal)
      {
	return;
      }
    const I n = m_size_block;

    has_inverse_diagonal = true;
    if ( NULL == m_idiagonal)
      {
	{
	  m_original_values = (pR)malloc(sizeof(R)*m_nc);
	  I n1=1;
	  dcopy(&m_nc,m_values,&n1,m_original_values,&n1);
	}
	
	this->m_idiagonal = (pR)malloc(sizeof(R)*m_nelm * n*n);
	this->m_diagonal  = (pR)malloc(sizeof(R)*m_nelm * n*n);
      }
    
    //    spy("roger.txt");
    for (I ielm=0;ielm<m_nelm;++ielm)
      {
	for (I i = ielm*n;i<(ielm+1)*n;++i)
	  {
	    //	    printf("ielm " ifmt " " ifmt " " ifmt "\n",ielm,m_begin[i],m_begin[i+1]);
	    for (I at = m_begin[i];at<m_begin[i+1];++at)
	      {
		//		printf("j " ifmt " ielm * n " ifmt "\n",m_index[at],ielm*n);		
		if (m_index[at] == ielm * n)
		  {
		    for (I j=0;j<n;++j)
		      {
			m_idiagonal[ielm * n*n + j * n + (i-ielm*n)] = m_values[at + j];
			m_diagonal[ielm * n*n + j * n + (i-ielm*n)] = m_values[at + j];
			m_values[at + j] = 0.0;
		      }
		    break;
		  }
	      }
	  }

      }
    
    I lcperm[1024];
    R tmp[1024];
    for (I ielm=0;ielm<m_nelm;++ielm)
      {
#if 0
	printf("BEFORE \n");
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		printf(" " rfmt "",m_idiagonal[ielm*n*n+j*n+i]);
	      }
	    printf("\n");
	  }
	if (ielm==3)
	exit(1);
#endif
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		tmp[j*n+i] = m_idiagonal[ielm*this->m_size_blockXsize_block+j*this->m_size_block+i];
	      }
	  }
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		m_idiagonal[ielm*n*n+j*n+i] = 0.0;
	      }
	  }

	for (I j=0;j<n;++j)
	  {
	    m_idiagonal[ielm*n*n+j*n+j] = 1.0;
	  }
	
	I info_lapack;
	dgesv(&this->m_size_block,
	      &this->m_size_block,
	      tmp,
	      &this->m_size_block,
	      lcperm,
	      &m_idiagonal[ielm*this->m_size_blockXsize_block],
	      &this->m_size_block,
	      &info_lapack);
#if 0
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		printf(" " rfmt "",m_idiagonal[ielm*n*n+j*n+i]);
	      }
	    printf("\n");
	  }
	printf("info_lapack %d\n",info_lapack);
#endif
	}		
  };

  void inverse_diagonal_gemv(pR y,cst_pI yoff)
  {
    R tmp[128];
    I n1=1;
    R r1=1.0;
    R r0=0.0;
    for (I ielm=0;ielm<m_nelm;++ielm)
      {
	for (I i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	
	Blas_dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_idiagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };

  void diagonal_gemv(pR y,cst_pI yoff)
  {
    R tmp[128];
    I n1=1;
    R r1=1.0;
    R r0=0.0;
    for (I ielm=0;ielm<m_nelm;++ielm)
      {
	for (I i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	Blas_dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_diagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };

  void update(I s,I n,pR x)
  {
    
    for (I i = 0;i<n;++i)
      {
	//	R y=0.0;
	for (I at = m_begin[s+i];at<m_begin[s+i+1];++at)
	  {
	    if (m_index[at] < s + i)
	      {
    //    printf("s+i " ifmt " " ifmt "\n",s+i,m_nelm*10);
#if 1
		x[s+i] -= m_values[at] * x[m_index[at]];
#endif
	      }
	  }
      }
  };

  void inverse_lower_gemv(pR y,cst_pI yoff)
  {
    R tmp[128];
    I n1=1;
    R r1=1.0;
    R r0=0.0;
    for (I ielm=0;ielm<m_nelm;++ielm)
      {

	update(ielm*m_size_block,
	       m_size_block,
	       y);
	
	for (I i=0;i<m_size_block;++i)
	  {
	    tmp[i] = y[yoff[0]*ielm+i];
	  }
	
	// -L1 * 
	// A1
	// L1 A2
	//	
	Blas_dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   &r1,
		   &this->m_idiagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   tmp,
		   &n1,
		   &r0,
		   &y[yoff[0]*ielm],
		   &n1);
      }
  };



  void extra_gemv(pR x,pR tmp)
  {    
    for (I i = 0;i<m_n;++i)
      {
	R y=0.0;
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    y += m_values[at] * x[m_index[at]];
	  }
	tmp[i]=y;
      }
    for (I i = 0;i<m_n;++i)
      {
	x[i] = -tmp[i];
      }    
  };

  void extra_gemv_upper(pR x,pR tmp)
  {    
    for (I i = 0;i<m_n;++i)
      {
	R y=0.0;
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	    if (m_index[at] > i)
	      {
		y += m_values[at] * x[m_index[at]];
	      }
	  }
	tmp[i]=y;
      }
    for (I i = 0;i<m_n;++i)
      {
	x[i] = -tmp[i];
      }    
  };


  DG_JACOBIAN(I 	nelm_,
	      I 	nfaceinelm_,
	      I 	size_block_,
	      cst_pI 	adj_,
	      cst_pI 	adjoff_)
  {
    this->m_size_block 			= size_block_;
    this->m_size_blockXsize_block 	= size_block_*size_block_;
    this->m_nelm			= nelm_;
    this->m_nfaceinelm			= nfaceinelm_;
    
    this->m_n 				= this->m_size_block * this->m_nelm;
    this->m_nc 				= this->m_size_block * this->m_size_block * (this->m_nfaceinelm+1) * this->m_nelm;
    this->m_begin 			= (pI)calloc(this->m_n+1,sizeof(I));
    this->m_index 			= (pI)malloc(this->m_nc*sizeof(I));
    this->m_values 			= (pR)malloc(this->m_nc*sizeof(R));
    
    for (I ielm=0;ielm<nelm_;++ielm)
      {
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_size_block;++k)
	  {
	    m_begin[ielm * m_size_block + k + 1] += m_size_block;
	  }
#if 1
	//
	// Extra diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		for (I k = 0;k <  this->m_size_block;++k)
		  {
		    m_begin[ielm * m_size_block + k + 1] += m_size_block;
		  }
	      }
	  }
#endif
      }
    
    for (I i=2;i<=this->m_n;++i)
      {
	m_begin[i]+=m_begin[i-1];
      }
#if 0
    fprintf(stdout," begin[" ifmt "] = " ifmt "\n",i,m_begin[i]);
    fprintf(stdout," " ifmt " \n",m_begin[this->m_n]);
    exit(1);
#endif
    for (I ielm=0;ielm<nelm_;++ielm)
      {
#if 1
	//
	// Before diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei-1 < ielm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }	
#endif	
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_size_block;++k)
	  {
	    for (I j = 0;j <  this->m_size_block;++j)
	      {		
		m_index[m_begin[ielm * m_size_block + k]] = ielm * m_size_block + j;
		m_begin[ielm * m_size_block + k]+=1;
	      }
	  }
#if 1
	//
	// After diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei - 1 > ielm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[ielm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[ielm * m_size_block + k]+=1;
			  }
		      }
		  }
	      }
	  }
#endif	
      }

    for (I i=this->m_n;i>0;--i)
      {
	m_begin[i] = m_begin[i-1];
      }
    m_begin[0] = 0;

    
    for (I i=0;i<this->m_n;++i)
      {
	qsort(&m_index[m_begin[i]],m_begin[i+1]-m_begin[i],sizeof(I),comp);
      }
    //    this->spy("dg.txt");

    {
      static constexpr bool use_fortran_indexing 	= false;
      static constexpr bool is_owner 			= true;
      this->m_sparsityPattern = new WLS::Sparse::SparsityPattern(use_fortran_indexing,
								 this->m_n,
								 this->m_n,
								 this->m_nc,
								 this->m_begin,
								 this->m_index,
								 is_owner);
    }
    
    this->m_matrix = new WLS::Sparse::Matrix<double>(this->m_sparsityPattern);
  };

  
  temp_jacvar operator*(const DG_VAR&x_)
  {
    return {*this,x_};
  };

  void addelm(I 	ielm_,
	      I 	jelm_,
	      R         s_,
	      pR 	block_x_,
	      I  	block_off_)
  {
    for (I i = 0;i<this->m_size_block;++i)
      {
	I bound = this->m_begin[ielm_*this->m_size_block + i + 1];
	for (I at = this->m_begin[ielm_*this->m_size_block + i];at<bound;++at)
	  {
	    if (this->m_index[at] == jelm_ * this->m_size_block)
	      {
		for (I k=0;k<this->m_size_block;++k)
		  {
		    this->m_values[at+k] = block_x_[block_off_*k+i] + s_ * this->m_values[at+k];
		  }
		break;
	      }
	  }
      }
  };

};
