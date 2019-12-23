#pragma once

struct DG_LINEAR_SOLVER_JACOBI
{
  I 		m_nblocks;
  I 		m_size_block;
  pSparse 	matrix;
  pR 		m_idiagonal;
  I m_size_blockXsize_block;
  ~DG_LINEAR_SOLVER_JACOBI()
  {
  };
  
  DG_LINEAR_SOLVER_JACOBI(DG_JACOBIAN&J)
  {
    m_idiagonal = (pR)malloc(J.m_n * J.size_block() * sizeof(R));
    matrix = Sparse_build(J.m_n,
			  J.m_n,
			  J.m_nc,
			  J.m_begin,
			  J.m_index,
			  J.m_values);
    this->m_nblocks = J.nelm();
    this->m_size_block = J.size_block();
    this->m_size_blockXsize_block = m_size_block*m_size_block;
  };

  
  void apply_preconditioner(DG_VAR& y)
  {
    R tmp[128];
    I n1=1;
    R r1=1.0;
    R r0=0.0;
    WLA::matrix_h&ymat = y.matrix();
    
    for (I ielm=0;ielm<m_nblocks;++ielm)
      {
	
	for (I i=0;i<m_size_block;++i)
	  {
	    //	    tmp[i] = y.m_values3.x[y.m_values3.ld*ielm+i];
	    tmp[i] = ymat.x[ymat.ld*ielm+i];
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
		   //		   &y.m_values3.x[y.m_values3.ld*ielm],
		   &ymat.x[ymat.ld*ielm],
		   &n1);
	//	y.m_values.x[y.m_values.ld*ielm] = tmp[0];
      }
  };
  
  
protected:  void residual(DG_JACOBIAN&	J,
		const DG_VAR&	x_,
		const DG_VAR&	b_,
		DG_VAR&		r_)
  {
    R mr1=-1.0;
    R r1=1.0;
    I n1=1;

    r_ = J * x_;    
    r_ *= -1.0;
    r_ += b_; 
  }

public:  void presolve(DG_JACOBIAN&	J)
  {
    const I n = m_size_block;    
    for (I ielm=0;ielm<m_nblocks;++ielm)
      {
	for (I i = ielm*n;i<(ielm+1)*n;++i)
	  {
	    for (I at = J.m_begin[i];at<J.m_begin[i+1];++at)
	      {
		if (J.m_index[at] == ielm * n)
		  {
		    for (I j=0;j<n;++j)
		      {
			m_idiagonal[ielm * n*n + j * n + (i-ielm*n)] = J.m_values[at + j];			
		      }
		    break;
		  }
	      }
	  }
      }


    I lcperm[1024];
    R tmp[1024];
    for (I ielm=0;ielm<m_nblocks;++ielm)
      {
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		tmp[j*n+i] = 0.0;
	      }
	  }
	
	for (I j=0;j<n;++j)
	  {
	    for (I i=0;i<n;++i)
	      {
		//	I i = j;
		tmp[j*n+i] = m_idiagonal[ielm*this->m_size_blockXsize_block+j*this->m_size_block+i];
	      }
	  }

	
	//	for (I j=0;j<n;++j)

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
	if (info_lapack!=0)
	  {
	    fprintf(stderr,"dgesv failed\n");
	    exit(1);
	  }
      }		

  };
  
public:  void solve(DG_JACOBIAN&	J,
		    DG_VAR&		F,
		    DG_VAR&		Residual,
		    DG_VAR&		E,
		    DG_VAR&		LinearResidual)
  {
    //
    //
    //
    LinearResidual.clear();
    R nrms[32];
    for (I i=0;i<800;++i)
      {
    
	//
	//
	//
	residual(J,E,Residual,LinearResidual);
	I nnn = LinearResidual.compute_nrms(nrms);
	printf(ifmt,i+1);
	for (I j=0;j<nnn;++j)
	  {
	    printf(" %e",nrms[j]);
	  }
	printf("\n");
	
	//
	//
	//
	apply_preconditioner(LinearResidual);
	
	//
	//
	//
	E += LinearResidual;

	if (nrms[0] < 1.0e-13)
	  {
	    printf("# convergence state reached.\n");
	    break;
	  }
	
      }
    
  };
  
};
