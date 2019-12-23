#pragma once

struct DG_LINEAR_SOLVER
{
  pSparse matrix;
  pExternPardiso pardiso{};

  ~DG_LINEAR_SOLVER()
  {
    if (nullptr != pardiso)
      {
	ExternPardiso_kill(pardiso);
      }
  };
  
  DG_LINEAR_SOLVER(DG_JACOBIAN&J)
  {
    pardiso = ExternPardiso_new();
  };
  
  void presolve(DG_JACOBIAN&J)
  {
    matrix = Sparse_build(J.m_n,
			  J.m_n,
			  J.m_nc,
			  J.m_begin,
			  J.m_index,
			  J.m_values);
    
    for (I i=0;i<=J.m_n;++i) J.m_begin[i]+=1;
    for (I i=0;i<J.m_nc;++i) J.m_index[i]+=1;    
    std::cout << "precompute: symbolic phase ... " << std::endl;
    std::cout << "            n  = " << J.m_n << std::endl;
    std::cout << "            nc = " << J.m_nc << std::endl;
    ExternPardiso_precompute(pardiso,
			     matrix);
    
    for (I i=0;i<=J.m_n;++i) J.m_begin[i]-=1;
    for (I i=0;i<J.m_nc;++i) J.m_index[i]-=1;
      
    std::cout << "precompute done. " << std::endl;         

  }

  void solve(DG_JACOBIAN&J,DG_VAR&F,DG_VAR&Residual,DG_VAR&E,DG_VAR&LR)
  {    
    matrix->format = 1;
    for (I i=0;i<=J.m_n;++i) J.m_begin[i]+=1;
    for (I i=0;i<J.m_nc;++i) J.m_index[i]+=1;    
    std::cout << "compute: numerical phase ..." << std::endl;
    ExternPardiso_compute(pardiso,
			  "N",
			  E.matrix().x,
			  Residual.matrix().x);
    for (I i=0;i<=J.m_n;++i) J.m_begin[i]-=1;
    for (I i=0;i<J.m_nc;++i) J.m_index[i]-=1;
  }
  
};
