#pragma once
struct DG_DATA
{
  
  //
  //
  //
  struct matrix_handle m_flux{};

  //
  // submatrix.
  //
  struct matrix_handle m_flux_part_corr{};
  
  //
  // submatrix.
  //
  struct matrix_handle m_flux_part_sol{};
  
  struct matrix_handle uvw_ldofs;

  struct vector_handle m_local_matrices;
  struct matrix_handle m_local_matrix_part_corr;  
  struct matrix_handle m_local_matrix_part_sol;

  struct vector_handle vec_neicorr;
  struct vector_handle vec_neisol;

  struct matrix_handle mat_belm;
  struct vector_handle vec_nrmelm[nfaceinelm];
  struct vector_handle vec_uface[nfaceinelm];

  struct vector_handle m_local_rhs;
  struct vector_handle hsol;

  pR m_udofs;
  pR m_local_matrices_memory;
  pR m_flux_memory{};
  R belm[dim*dim];
  R nrmelm[dim * nfaceinelm];
  R lcrhs[128];

  virtual ~DG_DATA()
  {
    if (m_flux_memory)
      {
	free(m_flux_memory);
	m_flux_memory = nullptr;
      }    
    if (m_udofs)
      {
	free(m_udofs);
	m_udofs = nullptr;
      }    
  };

  DG_DATA(I teta_n_,
	  I trial_n_,
	  I test_n_,
	  I teta_u_n_)
  {
    this->m_flux_memory = (pR)malloc(sizeof(R)*(trial_n_*test_n_+teta_n_*test_n_));
    matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
    matrix_handle_def(&this->m_flux_part_sol,test_n_,teta_n_,&this->m_flux_memory[trial_n_*test_n_],test_n_);    

    this->m_udofs = (pR)malloc(sizeof(R)*(teta_u_n_*dim));
    matrix_handle_def(&this->uvw_ldofs,teta_u_n_,dim,this->m_udofs,teta_u_n_);
    
    matrix_handle_def(&this->mat_belm,
		      2,2,belm,2);

    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&this->vec_nrmelm[localFaceIndex],dim,&nrmelm[dim * localFaceIndex],1);
      }

    I len  = trial_n_*test_n_+teta_n_*test_n_;
    m_local_matrices_memory = (pR)malloc(sizeof(R)*len);
    vector_handle_def(&this->m_local_matrices,len,m_local_matrices_memory,1);
    matrix_handle_def(&this->m_local_matrix_part_corr,test_n_,trial_n_,  m_local_matrices_memory, test_n_);
    matrix_handle_def(&this->m_local_matrix_part_sol, test_n_,teta_n_,   m_local_matrices_memory + trial_n_ * test_n_, test_n_);

    vector_handle_def(&this->m_local_rhs,trial_n_,lcrhs,1);
  };

};
