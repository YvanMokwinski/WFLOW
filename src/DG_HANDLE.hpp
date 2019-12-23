#pragma once
struct DG_HANDLE
{

  struct matrix_handle m_EVALU;
  struct matrix_handle m_UVWDOFS;
  struct matrix_handle m_BMAT;
  struct vector_handle m_BRHS;

  struct matrix_handle m_mat_tmpbrhs_uelm;
  struct matrix_handle m_brhs_uvw;
  struct matrix_handle m_mat_tmpbrhs_uface[nfaceinelm];
  
};
