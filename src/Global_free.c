#include "nsGLOBAL.h"

void Global_free(pGlobal const  self_)
{
#if 0
  if (G_->transport_usrptr)
    {
      G_->transport_usrptr = nsDG_transport_kill(G_->transport_usrptr);
    }
#endif

#if 0
  fprintf(stdout,"count_i   = "ifmt",count_r   = "ifmt"\n",mem_i,mem_r);
  fprintf(stdout,"mxcount_i = "ifmt",mxcount_r = "ifmt"\n",mem_mxni,mem_mxnr);
  if ( (mem_i>0) OR (mem_r>0) )
    {
      ns_memory_list(0);
    }
#endif
#if 0
  { size_t total = memc_i * sizeof(I);
    fprintf(stdout,"used integer memory :"ifmt" Go "ifmt" Mo "ifmt" Ko "ifmt" Oc\n",total/1024/1024/1024,(total/1024/1024)%1024,(total/1024)%1024,total%1024); }
  { size_t total = memc_r * sizeof(R);
    fprintf(stdout,"used double  memory :"ifmt" Go "ifmt" Mo "ifmt" Ko "ifmt" Oc\n",total/1024/1024/1024,(total/1024/1024)%1024,(total/1024)%1024,total%1024); }
#endif

  if (self_->blank)
    free(self_->blank);
  if (self_->iwork)
    free(self_->iwork);  
#if 0
  if (G_->params.linfo[__ens_linfo_vnsdg])
    {
      ns_var_free(&G_->var_dg);
    }
#endif
#if 0
  ns_var_free			(&G_->var_kappa);
  ns_var_free			(&G_->var_normal);
#endif

  
#if 0
  self_->linsysVectors 	= LinsysVectors_kill(self_->linsysVectors);
#endif
  
#if 0
  self_->linsysB       	= Linsys_kill(self_->linsysB);
#endif
  self_->sparseB	= Sparse_kill(self_->sparseB);
#if 0
  self_->linsysA       	= Linsys_kill(self_->linsysA);
#endif
  self_->sparseC	= Sparse_kill(self_->sparseC);
  self_->sparseS	= Sparse_kill(self_->sparseS);
  SparseStokes_free(self_->sparseStokes);


  if (self_->slip_iperm)
    {
      free(self_->slip_iperm);
    }
  extern void MKL_Free_Buffers(void);
  extern void MKL_Thread_Free_Buffers(void);
  MKL_Free_Buffers		();
  MKL_Thread_Free_Buffers	();
}
