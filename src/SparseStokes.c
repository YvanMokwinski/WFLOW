#include "SparseStokes.h"
#include <stdlib.h>

#ifndef NDEBUG
#endif
cst_pSparse		cst_SparseStokes_get_A	(cst_pSparseStokes const S_) {return S_->A;}
cst_pSparseBlock	cst_SparseStokes_get_UB	(cst_pSparseStokes const S_) {return &S_->block_UB;}
cst_pSparseBlock	cst_SparseStokes_get_UU	(cst_pSparseStokes const S_) {return &S_->block_UU;}
cst_pSparseBlock	cst_SparseStokes_get_BB	(cst_pSparseStokes const S_) {return &S_->block_BB;}
cst_pSparseBlock	cst_SparseStokes_get_BU	(cst_pSparseStokes const S_) {return &S_->block_BU;}
cst_pSparseBlock	cst_SparseStokes_get_SS	(cst_pSparseStokes const S_) {return &S_->block_SS;}
cst_pSparseBlock	cst_SparseStokes_get_SU	(cst_pSparseStokes const S_) {return &S_->block_SU;}
cst_pSparseBlock	cst_SparseStokes_get_US	(cst_pSparseStokes const S_) {return &S_->block_US;}

pSparse			SparseStokes_get_A		(pSparseStokes const S_) {return S_->A;}
pSparseBlock		SparseStokes_get_UB		(pSparseStokes const S_) {return &S_->block_UB;}
pSparseBlock		SparseStokes_get_UU		(pSparseStokes const S_) {return &S_->block_UU;}
pSparseBlock		SparseStokes_get_BB		(pSparseStokes const S_) {return &S_->block_BB;}
pSparseBlock		SparseStokes_get_BU		(pSparseStokes const S_) {return &S_->block_BU;}
pSparseBlock		SparseStokes_get_SS		(pSparseStokes const S_) {return &S_->block_SS;}
pSparseBlock		SparseStokes_get_SU		(pSparseStokes const S_) {return &S_->block_SU;}
pSparseBlock		SparseStokes_get_US		(pSparseStokes const S_) {return &S_->block_US;}


void SparseStokes_dirichlet_init(pSparseStokes 	const 	S_,
				 const I 	nedge_boundary_,
				 const I 	nshape_face_u_,
				 const I 	nshape_face_p_,
				 const I 	dec_ddlu_,
				 const I 	dec_ddlv_,
				 const I 	dec_ddlp_,
				 const I 	dec_ddls_)
{


  S_->dec_ddlu 		= dec_ddlu_;
  S_->ddlu_dirichlet 	= (pI)malloc(nedge_boundary_*(1+nshape_face_u_)*sizeof(I));

  S_->dec_ddlv 		= dec_ddlv_;
  S_->ddlv_dirichlet 	= (pI)malloc(nedge_boundary_*(1+nshape_face_u_)*sizeof(I));

  S_->dec_ddlp 		= dec_ddlp_;
  S_->ddlp_dirichlet 	= (pI)malloc(nedge_boundary_*(1+nshape_face_p_)*sizeof(I));

  S_->dec_ddls 		= dec_ddls_;

  if (S_->dec_ddls>0)
    S_->ddls_dirichlet	= (pI)malloc(nedge_boundary_*(1+nshape_face_u_)*sizeof(I));

}

static void ns_build_dirichlet(pR 		rhs_,
			       const I 	n_,
			       cst_pI 		dir_,
			       const I 	dir_dec_)
{
#if __ns_debug__
  __ns_debug(rhs_,ns_build_dirichlet);
  __ns_debug(dir_,ns_build_dirichlet);
#endif
  { I i; 
    for (i=0;i<n_;++i) 
      { 
	rhs_[dir_dec_+dir_[i]] = ((R)0.0);
      } } 
}


void SparseStokes_pressure_rhs(cst_pSparseStokes const 	S_,
			       cst_pR 			x_,
			       pR 			global_rhs_)
{
  if (!S_)
    {
      fprintf(stderr,"SparseStokes_pressure_rhs S_ is nullptr.");
      exit(1);
    }
  /* residu de pression (completement lineaire) */	  
  static const R regal1 = ((R)1.0);
#if 0
  printf("SparseStokes_pressure_rhs: CALL SparseBlock_gemv %p\n",&S_->block_BU);
#endif
  SparseBlock_gemv(&regal1,
		   &S_->block_BU,
		   x_,
		   &regal1,
		   global_rhs_);
#if 0
  printf("SparseStokes_pressure_rhs: CALL SparseBlock_gemv done\n");
#endif
#if 0
  /* 
     non inclu car le terme est compris dans le residu  (pas dans le jacobien, dans le residu de u (residu_u,residu_v), le terme en p est present
     ce qui n'est pas le cas de slip, puisque slip n'est pas present dans le residu (residu_u,residu_v)
  */
  SparseBlock_gemv(&regal1,&sparse_A_block_UB,x_,&regal1,global_rhs_);
#endif
}

void SparseStokes_slip_rhs(cst_pSparseStokes const	S_,
			   cst_pR 				x_,
			   pR 				global_rhs_)
{
  /* 
     on impose les conditions de dirichlet sur le slip 
     on rappel qu on impose le slip partout et qu on vient 
     ajouter une condition de dirichlet nulle sur le slip
     pour avoir du non slip
  */
  /* 
     ce block c'est le block identite pour les conditions de diricihlet sur le slip 
     on annule les conditions slip la ou elle ne sont pas imposees
  */
  static const R regal1 = ((R)1.0);
  printf("SparseStokes_slip_rhs: CALL SparseBlock_gemv\n");

  SparseBlock_gemv(&regal1,&S_->block_SS,x_,&regal1,global_rhs_);      
  SparseBlock_gemv(&regal1,&S_->block_SU,x_,&regal1,global_rhs_);
  SparseBlock_gemv(&regal1,&S_->block_US,x_,&regal1,global_rhs_);
}

void SparseStokes_dirichlet_rhs(cst_pSparseStokes const S_,pR global_rhs_)
{
  ns_build_dirichlet(global_rhs_,S_->nddlu_dirichlet,S_->ddlu_dirichlet,S_->dec_ddlu);
  ns_build_dirichlet(global_rhs_,S_->nddlv_dirichlet,S_->ddlv_dirichlet,S_->dec_ddlv);
  ns_build_dirichlet(global_rhs_,S_->nddlp_dirichlet,S_->ddlp_dirichlet,S_->dec_ddlp);
}


void SparseStokes_free(pSparseStokes const S_)
{
  if (S_)
    {
      if (S_->ddlp_dirichlet)
	free(S_->ddlp_dirichlet); 
      if (S_->ddlu_dirichlet) 
	free(S_->ddlu_dirichlet); 
      if (S_->ddlv_dirichlet) 
	free(S_->ddlv_dirichlet); 
      if (S_->ddls_dirichlet) 
	free(S_->ddls_dirichlet);   
      Sparse_kill	(SparseStokes_get_A(S_));
      SparseBlock_free	(SparseStokes_get_BB(S_));
      SparseBlock_free	(SparseStokes_get_BU(S_));
      SparseBlock_free	(SparseStokes_get_UB(S_));
      SparseBlock_free	(SparseStokes_get_SU(S_));
      SparseBlock_free	(SparseStokes_get_US(S_));
      SparseBlock_free	(SparseStokes_get_UU(S_));
      SparseBlock_free	(SparseStokes_get_SS(S_));
    }
}

void SparseStokes_clr(pSparseStokes const S_)
{
  Sparse_clr		(SparseStokes_get_A(S_));
  SparseBlock_clr	(SparseStokes_get_BB(S_));
  SparseBlock_clr	(SparseStokes_get_BU(S_));
  SparseBlock_clr	(SparseStokes_get_UB(S_));
  SparseBlock_clr	(SparseStokes_get_SU(S_));
  SparseBlock_clr	(SparseStokes_get_US(S_));
  SparseBlock_clr	(SparseStokes_get_UU(S_));
  SparseBlock_clr	(SparseStokes_get_SS(S_));
}

void SparseStokes_clear(pSparseStokes const  S_)
{
  I j;
  for (j=0;j<3;++j)
    {
      I i;
      for (i=0;i<3;++i)
	{
	  if ( (S_->block[j*3+i])&&(! S_->linear[j*3+i]) )
	    {
	      SparseBlock_clear(S_->block[j*3+i]);	      
	    }
	}
    }
}



#if 0


#ifndef NDEBUG
cst_pSparse		cst_SparseStokes_get_A	(cst_pSparseStokes const S_) {return S_->A;}
cst_pSparseBlock 	cst_SparseStokes_get_UB	(cst_pSparseStokes const S_) {return &S_->block_UB;}
cst_pSparseBlock 	cst_SparseStokes_get_UU	(cst_pSparseStokes const S_) {return &S_->block_UU;}
cst_pSparseBlock 	cst_SparseStokes_get_BB	(cst_pSparseStokes const S_) {return &S_->block_BB;}
cst_pSparseBlock 	cst_SparseStokes_get_BU	(cst_pSparseStokes const S_) {return &S_->block_BU;}
cst_pSparseBlock 	cst_SparseStokes_get_SS	(cst_pSparseStokes const S_) {return &S_->block_SS;}
cst_pSparseBlock 	cst_SparseStokes_get_SU	(cst_pSparseStokes const S_) {return &S_->block_SU;}
cst_pSparseBlock 	cst_SparseStokes_get_US	(cst_pSparseStokes const S_) {return &S_->block_US;}

pSparse			SparseStokes_get_A		(pSparseStokes const S_) {return S_->A;}
pSparseBlock 		SparseStokes_get_UB		(pSparseStokes const S_) {return &S_->block_UB;}
pSparseBlock 		SparseStokes_get_UU		(pSparseStokes const S_) {return &S_->block_UU;}
pSparseBlock 		SparseStokes_get_BB		(pSparseStokes const S_) {return &S_->block_BB;}
pSparseBlock 		SparseStokes_get_BU		(pSparseStokes const S_) {return &S_->block_BU;}
pSparseBlock 		SparseStokes_get_SS		(pSparseStokes const S_) {return &S_->block_SS;}
pSparseBlock 		SparseStokes_get_SU		(pSparseStokes const S_) {return &S_->block_SU;}
pSparseBlock 		SparseStokes_get_US		(pSparseStokes const S_) {return &S_->block_US;}
#endif

void SparseStokes_dirichlet_init(pSparseStokes const 	S_,
				 const I 		nboundary_edge_,
				 const I 		nshape_face_u_,
				 const I 		nshape_face_p_,
				 const I 		dec_ddlu_,
				 const I 		dec_ddlv_,
				 const I 		dec_ddlp_,
				 const I 		dec_ddls_)
{


  S_->dec_ddlu 		= dec_ddlu_;
  S_->ddlu_dirichlet 	= (pI)malloc(nboundary_edge_*(1+nshape_face_u_)*sizeof(I));

  S_->dec_ddlv 		= dec_ddlv_;
  S_->ddlv_dirichlet 	= (pI)malloc(nboundary_edge_*(1+nshape_face_u_)*sizeof(I));

  S_->dec_ddlp 		= dec_ddlp_;
  S_->ddlp_dirichlet 	= (pI)malloc(nboundary_edge_*(1+nshape_face_p_)*sizeof(I));

  S_->dec_ddls 		= dec_ddls_;

  if (S_->dec_ddls>0)
    S_->ddls_dirichlet	= (pI)malloc(nboundary_edge_*(1+nshape_face_u_)*sizeof(I));

}

static void ns_build_dirichlet(pR 		rhs_,
			       const I 	n_,
			       cst_pI 		dir_,
			       const I 	dir_dec_)
{
  { I i; 
    for (i=0;i<n_;++i) 
      { 
	rhs_[dir_dec_+dir_[i]] = ((R)0.0);
      } } 
}


void SparseStokes_pressure_rhs(cst_pSparseStokes const 	S_,
				     cst_pR 				x_,
				     pR 				global_rhs_)
{
  
  
  /* residu de pression (completement lineaire) */	  
  static const R regal1 = ((R)1.0);
#if 0
  printf("SparseStokes_pressure_rhs: CALL SparseBlock_gemv %d\n",__LINE__);
#endif
  SparseBlock_gemv(&regal1,
		       &S_->block_BU,
		       x_,
		       &regal1,
		       global_rhs_);
#if 0
  printf("SparseStokes_pressure_rhs: CALL SparseBlock_gemv\n");
#endif
#if 0
  /* 
     non inclu car le terme est compris dans le residu  (pas dans le jacobien, dans le residu de u (residu_u,residu_v), le terme en p est present
     ce qui n'est pas le cas de slip, puisque slip n'est pas present dans le residu (residu_u,residu_v)
  */
  SparseBlock_gemv(&regal1,&sparse_A_block_UB,x_,&regal1,global_rhs_);
#endif
}

void SparseStokes_slip_rhs(cst_pSparseStokes const	S_,
			       cst_pR 				x_,
			       pR 				global_rhs_)
{
      /* 
	 on impose les conditions de dirichlet sur le slip 
	 on rappel qu on impose le slip partout et qu on vient 
	 ajouter une condition de dirichlet nulle sur le slip
	 pour avoir du non slip
      */
      /* 
	 ce block c'est le block identite pour les conditions de diricihlet sur le slip 
	 on annule les conditions slip la ou elle ne sont pas imposees
      */
  static const R regal1 = ((R)1.0);
  SparseBlock_gemv(&regal1,&S_->block_SS,x_,&regal1,global_rhs_);      
  SparseBlock_gemv(&regal1,&S_->block_SU,x_,&regal1,global_rhs_);
  SparseBlock_gemv(&regal1,&S_->block_US,x_,&regal1,global_rhs_);

}

void SparseStokes_dirichlet_rhs(cst_pSparseStokes const S_,pR global_rhs_)
{
  ns_build_dirichlet(global_rhs_,S_->nddlu_dirichlet,S_->ddlu_dirichlet,S_->dec_ddlu);
  ns_build_dirichlet(global_rhs_,S_->nddlv_dirichlet,S_->ddlv_dirichlet,S_->dec_ddlv);
  ns_build_dirichlet(global_rhs_,S_->nddlp_dirichlet,S_->ddlp_dirichlet,S_->dec_ddlp);

}


void SparseStokes_free(pSparseStokes const S_)
{
  if (S_)
    {
      if (S_->ddlp_dirichlet)
	free(S_->ddlp_dirichlet); 
      if (S_->ddlu_dirichlet) 
	free(S_->ddlu_dirichlet); 
      if (S_->ddlv_dirichlet) 
	free(S_->ddlv_dirichlet); 
      if (S_->ddls_dirichlet) 
	free(S_->ddls_dirichlet);   

      Sparse_free	(SparseStokes_get_A(S_));
      SparseBlock_free	(SparseStokes_get_BB(S_));
      SparseBlock_free	(SparseStokes_get_BU(S_));
      SparseBlock_free	(SparseStokes_get_UB(S_));
      SparseBlock_free	(SparseStokes_get_SU(S_));
      SparseBlock_free	(SparseStokes_get_US(S_));
      SparseBlock_free	(SparseStokes_get_UU(S_));
      SparseBlock_free	(SparseStokes_get_SS(S_));

    }
}

void SparseStokes_clr(pSparseStokes const S_)
{
  Sparse_clr		(SparseStokes_get_A(S_));
  SparseBlock_clr	(SparseStokes_get_BB(S_));
  SparseBlock_clr	(SparseStokes_get_BU(S_));
  SparseBlock_clr	(SparseStokes_get_UB(S_));
  SparseBlock_clr	(SparseStokes_get_SU(S_));
  SparseBlock_clr	(SparseStokes_get_US(S_));
  SparseBlock_clr	(SparseStokes_get_UU(S_));
  SparseBlock_clr	(SparseStokes_get_SS(S_));
}

void SparseStokes_clear(pSparseStokes const  S_)
{
  I j;
  for (j=0;j<3;++j)
    {
      I i;
      for (i=0;i<3;++i)
	{
	  if ( (S_->block[j*3+i])&&(! S_->linear[j*3+i]) )
	    {
	      SparseBlock_clear(S_->block[j*3+i]);	      
	    }
	}
    }
}
#endif
