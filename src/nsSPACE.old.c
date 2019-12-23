#include "nsSPACE.h"
#include "ns_var.h"
#include "ns_constantes.h"
#include "ns_sys.h"
#include "SparseStokes.h"
#include <stdarg.h>
#ifndef NDEBUG
#if 0
cst_pR 	nsSPACE_mass		(cst_nsSPACE const S_,cst_ensBASIS shape_) 	{ return S_->MASS[shape_];}
cst_pR 	nsSPACE_convx		(cst_nsSPACE const S_,cst_ensBASIS shape_) 	{ return S_->CONV_DR[shape_];}
cst_pR 	nsSPACE_convy		(cst_nsSPACE const S_,cst_ensBASIS shape_) 	{ return S_->CONV_DS[shape_];}
#endif
cst_eElement 	nsSPACE_elm		(cst_nsSPACE const S_) 				{ return ns_mesh_elm(S_->mesh); }
const ns_mesh * nsSPACE_mesh		(cst_nsSPACE const S_) 				{ return S_->mesh; }
I 		nsSPACE_nddl		(cst_nsSPACE const S_) 				{ return S_->nddl; }
I 		nsSPACE_nddlelm		(cst_nsSPACE const S_) 				{ return S_->nddlelm; }
I 		nsSPACE_nddledge	(cst_nsSPACE const S_) 				{ return S_->nddledge; }
ensBASIS 	nsSPACE_shape		(cst_nsSPACE const S_) 				{ return S_->shape; }
cst_pR 	nsSPACE_ddlcoo		(cst_nsSPACE const S_) 				{ return S_->ddlcoo; }
cst_pI 	nsSPACE_ddlcnc		(cst_nsSPACE const S_) 				{ return S_->ddlcnc; }
cst_pI 	nsSPACE_ddlcod		(cst_nsSPACE const S_) 				{ return S_->ddlcod; }
#endif


pSparse nsSPACE_get_endomorphism(cst_nsSPACE const S_)
{

  pSparse back  = Sparse_buildReference(Sparse_get_n(S_->endomorphism),
					Sparse_get_m(S_->endomorphism),
					Sparse_get_nc(S_->endomorphism),
					Sparse_get_b(S_->endomorphism),
					Sparse_get_i(S_->endomorphism),
					NULL);
  return back;
#if 0
  nsSPARSE_ST * F = malloc(sizeof(nsSPARSE_ST));
  memcpy(F,&S_->endomorphism,sizeof(nsSPARSE_ST));
  F->own_x 	= malloc(sizeof(R)*F->nc);
  F->x 		= F->own_x;
  F->own_i 	= NULL;
  F->own_b 	= NULL;
  return F;
#endif
}

pSparse ns_mesh_get_endomorphism	(ns_mesh * 			mesh_,
					 cst_ensBASIS 	shape_)
{
  return nsSPACE_get_endomorphism( ( mesh_->spaces[shape_] ) ? mesh_->spaces[shape_] : nsSPACE_new(mesh_,shape_) );
}

I 		nsSPACE_iddl2ddlP6	(cst_nsSPACE 		SF_,
					 const I 		iddl_)
{
  if (iddl_<SF_->mesh->nvertex)
    {
      return iddl_;
    }
  else if (iddl_<SF_->mesh->nvertex+SF_->mesh->nedge*SF_->nddledge)
    {
      return SF_->mesh->nvertex + (iddl_-SF_->mesh->nvertex)/SF_->nddledge;
    }
  else
    {
#if 0
      printf("iddl_ = "ifmt"\n",iddl_);
      printf("elm = "ifmt"\n",(iddl_-SF_->mesh->nvertex-SF_->mesh->nedge*SF_->nddledge)/(SF_->nddlelm-SF_->nddledge*3-3));
#endif
      return SF_->mesh->nvertex + SF_->mesh->nedge + (iddl_-SF_->mesh->nvertex-SF_->mesh->nedge*SF_->nddledge)/(SF_->nddlelm-SF_->nddledge*3-3);
    }
}


static I nsSPACE_galerkin_line_fups( cst_nsSPACE 	SF_,
					 cst_nsSPACE 	SU_,
					 cst_nsSPACE 	SP_,
					 cst_nsSPACE 	SS_,
					 const I 		iddl_,
					 pI 			ddl_,
					 pI			blank_,
					 pI 			symbolic_cnc_,
					 pI 			slip_iperm_,
					 cst_pI 		decu_,
					 cst_pI		total_nddl_)
{

  const ns_mesh * mesh 	= nsSPACE_mesh(SU_);
  const I nu 	= nsSPACE_nddlelm(SU_);
  const I nf	= (SF_)?SF_->nddlelm:0;
  const I nddlelm 	= nf + 2*nu + SP_->nddlelm;
  I N2 		= 0;
  I N		= 0;  
  const I jddl 		= nsSPACE_iddl2ddlP6(SU_,iddl_);
  { I i;
    const I n = (jddl>=mesh->nvertex+mesh->nedge) ? 1  : mesh->bpatch[jddl+1]-mesh->bpatch[jddl];
    for (i=0;i<n;++i)
      {	
	const I jelm = (jddl>=mesh->nvertex+mesh->nedge) ? jddl-(mesh->nvertex+mesh->nedge) : mesh->patch[mesh->bpatch[jddl]+i];	
	
	if (SF_)
	  {
	    nsSPACE_cncelm_dec(SF_,&jelm,&symbolic_cnc_[0],0);
	  }
	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nf],decu_[0]);
	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nf+nu],decu_[0]+SU_->nddl);
	nsSPACE_cncelm_dec(SP_,&jelm,&symbolic_cnc_[nf+nu*2],decu_[0]+SU_->nddl*2);	
	{ I j;
	  for (j=0;j<nddlelm;++j)
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] = ++N;
		  ddl_[N2+N-1]=symbolic_cnc_[j];		  
		}
	    } }
	
	if (SS_)
	  {
	    { I j;
	      for (j=0;j<nu;++j)
		{
		  if (SU_->ddlcod[symbolic_cnc_[nf+j]-decu_[0]]!=mesh->noboundary_vcod)
		    {
		      const I k  = symbolic_cnc_[nf+j]-decu_[0];
		      const I k1 = slip_iperm_[k]-1+total_nddl_[0];
		      if (NOT blank_[k1])
			{
			  blank_[k1] 	= ++N;
			  ddl_[N-1] 	= k1;
			}
		    }
		} }
	  }      	
      } }
  
  
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[N2+i]]=(I)0;
      } }
  if (nf>0)
    {
      { I i;    
	for (i=0;i<N2;++i)
	  {
	    ddl_[N2+i]+=SF_->nddl;
	  } }
    }

  return N+N2;
}



static I nsSPACE_galerkin_line_p( cst_nsSPACE 	SU_,
				       cst_nsSPACE 	SP_,
				       const I 		iddl_,
				       pI 			ddl_,
				       pI			blank_,
				       pI 			symbolic_cnc_,
				       cst_pI			decu_)
{
  const ns_mesh * mesh = nsSPACE_mesh(SU_);
  const I nu = nsSPACE_nddlelm(SU_);
  const I nddlelm =  nu*2+SP_->nddlelm;
  I N=0;


  const I jddl 		= nsSPACE_iddl2ddlP6(SU_,iddl_);
  { I i;
    const I n 	= (jddl>=mesh->nvertex+mesh->nedge) ? 1  : mesh->bpatch[jddl+1]-mesh->bpatch[jddl];
    for (i=0;i<n;++i)
      {	
	const I jelm = (jddl>=mesh->nvertex+mesh->nedge) ? jddl-(mesh->nvertex+mesh->nedge) : mesh->patch[mesh->bpatch[jddl]+i];	
	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[0],decu_[0]);
	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nu],decu_[0]+SU_->nddl);
	nsSPACE_cncelm_dec(SP_,&jelm,&symbolic_cnc_[nu*2],decu_[0]+SU_->nddl*2);
	{ I j;
	  for (j=0;j<nddlelm;++j)
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] = ++N;
		  ddl_[N-1]=symbolic_cnc_[j];
		}
	    } }
      } }
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[i]]=(I)0;
      } }
  return N;
}






static I nsSPACE_galerkin_line( cst_nsSPACE	S_,
				     const I 	iddl_,
				     pI 		ddl_,
				     pI		blank_,
				     pI 		symbolic_cnc_)
{
  const I nddlelm 		= nsSPACE_nddlelm(S_);
  const ns_mesh * mesh 		= nsSPACE_mesh(S_);
  I N=0;
  /**/
  const I jddl 		= nsSPACE_iddl2ddlP6(S_,iddl_);
#ifndef NDEBUG
  DebugVerif(jddl<S_->nddl);
  DebugVerif(jddl+1>0);
#endif
  { I i;
    const I n 	= (jddl>=mesh->nvertex+mesh->nedge) ? 1  : mesh->bpatch[jddl+1]-mesh->bpatch[jddl];
    for (i=0;i<n;++i)
      {	
	const I k = (jddl>=mesh->nvertex+mesh->nedge) ? jddl-(mesh->nvertex+mesh->nedge) : mesh->patch[mesh->bpatch[jddl]+i];	
#ifndef NDEBUG
	DebugVerif(k+1>0);
#endif
#ifndef NDEBUG
	DebugVerif(k<S_->mesh->nelm);
#endif
	nsSPACE_cncelm(S_,
			&k,
			symbolic_cnc_);
	/**/
	{ I j;
	  for (j=0;j<nddlelm;++j)
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] = ++N;
		  ddl_[N-1]   = symbolic_cnc_[j];
		}
	    } }
	/**/
      } }
  /**/
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[i]]=(I)0;
      } }
  /**/
  return N;
}





static I nsSPACE_galerkin_line_s( cst_nsSPACE SU_,
				       cst_nsSPACE SS_,
				       const I 	iddl_,
				       pI 		ddl_,
				       pI		blank_,
				       pI 		symbolic_cnc_,
				       cst_pI 		decu_)
{
  const I nu = nsSPACE_nddlelm(SU_);
  const I nddlu = nsSPACE_nddl(SU_);
  const ns_mesh * mesh = nsSPACE_mesh(SU_);
  I N=0;


  I jddl 		= nsSPACE_iddl2ddlP6(SU_,iddl_);
  { I i;
    const I n 	= (jddl>=mesh->nvertex+mesh->nedge) ? 1  : mesh->bpatch[jddl+1]-mesh->bpatch[jddl];
    for (i=0;i<n;++i)
      {	
	const I jelm = (jddl>=mesh->nvertex+mesh->nedge) ? jddl-(mesh->nvertex+mesh->nedge) : mesh->patch[mesh->bpatch[jddl]+i];	

	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[0],decu_[0]);
	nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nu],decu_[0]+nddlu);	
	{ I j;
	  for (j=0;j<nu;++j)
	    {
	      if (SU_->ddlcod[symbolic_cnc_[j]-decu_[0]]!=mesh->noboundary_vcod)
		{
		  if (NOT blank_[symbolic_cnc_[j]])
		    {
		      blank_[symbolic_cnc_[j]] 	= ++N;
		      ddl_[N-1]   		= symbolic_cnc_[j];
		    }
		}
	    } }
	
	{ I j;
	  for (j=0;j<nu;++j)
	    {
	      if (SU_->ddlcod[symbolic_cnc_[nu+j]-decu_[0]-nddlu]!=mesh->noboundary_vcod)
		{
		  if (NOT blank_[symbolic_cnc_[nu+j]])
		    {
		      
		      blank_[symbolic_cnc_[nu+j]] = ++N;
		      ddl_[N-1]   = symbolic_cnc_[nu+j];
		    }
		}
	    } }
      } }
  
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[i]]=(I)0;
      } }
  /*
    on rajoute un element diagonal pour pouvoir annuler la contrainte
  */
#if 0
  ddl_[N] = iddl_;
  ++N;
#endif
  return N;
}


static void nsSPACE_galerkin_selectdiscbis	(cst_nsSPACE 	SF_,
						 cst_nsSPACE 	SU_,
						 const I 	jelm,
						 const I 	n_,
						 pI 		iwork_n_,
						 pI 		blank_,
						 pI 		iwork_,
						 pI 		symbolic_cnc_,
						 cst_pI 	decu_)
{

  const I nf = nsSPACE_nddlelm(SF_);
  const I nu = nsSPACE_nddlelm(SU_);
  const I nddlu = nsSPACE_nddl(SU_);
  { I q;
    for (q=0;q<n_;++q)
      {
	symbolic_cnc_[q] = jelm*nf+q;
      } }

  nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nf],decu_[0]);
  nsSPACE_cncelm_dec(SU_,&jelm,&symbolic_cnc_[nu+nf],decu_[0]+nddlu);	
  { I q;
    for (q=0;q<n_;++q)
      {
	const I jddl = symbolic_cnc_[q];
	if (NOT blank_[jddl])
	  {
	    iwork_[iwork_n_[0]] = jddl;	      
	    iwork_n_[0]		+=1;
	    blank_[jddl]	= iwork_n_[0];
	  }
      } }
  { I p;
    for (p=0;p<iwork_n_[0];++p)
      {
	blank_[iwork_[p]]=0;
      } }
}

static void nsSPACE_galerkin_selectdisc	(cst_nsSPACE 	SF_,
						 cst_nsSPACE 	SU_,
						 cst_pI 		jelm_,
						 const I 		n_,
						 pI 			iwork_n_,
						 pI 			iwork_,						 
						 pI 			blank_,						 
						 pI 			symbolic_cnc_,
						 cst_pI 		decu_)
{
  const ns_mesh * mesh = nsSPACE_mesh(SF_);
  const I nf = nsSPACE_nddlelm(SF_);
  { I q;
    for (q=0;q<3;++q)
      {
	if (mesh->adj[jelm_[0]*3+q])
	  {
	    const I jelm = mesh->adj[jelm_[0]*3+q]-1;
	    { I p;
	      for (p=0;p<nf;++p)
		{
		  iwork_[iwork_n_[0]] = jelm*nf+p;	      
		  iwork_n_[0]+=1;
		  blank_[jelm*nf+p]=iwork_n_[0];
		} }
	  }
      } }
  nsSPACE_galerkin_selectdiscbis(SF_,SU_,jelm_[0],
				  n_,
				  iwork_n_,
				  blank_,
				  iwork_,
				  symbolic_cnc_,
				  decu_);
}





void nsSPACE_divergence(cst_nsSPACE 	SU_,
			cst_nsSPACE 	SP_,
			cst_pI 		axi_,
			pSparseStokes const	S_,
			cst_pI		decu_)
{ 
#if 0
  static const I dimension 	= 2;
  const I 	nddlu 		= nsSPACE_nddl(SU_);
  const I 	nu 		= nsSPACE_nddlelm(SU_);
  const I 	np 		= nsSPACE_nddlelm(SP_);
  cst_ensBASIS shape_u 	= nsSPACE_shape(SU_);
  cst_ensBASIS shape_p 	= nsSPACE_shape(SP_);
  const ns_mesh * mesh  	= nsSPACE_mesh(SU_);
  const I 	nuxdim		= nu*dimension;

  pI blank 			= calloc(S_->A.n,sizeof(I));


  R 	dsu_p[32*32];
  R 	dru_p[32*32];
  R 	locmat_BU[64*32];
  R 	locmat_UB[64*32];
  R 	build_dsu_p[32*32*3];
  R 	build_dru_p[32*32*3];
  I 	locnumeru[64];
  I 	locnumerp[32];
  R 	belm[4],sbelm[4],alpha2;
  const I nuxnp 		= nu*np;  
  nsCUBA_TRIA_ST * cub 	= nsCUBA_TRIA_gauss_new(10);
  if (axi_[0]>0)
    {
      nsCUBA_TRIA_matrix_elm_dr_f	(cub,
						 shape_u,
						 shape_p,
						 build_dru_p,
						 &nuxnp);      
      nsCUBA_TRIA_matrix_elm_ds_f	(cub,
						 shape_u,
						 shape_p,
						 build_dsu_p,
						 &nuxnp);      
    }
  else
    {
      nsCUBA_TRIA_matrix_dr_f		(cub,
						 shape_u,
						 shape_p,
						 dru_p,
						 &np);
      nsCUBA_TRIA_matrix_ds_f		(cub,
						 shape_u,
						 shape_p,
						 dsu_p,
						 &np);
    }
  R cooelm[32];
  { I ielm;    
    for (ielm=0;ielm<mesh->nelm;++ielm)
      { 	
	nsSPACE_cncelm_dec(SU_,&ielm,locnumeru,decu_[0]);
	nsSPACE_cncelm_dec(SU_,&ielm,&locnumeru[nu],decu_[0]+nddlu);
	nsSPACE_cncelm_dec(SP_,&ielm,locnumerp,decu_[0]+nddlu*dimension);		
	ns_mesh_cooelm(mesh,&ielm,cooelm);
	alpha2 		= mesh->jacelm[ielm];
	belm[3]   	= mesh->trelm[4*ielm+3];
	belm[1]   	= mesh->trelm[4*ielm+1];
	belm[2]   	= mesh->trelm[4*ielm+2];
	belm[0]   	= mesh->trelm[4*ielm+0];
	sbelm[3]   	= belm[3]*alpha2;
	sbelm[1]   	= belm[1]*alpha2;
	sbelm[2]   	= belm[2]*alpha2;
	sbelm[0]   	= belm[0]*alpha2;				
	static const I negal1 = ((I)1);
	static const I negal3 = ((I)3);
	static const R regal0 = ((R)0.0);
	static const R regal1 = ((R)1.0);
	if (axi_[0]>0)
	  {
	    nsblas_dgemv(transN,&nuxnp,&negal3,&regal1,build_dru_p,&nuxnp,&cooelm[(axi_[0]-1)*3],&negal1,&regal0,dru_p,&negal1);
	    nsblas_dgemv(transN,&nuxnp,&negal3,&regal1,build_dsu_p,&nuxnp,&cooelm[(axi_[0]-1)*3],&negal1,&regal0,dsu_p,&negal1);
	  }	
	{ I i;
	  for (i=0;i<np;++i)
	    {
	      I j;
	      for (j=0;j<nu;++j)
		{
		  locmat_BU[j*np+i]    	= sbelm[0] * dru_p[j*np+i]+sbelm[2] * dsu_p[j*np+i];
		  locmat_UB[i*nuxdim+j] = -locmat_BU[j*np+i];		
		}
	      for (j=nu;j<nuxdim;++j)
		{
		  locmat_BU[j*np+i] 	= sbelm[1] * dru_p[(j-nu)*np+i]+sbelm[3] * dsu_p[(j-nu)*np+i];		
		  locmat_UB[i*nuxdim+j] = -locmat_BU[j*np+i];		
		}
	    } }   
	nsSPARSE_BLOCK_ass_blank(&S_->block_BU,&np,locnumerp,&nuxdim,locnumeru,locmat_BU,&np,blank);
	nsSPARSE_BLOCK_ass_blank(&S_->block_UB,&nuxdim,locnumeru,&np,locnumerp,locmat_UB,&nuxdim,blank);
      } }
  free(blank);
#endif

}





pSparse nsSPACE_galerkin_stokes_def(cst_nsSPACE 	SF_,
				    cst_nsSPACE 	SU_,
				    cst_nsSPACE 	SP_,
				    cst_nsSPACE 	SS_,
				    const L  	F_continuity_,
				    pI		slip_iperm_,
				    pI		nddl_slip_)
{
#if 0
  memset(A_,0,sizeof(nsSPARSE_ST));
#endif
  const ns_mesh * mesh 	= nsSPACE_mesh(SU_);
  I An 		= ((I)0);
  I Am 		= ((I)0);
  I Ancoeff 	= ((I)0);
  pI 	Ab 		= NULL;
  pI 	Ai 		= NULL;
  const I 	nddlu 		= nsSPACE_nddl(SU_);
  const I 	nddlp 		= nsSPACE_nddl(SP_);
  const I 	nddlf 		= (SF_)?nsSPACE_nddl(SF_):0;
  const I 	nf 		= (SF_)?nsSPACE_nddlelm(SF_):0;
  const I 	nu 		= nsSPACE_nddlelm(SU_);
  const I 	np 		= nsSPACE_nddlelm(SP_);

  if (slip_iperm_)
    {
      { I iddl,jddl=0;  
	for (iddl=0;iddl<nddlu;++iddl)
	  {      
	    if (SU_->ddlcod[iddl]!=mesh->noboundary_vcod)
	      {
		slip_iperm_[iddl] = jddl+1;
		++jddl;
	      }
	  } }
    }

  I 	nline,decu=0,decv=0,decp=0;

  I 	symbolic_cnc	[512];
  
  decu      	= nddlf;
  decv      	= decu+nddlu;
  decp      	= decv+nddlu;
  An 		= decp+nddlp;  

  const I total_nddlelm 	= nf+2*nu+np;
  const I An_without_slip 	= An;

  if (SS_)
    {
      An 	+= mesh->nedge_boundary * ensBASIS_n_onface(SS_->shape,mesh->elm);
    }

  Am		= An;
  Ab 		= malloc((An+1)*sizeof(I));
  if (SS_)
    {
      
      Ai = malloc(total_nddlelm*total_nddlelm*mesh->nelm+2*mesh->nedge_boundary*(12+1)*sizeof(I));     /*+1 pour l'element diagonal du slip */
    }
  else
    {
      Ai = malloc((total_nddlelm+1)*(total_nddlelm+1)*mesh->nelm*sizeof(I));
      printf("size "ifmt"\n",(total_nddlelm+1)*(total_nddlelm+1)*mesh->nelm);
    }


  pI  blank = calloc(An,sizeof(I));

  Ab[0] 	= 0;
  Ancoeff 	= 0;
  if (SF_)
    {
      if (F_continuity_)
	{
	  { I iddl;
	    for (iddl=0;iddl<nddlf;++iddl)
	      {
		nline = nsSPACE_galerkin_line(SF_,
					       iddl,
					       &Ai[Ancoeff],
					       blank,
					       symbolic_cnc);	       
		Ancoeff +=nline;
		Ab[iddl+1] = Ab[iddl] + nline;
	      } }
	}
      else
	{
	  pI iwork = malloc(sizeof(I)*An);
	  { I jelm;
	    for (jelm=0;jelm<mesh->nelm;++jelm)
	      {
		{ I j;
		  for (j=0;j<nf;++j)
		    {
		      const I iddl = jelm*nf+j;
		      I nselection = 0;	      
		      nsSPACE_galerkin_selectdisc(SF_,SU_,&jelm,total_nddlelm-np,&nselection,blank,iwork,symbolic_cnc,&decu);      	      
		      Ab[iddl+1] = Ab[iddl]+nselection;
		      memcpy(&Ai[Ancoeff],symbolic_cnc,sizeof(I)*nselection);
		      Ancoeff += nselection;
		    } }
	      } }
	}
    }

  { I iddl;	
    for (iddl=0;iddl<nddlu;++iddl)
      {
	const I nline = nsSPACE_galerkin_line_fups(SF_,SU_,SP_,SS_,iddl,&Ai[Ancoeff],blank,symbolic_cnc,slip_iperm_,&decu,&An_without_slip);
	Ancoeff +=nline;
	Ab[decu+iddl+1] = Ab[decu+iddl] + nline;
      } }  
  { I iddl;
    for (iddl=0;iddl<nddlu;++iddl)
      {
	const I nline = nsSPACE_galerkin_line_fups(SF_,SU_,SP_,SS_,iddl,&Ai[Ancoeff],blank,symbolic_cnc,slip_iperm_,&decu,&An_without_slip);
	Ancoeff +=nline;
	Ab[decv+iddl+1] = Ab[decv+iddl] + nline;
      } }
  { I iddl;
    for (iddl=0;iddl<nddlp;++iddl)
      {
	const I nline = nsSPACE_galerkin_line_p(SU_,SP_,iddl,&Ai[Ancoeff],blank,symbolic_cnc,&decu);
	Ancoeff +=nline;
	Ab[decp+iddl+1] = Ab[decp+iddl] + nline;
      } }
  
  
  if (SS_)
    {
      { I iddl,jddl = 0;
	for (iddl=0;iddl<nddlu;++iddl)
	  {
	    if (SU_->ddlcod[iddl]!=mesh->noboundary_vcod)
	      {
		nline = nsSPACE_galerkin_line_s(SU_,SS_,iddl,&Ai[Ancoeff],blank,symbolic_cnc,&decu);
		Ancoeff +=nline;
		/* on ajoute un element diagonal pour pouvoir annuler la condition */
		
		/* s'il n'y est pas normalement ca doit gueuler quand on vuet appliquer du dirichlet nul sur le slip*/
		Ai[Ancoeff]=An_without_slip+jddl;/* la notation est C pas fortran, la version fortran vient apres*/
		++Ancoeff;
		++nline;
		Ab[An_without_slip+jddl+1] = Ab[An_without_slip+jddl] + nline;
		++jddl;
	      }
	  } }      
      nddl_slip_[0] = An-An_without_slip;
    }
  else
    {
      nddl_slip_[0] = 0;
    }


  printf("Ancoeff = "ifmt"\n",Ancoeff);
  Ai    	= realloc(Ai,Ancoeff*sizeof(I));  
  free(blank);

  pSparse A  = Sparse_build(An,Am,Ancoeff,Ab,Ai,NULL);

#if 0
  A_->own_x    	= malloc(Ancoeff*sizeof(R));
  A_->own_b    	= Ab;
  A_->b 	= A_->own_b;
  A_->x 	= A_->own_x;
  A_->i 	= A_->own_i;
  A_->n    	= An;
  A_->m    	= Am;  
  A_->nc   	= Ancoeff;
#endif
  Sparse_sort(A);
  return A;
}


void SparseStokes_def(pSparseStokes const 	S_,
			 cst_nsSPACE 		SU_,
			 cst_nsSPACE 		SP_,
			 cst_nsSPACE 		SS_,
			 cst_pI 			axi_,
			 pI				slip_iperm_,
			 pI				nddl_slip_)
{ 
  memset(S_,0,sizeof(SparseStokes));  
  S_->A  = nsSPACE_galerkin_stokes_def(NULL,
				       SU_,
				       SP_,
				       SS_,
				       __emnsNO,
				       slip_iperm_,
				       nddl_slip_);  

  Sparse_fortran_indexation(S_->A);  

  static const I dim 	= 2;
  const I nddlf 	= 0;
  const I dec_ddlu 	= nddlf;
  const I dec_ddlv 	= dec_ddlu+SU_->nddl;
  const I dec_ddlp 	= dec_ddlv+SU_->nddl;
  const I total	= dec_ddlp+SP_->nddl;  
  const I nddluxdim	= SU_->nddl*dim;  
  { const I start_i_block_UU 	= dec_ddlu;
    const I start_j_block_UU 	= dec_ddlu;
    const I start_i_block_BU 	= dec_ddlp;
    const I start_j_block_BU 	= dec_ddlu;
    const I start_i_block_UB 	= dec_ddlu;
    const I start_j_block_UB 	= dec_ddlp;
    const I start_i_block_BB 	= dec_ddlp;
    const I start_j_block_BB 	= dec_ddlp;      
    SparseBlock_def	(&S_->block_UU,
			 S_->A,
			 &start_i_block_UU,
			 &nddluxdim,
			 &start_j_block_UU,
			 &nddluxdim); 
    SparseBlock_def	(&S_->block_BU,
			 S_->A,
			 &start_i_block_BU,
			 &SP_->nddl,
			 &start_j_block_BU,
			 &nddluxdim);            
    SparseBlock_def	(&S_->block_UB,
			 S_->A,
			 &start_i_block_UB,
			 &nddluxdim,
			 &start_j_block_UB,
			 &SP_->nddl);
    SparseBlock_def	(&S_->block_BB,
			 S_->A,
			 &start_i_block_BB,
			 &SP_->nddl,
			 &start_j_block_BB,
			 &SP_->nddl);  }
  if (SS_)
    {
      const I start_i_block_SS = total;
      const I start_i_block_SU = total;
      const I start_j_block_SU = dec_ddlu;	
      SparseBlock_def	(&S_->block_SS,
			 S_->A,
				 &start_i_block_SS,
				 nddl_slip_,
				 &start_i_block_SS,
				 nddl_slip_); 
      SparseBlock_def	(&S_->block_SU,
			 S_->A,
			 &start_i_block_SU,
				 nddl_slip_,
				 &start_j_block_SU,
				 &nddluxdim);    
      SparseBlock_def	(&S_->block_US,
			 S_->A,
				 &start_j_block_SU,
				 &nddluxdim,
				 &start_i_block_SU,
				 nddl_slip_); 
    }
  
  nsSPACE_divergence(SU_,
		      SP_,
		      axi_,
		      S_,
		      &dec_ddlu);

  if (SS_)
    {
      
    }
  /*
   */
}


#if 0
Err nsSPARSE_stokes_jacobian_transport	(ns_sparse_stokes*		S_,
						 cst_nsSPACE 		SF_,
						 cst_nsSPACE 		SU_,
						 const I 			elm_,
						 pR 				jacobian_,
						 cst_pI 			jacobianoff_)
{ 
  if (SF_)
    {
      {	I i;
	for (i=0;i<q_n;++i)
	  {
	    { I i;
	      for (i=0;i<n;++i)
		{
		  phi[i] = ref_phi[i] + 0.01*(ug[i] * dxfg[i] + ug[i] * dyfg[i]);
		} }
	    nsblas_dger(jacobian_,&a,phi,psi);	    
	    nsblas_dger(jacobian_,&b,phi,dxpsi);	    
	    nsblas_dger(jacobian_,&b,phi,dypsi);	    
	  } }
    }
}
#endif

#if 0
Err nsSPARSE_STOKES_def_eulerian(nsSPARSE_STOKES_ST*const 	S_,
				      cst_nsSPACE const 		SF_,
				      cst_nsSPACE const 		SU_,
				      cst_nsSPACE const 		SP_,
				      cst_nsSPACE const 		SS_,
				      const L  		F_continuity_,
				      cst_pI 			axi_,
				      pI			slip_iperm_,
				      pI			nddl_slip_)
{ 

  memset(S_,0,sizeof(nsSPARSE_STOKES_ST));  

  Err err =  nsSPACE_galerkin_stokes_def(&S_->A,
					      SF_,
					      SU_,
					      SP_,
					      SS_,
					      F_continuity_,
					      slip_iperm_,
					      nddl_slip_);  
  if (err)
    {
      ns_errmsg("nsSPARSE_STOKES_def_eulerian:nsSPACE_galerkin_stokes_def failed");
      return err;
    }
  nsSPARSE_fortran_indexation(&S_->A);  
  static const I dim = 2;
  const I nelm	= ns_mesh_nelm(nsSPACE_mesh(SU_));
  const I nddlf 	= (SF_) ? ( (F_continuity_) ? SF_->nddl : nelm * SF_->nddlelm ) : 0;
  const I dec_ddlu 	= nddlf;
  const I dec_ddlv 	= dec_ddlu+SU_->nddl;
  const I dec_ddlp 	= dec_ddlv+SU_->nddl;
  const I total	= dec_ddlp+SP_->nddl;  
  const I nddluxdim	= SU_->nddl*dim;  
  { const I start_i_block_UU 	= dec_ddlu;
    const I start_j_block_UU 	= dec_ddlu;
    const I start_i_block_BU 	= dec_ddlp;
    const I start_j_block_BU 	= dec_ddlu;
    const I start_i_block_UB 	= dec_ddlu;
    const I start_j_block_UB 	= dec_ddlp;
    const I start_i_block_BB 	= dec_ddlp;
    const I start_j_block_BB 	= dec_ddlp;      

    SparseBlock_def	(&S_->block_UU,
			 &S_->A,
			 &start_i_block_UU,
			 &nddluxdim,
			 &start_j_block_UU,
			 &nddluxdim); 

    SparseBlock_def	(&S_->block_BU,
			 &S_->A,
			 &start_i_block_BU,
			 &SP_->nddl,
			 &start_j_block_BU,
			 &nddluxdim);            

    SparseBlock_def	(&S_->block_UB,
			 &S_->A,
			 &start_i_block_UB,
			 &nddluxdim,
			 &start_j_block_UB,
			 &SP_->nddl);

    SparseBlock_def	(&S_->block_BB,
			 &S_->A,
			 &start_i_block_BB,
			 &SP_->nddl,
			 &start_j_block_BB,
			 &SP_->nddl);  }
  if (SS_)
    {
      const I start_i_block_SS = total;
      const I start_i_block_SU = total;
      const I start_j_block_SU = dec_ddlu;	
      SparseBlock_def	(&S_->block_SS,
				 &S_->A,
				 &start_i_block_SS,
				 nddl_slip_,
				 &start_i_block_SS,
				 nddl_slip_); 
      SparseBlock_def	(&S_->block_SU,
				 &S_->A,
				 &start_i_block_SU,
				 nddl_slip_,
				 &start_j_block_SU,
				 &nddluxdim);    
      SparseBlock_def	(&S_->block_US,
				 &S_->A,
				 &start_j_block_SU,
				 &nddluxdim,
				 &start_i_block_SU,
				 nddl_slip_); 
    }
  nsSPACE_divergence(SU_,
		      SP_,
		      axi_,
		      S_,
		      &dec_ddlu);
  return err;
}
#endif








pSparse nsSPACE_galerkin_new(cst_nsSPACE	S_)
{
  const ns_mesh * mesh_ = nsSPACE_mesh(S_);


  const I n 	= S_->nddl;
  const I m 	= S_->nddl;

  pI own_b 	= malloc((n+1)*sizeof(I));
  pI own_i	= malloc(S_->nddlelm*S_->nddlelm*mesh_->nelm*sizeof(I));
  own_b[0] 	= 0;
  I nc 		= 0;
#if 0
  F_->n		= S_->nddl;
  F_->m    	= S_->nddl;
  F_->own_b	= malloc((F_->n+1)*sizeof(I));
  F_->own_i	= malloc(S_->nddlelm*S_->nddlelm*mesh_->nelm*sizeof(I));
  F_->own_b[0] 	= 0;
  F_->nc 	= 0;
#endif

  I symbolic_cnc[128];
  I	nline;
  pI blank = calloc(S_->nddl,sizeof(I));
  { I iddl;
    for (iddl=0;iddl<S_->nddl;++iddl)
      {
	nline = nsSPACE_galerkin_line	(S_,
					 iddl,
					 &own_i[nc],
					 blank,
					 symbolic_cnc);

	nc 		+=nline;
	own_b[iddl+1] 	= own_b[iddl] + nline;
      } }
  free(blank);
  own_i	= realloc(own_i,nc*sizeof(I));  

  pSparse back = Sparse_build(n,
			      m,
			      nc,
			      own_b,
			      own_i,
			      NULL);

  Sparse_sort(back);
  Sparse_fortran_indexation(back);  
  return back;
}












#if 0
void nsSPACE_galerkin_add_mass_ns_var(cst_nsSPACE 	S_,
				       nsSPARSE_ST*		A_,
				       const ns_var* 		v_,
				       nsCUBA_TRIA_ST * 	c,
				       cst_pI 			rwork_n_,
				       pR			rwork_,
				       pI			iwork_,
				       pI			blank_ddl_)
{
  const I req_size = S_->nddlelm*S_->nddlelm+c->q_n+v_->nddlelm;
  if (rwork_n_[0]<req_size)
    {
      ns_err("nsSPACE_galerkin_add_mass_ns_var:not enough space");      
    }
  pR mat 		= rwork_;
  pR v_wei 		= &rwork_[S_->nddlelm*S_->nddlelm];
  pR v_gauss 		= &rwork_[S_->nddlelm*S_->nddlelm+v_->nddlelm];

  const ns_mesh * mesh = nsSPACE_mesh(S_);
  { I ielm;
    for (ielm=0;ielm<mesh->nelm;++ielm)
      {	

	{ R cooelm[6];
	  ns_mesh_cooelm		(mesh,
					 &ielm,
					 cooelm);
#if 0
	  pas la peine ici
	  nsCUBA_TRIA_setelm	(c,
					 v_->lagr_id,
					 cooelm,
					 mesh->trelm);
	  nsCUBA_TRIA_setelm	(c,
					 S_->shape,
					 cooelm,
					 mesh->trelm); 
#endif
	} 

	{ ns_var_welm			(v_,
					 &ielm,
					 v_wei,
					 &v_->nddlelm);
	  nsCUBA_TRIA_wei	(c,
					 v_->lagr_id,
					 v_->ortho,
					 &negal1,
					 &mesh->jacelm[ielm],
					 v_wei,
					 &v_->nddlelm,
					 &regal0,
					 v_gauss,
					 &negal1); }
	  
	
	{ memset(mat,0,sizeof(R)*S_->nddlelm*S_->nddlelm);	
	  { const I q_n = c->q_n;
	    I j;
	    for (j=0;j<q_n;++j)
	      {
		nsblas_dger		(&S_->nddlelm,
					 &S_->nddlelm,
					 &v_gauss[j],
					 c->wqelm_lagr[S_->shape],
					 &negal1,
					 c->qelm_lagr[S_->shape],	       
					 &q_n,
					 mat,
					 &S_->nddlelm);
	      } } }
	
	nsSPARSE_assmatelm	(A_,
				 &S_->nddlelm,
				 &S_->ddlcnc[ielm*S_->nddlelm],
				 mat,
				 &S_->nddlelm,
				 iwork_,
				 blank_ddl_);
	
      } }   
}
#endif



cst_nsSPACE ns_mesh_get_space(const ns_mesh*mesh_,
				 cst_ensBASIS shape_)
{
  return mesh_->spaces[shape_];
}

cst_nsSPACE ns_mesh_addspace	(ns_mesh*mesh_,
				 cst_ensBASIS shape_)
{
  if (mesh_->spaces[shape_])
    return mesh_->spaces[shape_];
  else
    {
      mesh_->spaces[shape_] = nsSPACE_new(mesh_,shape_);
      return mesh_->spaces[shape_];
    }
}

void  nsSPACE_get_ddlcoo(cst_nsSPACE	s_,
			 cst_pI 		nddl_,
			 cst_pI 		ddl_,
			 pR 			pos_,
			 cst_pI 		posoff_)
{
  const I dim = ns_mesh_dimension(s_->mesh);
  { I i;
    for (i=0;i<nddl_[0];++i)
      {
	{ I j;
	  for (j=0;j<dim;++j)
	    {
	      pos_[i*posoff_[0]+j] = s_->ddlcoo[dim*ddl_[i]+j];
	    } }
      } }
}



void nsSPACE_write_medit(cst_nsSPACE s_,const char * name_,...)
{
  { 
    char ctmp[512];

    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    FILE * fil = fopen(ctmp,"w");
    const I nddl 		= nsSPACE_nddl(s_);
    cst_pR ddlcoo 		= nsSPACE_ddlcoo(s_);
    cst_pI ddlcod 		= nsSPACE_ddlcod(s_);
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n"ifmt"\n",nddl);
    {I i;for (i=0;i<nddl;++i){ fprintf(fil,""rfmt" "rfmt" "ifmt"\n",ddlcoo[i*2+0],ddlcoo[i*2+1],ddlcod[i]); } }
    cst_pI 	tria_cnc 	= ensBASIS_cnctreilli(nsSPACE_shape(s_),nsSPACE_elm(s_));
    const I tria_n 		= ensBASIS_nelmtreilli(nsSPACE_shape(s_),nsSPACE_elm(s_));
    const I nelm 		= ns_mesh_nelm(nsSPACE_mesh(s_));

    fprintf(fil,"Triangles\n"ifmt"\n",nelm*tria_n); 
    I cncelm[128];
    { I i;
      for (i=0;i<nelm;++i)
	{
	  nsSPACE_cncelm(s_,&i,cncelm);
#if 0
	  I cncelm[128];
	  { I j;    
	    for (j=0;j<s_->nddlelm;++j)
	      {
		cncelm[j] = s_->mesh->cnc[i*s_->ddlcncoff+j];		
	      } }
	  printf(" %lld %lld %lld %lld %lld %lld %lld\n",s_->ddlcncoff,cncelm[0],cncelm[1],cncelm[2],cncelm[3],cncelm[4],cncelm[5]);
#endif
#if 0
	  if (s_->nddlelm==3)
	    {
	      { I j;
		for (j=0;j<tria_n;++j)
		  {
		    fprintf(fil,""ifmt" "ifmt" "ifmt" "ifmt"\n",cncelm[0]+1,cncelm[1]+1,cncelm[2]+1,((I)0));
		  } }
	    }
	  else
	    {
	    }
#endif
	  { I j;
	    for (j=0;j<tria_n;++j)
	      {
		fprintf(fil,""ifmt" "ifmt" "ifmt" "ifmt"\n",cncelm[tria_cnc[j*3+0]]+1,cncelm[tria_cnc[j*3+1]]+1,cncelm[tria_cnc[j*3+2]]+1,((I)0));
	      } }
	  
	} } 

    fprintf(fil,"End\n");						
    fclose(fil); }  
}



void nsSPACE_dg_write_medit(cst_nsSPACE s_,const char * name_,...)
{
  char ctmp[512];
  const ns_mesh*mesh= nsSPACE_mesh(s_);
  
  { char ctmp2[512];
    va_list args;
    va_start (args,name_);
    vsprintf(ctmp2,name_,args);
    va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }
  cst_ensBASIS 	shape_ 	= nsSPACE_shape(s_);
  cst_eElement 	elm_ 	= nsSPACE_elm(s_);
  R glcoo[128];
  const I subnelm 	= ensBASIS_nelmtreilli(shape_,elm_);
  const I n 	= ensBASIS_n(shape_,elm_);
  cst_pR lcoo 		= ensBASIS_cootreilli(shape_,elm_);
  const I nelm 	= ns_mesh_nelm(mesh);
  FILE * fich_ = fopen(ctmp,"w");
  fprintf(fich_,"MeshVersionFormatted 1\nDimension\n2\n");
  fprintf(fich_,"\nVertices\n" ifmt "\n",subnelm*nelm*3);
  cst_pI split = ensBASIS_cnctreilli(shape_,elm_);
  R elmcoo[6];
  { I j; 
    for (j=0;j<nelm;++j)
      {	
	ns_mesh_cooelm(mesh,&j,elmcoo);
	{ I k;
	  for (k=0;k<n;++k)
	    {
	      glcoo[k] 		= (((R)1.0)-lcoo[k]-lcoo[n+k])*elmcoo[0]+lcoo[k]*elmcoo[1]+lcoo[k+n]*elmcoo[2];
	      glcoo[n+k] 	= (((R)1.0)-lcoo[k]-lcoo[n+k])*elmcoo[3]+lcoo[k]*elmcoo[4]+lcoo[k+n]*elmcoo[5];
	    } }       
	{ I k;
	  for (k=0;k<subnelm;++k)
	    {
	      fprintf(fich_,"%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n",
		      glcoo[split[3*k+0]],
		      glcoo[n+split[3*k+0]],
		      glcoo[split[3*k+1]],
		      glcoo[n+split[3*k+1]],
		      glcoo[split[3*k+2]],
		      glcoo[n+split[3*k+2]]);		      
	    } }	
      } }
  fprintf(fich_,"\nTriangles\n" ifmt "\n",subnelm*nelm);      
  { I j;
    for (j=0;j<nelm*subnelm;++j)
      {	
	fprintf(fich_,"" ifmt " " ifmt " " ifmt " 0\n",j*3+1,j*3+2,j*3+3);	    
       } }
  fprintf(fich_,"\nEnd\n");
  fclose(fich_);

}


void 		nsSPACE_free	(nsSPACE_ST*const s_)
{
  if (s_)
    {
      if (s_->own_ddlcod) free(s_->own_ddlcod);
      if (s_->own_ddlcoo) free(s_->own_ddlcoo);
      if (s_->own_ddlcnc) free(s_->own_ddlcnc);
#if 0
      if (s_->own_cnc) free(s_->own_cnc);
      if (s_->own_bpatch) free(s_->own_patch);
      if (s_->own_patch) free(s_->own_bpatch);
#endif
      s_->endomorphism = Sparse_kill(s_->endomorphism);
      memset(s_,0,sizeof(nsSPACE_ST));
    }
}

nsSPACE_ST*	nsSPACE_kill	(nsSPACE_ST*const s_)
{
  if (s_)
    {
      nsSPACE_free(s_);
      free(s_);
    }
  return NULL;
}

void 		nsSPACE_def 	(nsSPACE_ST*		space_,  
				 const ns_mesh * 	mesh_,
				 cst_ensBASIS 	shape_)
{
#ifndef NDEBUG
  DebugVerif(space_);
  DebugVerif(mesh_);
  DebugVerif(shape_);
#endif

  memset(space_,0,sizeof(nsSPACE_ST));

  space_->mesh		= mesh_;
  space_->shape		= shape_;
  space_->nddl 		= ns_mesh_nddlspace(mesh_,shape_);
  space_->nddlelm	= ensBASIS_n(shape_,mesh_->elm);
  space_->nddledge 	= ensBASIS_n_onface(space_->shape,mesh_->elm);

  switch(space_->shape)
    {
    case __ensBASIS_LAGRANGEBUBBLE_0:
      {
	space_->own_ddlcnc 	= malloc(sizeof(R)*mesh_->nelm*2);
	space_->ddlcnc 		= space_->own_ddlcnc;
	space_->ddlcncoff 	= 2;
	{ I i;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      space_->own_ddlcnc[2*i+0] = i;
	      space_->own_ddlcnc[2*i+1] = mesh_->nelm+i;
	    } }
	break;
      }
    case __ensBASIS_L2ORTHONORMAL_0:
    case __ensBASIS_LAGRANGE_0:
      {
	space_->own_ddlcnc 	= malloc(sizeof(R)*mesh_->nelm);
	space_->ddlcnc 		= space_->own_ddlcnc;
	space_->ddlcncoff 	= 1;
	{ I i;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      space_->own_ddlcnc[i] = i;
	    } }
	break;
      }
    case __ensBASIS_LAGRANGE_3:
    case __ensBASIS_L2ORTHONORMAL_3:
    case __ensBASIS_LAGRANGEBUBBLE_2:
      {
	space_->own_ddlcnc 	= malloc(sizeof(I)*mesh_->nelm*space_->nddlelm);
	space_->ddlcnc 		= space_->own_ddlcnc;
	space_->ddlcncoff 	= space_->nddlelm;
	{ I i;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      { I j;
		for (j=0;j<3;++j)
		  {
		    space_->own_ddlcnc[space_->nddlelm*i+j] = mesh_->cnc[6*i+j];
		  } }
	    } }
	{ I i;
	  for (i=0;i<mesh_->nedge;++i)
	    {
	      { I nadj,ielm,jelm,ilocedge,jlocedge;
		ns_mesh_edgeinfo(mesh_,
				 &i,
				 &nadj,
				 &ielm,
				 &ilocedge,
				 &jelm,
				 &jlocedge);
		{ I j;
		  for (j=0;j<space_->nddledge;++j)
		    {
		      space_->own_ddlcnc[space_->nddlelm*ielm+3+space_->nddledge*ilocedge+j] = mesh_->nvertex + space_->nddledge*i + j+1;
		    } }
		if (nadj>1)
		  {
		    { I j;
		      for (j=0;j<space_->nddledge;++j)
			{
			  space_->own_ddlcnc[space_->nddlelm*jelm+3+space_->nddledge*jlocedge+space_->nddledge-1-j] = mesh_->nvertex + space_->nddledge*i + j+1;
			} }
		  } 				
	      }
	    } }
	{ I i;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      const I n = space_->nddlelm -(3+space_->nddledge*3);
	      { I j;
		for (j=0;j<n;++j)
		  {
		    space_->own_ddlcnc[space_->nddlelm*i+3+space_->nddledge*3+j] = mesh_->nvertex+space_->nddledge*mesh_->nedge+n*i+j+1;
		  } }
	    } }
	break;
      }

    case __ensBASIS_LAGRANGEBUBBLE_1:
      {
	break;
      }

    case __ensBASIS_L2ORTHONORMAL_1:
    case __ensBASIS_LAGRANGE_1:
    case __ensBASIS_L2ORTHONORMAL_2:
    case __ensBASIS_LAGRANGE_2:
      {
	space_->ddlcncoff  	= 6;
	space_->ddlcnc 		= mesh_->cnc;
	break;
      }  
    case __ensBASIS_n:
    case __ensBASIS_error:
      {
	break;
      }  
    }

  switch(space_->shape)
    {
    case __ensBASIS_LAGRANGE_3:
    case __ensBASIS_L2ORTHONORMAL_3:
    case __ensBASIS_LAGRANGEBUBBLE_2:
      {
	space_->own_ddlcoo 	= malloc(sizeof(R)*space_->nddl*2);
	space_->ddlcoo 		= space_->own_ddlcoo;
	{ I i;
	  R 	cooelm	[32];
	  R 	gl	[64*2];
	  I 	cnc	[64];
	  static const I negal3=3;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      ns_mesh_cooelm			(mesh_,&i,cooelm);
	      ensBASIS_lc2gltreilli		(shape_,mesh_->elm,cooelm,&negal3,gl,&space_->nddlelm);
	      nsSPACE_cncelm			(space_,&i,cnc);
	      { I j;
		for (j=0;j<space_->nddlelm;++j)
		  {
		    space_->own_ddlcoo[2*cnc[j]+0]	= gl[j];
		    space_->own_ddlcoo[2*cnc[j]+1]	= gl[space_->nddlelm+j];
		  } }
	    } }
	break;
      }  
    case __ensBASIS_LAGRANGE_0:
    case __ensBASIS_L2ORTHONORMAL_0:
      {
	space_->own_ddlcoo 	= malloc(sizeof(R)*mesh_->nelm*2);
	space_->ddlcoo 		= space_->own_ddlcoo;
	{ I i;
	  R 	cooelm	[32];
	  static const R a = ((R)1.0)/((R)3.0);
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      ns_mesh_cooelm(mesh_,&i,cooelm);
	      space_->own_ddlcoo[2*i+0]	= a*cooelm[0]+a*cooelm[1]+a*cooelm[2];
	      space_->own_ddlcoo[2*i+1]	= a*cooelm[3]+a*cooelm[4]+a*cooelm[5];
	    } }
	break;
      }  
    case __ensBASIS_LAGRANGE_1:
    case __ensBASIS_L2ORTHONORMAL_1:
    case __ensBASIS_LAGRANGE_2:
    case __ensBASIS_L2ORTHONORMAL_2:
      {
	space_->ddlcoo = mesh_->coo;
	break;
      }  
    case __ensBASIS_LAGRANGEBUBBLE_0:
    case __ensBASIS_LAGRANGEBUBBLE_1:
    case __ensBASIS_n:
    case __ensBASIS_error:
      {
	break;
      }  
    }
  
  
  
  
  
  switch(space_->shape)
    {
    case __ensBASIS_LAGRANGE_1:
    case __ensBASIS_L2ORTHONORMAL_1:
    case __ensBASIS_LAGRANGE_2:
    case __ensBASIS_L2ORTHONORMAL_2:
      {
	space_->ddlcod = mesh_->cod;
	break;
      }  
    case __ensBASIS_n:
    case __ensBASIS_error:
      {
	break;
      }  
    default:
      {
	space_->own_ddlcod 	= malloc(sizeof(I)*space_->nddl);
	space_->ddlcod 		= space_->own_ddlcod;
	break;
      }
    }

  switch(space_->shape)
    {
    case __ensBASIS_LAGRANGE_1:
    case __ensBASIS_L2ORTHONORMAL_1:
    case __ensBASIS_LAGRANGE_2:
    case __ensBASIS_L2ORTHONORMAL_2:
      {
	break;
      }  
    case __ensBASIS_LAGRANGE_0:
    case __ensBASIS_L2ORTHONORMAL_0:
      {
	{ I i;
	  for (i=0;i<space_->nddl;++i)
	    {
	      space_->own_ddlcod[i] = mesh_->noboundary_vcod;
	    } }	
	break;
      }  
    case __ensBASIS_LAGRANGEBUBBLE_2:
      {
	{ I i;
	  for (i=0;i<mesh_->nvertex+mesh_->nedge;++i)
	    {
	      space_->own_ddlcod[i]	= mesh_->cod[i];
	    } }
	{ I i;
	  for (i=0;i<mesh_->nelm;++i)
	    {
	      space_->own_ddlcod[mesh_->nvertex+mesh_->nedge+i]	= mesh_->noboundary_vcod;
	    } }

	break;
      }
    case __ensBASIS_LAGRANGE_3:
    case __ensBASIS_L2ORTHONORMAL_3:
      {
	{ I i;
	  for (i=0;i<mesh_->nvertex;++i)
	    {
	      space_->own_ddlcod[i] = mesh_->cod[i];
	    } }
	{ I i;
	  for (i=0;i<mesh_->nedge;++i)
	    {
	      { I j;
		for (j=0;j<2;++j)
		  {
		    space_->own_ddlcod[mesh_->nvertex+2*i+j] = mesh_->cod[mesh_->nvertex+i];
		  } }
	    } }
	{ I i;
	  for (i=mesh_->nvertex+2*mesh_->nedge;i<space_->nddl;++i)
	    {
	      space_->own_ddlcod[i] = mesh_->noboundary_vcod;
	    } }	
	break;
      }
    case __ensBASIS_LAGRANGEBUBBLE_0:
    case __ensBASIS_LAGRANGEBUBBLE_1:
    case __ensBASIS_error:
    case __ensBASIS_n:
      {
	break;
      }  
    }


#if 0
  I degree=ens_lagr_tria_degree(shape_);

  const I degree 		= 1 + ens_lagr_tria_degree(shape_);
  ns_mesh_addcubature		(mesh_,(degree<2) ? 2 : degree);
  ns_mesh_addcubature_interval	(mesh_,(degree<2) ? 2 : degree);
  nsSPACE_galerkin_def(space_,
			&space_->endomorphism);      
#endif
#if 0
  { ens_lagr_tria j 		= __ens_lagr_tria_error;
    const I n 		= space_->nddlelm;
    nsCUBA_TRIA_ST * c 	= nsCUBA_TRIA_gauss_new(10);
    for (++j;j<__ens_lagr_tria_n;++j)
      {
	const I nj		= ens_lagr_tria_n(j);
	S_->MASS 	[j] 	= calloc(nj*n*n,sizeof(R));
	S_->CONV_DR 	[j] 	= calloc(nj*n*n,sizeof(R));
	S_->CONV_DS 	[j] 	= calloc(nj*n*n,sizeof(R));
	{ I k;
	  for (k=0;k<nj;++k)
	    {
	      { I i;
		for (i=0;i<c->q_n;++i)
		  {
		    nsblas_dger(jacobian_,
				c->qelm_lagr[j] + c->q_n * k + i,
				c->wqelm_lagr[space_->shape] + i*n,
				&n,
				c->qelm_lagr[space_->shape] + i,
				&c->q_n,
				S_->MASS[j] + n*n*k,
				&n);
		    nsblas_dger(jacobian_,
				c->qelm_lagr[j] + c->q_n * k + i,
				c->wqelm_lagr[space_->shape] + i*n,
				&n,
				c->qelm_lagr_ds[space_->shape] + i,
				&c->q_n,
				S_->CONV_DS[j] + n*n*k,
				&n);
		    nsblas_dger(jacobian_,
				c->qelm_lagr[j] + c->q_n * k + i,
				c->wqelm_lagr[space_->shape] + i*n,
				&n,
				c->qelm_lagr_dr[space_->shape] + i,
				&c->q_n,
				S_->CONV_DR[j] + n*n*k,
				&n);



		  } }
	    } }
	c = nsCUBA_TRIA_kill(c);
      } }
#endif
}

nsSPACE_ST* 	 nsSPACE_new 	(const ns_mesh * 	mesh_,
				 cst_ensBASIS 	shape_)
{
  nsSPACE_ST*space =malloc(sizeof(nsSPACE_ST));
  nsSPACE_def 	(space,  
		 mesh_,
		 shape_);
  return space;
}



void 		nsSPACE_get_ddlelm		(cst_nsSPACE  S_,cst_pI ielm_,cst_pR x_,pR ddlelm_)
{
  { I i;
    for (i=0;i<S_->nddlelm;++i)
      {
	ddlelm_[i] = x_[S_->ddlcnc[ielm_[0]*S_->ddlcncoff+i]-1];
      } }
}


void 		nsSPACE_cncelm_dec	(cst_nsSPACE	space_,
					 cst_pI 		ielm_,
					 pI   			ddlelm_,
					 const I 		dec_)
{
  { I i;    
    for (i=0;i<space_->nddlelm;++i)
      {
	ddlelm_[i] = space_->ddlcnc[ielm_[0]*space_->ddlcncoff+i]-1 + dec_;
      } }
}

void 		nsSPACE_cncelm	(cst_nsSPACE	space_,
				 cst_pI 		ielm_,
				 pI   			ddlelm_)
{
  { I i;    
    for (i=0;i<space_->nddlelm;++i)
      {
	ddlelm_[i] = space_->ddlcnc[ielm_[0]*space_->ddlcncoff+i]-1;
      } }
}
#if 0
void nsSPACE_codelm	(cst_nsSPACE	space_,
			 cst_pI 		ielm_,
			 pI   			codelm_)
{
  const ns_mesh * mesh  = nsSPACE_mesh(space_);

  const I refcod[7] = {mesh->cod[mesh->cnc[ielm_[0]*6+0]-1],
			   mesh->cod[mesh->cnc[ielm_[0]*6+1]-1],
			   mesh->cod[mesh->cnc[ielm_[0]*6+2]-1],
			   mesh->cod[mesh->cnc[ielm_[0]*6+3]-1],
			   mesh->cod[mesh->cnc[ielm_[0]*6+4]-1],
			   mesh->cod[mesh->cnc[ielm_[0]*6+5]-1],
			   mesh->noboundary_vcod};  

  switch(space_->shape)
    {
    case __ens_lagr_tria_error:
    case __ens_lagr_tria_n:
      {
	ns_errmsg("nsSPACE_codelm:switch failed on __ens_lagr_tria");
	break;
      }
    case __ens_lagr_tria_01:
      {
	codelm_[0]  = refcod[6];
	break;
      }  
    case __ens_lagr_tria_03:
      {
	codelm_[0] =  refcod[0];
	codelm_[1] =  refcod[1];
	codelm_[2] =  refcod[2];
	break;
      }  
    case __ens_lagr_tria_06:
      {
	  codelm_[0] = refcod[0];
	  codelm_[1] = refcod[1];
	  codelm_[2] = refcod[2];
	  codelm_[3] = refcod[3];
	  codelm_[4] = refcod[4];
	  codelm_[5] = refcod[5];
	break;
      }  
    case __ens_lagr_tria_07:
      {
	codelm_[0] = refcod[0];
	codelm_[1] = refcod[1];
	codelm_[2] = refcod[2];
	codelm_[3] = refcod[3];
	codelm_[4] = refcod[4];
	codelm_[5] = refcod[5];
	codelm_[6] = refcod[6];
	break;
      }
    case __ens_lagr_tria_10:
      { 
	codelm_[0] = refcod[0];
	codelm_[1] = refcod[1];
	codelm_[2] = refcod[2];

	codelm_[3] = refcod[3];
	codelm_[4] = refcod[3];

	codelm_[5] = refcod[4];
	codelm_[6] = refcod[4];

	codelm_[7] = refcod[5];
	codelm_[8] = refcod[5];

	codelm_[9] = refcod[6];
	break;
      }
    case __ens_lagr_tria_15:
      {
	codelm_[0] = refcod[0];
	codelm_[1] = refcod[1];
	codelm_[2] = refcod[2];


	codelm_[3] = refcod[3];
	codelm_[4] = refcod[3];
	codelm_[5] = refcod[3];

	codelm_[6] = refcod[4];
	codelm_[7] = refcod[4];
	codelm_[8] = refcod[4];

	codelm_[9] = refcod[5];
	codelm_[10] = refcod[5];
	codelm_[11] = refcod[5];

	codelm_[12] = refcod[6];
	codelm_[13] = refcod[6];
	codelm_[14] = refcod[6];
	break;
      }
    case __ens_lagr_tria_21:
      {
	codelm_[0] = refcod[0];
	codelm_[1] = refcod[1];
	codelm_[2] = refcod[2];


	codelm_[3] = refcod[3];
	codelm_[4] = refcod[3];
	codelm_[5] = refcod[3];
	codelm_[6] = refcod[3];

	codelm_[7] = refcod[4];
	codelm_[8] = refcod[4];
	codelm_[9] = refcod[4];
	codelm_[10] = refcod[4];

	codelm_[11] = refcod[5];
	codelm_[12] = refcod[5];
	codelm_[13] = refcod[5];
	codelm_[14] = refcod[5];


	codelm_[15] = refcod[6];
	codelm_[16] = refcod[6];
	codelm_[17] = refcod[6];
	codelm_[18] = refcod[6];
	codelm_[19] = refcod[6];
	codelm_[20] = refcod[6];

	break;
      }
    case __ens_lagr_tria_28:
      {
	codelm_[0] = refcod[0];
	codelm_[1] = refcod[1];
	codelm_[2] = refcod[2];

	codelm_[3] = refcod[3];
	codelm_[4] = refcod[3];
	codelm_[5] = refcod[3];
	codelm_[6] = refcod[3];
	codelm_[7] = refcod[3];

	codelm_[8] = refcod[4];
	codelm_[9] = refcod[4];
	codelm_[10] = refcod[4];
	codelm_[11] = refcod[4];
	codelm_[12] = refcod[4];

	codelm_[13] = refcod[5];
	codelm_[14] = refcod[5];
	codelm_[15] = refcod[5];
	codelm_[16] = refcod[5];
	codelm_[17] = refcod[5];

	codelm_[18] = refcod[6];
	codelm_[19] = refcod[6];
	codelm_[20] = refcod[6];
	codelm_[21] = refcod[6];
	codelm_[22] = refcod[6];
	codelm_[23] = refcod[6];
	codelm_[24] = refcod[6];
	codelm_[25] = refcod[6];
	codelm_[26] = refcod[6];
	codelm_[27] = refcod[6];
	break;
      } 
    }
}
#endif

void nsSPACE_elm2ddl(cst_nsSPACE 	space_,
		      cst_pI		select_n_,
		      cst_pI		select_,
		      pI       	ddl_,
		      pI   		ddl_n_,
		      pI		blank_)
{
  cst_pI   nddlelm_ 	= &space_->nddlelm;
  cst_pI   ddlcnc_ 	= space_->ddlcnc;
  cst_pI   ddlcncoff_ 	= &space_->nddlelm;
		      
  ddl_n_[0]=(I)0;
  { I i;
    for (i=0;i<select_n_[0];++i)
      {
	const I ielm = select_[i];
	{ I j;
	  for (j=0;j<nddlelm_[0];++j)
	    {
	      const I iddl 		= ddlcnc_[ddlcncoff_[0]*ielm+j]-1;
	      if (NOT blank_[iddl])
		{
		  blank_[iddl] 		= ++ddl_n_[0]; 
		  ddl_[ddl_n_[0]-1] 	= iddl;
		}
	    } }
      } }
  { I j;
    for (j=0;j<ddl_n_[0];++j)
      {
	blank_[ddl_[j]]=(I)0;	    
      } }
}


void nsSPACE_elm2ddl_hat(cst_nsSPACE 	space_,
			  cst_pI		ivertex_,
			  cst_pI		select_n_,
			  cst_pI		select_,
			  pI       		ddl_,
			  pI   		ddl_n_,
			  pR 			hat_,
			  pI			blank_)
{
  cst_pI   nddlelm_ 	= &space_->nddlelm;
  cst_pI   ddlcnc_ 	= space_->ddlcnc;
  cst_pI   ddlcncoff_ 	= &space_->nddlelm;
  cst_pR   lc2gl	= ensBASIS_lc2gltreilli_ptr(space_->shape,nsSPACE_elm(space_));
  ddl_n_[0]=(I)0;
  { I i;
    for (i=0;i<select_n_[0];++i)
      {
	const I ielm = select_[i];
	I barycentric=0;
	if (ivertex_[0]==ddlcnc_[ddlcncoff_[0]*ielm+1]-1)
	  barycentric=1;
	if (ivertex_[0]==ddlcnc_[ddlcncoff_[0]*ielm+2]-1)
	  barycentric=2;
	{ I j;
	  for (j=0;j<nddlelm_[0];++j)
	    {
	      const I iddl 		= ddlcnc_[ddlcncoff_[0]*ielm+j]-1;
	      if (NOT blank_[iddl])
		{
		  blank_[iddl] 		= ++ddl_n_[0]; 
		  ddl_[ddl_n_[0]-1] 	= iddl;
		  hat_[ddl_n_[0]-1]	= lc2gl[nddlelm_[0]*barycentric+j];
		}
	    } }
      } }
  { I j;
    for (j=0;j<ddl_n_[0];++j)
      {
	blank_[ddl_[j]]=(I)0;	    
      } }
}



void nsSPACE_patch_vertices2elm(cst_nsSPACE 	space_,
				 cst_pI 		level_,
				 cst_pI 		nvertices_,
				 cst_pI 		vertices_,
				 pI 			select_n_,
				 pI 			select_,
				 pI 			blank_)
{
  const ns_mesh * mesh = nsSPACE_mesh(space_);
  select_n_[0]=(I)0;
  { I ivertex;
    for (ivertex=0;ivertex<nvertices_[0];++ivertex)
      {
	const I N1=mesh->bpatch[vertices_[ivertex]+1]-mesh->bpatch[vertices_[ivertex]];
	{ I i;
	  for (i=0;i<N1;++i)
	    {
	      const I q = mesh->patch[mesh->bpatch[vertices_[ivertex]]+i];
	      if (NOT blank_[q])
		{
		  select_n_[0]			= select_n_[0]+1;
		  blank_[q]  			= select_n_[0];
		  select_[select_n_[0]-1] 	= q;
		}
	    } }
      } }
  { I j;
    for (j=0;j<select_n_[0];++j)
      {
	blank_[select_[j]]=(I)0;	    
      } }
}


void	nsSPACE_patch_level_vertices(cst_nsSPACE	space_,
				      cst_pI 		level_,					      
				      pI		select_n_,
				      pI		select_,
				      pI		blank_)
{
  const ns_mesh * mesh = nsSPACE_mesh(space_);

  cst_pI   	cnc_ 		= space_->ddlcnc;
  cst_pI   	cncoff_ 	= &space_->nddlelm;		      
  cst_pI	bpatch_		= mesh->bpatch;
  cst_pI 	patch_		= mesh->patch;
  { I lev;
    for (lev=1;lev<level_[0];++lev)
      {      
	const I N = select_n_[0];
	{ I s;
	  for (s=0;s<N;++s)
	    {
	      const I iv = select_[s];
	      const I N1 = bpatch_[iv+1]-bpatch_[iv];
	      { I i;
		for (i=0;i<N1;++i)
		  {
		    const I k  = patch_[bpatch_[iv]+i];
		    { I j;
		      for (j=0;j<3;++j)
			{
			  const I  q = cnc_[cncoff_[0]*k+j]-1;
			  if ( NOT blank_[q] )
			    {			      
			      select_[select_n_[0]++]=q;
			      blank_[q]=select_n_[0];
			    }
			} }
		  } }
	    } }
      } }
  { I j;
    for (j=0;j<select_n_[0];++j)
      {
	blank_[select_[j]]=(I)0;	    
      } }
}



void nsSPACE_infopatch(cst_nsSPACE	space_,
			cst_pI 	ivertex_,
			cst_pI 	level_,
			pI 		patchvertex_n_,
			pI 		patchvertex_,
			pI 		patchelm_n_,
			pI 		patchelm_,
			pI 		patchddl_n_,
			pI 		patchddl_,
			pI 		blank_)
{
  patchvertex_n_[0] 	= (I)1;
  blank_[ivertex_[0]]	= (I)1;
  patchvertex_[0]	= ivertex_[0];  
  nsSPACE_patch_level_vertices	(space_,
				 level_,
				 patchvertex_n_,
				 patchvertex_,
				 blank_);  
  patchelm_n_[0] = (I)0;
  nsSPACE_patch_vertices2elm	(space_,
				 level_,
				 patchvertex_n_,
				 patchvertex_,
				 patchelm_n_,
				 patchelm_,
				 blank_);
  patchddl_n_[0] = 0;
  nsSPACE_elm2ddl		(space_,
				 patchelm_n_,
				 patchelm_,
				 patchddl_,
				 patchddl_n_,
				 blank_);
}


void nsSPACE_infopatch_hat(cst_nsSPACE	space_,
			    cst_pI 	ivertex_,
			    pI 	patchvertex_n_,
			    pI 	patchvertex_,
			    pI 	patchelm_n_,
			    pI 	patchelm_,
			    pI 	patchddl_n_,
			    pI 	patchddl_,
			    pR 	hatddl_,
			    pI 	blank_)
{
  I level_[1] = {1};
  patchvertex_n_[0] 	= (I)1;
  blank_[ivertex_[0]]	= (I)1;
  patchvertex_[0]	= ivertex_[0];  
  nsSPACE_patch_level_vertices	(space_,
				 level_,
				 patchvertex_n_,
				 patchvertex_,
				 blank_);  
  patchelm_n_[0] = (I)0;
  nsSPACE_patch_vertices2elm	(space_,
				 level_,
				 patchvertex_n_,
				 patchvertex_,
				 patchelm_n_,
				 patchelm_,
				 blank_);
  patchddl_n_[0] = 0;
  nsSPACE_elm2ddl_hat		(space_,
				 ivertex_,
				 patchelm_n_,
				 patchelm_,
				 patchddl_,
				 patchddl_n_,
				 hatddl_,
				 blank_);
}




#if 0
nsSPARSE_ST*nsSPACE_smooth_matrix(cst_nsSPACE space_)
{
  cst_emk_smoothing_method 	method	= __emk_smoothing_method_L2_q;
  cst_emk_smoothing_op		op	= __emk_smoothing_op_nabla;
  I mat_level		= 1;
  I weighted		= (I)NO;
  I zz			= (I)NO;
  I apply_metric 		= (I)NO;
  mkS_st 	shape_V;
  mkS_definit(&shape_V,mkZ_elm(zone),__emkS_FAMILY_canonic,kf+1,__emk_discontinuous,&err);  
  mkSP dx_matrix = mok_patch_matrix_new_(method,
					 op,
					 &mat_level,
					 "Quadrature",
					 &zz,
					 &weighted,
					 &apply_metric,			  
					 zone,
					 space_f,
					 space_f,
					 &shape_V,  
					 NULL,
					 &err);

#if 0
  space_->nddlelm	= ens_lagr_tria_n(shape_);
  space_->nddl		= ns_mesh_nddlspace(mesh_,shape_);
  switch(lagr_)
    {
    case __ens_lagr_tria_01:
      {
	break;
      }  
    case __ens_lagr_tria_03:
      {
	space_->cnc 	= mesh_->cnc;
	space_->coo 	= mesh_->coo;
	space_->vcod 	= mesh_->cod;
	space_->bpatch 	= mesh_->bpatch;
	space_->patch 	= mesh_->patch;
	break;
      }  
    case __ens_lagr_tria_06:
      {
	space_->cnc 	= mesh_->cnc;
	space_->coo 	= mesh_->coo;
	space_->vcod 	= mesh_->cod;
	space_->bpatch 	= mesh_->bpatch;
	space_->patch 	= mesh_->patch;
	break;
      }  
    case __ens_lagr_tria_07:
    case __ens_lagr_tria_10:
    case __ens_lagr_tria_15:
    case __ens_lagr_tria_21:
    case __ens_lagr_tria_28:
      {
	space_->own_cnc 	= malloc(sizeof(I)*space_->nddlelm);
	space_->own_coo 	= malloc(sizeof(R)*space_->nddl*2);
	space_->own_vcod 	= malloc(sizeof(I)*space_->nddl);
	space_->own_bpatch 	= malloc(sizeof(I)*(space_->nddl+1));
	space_->own_patch 	= malloc(sizeof(I)*space_->nddlelm);
	space_->cnc 		= space_->own_cnc;
	space_->coo 		= space_->own_coo;
	space_->vcod 		= space_->own_vcod;
	{ I i;
	  for (i=0;i<mesh_->nvertex;++i)
	    {
	      space_->own_coo[2*i+0]	= mesh_->coo[2*i+0];
	      space_->own_coo[2*i+1]	= mesh_->coo[2*i+1];
	    } }

	break;
      }  
    }
#endif
#endif






#if 0
void nsSPACE_residu_initial(pR 	global_rhs_,
			     cst_pR 	x_,
			     cst_pR 	xi_,
			     cst_pR 	xii_,
			     ns_var * 	var_tension_)
{
#ifndef NDEBUG
  wrong_param(NOT global_rhs_,1,ns_build_residu_gauss);
  wrong_param(NOT x_,2,ns_build_residu_gauss);
  wrong_param(NOT xi_,3,ns_build_residu_gauss);
  wrong_param(NOT xii_,4,ns_build_residu_gauss);
  if (ns_linfo[__ens_linfo_tension])
    {
      wrong_param(NOT var_tension_,5,ns_build_residu_gauss);      
    }
#endif
  /* residu de pression (completement lineaire) */	  
  nsSPARSE_BLOCK_gemv(&regal1,
		       &sparse_A_block_BU,
		       x_,
		       &regal1,
		       global_rhs_);
  
  /*
    on calcule le residu element par element
  */
  { I ielm;
    for (ielm=0;ielm<mesh->nelm;++ielm)
      {       
	/* INFO ELEMENT */
	ns_setelm		(x_,xi_,xii_);            
	/* CLEAR LOCAL SYSTEM */
	clr_dvect		(total_nddlelmxtotal_nddlelm,locmat);
	clr_dvect		(_total_nddlelm,locresidu);                  
	/* COMPUTE GAUSS  */
	ns_gauss_all();                  
	/* COMPUTE RESIDU F : masse */
	gauss_residu_f();
	gauss_residu_u();
	gauss_residu_v();
	if ( (ns_linfo[__ens_linfo_transport]) AND (NOT ns_linfo[__ens_linfo_transport_galerkin]) )
	  {
	    printf("a corriger\n");exit(1);
	    integrator_freesurface_dg_system();
	  }      
	ns_get_ddlcnc(ielm,locnumer);
	{ I i;
	  for (i=0;i<_total_nddlelm;++i)
	    {
	      global_rhs_[locnumer[i]] += locresidu[i];
	    } }    
      } }

  if (ns_linfo[__ens_linfo_slip])
    {
      nsSPARSE_BLOCK_gemv(&regal1,&sparse_A_block_SS,x_,&regal1,global_rhs_);      
      nsSPARSE_BLOCK_gemv(&regal1,&sparse_A_block_SU,x_,&regal1,global_rhs_);
      nsSPARSE_BLOCK_gemv(&regal1,&sparse_A_block_US,x_,&regal1,global_rhs_);
    }   
  ns_build_dirichlet	(global_rhs_,nddlf_dirichlet,ddlf_dirichlet,((I)0));
  ns_build_dirichlet	(global_rhs_,nddlu_dirichlet,ddlu_dirichlet,dec_ddlu);
  ns_build_dirichlet	(global_rhs_,nddlv_dirichlet,ddlv_dirichlet,dec_ddlu+nddlu);
  ns_build_dirichlet	(global_rhs_,nddlp_dirichlet,ddlp_dirichlet,dec_ddlp);
}
#endif

#if 0
static cst_pR nsSPACE_K_M	(cst_nsSPACE 	S_,
				 const ens_lagr_tria 	coeff_shape_)
{
  return S_->MASS[coeff_shape_];
}

static cst_pR nsSPACE_K_CDR	(cst_nsSPACE 	S_,
				 const ens_lagr_tria 	coeff_shape_)
{
  return S_->CONV_DR[coeff_shape_];
}

static cst_pR nsSPACE_K_CDS	(cst_nsSPACE 	S_,
				 const ens_lagr_tria 	coeff_shape_)
{
  return S_->CONV_DS[coeff_shape_];
}

static void nsSPACE_K_CDR_matrix(cst_nsSPACE 	space_,
				  const ens_lagr_tria 	coeff_shape_,
				  cst_pR 		xcoeff_,
				  cst_pR		coeff_,
				  cst_pI		coeffoff_,
				  cst_pR		xmatrix_,
				  pR 			matrix_,
				  cst_pI 		matrixoff_)
{
  const I n = ens_lagr_tria_n(coeff_shape_);
  const I N = space_->nddlelm*space_->nddlelm;
  nsblas_dgemm(transN,
	       transN,
	       &N,
	       &negal1,
	       &n,
	       xcoeff_,
	       nsSPACE_K_CDR(space_,coeff_shape_),
	       &N,
	       coeff_,
	       coeffoff_,
	       xmatrix_,
	       matrix_,
	       matrixoff_);
}

static void nsSPACE_K_CDS_matrix(cst_nsSPACE 		space_,
			   const ens_lagr_tria 		coeff_shape_,
			   cst_pR 			xcoeff_,
			   cst_pR			coeff_,
			   cst_pI			coeffoff_,
			   cst_pR			xmatrix_,
			   pR 				matrix_,
			   cst_pI 			matrixoff_)
{
  const I n = ens_lagr_tria_n(coeff_shape_);
  const I N = space_->nddlelm*space_->nddlelm;
  nsblas_dgemm(transN,
	       transN,
	       &N,
	       &negal1,
	       &n,
	       xcoeff_,
	       nsSPACE_K_CDS(space_,coeff_shape_),
	       &N,
	       coeff_,
	       coeffoff_,
	       xmatrix_,
	       matrix_,
	       matrixoff_);
}



void nsSPACE_K_M_matrix(cst_nsSPACE 		space_,
			 const ens_lagr_tria 		coeff_shape_,
			 cst_pR 			xcoeff_,
			 cst_pR			coeff_,
			 cst_pI			coeffoff_,
			 cst_pR			xmatrix_,
			 pR 				matrix_,
			 cst_pI 			matrixoff_)
{
  const I n = ens_lagr_tria_n(coeff_shape_);
  const I N = space_->nddlelm*space_->nddlelm;
  nsblas_dgemm(transN,
	       transN,
	       &N,
	       &negal1,
	       &n,
	       xcoeff_,
	       nsSPACE_K_M(space_,coeff_shape_),
	       &N,
	       coeff_,
	       coeffoff_,
	       xmatrix_,
	       matrix_,
	       matrixoff_);
}

void nsSPACE_K_MASS_matrix(cst_nsSPACE 		space_,
			    const ens_lagr_tria         shape_density_,
			    cst_pR			jacelm_,
			    cst_pR 			xcoeff_,
			    cst_pR			coeff_,
			    cst_pI			coeffoff_,
			    pR 			matrix_,
			    cst_pI 			matrixoff_)
{
  const R a = jacelm_[0] * xcoeff_[0];
  nsSPACE_K_M_matrix		(space_,
				 shape_density_,
				 &a,
				 coeff_,
				 coeffoff_,
				 &regal1,
				 matrix_,
				 matrixoff_);
}

void nsSPACE_K_CONVECTION_matrix(cst_nsSPACE 		space_,
				  const ens_lagr_tria 		shape_velocity_,
				  cst_pR			belm_,
				  cst_pR 			xcoeff_,
				  cst_pR			coeff_,
				  cst_pI			coeffoff_,
				  pR 				matrix_,
				  cst_pI 			matrixoff_)
{
  const R belm [4] = {xcoeff_[0]*belm_[0],xcoeff_[0]*belm_[1],xcoeff_[0]*belm_[2],xcoeff_[0]*belm_[3]};
  nsSPACE_K_CDR_matrix		(space_,
				 shape_velocity_,
				 belm,
				 coeff_,
				 coeffoff_,
				 &regal1,
				 matrix_,
				 matrixoff_);
  nsSPACE_K_CDS_matrix		(space_,
				 shape_velocity_,
				 &belm[2],
				 coeff_,
				 coeffoff_,
				 &regal1,
				 matrix_,
				 matrixoff_);
  nsSPACE_K_CDR_matrix		(space_,
				 shape_velocity_,
				 &belm[1],
				 &coeff_[coeffoff_[0]],
				 coeffoff_,
				 &regal1,
				 matrix_,
				 matrixoff_);
  nsSPACE_K_CDS_matrix		(space_,
				 shape_velocity_,
				 &belm[3],
				 &coeff_[coeffoff_[0]],
				 coeffoff_,
				 &regal1,
				 matrix_,
				 matrixoff_);
}
#endif
