#include "nsSPACE_PRODUCT.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
I 		nsSPACE_PRODUCT_get_nddl	(cst_nsSPACE_PRODUCT const self_){ return self_->nddl; }
I 		nsSPACE_PRODUCT_get_nddlelm	(cst_nsSPACE_PRODUCT const self_){ return self_->nddlelm; }
I 		nsSPACE_PRODUCT_get_nspace	(cst_nsSPACE_PRODUCT const self_){ return self_->nspace; }
pSpaceReadOnly 	nsSPACE_PRODUCT_get_space	(cst_nsSPACE_PRODUCT const self_,const I i_){ return self_->spaces[i_]; }
const ns_mesh* 	nsSPACE_PRODUCT_get_mesh	(cst_nsSPACE_PRODUCT const self_){ return nsSPACE_mesh(self_->spaces[0]); }
void 		nsSPACE_PRODUCT_free		(nsSPACE_PRODUCT 	self_) 
{ 
  if (self_)
    {
      if (self_->dep)
	free(self_->dep);
      memset(self_,0,sizeof(nsSPACE_PRODUCT_ST));
    }
}

void 		nsSPACE_PRODUCT_def	(nsSPACE_PRODUCT 	self_,
					 cst_pI		dep_,
					 const I 		nspace_,
					 ...)
{
  memset(self_,0,sizeof(nsSPACE_PRODUCT_ST));
  va_list args;
  va_start (args,nspace_);
  { I i;
    for (i=0;i<nspace_;++i)
      {
	self_->spaces[i] = va_arg(args,pSpace);
	self_->nddl 	+= nsSPACE_nddl		(self_->spaces[i]);
	self_->nddlelm 	+= nsSPACE_nddlelm	(self_->spaces[i]);
      } }
  va_end(args);  
  self_->nspace = nspace_;
  self_->dep 	= (pI)calloc(nspace_*nspace_,sizeof(I));
  memcpy(self_->dep,dep_,nspace_*nspace_);
}


static I nsSPACE_PRODUCT_galerkin_line( cst_nsSPACE_PRODUCT	self_,
					    const I 	ispace_,
					    const I 	iddl_,
					    pI 		ddl_,
					    pI		blank_,
					    pI 		iwork_)
{
  pSpaceReadOnly spacei 		= nsSPACE_PRODUCT_get_space(self_,ispace_);
  const ns_mesh * mesh 		= nsSPACE_PRODUCT_get_mesh(self_);
  I N=0;
  const I numGridNodes = ns_mesh_get_numNodes(mesh);
  const I numGridEdges = ns_mesh_get_numEdges(mesh);
  const I jddl 		= nsSPACE_iddl2ddlP6(spacei,iddl_);
  { I i;
    const I n 	= (jddl>=numGridNodes+numGridEdges) ? 1  : mesh->bpatch[jddl+1]-mesh->bpatch[jddl];
    for (i=0;i<n;++i)
      {	
	const I k = (jddl>=numGridNodes+numGridEdges) ? jddl-(numGridNodes+numGridEdges) : mesh->patch[mesh->bpatch[jddl]+i];	
	{ I jspace=0;
	  for (jspace=0;jspace<self_->nspace;++jspace)
	    {
	      pSpaceReadOnly spacej 	= nsSPACE_PRODUCT_get_space(self_,jspace);
	      const I nddlelmj 	= nsSPACE_nddlelm(spacei);
	      nsSPACE_cncelm	(spacej,
				 &k,
				 iwork_);	      
	      /**/
	      { I j;
		for (j=0;j<nddlelmj;++j)
		  {
		    if (NOT blank_[iwork_[j]])
		      {
			blank_[iwork_[j]] = ++N;
			ddl_[N-1]   = iwork_[j];
		      }
		  } }
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


pSparse nsSPACE_PRODUCT_galerkin_def(cst_nsSPACE_PRODUCT	self_)
{
  const I nspace 	= nsSPACE_PRODUCT_get_nspace(self_);
  const I nddl 		= nsSPACE_PRODUCT_get_nddl(self_);
  const I nddlelm 	= nsSPACE_PRODUCT_get_nddlelm(self_);
  
#if 0
  Sparse_clr(F_);
#endif
#if 0
  memset(F_,0,sizeof(nsSPARSE_ST));
#endif

  const I n			= nddl;
  const I m    		= nddl;
  pI own_b		= (pI)malloc((n+1)*sizeof(I));

  pI own_i;
  { const ns_mesh * mesh= nsSPACE_PRODUCT_get_mesh(self_);
    const I nelm 	= ns_mesh_nelm(mesh);
    own_i		= (pI)malloc(nddlelm*nddlelm*nelm*sizeof(I)); }

  own_b[0] 		= 0;
  I nc 		= 0;
  I	nline;

  pI blank 		= (pI)calloc(nddl,sizeof(I));
  pI iwork 		= (pI)calloc(nddl,sizeof(I));

  { I ispace;
    for (ispace=0;ispace<nspace;++ispace)
      {


	{ pSpaceReadOnly spacei	= nsSPACE_PRODUCT_get_space(self_,ispace);
	  const I nddli 	= nsSPACE_nddl(spacei);
	  I iddl;
	  for (iddl=0;iddl<nddli;++iddl)
	    {
	      nline = nsSPACE_PRODUCT_galerkin_line	(self_,
							 ispace,
							 iddl,
							 &own_i[nc],
							 blank,
							 iwork);
	      nc 			+=nline;
	      own_b[iddl+1] 	= own_b[iddl] + nline;
	    } }


      } }


  free(iwork);
  free(blank);
  own_i	= (pI)realloc(own_i,nc*sizeof(I));  

  pSparse back = Sparse_build(n,m,nc,own_b,own_i,NULL);
  Sparse_sort(back);
  Sparse_fortran_indexation(back);
  return back;
}
