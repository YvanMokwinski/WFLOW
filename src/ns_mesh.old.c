#include "ns_mesh.h"
#include "Medit.h"
#include "ns_sys.h"
#include "eFileSuffix.h"

#include <math.h>
#include <stdarg.h>
#ifndef NDEBUG
#define wrong_param_ns_mesh(_cond,_i,_s) if ( (_cond) ) { fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); exit(1);} ((void)0)
#else
#define wrong_param_ns_mesh(_cond,_i,_s) ((void)0)
#endif

#ifndef NDEBUG
eDim		ns_mesh_dimension	(const ns_mesh*		mesh_){return mesh_->dimension;}
L 		ns_mesh_istransient	(const 	ns_mesh*mesh_) 	
{ 
  DebugVerif(mesh_); 
  return ( (mesh_->time_x) ? __emnsYES : __emnsNO ); 
}

eElement	ns_mesh_elm		(const ns_mesh*		mesh_){return mesh_->elm;}

I 		ns_mesh_nelm		(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->nelm; }
I 		ns_mesh_nvertex		(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->nvertex; }
I 		ns_mesh_nedge_boundary	(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->nedge_boundary; }
I 		ns_mesh_nedge		(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->nedge; }
cst_pR 	ns_mesh_trelm		(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->trelm; }
cst_pR 	ns_mesh_jacelm		(const 	ns_mesh*mesh_)	{ DebugVerif(mesh_); return mesh_->jacelm; }
I 		ns_mesh_nddlspace	(const 	ns_mesh * 	mesh_,
					 const ensBASIS 	shape_)
{
  DebugVerif(mesh_);
  return mesh_->nddl_lagr[shape_];
}

#endif
  void ns_mesh_get_trelm(const ns_mesh*mesh_,cst_pI jelm_,pR belm_)
  {
    belm_[3]   	= mesh_->trelm[4*jelm_[0]+3];
    belm_[1]   	= mesh_->trelm[4*jelm_[0]+1];
    belm_[2]   	= mesh_->trelm[4*jelm_[0]+2];
    belm_[0]   	= mesh_->trelm[4*jelm_[0]+0];
  }


void ns_mesh_vertex_patchelm(const ns_mesh * 	mesh_,
			     cst_pI 		ivertex_,
			     pI 		patchelm_n_,
			     pI 		patchelm_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT ivertex_,2,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT patchelm_n_,3,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT patchelm_,4,ns_mesh_patchvertex);
  wrong_param_ns_mesh(ivertex_[0]+1<1,2,ns_mesh_patchvertex);
  wrong_param_ns_mesh(ivertex_[0]>=mesh_->nvertex,2,ns_mesh_patchvertex);
#endif
  { I i,j=0;
    for (i=mesh_->bpatch[ivertex_[0]];i<mesh_->bpatch[ivertex_[0]];++i)
      {
	patchelm_[j++] = mesh_->patch[i];
      } 
    patchelm_n_[0] = j; }
}


#if 0
void ns_mesh_patchvertex(const ns_mesh * 	mesh_,
			 cst_pI 		ivertex_,
			 pI 			patchelm_n_,
			 pI 			patchelm_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT ivertex_,2,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT patchelm_n_,3,ns_mesh_patchvertex);
  wrong_param_ns_mesh(NOT patchelm_,4,ns_mesh_patchvertex);
  wrong_param_ns_mesh(ivertex_[0]+1<1,2,ns_mesh_patchvertex);
  wrong_param_ns_mesh(ivertex_[0]>=mesh_->nvertex,2,ns_mesh_patchvertex);
#endif
  { I i,j=0;
    for (i=mesh_->bpatch[ivertex_[0]];i<mesh_->bpatch[ivertex_[0]];++i)
      {
	patchelm_[j++] = mesh_->patch[i]-1;
      } 
    patchelm_n_[0] = j; }
}
#endif

void  ns_mesh_get_coovertex(const ns_mesh*mesh_,cst_pI ivertex_,pR xc_)
{
  const I dimension = ns_mesh_dimension(mesh_);
  I j;
  for (j=0;j<dimension;++j)
      xc_[j] = mesh_->coo[dimension*ivertex_[0]+j];
}

void ns_mesh_cncelm(const ns_mesh*	mesh_,
		    cst_pI 		ielm_,
		    pI         	cncelm_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_cncelm); 
  wrong_param_ns_mesh(NOT ielm_,2,ns_mesh_cncelm); 
  wrong_param_ns_mesh(NOT cncelm_,3,ns_mesh_cncelm); 
  wrong_param_ns_mesh(ielm_[0]+1<1,2,ns_mesh_cncelm); 
  wrong_param_ns_mesh(ielm_[0]+1>mesh_->nelm,2,ns_mesh_cncelm); 
#endif
  cncelm_[0] = mesh_->cnc[ielm_[0]*6+0]-1;
  cncelm_[1] = mesh_->cnc[ielm_[0]*6+1]-1;
  cncelm_[2] = mesh_->cnc[ielm_[0]*6+2]-1;
  cncelm_[3] = mesh_->cnc[ielm_[0]*6+3]-1;
  cncelm_[4] = mesh_->cnc[ielm_[0]*6+4]-1;
  cncelm_[5] = mesh_->cnc[ielm_[0]*6+5]-1;
}


void ns_mesh_cncelm0(const ns_mesh*	mesh_,
		     cst_pI 		ielm_,
		     pI         	cncelm_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_cncelm); 
  wrong_param_ns_mesh(NOT ielm_,2,ns_mesh_cncelm); 
  wrong_param_ns_mesh(NOT cncelm_,3,ns_mesh_cncelm); 
  wrong_param_ns_mesh(ielm_[0]+1<1,2,ns_mesh_cncelm); 
  wrong_param_ns_mesh(ielm_[0]+1>mesh_->nelm,2,ns_mesh_cncelm); 
#endif
  cncelm_[0] = mesh_->cnc[ielm_[0]*6+0]-1;
  cncelm_[1] = mesh_->cnc[ielm_[0]*6+1]-1;
  cncelm_[2] = mesh_->cnc[ielm_[0]*6+2]-1;
}

void ns_mesh_patch_cooelm(const ns_mesh*	mesh_,
			  cst_pI 		nelm_,
			  cst_pI            	elm_,
			  pR         		cooelms_)
{
  I cncelm[3];
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_cooelm); 
  wrong_param_ns_mesh(NOT nelm_,2,ns_mesh_cooelm); 
  wrong_param_ns_mesh(NOT elm_,3,ns_mesh_cooelm); 
  wrong_param_ns_mesh(NOT cooelms_,4,ns_mesh_cooelm); 
#endif
  { I i;
    for (i=0;i<nelm_[0];++i)
      {
	ns_mesh_cncelm0(mesh_,elm_+i,cncelm);
	cooelms_[i*3+0] 		= mesh_->coo[cncelm[0]*2+0];
	cooelms_[i*3+1] 		= mesh_->coo[cncelm[1]*2+0];
	cooelms_[i*3+2] 		= mesh_->coo[cncelm[2]*2+0];       
	cooelms_[(i+nelm_[0])*3+0] 	= mesh_->coo[cncelm[0]*2+1];
	cooelms_[(i+nelm_[0])*3+1] 	= mesh_->coo[cncelm[1]*2+1];
	cooelms_[(i+nelm_[0])*3+2] 	= mesh_->coo[cncelm[2]*2+1];
      } }
}


void ns_mesh_cooelm(const ns_mesh*	mesh_,
		    cst_pI 		ielm_,
		    pR         	cooelm_)
{
  I cncelm[3];
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_cooelm); 
  wrong_param_ns_mesh(NOT ielm_,2,ns_mesh_cooelm); 
  wrong_param_ns_mesh(NOT cooelm_,3,ns_mesh_cooelm); 
  wrong_param_ns_mesh(ielm_[0]+1<1,2,ns_mesh_cooelm); 
  wrong_param_ns_mesh(ielm_[0]+1>mesh_->nelm,2,ns_mesh_cooelm); 
#endif

  ns_mesh_cncelm0(mesh_,ielm_,cncelm);

  cooelm_[0] = mesh_->coo[cncelm[0]*2+0];
  cooelm_[1] = mesh_->coo[cncelm[1]*2+0];
  cooelm_[2] = mesh_->coo[cncelm[2]*2+0];

  cooelm_[3] = mesh_->coo[cncelm[0]*2+1];
  cooelm_[4] = mesh_->coo[cncelm[1]*2+1];
  cooelm_[5] = mesh_->coo[cncelm[2]*2+1];
}


void ns_mesh_cooelms(const ns_mesh*	mesh_,
		     cst_pI 		start_ielm_,
		     cst_pI 		stop_ielm_,
		     pR         	cooelm_,
		     cst_pI            cooelmoff_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_cooelms); 
  wrong_param_ns_mesh(NOT start_ielm_,2,ns_mesh_cooelms); 
  wrong_param_ns_mesh(NOT stop_ielm_,3,ns_mesh_cooelms); 
  wrong_param_ns_mesh(start_ielm_[0]+1<1,2,ns_mesh_cooelms); 
  wrong_param_ns_mesh(start_ielm_[0]+1>mesh_->nelm,2,ns_mesh_cooelms); 
  wrong_param_ns_mesh(stop_ielm_[0]+1<1,3,ns_mesh_cooelm); 
  wrong_param_ns_mesh(stop_ielm_[0]+1>mesh_->nelm,3,ns_mesh_cooelms); 
  wrong_param_ns_mesh(NOT cooelm_,4,ns_mesh_cooelms); 
  wrong_param_ns_mesh(NOT cooelmoff_,5,ns_mesh_cooelms); 
  wrong_param_ns_mesh(cooelmoff_[0]<4,5,ns_mesh_cooelms); 
#endif
  const I n = stop_ielm_[0]-start_ielm_[0];
  cst_pR coo 	= mesh_->coo;
  cst_pI cnc 	= mesh_->cnc;
  { I i;
    for (i=0;i<n;++i)
      {
	const I j 	= start_ielm_[0]+i;
	cst_pI cncelm 	= &cnc[j*6+0];
	cooelm_[cooelmoff_[0]*i+0] = coo[cncelm[0]*2+0];
	cooelm_[cooelmoff_[0]*i+1] = coo[cncelm[1]*2+0];
	cooelm_[cooelmoff_[0]*i+2] = coo[cncelm[2]*2+0];	
	cooelm_[cooelmoff_[0]*i+3] = coo[cncelm[0]*2+1];
	cooelm_[cooelmoff_[0]*i+4] = coo[cncelm[1]*2+1];
	cooelm_[cooelmoff_[0]*i+5] = coo[cncelm[2]*2+1];	
      } }
}


void ns_mesh_trelms(const ns_mesh*	mesh_,
		    cst_pI 		start_ielm_,
		    cst_pI 		stop_ielm_,
		    pR         	trelm_,
		    cst_pI            	trelmoff_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_trelms); 
  wrong_param_ns_mesh(NOT start_ielm_,2,ns_mesh_trelms); 
  wrong_param_ns_mesh(NOT stop_ielm_,3,ns_mesh_trelms); 
  wrong_param_ns_mesh(start_ielm_[0]+1<1,2,ns_mesh_trelms); 
  wrong_param_ns_mesh(start_ielm_[0]+1>mesh_->nelm,2,ns_mesh_trelms); 
  wrong_param_ns_mesh(stop_ielm_[0]+1<1,3,ns_mesh_trelm); 
  wrong_param_ns_mesh(stop_ielm_[0]+1>mesh_->nelm,3,ns_mesh_trelms); 
  wrong_param_ns_mesh(NOT trelm_,4,ns_mesh_trelms); 
  wrong_param_ns_mesh(NOT trelmoff_,5,ns_mesh_trelms); 
  wrong_param_ns_mesh(trelmoff_[0]<4,5,ns_mesh_trelms); 
#endif
  
  const I n = stop_ielm_[0]-start_ielm_[0];
  cst_pR trelm = mesh_->trelm;
  /*  cst_pI cnc 	= mesh_->cnc;*/
  { I i;
    for (i=0;i<n;++i)
      {
	const I j 			= start_ielm_[0]+i;
	/*	cst_pI cncelm 			= &cnc[j*6+0];*/
	trelm_[trelmoff_[0]*i+0] 	= trelm[j*4+0];
	trelm_[trelmoff_[0]*i+1] 	= trelm[j*4+1];
	trelm_[trelmoff_[0]*i+2] 	= trelm[j*4+2];	
	trelm_[trelmoff_[0]*i+3] 	= trelm[j*4+3];
      } }
}


void ns_mesh_get_ddlcnc(const ns_mesh*	mesh_,
			const L    is_discont_,
			const I 	n_,
			const I 	ielm_,
			pI 		ddl_,
			const I 	dec)
{
  if (NOT is_discont_)
    {
      if (n_==1)
	{
	  ddl_[0]  = dec + ielm_;
	}
      else if (n_==3)
	{ 
	  ddl_[0] = dec + mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec + mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec + mesh_->cnc[ielm_*6+2]-1;
	}
      else if (n_==6)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +mesh_->cnc[ielm_*6+3]-1;
	  ddl_[4] = dec +mesh_->cnc[ielm_*6+4]-1;
	  ddl_[5] = dec +mesh_->cnc[ielm_*6+5]-1;
	}
      else if (n_==7)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +mesh_->cnc[ielm_*6+3]-1;
	  ddl_[4] = dec +mesh_->cnc[ielm_*6+4]-1;
	  ddl_[5] = dec +mesh_->cnc[ielm_*6+5]-1;
	  ddl_[6] = dec +mesh_->nvertex+mesh_->nedge+ielm_;
	}
      else if (n_==10)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +2*(mesh_->cnc[ielm_*6+3]-1);
	  ddl_[4] = dec +2*(mesh_->cnc[ielm_*6+3]-1)+1;
	  ddl_[5] = dec +2*(mesh_->cnc[ielm_*6+4]-1);
	  ddl_[6] = dec +2*(mesh_->cnc[ielm_*6+4]-1)+1;
	  ddl_[7] = dec +2*(mesh_->cnc[ielm_*6+5]-1);
	  ddl_[8] = dec +2*(mesh_->cnc[ielm_*6+5]-1)+1;
	  ddl_[9] = dec +mesh_->nvertex+2*mesh_->nedge+ielm_;
	}
      else if (n_==15)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +3*(mesh_->cnc[ielm_*6+3]-1);
	  ddl_[4] = dec +3*(mesh_->cnc[ielm_*6+3]-1)+1;
	  ddl_[5] = dec +3*(mesh_->cnc[ielm_*6+3]-1)+2;
	  ddl_[6] = dec +3*(mesh_->cnc[ielm_*6+4]-1);
	  ddl_[7] = dec +3*(mesh_->cnc[ielm_*6+4]-1)+1;
	  ddl_[8] = dec +3*(mesh_->cnc[ielm_*6+4]-1)+2;
	  ddl_[9] = dec +3*(mesh_->cnc[ielm_*6+5]-1);
	  ddl_[10] = dec +3*(mesh_->cnc[ielm_*6+5]-1)+1;
	  ddl_[11] = dec +3*(mesh_->cnc[ielm_*6+5]-1)+2;
	  ddl_[12] = dec +mesh_->nvertex+3*mesh_->nedge+ielm_;
	  ddl_[13] = dec +mesh_->nvertex+3*mesh_->nedge+ielm_+1;
	  ddl_[14] = dec +mesh_->nvertex+3*mesh_->nedge+ielm_+2;
	}
      else if (n_==21)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +4*(mesh_->cnc[ielm_*6+3]-1);
	  ddl_[4] = dec +4*(mesh_->cnc[ielm_*6+3]-1)+1;
	  ddl_[5] = dec +4*(mesh_->cnc[ielm_*6+3]-1)+2;
	  ddl_[6] = dec +4*(mesh_->cnc[ielm_*6+3]-1)+3;
	  ddl_[7] = dec +4*(mesh_->cnc[ielm_*6+4]-1);
	  ddl_[8] = dec +4*(mesh_->cnc[ielm_*6+4]-1)+1;
	  ddl_[9] = dec +4*(mesh_->cnc[ielm_*6+4]-1)+2;
	  ddl_[10] = dec +4*(mesh_->cnc[ielm_*6+4]-1)+3;
	  ddl_[11] = dec +4*(mesh_->cnc[ielm_*6+5]-1);
	  ddl_[12] = dec +4*(mesh_->cnc[ielm_*6+5]-1)+1;
	  ddl_[13] = dec +4*(mesh_->cnc[ielm_*6+5]-1)+2;
	  ddl_[14] = dec +4*(mesh_->cnc[ielm_*6+5]-1)+3;
	  ddl_[15] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_;
	  ddl_[16] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+1;
	  ddl_[17] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+2;
	  ddl_[18] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+3;
	  ddl_[19] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+4;
	  ddl_[20] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+5;
	}
      else if (n_==28)
	{ 
	  ddl_[0] = dec +mesh_->cnc[ielm_*6+0]-1;
	  ddl_[1] = dec +mesh_->cnc[ielm_*6+1]-1;
	  ddl_[2] = dec +mesh_->cnc[ielm_*6+2]-1;
	  ddl_[3] = dec +5*(mesh_->cnc[ielm_*6+3]-1);
	  ddl_[4] = dec +5*(mesh_->cnc[ielm_*6+3]-1)+1;
	  ddl_[5] = dec +5*(mesh_->cnc[ielm_*6+3]-1)+2;
	  ddl_[6] = dec +5*(mesh_->cnc[ielm_*6+3]-1)+3;
	  ddl_[7] = dec +5*(mesh_->cnc[ielm_*6+3]-1)+4;
	  ddl_[8] = dec +5*(mesh_->cnc[ielm_*6+4]-1);
	  ddl_[9] = dec +5*(mesh_->cnc[ielm_*6+4]-1)+1;
	  ddl_[10] = dec +5*(mesh_->cnc[ielm_*6+4]-1)+2;
	  ddl_[11] = dec +5*(mesh_->cnc[ielm_*6+4]-1)+3;
	  ddl_[12] = dec +5*(mesh_->cnc[ielm_*6+4]-1)+4;
	  ddl_[13] = dec +5*(mesh_->cnc[ielm_*6+5]-1);
	  ddl_[14] = dec +5*(mesh_->cnc[ielm_*6+5]-1)+1;
	  ddl_[15] = dec +5*(mesh_->cnc[ielm_*6+5]-1)+2;
	  ddl_[16] = dec +5*(mesh_->cnc[ielm_*6+5]-1)+3;
	  ddl_[17] = dec +5*(mesh_->cnc[ielm_*6+5]-1)+4;
	  ddl_[18] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_;
	  ddl_[19] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+1;
	  ddl_[20] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+2;
	  ddl_[21] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+3;
	  ddl_[22] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+4;
	  ddl_[23] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+5;
	  ddl_[24] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+6;
	  ddl_[25] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+7;
	  ddl_[26] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+8;
	  ddl_[27] = dec +mesh_->nvertex+4*mesh_->nedge+ielm_+9;
	}
#if 0
      else
	{
	  ns_err("n_>28");
	}
#endif
    }
  else
    {
      if (n_==1)
	{
	  ddl_[0]  =  ielm_;
	}
      else if (n_==3)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	}
      else if (n_==6)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	}
      else if (n_==7)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	  ddl_[6] = ielm_*n_+6;
	}
      else if (n_==10)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	  ddl_[6] = ielm_*n_+6;
	  ddl_[7] = ielm_*n_+7;
	  ddl_[8] = ielm_*n_+8;
	  ddl_[9] = ielm_*n_+9;
	}
      else if (n_==15)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	  ddl_[6] = ielm_*n_+6;
	  ddl_[7] = ielm_*n_+7;
	  ddl_[8] = ielm_*n_+8;
	  ddl_[9] = ielm_*n_+9;
	  ddl_[10] = ielm_*n_+10;
	  ddl_[11] = ielm_*n_+11;
	  ddl_[12] = ielm_*n_+12;
	  ddl_[13] = ielm_*n_+13;
	  ddl_[14] = ielm_*n_+14;
	}
      else if (n_==21)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	  ddl_[6] = ielm_*n_+6;
	  ddl_[7] = ielm_*n_+7;
	  ddl_[8] = ielm_*n_+8;
	  ddl_[9] = ielm_*n_+9;
	  ddl_[10] = ielm_*n_+10;
	  ddl_[11] = ielm_*n_+11;
	  ddl_[12] = ielm_*n_+12;
	  ddl_[13] = ielm_*n_+13;
	  ddl_[14] = ielm_*n_+14;
	  ddl_[15] = ielm_*n_+15;
	  ddl_[16] = ielm_*n_+16;
	  ddl_[17] = ielm_*n_+17;
	  ddl_[18] = ielm_*n_+18;
	  ddl_[19] = ielm_*n_+19;
	  ddl_[20] = ielm_*n_+20;
	}
      else if (n_==28)
	{ 
	  ddl_[0] = ielm_*n_+0;
	  ddl_[1] = ielm_*n_+1;
	  ddl_[2] = ielm_*n_+2;
	  ddl_[3] = ielm_*n_+3;
	  ddl_[4] = ielm_*n_+4;
	  ddl_[5] = ielm_*n_+5;
	  ddl_[6] = ielm_*n_+6;
	  ddl_[7] = ielm_*n_+7;
	  ddl_[8] = ielm_*n_+8;
	  ddl_[9] = ielm_*n_+9;
	  ddl_[10] = ielm_*n_+10;
	  ddl_[11] = ielm_*n_+11;
	  ddl_[12] = ielm_*n_+12;
	  ddl_[13] = ielm_*n_+13;
	  ddl_[14] = ielm_*n_+14;
	  ddl_[15] = ielm_*n_+15;
	  ddl_[16] = ielm_*n_+16;
	  ddl_[17] = ielm_*n_+17;
	  ddl_[18] = ielm_*n_+18;
	  ddl_[19] = ielm_*n_+19;
	  ddl_[20] = ielm_*n_+20;
	  ddl_[21] = ielm_*n_+21;
	  ddl_[22] = ielm_*n_+22;
	  ddl_[23] = ielm_*n_+23;
	  ddl_[24] = ielm_*n_+24;
	  ddl_[25] = ielm_*n_+25;
	  ddl_[26] = ielm_*n_+26;
	  ddl_[27] = ielm_*n_+27;
	}
    }
}




Err ns_mesh_print(const ns_mesh * 	mesh_,
		      const char * 	name_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT mesh_,1,ns_mesh_print); 
  wrong_param_ns_mesh(NOT name_,2,ns_mesh_print); 
#endif

  { char ctmp[512];
    sprintf(ctmp,"%s.mesh",name_);
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n"ifmt"\n",mesh_->nvertex);
    { I i;
      for (i=0;i<mesh_->nvertex;++i)
	{ 
	  fprintf(fil,
		  ""rfmt" "rfmt" "ifmt"\n",
		  mesh_->coo[i*2+0],
		  mesh_->coo[i*2+1],
		  mesh_->cod[i]); 
	} }
    fprintf(fil,"Triangles\n"ifmt"\n",mesh_->nelm); 
    { I i;
      for (i=0;i<mesh_->nelm;++i)
	{
	  fprintf(fil,
		  ""ifmt" "ifmt" "ifmt" "ifmt"\n",
		  mesh_->cnc[i*6+0],
		  mesh_->cnc[i*6+1],
		  mesh_->cnc[i*6+2],
		  ((I)0));
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  

  if (mesh_->time_x)
    {
      { char name2[256];
	sprintf(name2,"%s.nstime.mesh",name_);   
	FILE * fich = fopen(name2,"w");
	if (fich)
	  {
#if 0
#ifndef NDEBUG
	    ns_printf("read mesh time vector from '%s'\n",name2);
#endif
#endif
	    fprintf(fich,
		    ""ifmt"\n",
		    mesh_->time_n);
	    { I i;
	      for (i=0;i<mesh_->time_n;++i)
		{
		  fprintf(fich,""rfmt"\n",
			  mesh_->time_x[i]);
		} }
	    fclose(fich);
	  } }
    }  
  return __eErr_no;
}



void ns_mesh_edgejac(const ns_mesh*	mesh_,cst_pI iedge_,pR len_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT iedge_,1,ns_mesh_edgejac);
  wrong_param_ns_mesh(NOT len_,2,ns_mesh_edgejac);
  wrong_param_ns_mesh(iedge_[0]+1<1,1,ns_mesh_edgejac);
  wrong_param_ns_mesh(iedge_[0]>mesh_->nedge,1,ns_mesh_edgejac);
#endif
  const I ielm_[1] = {mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]]/*-1*/};
  if (mesh_->cnc[ielm_[0]*6+3+0]==mesh_->nvertex+iedge_[0]+1)
    {
      const I i = mesh_->cnc[ielm_[0]*6+0]-1;
      const I j = mesh_->cnc[ielm_[0]*6+1]-1;
      len_[0] = nsSQRT((mesh_->coo[2*i+0]-mesh_->coo[2*j+0])*(mesh_->coo[2*i+0]-mesh_->coo[2*j+0])+(mesh_->coo[2*i+1]-mesh_->coo[2*j+1])*(mesh_->coo[2*i+1]-mesh_->coo[2*j+1]))*((R)0.5);
      return;
    }
  if (mesh_->cnc[ielm_[0]*6+3+1]==mesh_->nvertex+iedge_[0]+1)
    {
      const I i = mesh_->cnc[ielm_[0]*6+1]-1;
      const I j = mesh_->cnc[ielm_[0]*6+2]-1;
      len_[0] = nsSQRT((mesh_->coo[2*i+0]-mesh_->coo[2*j+0])*(mesh_->coo[2*i+0]-mesh_->coo[2*j+0])+(mesh_->coo[2*i+1]-mesh_->coo[2*j+1])*(mesh_->coo[2*i+1]-mesh_->coo[2*j+1]))*((R)0.5);
      return;
    }
  if (mesh_->cnc[ielm_[0]*6+3+2]==mesh_->nvertex+iedge_[0]+1)
    {
      const I i = mesh_->cnc[ielm_[0]*6+2]-1;
      const I j = mesh_->cnc[ielm_[0]*6+0]-1;
      len_[0] = nsSQRT((mesh_->coo[2*i+0]-mesh_->coo[2*j+0])*(mesh_->coo[2*i+0]-mesh_->coo[2*j+0])+(mesh_->coo[2*i+1]-mesh_->coo[2*j+1])*(mesh_->coo[2*i+1]-mesh_->coo[2*j+1]))*((R)0.5);
      return;    
    }
}



void ns_mesh_normal(const ns_mesh * 	mesh_,
		    cst_pI 		iedge_,
		    pR 		normal_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT iedge_,1,ns_mesh_normal);
  wrong_param_ns_mesh(NOT normal_,2,ns_mesh_normal);
  wrong_param_ns_mesh(iedge_[0]+1<1,1,ns_mesh_normal);
  wrong_param_ns_mesh(iedge_[0]>mesh_->nedge,1,ns_mesh_normal);
#endif
  const I ielm_[1] = {mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]]/*-1*/};
  I k;
  for (k=0;k<3;++k)
    if (mesh_->cnc[ielm_[0]*6+3+k]==mesh_->nvertex+iedge_[0]+1)
      break;
#ifndef NDEBUG
  if (k>=3)
    {
      fprintf(stderr,"ns_mesh_normal:severe error "ifmt"\n",iedge_[0]);
      exit(1);
    }
#endif
  const I i  = mesh_->cnc[ielm_[0]*6+k]-1;
  const I j  = mesh_->cnc[ielm_[0]*6+(k+1)%3]-1;
  normal_[1] 	 = mesh_->coo[2*i+0]-mesh_->coo[2*j+0];
  normal_[0] 	 = mesh_->coo[2*j+1]-mesh_->coo[2*i+1];
  const R a = ((R)1.0)/nsSQRT(normal_[0]*normal_[0]+normal_[1]*normal_[1]);
  normal_[0] 	*= a;
  normal_[1] 	*= a;
}



void ns_mesh_edge(const ns_mesh * mesh_,
		  cst_pI 	iedge_,
		  pI 		nadj_,
		  pI 		ielm_,
		  pI           ilocedge_,
		  pI 		jelm_,
		  pI           jlocedge_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT iedge_,1,ns_mesh_edge);
  wrong_param_ns_mesh(iedge_[0]+1<1,1,ns_mesh_edge);
  wrong_param_ns_mesh(iedge_[0]>mesh_->nedge,1,ns_mesh_edge);
  wrong_param_ns_mesh(NOT nadj_,2,ns_mesh_edge);
  wrong_param_ns_mesh(NOT ielm_,3,ns_mesh_edge);
  wrong_param_ns_mesh(NOT ilocedge_,4,ns_mesh_edge);
  wrong_param_ns_mesh(NOT jelm_,5,ns_mesh_edge);
  wrong_param_ns_mesh(NOT jlocedge_,6,ns_mesh_edge);
#endif

  nadj_[0] 	= mesh_->bpatch[mesh_->nvertex+iedge_[0]+1]-mesh_->bpatch[mesh_->nvertex+iedge_[0]];
  ielm_[0] 	= mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]]-1;
  ilocedge_[0] 	= 0;
  if (mesh_->cnc[ielm_[0]*6+3+0]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 0;
  else if (mesh_->cnc[ielm_[0]*6+3+1]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 1;
  else if (mesh_->cnc[ielm_[0]*6+3+2]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 2;
  if (nadj_[0]==2)
    {
      jelm_[0] = mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]+1]-1;
      if (mesh_->cnc[jelm_[0]*6+3+0]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 0;
      else if (mesh_->cnc[jelm_[0]*6+3+1]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 1;
      else if (mesh_->cnc[jelm_[0]*6+3+2]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 2;
    } 
  else if (nadj_[0]!=1)
    {
      fprintf(stderr,"ns_mesh_edge:erreur nadj = "ifmt"\n",nadj_[0]);
    }
}


void ns_mesh_edgeinfo(const ns_mesh * mesh_,
		  cst_pI 	iedge_,
		  pI 		nadj_,
		  pI 		ielm_,
		  pI           ilocedge_,
		  pI 		jelm_,
		  pI           jlocedge_)
{
#ifndef NDEBUG
  wrong_param_ns_mesh(NOT iedge_,1,ns_mesh_edge);
  wrong_param_ns_mesh(iedge_[0]+1<1,1,ns_mesh_edge);
  wrong_param_ns_mesh(iedge_[0]>mesh_->nedge,1,ns_mesh_edge);
  wrong_param_ns_mesh(NOT nadj_,2,ns_mesh_edge);
  wrong_param_ns_mesh(NOT ielm_,3,ns_mesh_edge);
  wrong_param_ns_mesh(NOT ilocedge_,4,ns_mesh_edge);
  wrong_param_ns_mesh(NOT jelm_,5,ns_mesh_edge);
  wrong_param_ns_mesh(NOT jlocedge_,6,ns_mesh_edge);
#endif

  nadj_[0] 	= mesh_->bpatch[mesh_->nvertex+iedge_[0]+1]-mesh_->bpatch[mesh_->nvertex+iedge_[0]];
  ielm_[0] 	= mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]];
  ilocedge_[0] 	= 0;
  if (mesh_->cnc[ielm_[0]*6+3+0]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 0;
  else if (mesh_->cnc[ielm_[0]*6+3+1]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 1;
  else if (mesh_->cnc[ielm_[0]*6+3+2]==mesh_->nvertex+iedge_[0]+1)
    ilocedge_[0] = 2;
  if (nadj_[0]==2)
    {
      jelm_[0] = mesh_->patch[mesh_->bpatch[mesh_->nvertex+iedge_[0]]+1];
      if (mesh_->cnc[jelm_[0]*6+3+0]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 0;
      else if (mesh_->cnc[jelm_[0]*6+3+1]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 1;
      else if (mesh_->cnc[jelm_[0]*6+3+2]==mesh_->nvertex+iedge_[0]+1)
	jlocedge_[0] = 2;
    } 
  else if (nadj_[0]!=1)
    {
      fprintf(stderr,"ns_mesh_edge:erreur nadj = "ifmt"\n",nadj_[0]);
    }
}


Err ns_mesh_read(ns_mesh * 	mesh_,
		 const char * 	name_,...)
{  
  STR name;
  { va_list args;
    va_start (args,name_);
    vsprintf(name,name_,args);
    va_end(args); }    

  pMedit medit_reader = Medit_new(name);
  if (!medit_reader)
    {
      fprintf(stderr,"ns_mesh_read:can not read '%s'\n",name);
      return __eErr_file;
    }

  const I medit_dimension = Medit_get_eDim(medit_reader);
  if (((I)2)!=medit_dimension)
    {
      fprintf(stderr,"ns_mesh_read:expect dimension 2\n");
      return __eErr_file;
    }
  memset(mesh_,0,sizeof(ns_mesh));
  mesh_->dimension	=  __eDim_2;
  mesh_->nvertex	=  Medit_get_nbVertices(medit_reader);
  mesh_->nelm  		=  Medit_get_nbFaces(medit_reader,__eFace_TRIANGLE);  
  mesh_->elm 		=  __eElement_TRIA;

#if 0
  int ii=0;
  if (mesh_->nelm==0)
    {
      ii=1;
      mesh_->nelm =  (I)GmfStatKwd(inm,GmfTrianglesP2);  
      printf("pppppp "ifmt"\n",mesh_->nelm);
    }
#endif

  mesh_->own_cnc  	= calloc(mesh_->nelm*((I)6),sizeof(I));
  mesh_->own_adj  	= malloc(sizeof(I)*mesh_->nelm*((I)3));
  mesh_->cnc 		= mesh_->own_cnc;
  mesh_->adj		= mesh_->own_adj;
  
  const I cncoff = 6;
  const I adjoff = 3;

  unsigned char * medit_voy = Medit_get_adjcncs		(medit_reader,
							 __emnsNO,/*const L 	CIndexation_*/
							 mesh_->own_cnc,
							 &cncoff,
							 NULL,
							 NULL,
							 mesh_->own_adj,
							 &adjoff);

  /*\brief connectivity*/
  { I N = mesh_->nvertex;
    I i;
    for (i=0;i<mesh_->nelm;++i)
      {
	{ I j;
	  for (j=0;j<3;++j)
	    {
	      if (NOT mesh_->own_cnc[6*i+(j+3)])
		{
		  if (!mesh_->own_adj[3*i+j])
		    {
		      mesh_->own_cnc[6*i+(j+3)] = ++N;
		    }
		  else
		    {
		      if (NOT mesh_->own_cnc[6*(mesh_->own_adj[3*i+j]-1)+((medit_voy[1+3*i+j]+1)%3+3)])
			{
			  mesh_->own_cnc[6*i+(j+3)] = ++N;
			  mesh_->own_cnc[6*(mesh_->own_adj[3*i+j]-1)+((medit_voy[1+3*i+j]+1)%3+3)] 	= N;		      
			}
		      else
			{
			  mesh_->own_cnc[6*i+(j+3)] = mesh_->own_cnc[6*(mesh_->own_adj[3*i+j]-1)+((medit_voy[1+3*i+j]+1)%3+3)];
			}
		    }
		}
	    } }
      }       
    mesh_->nedge = N-mesh_->nvertex;  }
  
  free(medit_voy); 
  mesh_->nedge_boundary = 0;    
  { I i;
    for (i=0;i<mesh_->nelm;++i)
      {
	{ I j;
	  for (j=0;j<3;++j)
	    {
	      if (NOT mesh_->own_adj[3*i+j])
		{
		  ++mesh_->nedge_boundary;
		}
	    } } 
      } }

  cst_eElement elm = __eElement_TRIA;

  mesh_->nddl_lagr[__ensBASIS_error] 		= (I)0;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGE_0] 	= mesh_->nelm;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGE_1] 	= mesh_->nvertex;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGE_2] 	= mesh_->nvertex+mesh_->nedge;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGE_3] 	= mesh_->nvertex+mesh_->nedge*2+mesh_->nelm;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGEBUBBLE_0]	= 2*mesh_->nelm;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGEBUBBLE_1]	= mesh_->nvertex+mesh_->nelm;
  mesh_->nddl_lagr[__ensBASIS_LAGRANGEBUBBLE_2]	= mesh_->nvertex+mesh_->nedge + mesh_->nelm;
  mesh_->nddl_lagr[__ensBASIS_L2ORTHONORMAL_0] 	= mesh_->nelm*ensBASIS_n(__ensBASIS_L2ORTHONORMAL_0,elm);
  mesh_->nddl_lagr[__ensBASIS_L2ORTHONORMAL_1] 	= mesh_->nelm*ensBASIS_n(__ensBASIS_L2ORTHONORMAL_1,elm);
  mesh_->nddl_lagr[__ensBASIS_L2ORTHONORMAL_2] 	= mesh_->nelm*ensBASIS_n(__ensBASIS_L2ORTHONORMAL_2,elm);
  mesh_->nddl_lagr[__ensBASIS_L2ORTHONORMAL_3] 	= mesh_->nelm*ensBASIS_n(__ensBASIS_L2ORTHONORMAL_3,elm);

  mesh_->nddlP6				= mesh_->nddl_lagr[__ensBASIS_LAGRANGE_2];
  mesh_->own_coo 			= malloc(sizeof(R)*2*mesh_->nddlP6);
  mesh_->own_cod 			= malloc(sizeof(I)*mesh_->nddlP6);
  mesh_->cod = mesh_->own_cod;
  mesh_->coo = mesh_->own_coo;
  
  const I notranspose 	= 0;
  const I cooff 	= 2;
  const I codoff 	= 1;
  Medit_get_cncv		(medit_reader,
				 &notranspose,
				 mesh_->own_coo,
				 &cooff,
				 mesh_->own_cod,
				 &codoff);

  mesh_->own_jacelm	= malloc(sizeof(R)*mesh_->nelm);
  mesh_->own_trelm	= malloc(sizeof(R)*mesh_->nelm*4);
  mesh_->jacelm=mesh_->own_jacelm;
  mesh_->trelm=mesh_->own_trelm;
  { I i;
    R tmp[3],tmpcooelm[6];
    for (i=0;i<mesh_->nelm;++i)
      {
	tmpcooelm[0] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)];
	tmpcooelm[1] = mesh_->coo[2*(mesh_->cnc[6*i+1]-1)];
	tmpcooelm[2] = mesh_->coo[2*(mesh_->cnc[6*i+2]-1)];
	tmpcooelm[3] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)+1];
	tmpcooelm[4] = mesh_->coo[2*(mesh_->cnc[6*i+1]-1)+1];
	tmpcooelm[5] = mesh_->coo[2*(mesh_->cnc[6*i+2]-1)+1];
	tmp[0]    = nsSQRT( (tmpcooelm[0]-tmpcooelm[1])*(tmpcooelm[0]-tmpcooelm[1]) + (tmpcooelm[3]-tmpcooelm[4])*(tmpcooelm[3]-tmpcooelm[4]) ); 
	tmp[1]    = nsSQRT( (tmpcooelm[1]-tmpcooelm[2])*(tmpcooelm[1]-tmpcooelm[2]) + (tmpcooelm[4]-tmpcooelm[5])*(tmpcooelm[4]-tmpcooelm[5]) ); 
	tmp[2]    = nsSQRT( (tmpcooelm[2]-tmpcooelm[0])*(tmpcooelm[2]-tmpcooelm[0]) + (tmpcooelm[5]-tmpcooelm[3])*(tmpcooelm[5]-tmpcooelm[3]) ); 	    
	if (tmp[1]<tmp[2]) { const R x = tmp[1]; tmp[1] = tmp[2]; tmp[2] = x; }
	if (tmp[0]<tmp[1]) { const R x = tmp[0]; tmp[0] = tmp[1]; tmp[1] = x; }
	if (tmp[1]<tmp[2]) { const R x = tmp[1]; tmp[1] = tmp[2]; tmp[2] = x; }
	mesh_->own_jacelm[i] = nsSQRT((tmp[0]+tmp[1]+tmp[2])*(tmp[2]-(tmp[0]-tmp[1]))*(tmp[2]+(tmp[0]-tmp[1]) )*(tmp[0]+(tmp[1]-tmp[2]) ))*((R)0.5);  
      } }


  { I i;
    R tmpcooelm[6];
    for (i=0;i<mesh_->nelm;++i)
      {
	tmpcooelm[0] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+0]-1)];
	tmpcooelm[1] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+1]-1)];
	tmpcooelm[2] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+2]-1)];
	tmpcooelm[3] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+0]-1)+1];
	tmpcooelm[4] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+1]-1)+1];
	tmpcooelm[5] = mesh_->own_coo[2*(mesh_->own_cnc[6*i+2]-1)+1];
	const R a 	= mesh_->own_jacelm[i];
	const R b 	= ((R)1.0)/a;
	mesh_->own_trelm[4*i+3]   	= (tmpcooelm[1]-tmpcooelm[0])*b;
	mesh_->own_trelm[4*i+1]   	= (tmpcooelm[0]-tmpcooelm[2])*b;
	mesh_->own_trelm[4*i+2]   	= (tmpcooelm[3]-tmpcooelm[4])*b;
	mesh_->own_trelm[4*i+0]   	= (tmpcooelm[5]-tmpcooelm[3])*b;
      } }

  { I k;
    for (k=0;k<mesh_->nelm;++k)
      {
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+0]-1)+0] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+0]-1)+0]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+1]-1)+0])*((R)0.5);
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+0]-1)+1] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+0]-1)+1]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+1]-1)+1])*((R)0.5);
	
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+1]-1)+0] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+1]-1)+0]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+2]-1)+0])*((R)0.5);
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+1]-1)+1] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+1]-1)+1]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+2]-1)+1])*((R)0.5);
	
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+2]-1)+0] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+2]-1)+0]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+0]-1)+0])*((R)0.5);
	mesh_->own_coo[2*(mesh_->own_cnc[6*k+3+2]-1)+1] = (mesh_->own_coo[2*(mesh_->own_cnc[6*k+2]-1)+1]+mesh_->own_coo[2*(mesh_->own_cnc[6*k+0]-1)+1])*((R)0.5);
	
	mesh_->own_cod[mesh_->own_cnc[6*k+3+0]-1] 	    = MAX(mesh_->own_cod[mesh_->own_cnc[6*k+0]-1],mesh_->own_cod[mesh_->own_cnc[6*k+1]-1]);
	mesh_->own_cod[mesh_->own_cnc[6*k+3+1]-1] 	    = MAX(mesh_->own_cod[mesh_->own_cnc[6*k+1]-1],mesh_->own_cod[mesh_->own_cnc[6*k+2]-1]);
	mesh_->own_cod[mesh_->own_cnc[6*k+3+2]-1] 	    = MAX(mesh_->own_cod[mesh_->own_cnc[6*k+2]-1],mesh_->own_cod[mesh_->own_cnc[6*k+0]-1]);
      } }    

  mesh_->own_bpatch		= calloc(mesh_->nddlP6+1,sizeof(I));
  mesh_->own_patch		= malloc(sizeof(I)*mesh_->nelm*6);
  mesh_->bpatch=mesh_->own_bpatch;
  mesh_->patch=mesh_->own_patch;
  { I k,j;
    for (k=0; k<mesh_->nelm; k++) 
      {
	mesh_->own_bpatch[mesh_->own_cnc[k*6+0]]+=1;
	mesh_->own_bpatch[mesh_->own_cnc[k*6+1]]+=1;
	mesh_->own_bpatch[mesh_->own_cnc[k*6+2]]+=1;
	mesh_->own_bpatch[mesh_->own_cnc[k*6+3]]+=1;
	mesh_->own_bpatch[mesh_->own_cnc[k*6+4]]+=1;
	mesh_->own_bpatch[mesh_->own_cnc[k*6+5]]+=1;
      }        
    for (k=2;k<=mesh_->nddlP6;++k)
      mesh_->own_bpatch[k]+=mesh_->own_bpatch[k-1];
    for (k=0;k<mesh_->nelm;++k)
      for (j=0;j<6;++j)
	{
	  mesh_->own_patch[mesh_->own_bpatch[mesh_->own_cnc[k*6+j]-1]]=k;
	  mesh_->own_bpatch[mesh_->own_cnc[k*6+j]-1]+=1;
	}
    for (k=0;k<=mesh_->nddlP6;++k)
      mesh_->own_bpatch[k]=0;
    for (k=0;k<mesh_->nelm;++k)
      for (j=0;j<6;++j)
	mesh_->own_bpatch[mesh_->own_cnc[k*6+j]]+=1;
    for (k=2;k<=mesh_->nddlP6;++k)
      mesh_->own_bpatch[k]+=mesh_->own_bpatch[k-1];
  }    

  mesh_->own_jacedge = malloc(sizeof(R)*mesh_->nedge);
  mesh_->jacedge=mesh_->own_jacedge;
  {
    I iedge;
    for (iedge = 0;iedge<mesh_->nedge;++iedge)
      {
	ns_mesh_edgejac(mesh_,&iedge,&mesh_->own_jacedge[iedge]);
      }
  }
  mesh_->own_normaledge = malloc(sizeof(R)*mesh_->nedge*2);
  mesh_->normaledge=mesh_->own_normaledge;
  {
    /*
      on calcule la normale vers l'exterieur du premier element du patch de l'arete
    */
    I iedge;
    for (iedge = 0;iedge<mesh_->nedge;++iedge)
      {
	ns_mesh_normal(mesh_,&iedge,&mesh_->own_normaledge[2*iedge]);	
      }
  }
  mesh_->noboundary_vcod = 100;

  /*
    read time 
  */
  { 
    STR name2,basenam;
    cst_eFileSuffix fileSuffix = eFileSuffix_get_suffixAndBasename(name,
								   basenam);
    
    sprintf(name2,"%s.nstime.mesh",basenam);   
    FILE * fich = fopen(name2,"r");
    if (fich)
      {
#if 0
#ifndef NDEBUG
	ns_printf("read mesh time vector from '%s'\n",name2);
#endif
#endif
	I nn=0;
	int ii = fscanf(fich,"%"nsFORMAT_INTEGER" %"nsFORMAT_INTEGER"",&mesh_->time_n,&nn);
	if (ii)
	  {

	  }
	mesh_->time_x = malloc(sizeof(R)*mesh_->time_n);
	{ I i;
	  for (i=0;i<mesh_->time_n;++i)
	    {
	      int ii = fscanf(fich,"%lf",&mesh_->time_x[i]);
	      if (ii)
		{}
	    } }
	fclose(fich);
#if 0
	ns_mesh_extrude_medit	(mesh_,
				 "extrude",
				 mesh_->time_n,
				 mesh_->time_n,
				 mesh_->time_x,
				 1.0);
#endif

      } }
  return __eErr_no;

}




#if 0
const nsCUBA_TRIA_INTERVAL_ST * ns_mesh_getcubature_interval(const ns_mesh * 	mesh_,
								   const I 		degree_)
{
  if (mesh_->instances_cubature_triangle_interval[degree_])
    return mesh_->instances_cubature_triangle_interval[degree_];
  else
    {
      ns_errmsg("ns_mesh_getcubature_interval:mesh_ do not have gauss cubature on triangle_interval with degree_="ifmt"",degree_);
      return NULL;
    }
}


const nsCUBA_TRIA_INTERVAL_ST * ns_mesh_addcubature_interval(ns_mesh * mesh_,const I degree_)
{
  if (NOT mesh_->instances_cubature_triangle_interval[degree_])
    mesh_->instances_cubature_triangle_interval[degree_]=nsCUBA_TRIA_INTERVAL_gauss_new(&degree_);
  return mesh_->instances_cubature_triangle_interval[degree_];
}



const nsCUBA_TRIA_ST * ns_mesh_getcubature(const ns_mesh * 	mesh_,
						 const I 		degree_)
{
  if (mesh_->instances_cubature_triangle[degree_])
    return mesh_->instances_cubature_triangle[degree_];
  else
    {
      ns_errmsg("ns_mesh_getcubature:mesh_ do not have gauss cubature on triangle with degree_="ifmt"",degree_);
      return NULL;
    }
}


const nsCUBA_TRIA_ST * ns_mesh_addcubature(ns_mesh * mesh_,const I degree_)
{
  if (NOT mesh_->instances_cubature_triangle[degree_])
    mesh_->instances_cubature_triangle[degree_]=nsCUBA_TRIA_gauss_new(degree_);
  return mesh_->instances_cubature_triangle[degree_];
}
#endif




#if 0
Err ns_mesh_from_mkZ(ns_mesh * 	mesh_,
			 const char * 	name_)
{  
  int ver;
  int medit_dim;
  char name[512];
  sprintf(name,"%s",name_);
  memset(mesh_,0,sizeof(ns_mesh));
  mesh_->dimension	= (I)zone_->dimension;
  mesh_->nvertex	= mkZ_nvertex(zone_);
  mesh_->nelm  		= mkZ_nelm(zone_);

  cst_mkFE spaceP2 	= mkZ_addspace_quadratic(zone_);
  mesh_->cnc  		= cst_pIM_x(cst_mkFE_cnc(spaceP2));
  mesh_->adj  		= cst_pIM_x(cst_mkZ_adj(zone));
  mesh_->coo 		= malloc(sizeof(R)*medit_dim*mesh_->nddlP6);
  mesh_->cod 		= malloc(sizeof(I)*mesh_->nddlP6);
  mesh_->nedge 		= N-mesh_->nvertex;  
  mesh_->nedge_boundary = 0;    
  mesh_->nddl_lagr[__ens_lagr_tria_error] = (I)0;
  mesh_->nddl_lagr[__ens_lagr_tria_01] 	= mesh_->nelm;
  mesh_->nddl_lagr[__ens_lagr_tria_03] 	= mesh_->nvertex;
  mesh_->nddl_lagr[__ens_lagr_tria_06] 	= mesh_->nvertex+mesh_->nedge;
  mesh_->nddl_lagr[__ens_lagr_tria_07] 	= mesh_->nvertex+mesh_->nedge + mesh_->nelm;
  mesh_->nddl_lagr[__ens_lagr_tria_10] 	= mesh_->nvertex+mesh_->nedge*2+mesh_->nelm;
  mesh_->nddl_lagr[__ens_lagr_tria_15] 	= mesh_->nvertex+mesh_->nedge*3+mesh_->nelm*3;
  mesh_->nddl_lagr[__ens_lagr_tria_21] 	= mesh_->nvertex+mesh_->nedge*4+mesh_->nelm*6;
  mesh_->nddl_lagr[__ens_lagr_tria_28]	= mesh_->nvertex+mesh_->nedge*5+mesh_->nelm*10;
  mesh_->nddlP6 			= mesh_->nddl_lagr[__ens_lagr_tria_06];
  mesh_->jacelm		= malloc(sizeof(R)*mesh_->nelm);
  mesh_->trelm		= malloc(sizeof(R)*mesh_->nelm*4);
  { I i;
    R tmpcooelm[6];
    for (i=0;i<mesh_->nelm;++i)
      {
	tmpcooelm[0] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)];
	tmpcooelm[1] = mesh_->coo[2*(mesh_->cnc[6*i+1]-1)];
	tmpcooelm[2] = mesh_->coo[2*(mesh_->cnc[6*i+2]-1)];
	tmpcooelm[3] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)+1];
	tmpcooelm[4] = mesh_->coo[2*(mesh_->cnc[6*i+1]-1)+1];
	tmpcooelm[5] = mesh_->coo[2*(mesh_->cnc[6*i+2]-1)+1];
	const R a 	= mesh_->jacelm[i];
	const R b 	= ((R)1.0)/a;
	mesh_->trelm[4*i+3]   	= (tmpcooelm[1]-tmpcooelm[0])*b;
	mesh_->trelm[4*i+1]   	= (tmpcooelm[0]-tmpcooelm[2])*b;
	mesh_->trelm[4*i+2]   	= (tmpcooelm[3]-tmpcooelm[4])*b;
	mesh_->trelm[4*i+0]   	= (tmpcooelm[5]-tmpcooelm[3])*b;
      } }
  mesh_->jacedge = malloc(sizeof(R)*mesh_->nedge);
  {
    I iedge;
    for (iedge = 0;iedge<mesh_->nedge;++iedge)
      {
	ns_mesh_edgejac(mesh_,&iedge,&mesh_->jacedge[iedge]);
      }
  }
  mesh_->normaledge = malloc(sizeof(R)*mesh_->nedge*2);
  {
    /*
      on calcule la normale vers l'exterieur du premier element du patch de l'arete
    */
    I iedge;
    for (iedge = 0;iedge<mesh_->nedge;++iedge)
      {
	ns_mesh_normal(mesh_,&iedge,&mesh_->normaledge[2*iedge]);	
      }
  }
  mesh_->noboundary_vcod = 100;
  return __eErr_no;
}
#endif








void ns_mesh_set_periodic(ns_mesh*mesh_)  
{
  I i,j;
  I count23=0,*id23;
  I count14=0,*id14;
  for (i=0;i<mesh_->nvertex;++i)    
  {
    R yi,yj;
    for (i=0;i<mesh_->nvertex;++i)    
      {
	if (mesh_->cod[i]==23)
	  {
	    yi = mesh_->coo[2*i+1];
	    for (j=0;j<mesh_->nvertex;++j)    
	      {
		if (mesh_->cod[j]==41)
		  {
		    yj = mesh_->coo[2*j+1];
		    if (nsFABS(yi-yj)<1.0e-13)
		      break;
		  }
	      }
#if 0
	    if (j>mesh_->nvertex)
	      ns_err("alllllllll");
#endif
	  }
      }
  }
  for (i=0;i<mesh_->nelm;++i)    
    {
      for (j=0;j<3;++j)    
	{
	  if (NOT mesh_->adj[3*i+j])
	    {
	      if (
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==23)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==2)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==3)
		  )
		{
		  if (
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==23)
		      OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==2)
		      OR
			(mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==3)
		      )
		    {	
		      ++count23;
		      break;
		    }
		}
	    }
	}
    }
  for (i=0;i<mesh_->nelm;++i)    
    {
      for (j=0;j<3;++j)    
	{
	  if (NOT mesh_->adj[3*i+j])
	    {
	      if (
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==41)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==1)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==4)
		  )
		{
		  if (
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==41)
		      OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==1)
		      OR
			(mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==4)
		      )
		    {	
		      ++count14;
			break;
		    }
		}
	    }
	}
    }
  I * jd14,*jd23;
  id14 = (I*)malloc(sizeof(I)*(count14+1));
  id23 = (I*)malloc(sizeof(I)*(count23+1));
  jd14 = (I*)malloc(sizeof(I)*(count14+1));
  jd23 = (I*)malloc(sizeof(I)*(count23+1));
#if 0
  ns_printf("modif mesh %ld %ld\n",count14,count23);
#endif
  count23=0;
  count14=0;  
  for (i=0;i<mesh_->nelm;++i)    
    {
      for (j=0;j<3;++j)    
	{
	  if (NOT mesh_->adj[3*i+j])
	    {
	      if (
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==23)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==2)
		  OR
		    (mesh_->cod[mesh_->cnc[6*i+j]-1]==3)
		  )
		{
		  if (
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==23)
			OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==2)
		      OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==3)
		      )
		    {	
		      id23[count23] = i;
		      jd23[count23] = j;
		      count23++;
		      break;
		    }
		}
	    }
	}
    }
  
  for (i=0;i<mesh_->nelm;++i)    
    {
      for (j=0;j<3;++j)    
	{
	  if (NOT mesh_->adj[3*i+j])
	    {
	      if (
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==41)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==1)
		  OR
		  (mesh_->cod[mesh_->cnc[6*i+j]-1]==4)
		  )
		{
		  if (
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==41)
		      OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==1)
			OR
		      (mesh_->cod[mesh_->cnc[6*i+(j+1)%3]-1]==4)
		      )
		    {	
		      id14[count14] = i;
		      jd14[count14] = j;
		      ++count14;
		      break;
		    }
		}
	    }
	}
    }
			  
#if 0			  
  if (count23!=count14) ns_err("prprprprpr "ifmt" "ifmt"",count23,count14);  
#endif
  for (i=0;i<count23;++i)    
    {
      mesh_->own_adj[3*id23[i]+jd23[i]] = id14[i]+1;
    }

  for (i=0;i<count14;++i)    
    {
      mesh_->own_adj[3*id14[i]+jd14[i]] = id23[i]+1;
    }
  for (i=0;i<mesh_->nvertex;++i)    
    {
      if ( (mesh_->cod[i]==41) OR (mesh_->cod[i]==23) )
	{
	  mesh_->own_cod[i] = mesh_->noboundary_vcod;
	}
    }
  free(id14);
  free(id23);
  free(jd14);
  free(jd23);
}





void ns_mesh_write_medit_with_codelm(const ns_mesh*s_,cst_pI codelm_,const char * name_,...)
{
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n"ifmt"\n",s_->nvertex);
    {I i;for (i=0;i<s_->nvertex;++i){ fprintf(fil,""rfmt" "rfmt" "ifmt"\n",s_->coo[i*2+0],s_->coo[i*2+1],s_->cod[i]); } }
    fprintf(fil,"Triangles\n"ifmt"\n",s_->nelm); 
    I cncelm[3];
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_cncelm0(s_,&i,cncelm);
	  fprintf(fil,""ifmt" "ifmt" "ifmt" "ifmt"\n",cncelm[0]+1,cncelm[1]+1,cncelm[2]+1,codelm_[i]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}

void ns_mesh_write_medit(const ns_mesh*s_,const char * name_,...)
{
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n"ifmt"\n",s_->nvertex);
    {I i;for (i=0;i<s_->nvertex;++i){ fprintf(fil,""rfmt" "rfmt" "ifmt"\n",s_->coo[i*2+0],s_->coo[i*2+1],s_->cod[i]); } }
    fprintf(fil,"Triangles\n"ifmt"\n",s_->nelm); 
    I cncelm[3];
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_cncelm0(s_,&i,cncelm);
	  fprintf(fil,""ifmt" "ifmt" "ifmt" "ifmt"\n",cncelm[0]+1,cncelm[1]+1,cncelm[2]+1,((I)0));
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void ns_mesh_extrude_scalar	(const ns_mesh * 	mesh_,
				 const char * 	name_,
				 const R 	coeff_,
				 const I	nframe__,
				 const I	nz_,
				 cst_pR 	z_,
				 cst_pR 	*scalars_)
{  
  const I 	nvertex_ 	= mesh_->nvertex;
  const I 	nelm_ 		= mesh_->nelm;
  cst_pI	cnc_		= mesh_->cnc;
  cst_pI	cod_		= mesh_->cod;
  cst_pR	coo_		= mesh_->coo;
				 
  I irot;
  I i,j;
  STR obasename;
  cst_eFileSuffix fileSuffix = eFileSuffix_get_suffixAndBasename(name_,
								 obasename);

  char bbname[256];
  I nframe = nframe__;
  for (irot=1;irot<nframe__;++irot)
    {
      if (z_[irot]==z_[0])
	break;
    }
  nframe = irot;  
  if (scalars_)
    { 
      sprintf(bbname,"%s.bb",obasename);
      FILE * out = fopen(bbname,"w");
      fprintf(out,"3 1 "ifmt" 2\n",nframe*nvertex_);  
      for (irot=0;irot<nframe;++irot)
	{
	  cst_pR scalar_ = (nz_==nframe__)?scalars_[irot]:scalars_[0];
	  for (i=0;i<nvertex_;++i)
	    {      
	      fprintf(out," "rfmt"\n",scalar_[i]);	  
	    } 
	}
      fclose(out);
    }

  sprintf(bbname,"%s.mesh",obasename);
  
  FILE * out 			= fopen(bbname,"w");
  fprintf(out,"MeshVersionFormatted 1\nDimension 3\nVertices "ifmt"\n",nframe*nvertex_);
#if 0
  for (i=0;i<nvertex_;++i)
    {
      fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",coo_[2*i+0],coo_[2*i+1],((R)0.0),cod_[i]);	  
    }  
#endif
  for (irot=0;irot<nframe;++irot)
    {
      for (i=0;i<nvertex_;++i)
	{
	  fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",coo_[2*i+0],coo_[2*i+1],(nz_==nframe__)?z_[irot]*coeff_:z_[0]*((R)irot),cod_[i]);	  
	}  
    }
  fprintf(out,"Tetrahedra\n"ifmt"\n",3*nelm_*(nframe-1));
  for (irot=1;irot<nframe;++irot)
    {
      {
	for (i=0;i<nelm_;++i)
	  {
	    I mx = cnc_[6*i+0];
	    for (j=1;j<3;++j)
	      if (mx<cnc_[6*i+j])
		{
		  mx = cnc_[6*i+j];
		}	    
	    {
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			(irot-1)*nvertex_+cnc_[6*i+0],
			(irot-1)*nvertex_+cnc_[6*i+1],
			(irot-1)*nvertex_+cnc_[6*i+2],
			nvertex_*irot+mx);
	      }
	    }
	  }
      }
      {
	for (i=0;i<nelm_;++i)
	  {
	    I mx = cnc_[6*i+0];
	    for (j=1;j<3;++j)
	      if (mx>cnc_[6*i+j])
		{
		  mx = cnc_[6*i+j];
		}
	    
	    {
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			nvertex_*irot+cnc_[6*i+2],
			nvertex_*irot+cnc_[6*i+1],
			nvertex_*irot+cnc_[6*i+0],
			(irot-1)*nvertex_+mx);
	      }
	    }
	  }
      }
      {
	for (i=0;i<nelm_;++i)
	  {
	    I j1,j2,j3,mx1 = cnc_[6*i+0];
	    j1=0;
	    j2=0;
	    j3=0;
	  for (j=1;j<3;++j)
	    if (mx1<cnc_[6*i+j])
	      {
		mx1 = cnc_[6*i+j];
		j1=j;
	      }
	  I mx2 = cnc_[6*i+0];
	  for (j=1;j<3;++j)
	    if (mx2>cnc_[6*i+j])
	      {
		mx2 = cnc_[6*i+j];
		j2=j;
	      }
	  if ((j1==0)AND(j2==1))
	    j3 = 2;
	  else if ((j1==1)AND(j2==2))
	    j3 = 0;
	  else if ((j1==0)AND(j2==2))
	    j3 = 1;
	  else if ((j1==1)AND(j2==0))
	    j3 = 2;
	  else if ((j1==2)AND(j2==1))
	    j3 = 0;
	  else if ((j1==2)AND(j2==0))
	    j3 = 1;
	  else 
	    {
	      fprintf(stderr,"coucou\n");
	      exit(1);
	    }
	  {
	    {
	      fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
		      (irot-1)*nvertex_+cnc_[6*i+j3],
		      (irot-1)*nvertex_+mx2,
		      nvertex_*irot+cnc_[6*i+j3],
		      nvertex_*irot+mx1);
	    }
	  }
	  }
      }
    }
  fprintf(out,"End");
  fclose(out);
}




void ns_mesh_extrude_medit	(const ns_mesh*	mesh_,
				 const char * 	name_,
				 const I 	nframe_,
				 const I 	nz_,
				 cst_pR 	z_,
				 const R 	coeff_)
{

  ns_mesh_extrude_scalar	(mesh_,
				 name_,
				 coeff_,
				 nframe_,
				 nz_,
				 z_,
				 NULL);
}




void ns_mesh_extrude_vector	(const ns_mesh*		mesh_,
				 const char * 		name_,
				 const R 		coeff_,
				 const I 		nframe__,
				 const I 		nz_,
				 cst_pR 		z_,
				 cst_pR 		*scalarsx_,
				 cst_pR 		*scalarsy_)
{  
  const I 	nvertex_ = mesh_->nvertex;
  const I 	nelm_ = mesh_->nelm;
  cst_pI	cnc_=mesh_->cnc;
  cst_pI	cod_=mesh_->cod;
  cst_pR	coo_=mesh_->coo;
  I irot;
  I i,j;
  STR obasename;
  cst_eFileSuffix fileSuffix = eFileSuffix_get_suffixAndBasename(name_,
								 obasename);
  char bbname[256];

  I nframe = nframe__;
  for (irot=1;irot<nframe__;++irot)
    {
      if (z_[irot]==z_[0])
	break;
    }
  nframe = irot;

  
  if (scalarsx_)
    { 
      sprintf(bbname,"%s.bb",obasename);
      FILE * out = fopen(bbname,"w");
      fprintf(out,"3 3 "ifmt" 2\n",nframe*nvertex_);  
      for (irot=0;irot<nframe;++irot)
	{
	  cst_pR scalarx_ = scalarsx_[irot];
	  cst_pR scalary_ = scalarsy_[irot];
	  for (i=0;i<nvertex_;++i)
	    {      
	      fprintf(out," "rfmt" "rfmt" "rfmt"\n",((R)0.0),scalarx_[i],scalary_[i]);	  
	    } 
	}
      fclose(out);
    }


  sprintf(bbname,"%s.mesh",obasename);

#if 0
  free(obasename);
#endif
  
  FILE * out 			= fopen(bbname,"w");
  fprintf(out,"MeshVersionFormatted 1\nDimension 3\nVertices "ifmt"\n",nframe*nvertex_);
#if 0
  for (i=0;i<nvertex_;++i)
    {
      fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",((R)0.0),coo_[2*i+0],coo_[2*i+1],cod_[i]);	  
    }  
#endif
  for (irot=0;irot<nframe;++irot)
    {
      for (i=0;i<nvertex_;++i)
	{
	  fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",coo_[2*i+0],(nz_==nframe__)?z_[irot]*coeff_:z_[0]*((R)irot),coo_[2*i+1],cod_[i]);	  
	}  
    }
  fprintf(out,"Tetrahedra\n"ifmt"\n",3*nelm_*(nframe-1));
  for (irot=1;irot<nframe;++irot)
    {
      {
	for (i=0;i<nelm_;++i)
	  {
	    I mx = cnc_[6*i+0];
	    for (j=1;j<3;++j)
	      if (mx<cnc_[6*i+j])
		{
		  mx = cnc_[6*i+j];
		}	    
	    {
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			(irot-1)*nvertex_+cnc_[6*i+0],
			(irot-1)*nvertex_+cnc_[6*i+1],
			(irot-1)*nvertex_+cnc_[6*i+2],
			nvertex_*irot+mx);
	      }
	    }
	  }
      }
      {
	for (i=0;i<nelm_;++i)
	  {
	    I mx = cnc_[6*i+0];
	    for (j=1;j<3;++j)
	      if (mx>cnc_[6*i+j])
		{
		  mx = cnc_[6*i+j];
		}
	    
	    {
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			nvertex_*irot+cnc_[6*i+2],
			nvertex_*irot+cnc_[6*i+1],
			nvertex_*irot+cnc_[6*i+0],
			(irot-1)*nvertex_+mx);
	      }
	    }
	  }
      }
      {
	for (i=0;i<nelm_;++i)
	  {
	    I j1,j2,j3,mx1 = cnc_[6*i+0];
	    j1=0;
	    j2=0;
	    j3=0;
	  for (j=1;j<3;++j)
	    if (mx1<cnc_[6*i+j])
	      {
		mx1 = cnc_[6*i+j];
		j1=j;
	      }
	  I mx2 = cnc_[6*i+0];
	  for (j=1;j<3;++j)
	    if (mx2>cnc_[6*i+j])
	      {
		mx2 = cnc_[6*i+j];
		j2=j;
	      }
	  if ((j1==0)AND(j2==1))
	    j3 = 2;
	  else if ((j1==1)AND(j2==2))
	    j3 = 0;
	  else if ((j1==0)AND(j2==2))
	    j3 = 1;
	  else if ((j1==1)AND(j2==0))
	    j3 = 2;
	  else if ((j1==2)AND(j2==1))
	    j3 = 0;
	  else if ((j1==2)AND(j2==0))
	    j3 = 1;
	  else 
	    {
	      fprintf(stderr,"coucou\n");
	      exit(1);
	    }
	  {
	    {
	      fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
		      (irot-1)*nvertex_+cnc_[6*i+j3],
		      (irot-1)*nvertex_+mx2,
		      nvertex_*irot+cnc_[6*i+j3],
		      nvertex_*irot+mx1);
	    }
	  }
	  }
      }
    }
  fprintf(out,"End");
  fclose(out);

}




void print_mesh_extrude_axi_permutationP2(const ns_mesh * 	mesh_,
					  cst_pI 		axe_,
					  pI 			perm_,
					  pI 			rperm_,
					  pI                   nvertex_axi_,
					  pI                   nvertex_noaxi_)
{  
  R value_axe 	= 0.0;
  nvertex_noaxi_[0] 	= 0;
  nvertex_axi_[0]	= 0;
  const I N = mesh_->nvertex+mesh_->nedge;
  { I i;
    for (i=0;i<N;++i)
      {
	R y = (axe_[0]==1)?mesh_->coo[2*i+1]:mesh_->coo[2*i];
	if (y<value_axe+1.0e-5)
	  {
	    perm_[i]=N-1-nvertex_axi_[0];
	    nvertex_axi_[0]++;
	  }
	else
	  {
	    perm_[i]=nvertex_noaxi_[0];
	    nvertex_noaxi_[0]++;
	  }
      } }
  { I i;
    for (i=0;i<N;++i)
      {
	rperm_[perm_[i]]=i;
      } }
}



void print_mesh_extrude_axi_permutation(const ns_mesh * 	mesh_,
					cst_pI 		axe_,
					pI 			perm_,
					pI 			rperm_,
					pI                     nvertex_axi_,
					pI                     nvertex_noaxi_)
{  
  R value_axe 	= 0.0;
  nvertex_noaxi_[0] 	= 0;
  nvertex_axi_[0]	= 0;
  { I i;
    for (i=0;i<mesh_->nvertex;++i)
      {
	R y = (axe_[0]==1)?mesh_->coo[2*i+1]:mesh_->coo[2*i];
	if (y<value_axe+1.0e-5)
	  {
	    perm_[i]=mesh_->nvertex-1-nvertex_axi_[0];
	    nvertex_axi_[0]++;
	  }
	else
	  {
	    perm_[i]=nvertex_noaxi_[0];
	    nvertex_noaxi_[0]++;
	  }
      } }
  { I i;
    for (i=0;i<mesh_->nvertex;++i)
      {
	rperm_[perm_[i]]=i;
      } }
}



void print_vector_extrude_axi(const ns_mesh * 	mesh_,
			      const char * 		basename_,
			      cst_pI 			nrot_,
			      cst_pI                   axe_,
			      pI                       rperm_,
			      cst_pI 			nvertex_axi_,
			      cst_pI                   nvertex_noaxi_,
			      cst_pR                   vectorx_,
			      cst_pR                   vectory_) 
{
  char ctmp[512];
  sprintf(ctmp,"%s.bb",basename_);
  FILE * out = fopen(ctmp,"w");
  fprintf(out,"3 3 "ifmt" 2\n",nvertex_axi_[0]+nrot_[0]*nvertex_noaxi_[0]);  
  { I i;
    for (i=0;i<nvertex_noaxi_[0];++i)
      {      
	if (axe_[0]==0)
	  {
	    fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*1.0);  		
	    fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*0.0);	  
	    fprintf(out,"" rfmt "\n",vectorx_[rperm_[i]]);	  
	  }
	else
	  {
	    fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*1.0);  		
	    fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*0.0);	  
	    fprintf(out,"" rfmt "\n",vectory_[rperm_[i]]);	  
	  }

      } }
  { I i;
    for (i=nvertex_noaxi_[0];i<mesh_->nvertex;++i)
      {
	if (axe_[0]==0)
	  {
	    fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*1.0);  		
	    fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*0.0);	  
	    fprintf(out,"" rfmt "\n",vectorx_[rperm_[i]]);	  
	  }
	else
	  {
	    fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*1.0);  		
	    fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*0.0);	  
	    fprintf(out,"" rfmt "\n",vectory_[rperm_[i]]);	  
	  }
      } }
  { I irot;
    for (irot=1;irot<nrot_[0];++irot)
      {
	{ I i;
	  for (i=0;i<nvertex_noaxi_[0];++i)
	    {
	      const R teta  = nsACOS(-1.0)*((R)2.0)/((R)nrot_[0])*((R)irot);
	      if (axe_[0]==0)
		{
		  fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*nsCOS(teta));  		
		  fprintf(out,"" rfmt " ",vectory_[rperm_[i]]*nsSIN(teta));	  
		  fprintf(out,"" rfmt "\n",vectorx_[rperm_[i]]);	  
		}
	      else
		{
		  fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*nsCOS(teta));  		
		  fprintf(out,"" rfmt " ",vectorx_[rperm_[i]]*nsSIN(teta));	  
		  fprintf(out,"" rfmt "\n",vectory_[rperm_[i]]);	  
		}
	    } }
      } }
  fclose(out);
}    


void print_scalar_extrude_axi(const ns_mesh * 	mesh_,
			      const char * 		basename_,
			      cst_pI 			nrot_,
			      pI                       rperm_,
			      cst_pI 			nvertex_axi_,
			      cst_pI                   nvertex_noaxi_,
			      cst_pR                   scalar_)
{
  char ctmp[512];
  sprintf(ctmp,"%s.bb",basename_);
  FILE * out = fopen(ctmp,"w");
  fprintf(out,"3 1 "ifmt" 2\n",nvertex_axi_[0]+nrot_[0]*nvertex_noaxi_[0]);  
  { I i;
    for (i=0;i<nvertex_noaxi_[0];++i)
      {      
	fprintf(out," "rfmt"\n",scalar_[rperm_[i]]);	  
      } }
  { I i;
    for (i=nvertex_noaxi_[0];i<mesh_->nvertex;++i)
      {
	fprintf(out," "rfmt"\n",scalar_[rperm_[i]]);	  
      } }
  { I irot;
    for (irot=1;irot<nrot_[0];++irot)
      {
	{ I i;
	  for (i=0;i<nvertex_noaxi_[0];++i)
	    {
	      fprintf(out," "rfmt"\n",scalar_[rperm_[i]]);	  
	    } }
      } }
  fclose(out);
}    

void print_mesh_extrude_axi(const ns_mesh * 	mesh_,
			    const char * 		basename_,
			    cst_pI 			nrot_,
			    cst_pI 			axe_,
			    pI 			perm_,
			    pI                         rperm_,
			    cst_pI 			nvertex_axi_,
			    cst_pI                     nvertex_noaxi_)
{  
  R value_axe = 0.0;
  I i;
  char ctmp[512];
  sprintf(ctmp,"%s.mesh",basename_);
  FILE * out 			= fopen(ctmp,"w");
  const I   nvertex_axi 	= nvertex_axi_[0];
  const I   nvertex_noaxi	= nvertex_noaxi_[0];
  fprintf(out,"MeshVersionFormatted 1\nDimension 3\nVertices "ifmt"\n",nvertex_axi+nrot_[0]*nvertex_noaxi);
  for (i=0;i<nvertex_noaxi;++i)
    {      
      I p = rperm_[i];
      if (axe_[0]==1)
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe),mesh_->coo[2*p],mesh_->cod[p]);	  
	}
      else
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	}
    } 
  for (i=nvertex_noaxi;i<mesh_->nvertex;++i)
    {      
      I p = rperm_[i];
      if (axe_[0]==1)
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe),mesh_->coo[2*p],mesh_->cod[p]);	  
	}
      else
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	}
    }
  I irot=0;
  for (irot=1;irot<nrot_[0];++irot)
    for (i=0;i<nvertex_noaxi;++i)
      {      
	I p 	= rperm_[i];
	R teta 	= nsACOS(-1.0)*((R)2.0)/((R)nrot_[0])*((R)irot);      
	if (axe_[0]==1)
	  {
	    fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe) * nsCOS(teta),(mesh_->coo[2*p+1]-value_axe) * nsSIN(teta),mesh_->coo[2*p],mesh_->cod[p]);	  
	  }
	else
	  {
	    fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe) * nsCOS(teta),(mesh_->coo[2*p]-value_axe) * nsSIN(teta),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	  }
      } 
  I myN = 0;
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		myN+=3;
	      }	    
	  }
      }
    }
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }	    
	    else if ( (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }	    
	    else if ( (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }
	  }
      } 
    }


  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) )
	      {
		myN+=2;
	      }	    
	    else if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) )
	      {
		myN+=2;

	      }	    
	    else if ( (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) )
	      {
		myN+=2;

	      }	    
	  }
      } 
    }
  fprintf(out,"Tetrahedra\n"ifmt"\n",myN);
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}
      {
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx = p1;
		if (mx<p2)
		  {
		    mx=p2; 
		  }      
		if (mx<p3)
		  {
		    mx=p3; 
		  }	      
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx-1]);
	      }	    
	  }
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx = p1;
		if (mx>p2)
		  {
		    mx=p2; 
		  }      
		if (mx>p3)
		  {
		    mx=p3; 
		  }	      
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx-1]);
	      }
	  }
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx1 = p1;
		I j1=0;
		if (mx1<p2)
		  {
		    j1=1;
		    mx1=p2; 
		  }      
		if (mx1<p3)
		  {
		    j1=2;
		    mx1=p3; 
		  }	      
		I mx2 = p1;
		I j2=0;
		if (mx2>p2)
	      {
		j2=1;
		mx2=p2; 
	      }      
	    if (mx2>p3)
	      {
		j2=2;
		mx2=p3; 
	      }	      
	    I j3=0;
	    if ((j1==0)AND(j2==1))
	      j3 = 2;
	    else if ((j1==1)AND(j2==2))
	      j3 = 0;
	    else if ((j1==0)AND(j2==2))
	      j3 = 1;
	    else if ((j1==1)AND(j2==0))
	      j3 = 2;
	    else if ((j1==2)AND(j2==1))
	      j3 = 0;
	    else if ((j1==2)AND(j2==0))
	      j3 = 1;
	    else 
	      {
		fprintf(stderr,"integrator_routine:coucou\n");
		exit(1);
	      }
	    fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
		    dec+(jrot-1)*nvertex_noaxi+1+perm_[mesh_->cnc[6*i+j3]-1],
		    dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1],
		    dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mesh_->cnc[6*i+j3]-1],
		    dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
	  }
	  }
      }
    }
   
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}
      {
	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1]);
	      }	    
	    else if ( (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1]);
	      }	    
	    else if ( (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1]);
	      }
	  }
      } 
    }


  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}


      {	for (i=0;i<mesh_->nelm;++i)
	  {
	    I p1 = mesh_->cnc[6*i+0];
	    I p2 = mesh_->cnc[6*i+1];
	    I p3 = mesh_->cnc[6*i+2];
	    if ( (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) )
	      {
		I mx1 = p2;
		I mx2 = p2;
		if (mx1<p3)
		  mx1=p3;
		if (mx2>p3)
		  mx2=p3;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	    else if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) )
	      {
		I mx1 = p1;
		I mx2 = p1;
		if (mx1<p2)
		  mx1=p2;
		if (mx2>p2)
		  mx2=p2;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	    else if ( (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) )
	      {
		I mx1 = p1;
		I mx2 = p1;
		if (mx1<p3)
		  mx1=p3;
		if (mx2>p3)
		  mx2=p3;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	  }
      } 
    }

    fprintf(out,"End");
    fclose(out);
    
}



void print_mesh_extrude_axiP2(const ns_mesh * 	mesh_,
			      const char * 		basename_,
			    cst_pI 			nrot_,
			    cst_pI 			axe_,
			    pI 			perm_,
			    pI                         rperm_,
			    cst_pI 			nvertex_axi_,
			    cst_pI                     nvertex_noaxi_)
{  
  R value_axe = 0.0;
  I i;
  char ctmp[512];
  sprintf(ctmp,"%s.mesh",basename_);
  const I nelm2 = mesh_->nelm*4;
  pI ddlcnc2 = calloc(nelm2*3,sizeof(I));
  {
    for (i=0;i<mesh_->nelm;++i)
      {
	ddlcnc2[i*12+0]	= mesh_->cnc[i*6+0];
	ddlcnc2[i*12+1]	= mesh_->cnc[i*6+3];
	ddlcnc2[i*12+2]	= mesh_->cnc[i*6+5];

	ddlcnc2[i*12+3]	= mesh_->cnc[i*6+3];
	ddlcnc2[i*12+4]	= mesh_->cnc[i*6+1];
	ddlcnc2[i*12+5]	= mesh_->cnc[i*6+4];

	ddlcnc2[i*12+6]	= mesh_->cnc[i*6+4];
	ddlcnc2[i*12+7]	= mesh_->cnc[i*6+2];
	ddlcnc2[i*12+8]	= mesh_->cnc[i*6+5];

	ddlcnc2[i*12+9]	= mesh_->cnc[i*6+3];
	ddlcnc2[i*12+10]= mesh_->cnc[i*6+4];
	ddlcnc2[i*12+11]= mesh_->cnc[i*6+5];
      }
  }
  FILE * out 			= fopen(ctmp,"w");
  const I   nvertex_axi 	= nvertex_axi_[0];
  const I   nvertex_noaxi	= nvertex_noaxi_[0];
  fprintf(out,"MeshVersionFormatted 1\nDimension 3\nVertices "ifmt"\n",nvertex_axi+nrot_[0]*nvertex_noaxi);
  for (i=0;i<nvertex_noaxi;++i)
    {      
      I p = rperm_[i];
      if (axe_[0]==1)
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe),mesh_->coo[2*p],mesh_->cod[p]);	  
	}
      else
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	}
    } 
  for (i=nvertex_noaxi;i<mesh_->nvertex+mesh_->nedge;++i)
    {      
      I p = rperm_[i];
      if (axe_[0]==1)
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe),mesh_->coo[2*p],mesh_->cod[p]);	  
	}
      else
	{
	  fprintf(out," "rfmt" 0.0 "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	}
    }
  I irot=0;
  for (irot=1;irot<nrot_[0];++irot)
    for (i=0;i<nvertex_noaxi;++i)
      {      
	I p 	= rperm_[i];
	R teta 	= nsACOS(-1.0)*((R)2.0)/((R)nrot_[0])*((R)irot);      
	if (axe_[0]==1)
	  {
	    fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",(mesh_->coo[2*p+1]-value_axe) * nsCOS(teta),(mesh_->coo[2*p+1]-value_axe) * nsSIN(teta),mesh_->coo[2*p],mesh_->cod[p]);	  
	  }
	else
	  {
	    fprintf(out," "rfmt" "rfmt" "rfmt" "ifmt"\n",(mesh_->coo[2*p]-value_axe) * nsCOS(teta),(mesh_->coo[2*p]-value_axe) * nsSIN(teta),mesh_->coo[2*p+1],mesh_->cod[p]);	  
	  }
      } 
  I myN = 0;
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		myN+=3;
	      }	    
	  }
      }
    }
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }	    
	    else if ( (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }	    
	    else if ( (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) )
	      {
		myN+=1;
	      }
	  }
      } 
    }


  for (irot=1;irot<=nrot_[0];++irot)
    {  
      {
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) )
	      {
		myN+=2;
	      }	    
	    else if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) )
	      {
		myN+=2;

	      }	    
	    else if ( (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) )
	      {
		myN+=2;

	      }	    
	  }
      } 
    }
  fprintf(out,"Tetrahedra\n"ifmt"\n",myN);
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nvertex+mesh_->nedge+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nedge+mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nedge+mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nedge+mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}
      {
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx = p1;
		if (mx<p2)
		  {
		    mx=p2; 
		  }      
		if (mx<p3)
		  {
		    mx=p3; 
		  }	      
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx-1]);
	      }	    
	  }
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx = p1;
		if (mx>p2)
		  {
		    mx=p2; 
		  }      
		if (mx>p3)
		  {
		    mx=p3; 
		  }	      
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx-1]);
	      }
	  }
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		I mx1 = p1;
		I j1=0;
		if (mx1<p2)
		  {
		    j1=1;
		    mx1=p2; 
		  }      
		if (mx1<p3)
		  {
		    j1=2;
		    mx1=p3; 
		  }	      
		I mx2 = p1;
		I j2=0;
		if (mx2>p2)
	      {
		j2=1;
		mx2=p2; 
	      }      
	    if (mx2>p3)
	      {
		j2=2;
		mx2=p3; 
	      }	      
	    I j3=0;
	    if ((j1==0)AND(j2==1))
	      j3 = 2;
	    else if ((j1==1)AND(j2==2))
	      j3 = 0;
	    else if ((j1==0)AND(j2==2))
	      j3 = 1;
	    else if ((j1==1)AND(j2==0))
	      j3 = 2;
	    else if ((j1==2)AND(j2==1))
	      j3 = 0;
	    else if ((j1==2)AND(j2==0))
	      j3 = 1;
	    else 
	      {
		fprintf(stderr,"integrator_routine:coucou\n");
		exit(1);
	      }
	    fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
		    dec+(jrot-1)*nvertex_noaxi+1+perm_[ddlcnc2[3*i+j3]-1],
		    dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1],
		    dec2+(jrot2-1)*nvertex_noaxi+1+perm_[ddlcnc2[3*i+j3]-1],
		    dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
	  }
	  }
      }
    }
   
  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nedge+mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nedge+mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nedge+mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nedge+mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}
      {
	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1]);
	      }	    
	    else if ( (perm_[p2-1]>=nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1]);
	      }	    
	    else if ( (perm_[p3-1]>=nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) )
	      {
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1]);
	      }
	  }
      } 
    }


  for (irot=1;irot<=nrot_[0];++irot)
    {  
      I dec,dec2;
      I jrot	= irot;
      I jrot2 	= irot;
      if (irot==nrot_[0])
	{      
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nedge+mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = 0;      
	}
      else if (irot==1)
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = 0;      
	  dec2  = mesh_->nedge+mesh_->nvertex;
	}
      else
	{
	  jrot	= 1;
	  jrot2 = 1;
	  dec   = mesh_->nedge+mesh_->nvertex+(irot-2)*nvertex_noaxi;
	  dec2  = mesh_->nedge+mesh_->nvertex+(irot-1)*nvertex_noaxi;      
	}


      {	for (i=0;i<nelm2;++i)
	  {
	    I p1 = ddlcnc2[3*i+0];
	    I p2 = ddlcnc2[3*i+1];
	    I p3 = ddlcnc2[3*i+2];
	    if ( (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]>=nvertex_noaxi) )
	      {
		I mx1 = p2;
		I mx2 = p2;
		if (mx1<p3)
		  mx1=p3;
		if (mx2>p3)
		  mx2=p3;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	    else if ( (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]<nvertex_noaxi) AND (perm_[p3-1]>=nvertex_noaxi) )
	      {
		I mx1 = p1;
		I mx2 = p1;
		if (mx1<p2)
		  mx1=p2;
		if (mx2>p2)
		  mx2=p2;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p2-1],
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			1+perm_[p3-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	    else if ( (perm_[p3-1]<nvertex_noaxi) AND (perm_[p1-1]<nvertex_noaxi) AND (perm_[p2-1]>=nvertex_noaxi) )
	      {
		I mx1 = p1;
		I mx2 = p1;
		if (mx1<p3)
		  mx1=p3;
		if (mx2>p3)
		  mx2=p3;
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p3-1],
			1+perm_[p2-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[p1-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[mx1-1]);
		fprintf(out,""ifmt" "ifmt" "ifmt" "ifmt" 0\n",
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p1-1],
			1+perm_[p2-1],
			dec2+(jrot2-1)*nvertex_noaxi+1+perm_[p3-1],
			dec+(jrot-1)*nvertex_noaxi+1+perm_[mx2-1]);
	      }	    
	  }
      } 
    }

    fprintf(out,"End");
    fclose(out);
    free(ddlcnc2);
}



#include "ns_sys.h"

#include "ExternABDPY.h"

L elm_contains(cst_pR cooelm,
		       cst_pR coo)
{
  R a,b,c,d;
  nsFLAG i,rrr;
 L contains = __emnsNO;
  for(i=0;i<3;i++)
    {
      a=cooelm[i]-coo[0];
      b=cooelm[3+i]-coo[1];
      c=cooelm[((i+1)%3)]-coo[0];
      d=cooelm[3+((i+1)%3)]-coo[1];
#if 1
      if (nsFABS(a)<1.0e-8)
	a=(R)0.0;
      if (nsFABS(b)<1.0e-8)
	b=(R)0.0;
      if (nsFABS(c)<1.0e-8)
	c=(R)0.0;
      if (nsFABS(d)<1.0e-8)
	d=(R)0.0;
#endif
      rrr = ABDPY_det2x2((double)a,(double)b,(double)c,(double)d);
      if ( (rrr==-1) OR (rrr>1) )
	break;      
    }/* end for */
  if (rrr==0)
    contains = __emnsYES;
  else  if (rrr==1)
    contains = __emnsYES;  
  else  if (rrr==-1)
    contains = __emnsNO;
  else 
    {

      fprintf(stderr,"nsmshelm_contains:predicate (ABDPY_det2x2) returns error code 2 or 3\n");
      contains = __emnsNO;
    }
  return contains;
}




#ifndef __HEADER_INTEGRATOR_ADT_H__
#define __HEADER_INTEGRATOR_ADT_H__

/*#include "mok_list.h"*/

typedef struct __integrator_adt2d_node_st
{
  I		noeud_utilise;
  I         donnee;
  nsOBJ         data;
  R	ptref[4];
  struct __integrator_adt2d_node_st * gauche;
  struct __integrator_adt2d_node_st * droite;
} integrator_adt2d_node_st,*integrator_adt2d_node;

typedef struct
{
  R translation[2];
  R homothetie[2];
  R ptdim2[4];
  R coin_min[4];
  R coin_max[4];
  integrator_adt2d_node racine;
} integrator_adt2d_st,* integrator_adt2d;
void integrator_adt2d_insert( integrator_adt2d,R*,  R*, const I  );
I integrator_adt2d_select_mok_list(integrator_adt2d P,
				R* pta,
				R* ptb,
				nsOBJ list);

integrator_adt2d 	integrator_adt2d_new		(R*,R*);
integrator_adt2d 	integrator_adt2d_kill		(integrator_adt2d);
integrator_adt2d_node 	integrator_adt2d_node_new_data	(R *pt,nsOBJ data_);

I 		integrator_adt2d_select(integrator_adt2d, R*,  R*, I *,I* );
void 		integrator_adt2d_node_insert(integrator_adt2d_node,
				     I, 
				     R*, 
				     R*, 
				     R*,
				     const I);

void integrator_adt2d_init_corner(integrator_adt2d);
 integrator_adt2d_node integrator_adt2d_node_new	(R*,const I );

I integrator_adt2d_select_mpslist(integrator_adt2d P,R* pta,R* ptb,nsOBJ list);
 I integrator_adt2d_node_delete_data	( 
				 integrator_adt2d_node,
				 I,
				 R*,
				 R*,
				 R*,
				 nsOBJ d_
				 );



integrator_adt2d_node integrator_adt2d_node_new_data(R *pt,nsOBJ data_);
  void integrator_adt2d_node_insert_data(integrator_adt2d_node N,I iniv,R*cmin,R*cmax,R*pt,nsOBJ d);
I integrator_adt2d_node_select_mpslist(integrator_adt2d_node N,I iniv,R*cmin,R*cmax,R*pt,nsOBJ list);

void  integrator_adt_localize_recover(void   * 	localize,
			       cst_pR		cooelm_,
			       pI 		select_,
			       pI		select_n_);

#endif

static I __rintegrator_adt_rectangles_chevauchent(R* rect1,R* rect2) 
{
  I idim;
  I ret = (I)YES;
  for( idim = 0; idim < 2; ++idim )
    {
      if( rect2[idim] > rect1[2+idim] ||
          rect1[idim] > rect2[2+idim] )
	{
	  ret = (I)NO;
	  break;
	}
   }
  return ret;
}

#if 0
I integrator_adt2d_node_select_mok_list(integrator_adt2d_node N,I iniv,R*cmin,R*cmax,R*pt,nsOBJ list)
{
  I ok;
  I c = 0;
  I icoord = iniv%4;
  I jcoord = iniv%2;
  R milieu = (cmin[icoord] + cmax[icoord])/2.0;
  
  if( N->noeud_utilise) 
    {
      ok = __rintegrator_adt_rectangles_chevauchent( pt, N->ptref );
      if (ok)
	{
	  mok_list_pushFront(list,N->data);
#if 0
	  select[select_n[0]] = N->donnee;
	  select_n[0]++;
#endif
	  c++;
	}
    }
  
  if( cmin[jcoord] <= pt[jcoord+2] && pt[jcoord] <= cmax[jcoord+2] )
    {
      if( N->gauche != NULL )
	{
	  R sauve_max = cmax[icoord];
	  cmax[icoord] = milieu;
	  c += integrator_adt2d_node_select_mok_list( N->gauche,iniv+1, cmin, cmax, pt, list );
	  cmax[icoord] = sauve_max;
	}
      if( N->droite != NULL )
	{
	  R sauve_min = cmin[icoord];
	  cmin[icoord] = milieu;
	  c += integrator_adt2d_node_select_mok_list( N->droite, iniv+1, cmin, cmax, pt, list );
	  cmin[icoord] = sauve_min;
	}
    }
  return c;
}
#endif






static integrator_adt2d_node integrator_adt2d_node_kill	(integrator_adt2d_node);
static I integrator_adt2d_node_select	(integrator_adt2d_node,
					 I,
					 R*,
					 R*,
					 R*,
					 I *,
					 I *);

static void integrator_adt2d_normalize_vertex( integrator_adt2d, 
				       R*,  
				       R*);





/**
   \brief Construire un ADT
   @param ncoo dimension (ou nombre de composante d'un point)
   @param pt0 point inferieur gauche de la boite englobante
   @param pt1 point superieur droit de la boite englobante
   @return la structure cree
 */
integrator_adt2d  integrator_adt2d_new(R* pt0,R* pt1)
{
  integrator_adt2d P = malloc(sizeof(integrator_adt2d_st));
  I idim;
  P->racine = NULL;
  for( idim = 0; idim < ((I)2); ++idim )
    {
      P->translation[idim] = - pt0[idim];
      if (nsFABS(pt1[idim] - pt0[idim]) > ((R)1.2e-16))
	  P->homothetie[idim] =	  ((R)1.0) /(pt1[idim] - pt0[idim]);
      else
	  P->homothetie[idim] = (R)1.0;
    }
  return P;
}




/**
   \brief Tuer un ADT
   @param P adresse de la structure cree
   @return NULL
 */
integrator_adt2d  integrator_adt2d_kill(integrator_adt2d  P)
{
  if (P)
    {
      P->racine = integrator_adt2d_node_kill(P->racine);
      free(P);
    }
  return NULL;
}

void integrator_adt2d_normalize_vertex(integrator_adt2d P,R * pt0,R * pt1)
{
  P->ptdim2[0] = (pt0[0] + P->translation[0]) * P->homothetie[0];
  P->ptdim2[2] = (pt1[0] + P->translation[0]) * P->homothetie[0];
  P->ptdim2[1] = (pt0[1] + P->translation[1]) * P->homothetie[1];
  P->ptdim2[3] = (pt1[1] + P->translation[1]) * P->homothetie[1];     
}

void integrator_adt2d_init_corner(integrator_adt2d P)
{
  
  P->coin_min[0] = (R)0.0;
  P->coin_max[0] = (R)1.0;
  P->coin_min[1] = (R)0.0;
  P->coin_max[1] = (R)1.0;
  P->coin_min[2] = (R)0.0;
  P->coin_max[2] = (R)1.0;
  P->coin_min[3] = (R)0.0;
  P->coin_max[3] = (R)1.0;
}

void integrator_adt2d_insert(integrator_adt2d P,R*pt0,R*pt1,const I donnee)
{
   integrator_adt2d_normalize_vertex( P, pt0, pt1 );
   if( P->racine == NULL )
      P->racine = integrator_adt2d_node_new( P->ptdim2, donnee );
   else
     {
       integrator_adt2d_init_corner(P);
       integrator_adt2d_node_insert(
			    P->racine,			    
			    (I)0,
			    P->coin_min,
			    P->coin_max,
			    P->ptdim2, 
			    donnee 
			    );
     }
}
#if 0
I integrator_adt2d_select_mok_list(integrator_adt2d P,
				R* pta,
				R* ptb,
				nsOBJ list)
{
  I c;
  (void)integrator_adt2d_normalize_vertex(P, pta, ptb );
  c = (I)0;
  if( P->racine != NULL )
    {
      (void)integrator_adt2d_init_corner(P);
      c = integrator_adt2d_node_select_mok_list( P->racine,(I)0, P->coin_min, P->coin_max, P->ptdim2, list );
    }
  return c;
}
#endif
I integrator_adt2d_select(integrator_adt2d P,R* pta,R* ptb,I * select,I * select_n)
{
  I c;
#if 0
  (void)integrator_adt2d_normalize_vertex(P, pta, ptb );
#else
  R ptdim2[4],coin_min[4],coin_max[4];
  ptdim2[0] = (pta[0] + P->translation[0]) * P->homothetie[0];
  ptdim2[2] = (ptb[0] + P->translation[0]) * P->homothetie[0];
  ptdim2[1] = (pta[1] + P->translation[1]) * P->homothetie[1];
  ptdim2[3] = (ptb[1] + P->translation[1]) * P->homothetie[1];     
#endif

  c = (I)0;
  if( P->racine != NULL )
    {
#if 0
      (void)integrator_adt2d_init_corner(P);
      c = integrator_adt2d_node_select( P->racine,(I)0, P->coin_min, P->coin_max, P->ptdim2, select,select_n );
#else
      coin_min[0] = (R)0.0;
      coin_max[0] = (R)1.0;
      coin_min[1] = (R)0.0;
  coin_max[1] = (R)1.0;
  coin_min[2] = (R)0.0;
  coin_max[2] = (R)1.0;
  coin_min[3] = (R)0.0;
  coin_max[3] = (R)1.0;
      c = integrator_adt2d_node_select( P->racine,(I)0, coin_min, coin_max, ptdim2, select,select_n );
#endif
    }
  return c;
}



integrator_adt2d_node integrator_adt2d_node_new(R *pt,const I  d)
{
  integrator_adt2d_node N = (integrator_adt2d_node)malloc(sizeof(integrator_adt2d_node_st));
  N->noeud_utilise = 1;
  N->donnee = d;
  N->ptref[0] = pt[0];
  N->ptref[1] = pt[1];
  N->ptref[2] = pt[2];
  N->ptref[3] = pt[3];
  N->gauche = NULL;
  N->droite = NULL;
  return N;
}


integrator_adt2d_node integrator_adt2d_node_new_data(R *pt,nsOBJ data_)
{
  integrator_adt2d_node N = (integrator_adt2d_node)malloc(sizeof(integrator_adt2d_node_st));
  N->noeud_utilise = 1;
  N->donnee = (I)0;
  N->data   = data_;
  N->ptref[0] = pt[0];
  N->ptref[1] = pt[1];
  N->ptref[2] = pt[2];
  N->ptref[3] = pt[3];
  N->gauche = NULL;
  N->droite = NULL;
  return N;
}

/**
 *
 *  Liberer la memoire utilisee par un noeud d'un arbre digital alterne et
 *  par ses descendants
 * 
 * 
 */
integrator_adt2d_node integrator_adt2d_node_kill(integrator_adt2d_node N)
{
  if (N)
    {
      N->gauche = integrator_adt2d_node_kill(N->gauche);
      N->droite = integrator_adt2d_node_kill(N->droite);
      free(N);
      N  = NULL;
    }
  return N;
}


void integrator_adt2d_node_insert(integrator_adt2d_node N,
			  I iniv,
			  R*cmin,
			  R*cmax,
			  R*pt,
			  const I d)
{
  if ( N->noeud_utilise )
    {
      I  icoord = iniv%4;
      R milieu = (cmin[icoord] + cmax[icoord])/2.0;      
      if( pt[icoord] <= milieu )
	{
	  if( N->gauche == NULL )
	    N->gauche = integrator_adt2d_node_new(pt, d );
	  else
	    {
	      cmax[icoord] = milieu;
	      integrator_adt2d_node_insert(N->gauche , ++iniv, cmin, cmax, pt, d );
	    }
	}
      else
	{
	  
	  if( N->droite == NULL )
	    N->droite = integrator_adt2d_node_new( pt, d );
	  else
	    {
	      cmin[icoord] = milieu;
	      integrator_adt2d_node_insert( N->droite, ++iniv, cmin, cmax, pt, d );
	    }
	}
    }
  else
    {
      N->noeud_utilise = 1;
      N->donnee = d;
      N->ptref[0] = pt[0];
      N->ptref[1] = pt[1];
      N->ptref[2] = pt[2];
      N->ptref[3] = pt[3];
    }
}





void integrator_adt2d_node_insert_data(integrator_adt2d_node N,I iniv,R*cmin,R*cmax,R*pt,nsOBJ d)
{
  if ( N->noeud_utilise )
    {
      I  icoord = iniv%4;
      R milieu = (cmin[icoord] + cmax[icoord])/2.0;      
      if( pt[icoord] <= milieu )
	{
	  if( N->gauche == NULL )
	    N->gauche = integrator_adt2d_node_new_data(pt, d );
	  else
	    {
	      cmax[icoord] = milieu;
	      integrator_adt2d_node_insert_data(N->gauche , ++iniv, cmin, cmax, pt, d );
	    }
	}
      else
	{
	  
	  if( N->droite == NULL )
	    N->droite = integrator_adt2d_node_new_data( pt, d );
	  else
	    {
	      cmin[icoord] = milieu;
	      integrator_adt2d_node_insert_data( N->droite, ++iniv, cmin, cmax, pt, d );
	    }
	}
    }
  else
    {
      N->noeud_utilise = 1;
      N->donnee = (I)0;
      N->data = d;
      N->ptref[0] = pt[0];
      N->ptref[1] = pt[1];
      N->ptref[2] = pt[2];
      N->ptref[3] = pt[3];
    }
}

I integrator_adt2d_node_select(integrator_adt2d_node N,I iniv,R*cmin,R*cmax,R*pt,I*select,I*select_n )
{
  I ok;
  I c = 0;
  I icoord = iniv%4;
  I jcoord = iniv%2;
  R milieu = (cmin[icoord] + cmax[icoord])/2.0;
  
  if( N->noeud_utilise) 
    {
      ok = __rintegrator_adt_rectangles_chevauchent( pt, N->ptref );
      if (ok)
	{
	  select[select_n[0]] = N->donnee;
	  select_n[0]++;
	  c++;
	}
    }
  
  if( cmin[jcoord] <= pt[jcoord+2] && pt[jcoord] <= cmax[jcoord+2] )
    {
      if( N->gauche != NULL )
	{
	  R sauve_max = cmax[icoord];
	  cmax[icoord] = milieu;
	  c += integrator_adt2d_node_select( N->gauche,iniv+1, cmin, cmax, pt, select,select_n  );
	  cmax[icoord] = sauve_max;
	}
      if( N->droite != NULL )
	{
	  R sauve_min = cmin[icoord];
	  cmin[icoord] = milieu;
	  c += integrator_adt2d_node_select( N->droite, iniv+1, cmin, cmax, pt, select,select_n  );
	  cmin[icoord] = sauve_min;
	}
    }
  return c;
}



I integrator_adt2d_node_delete_data(integrator_adt2d_node N,I iniv,R * cmin,R * cmax,R * pt,nsOBJ d)
{
  I icoord;
  R milieu;

  if( N->noeud_utilise && N->data == d )
    {
      N->data = NULL;
      N->noeud_utilise = 0;
      return 1;
    }
  
  icoord = iniv%4;
  milieu = (cmin[icoord] + cmax[icoord])/2.0;
  
  if( pt[icoord] <= milieu )
    {
      
      cmax[icoord] = milieu;
      return (N->gauche != NULL) ? integrator_adt2d_node_delete_data(N->gauche, ++iniv,
							 cmin, cmax, pt, d ) : 0;
    }
  else
    {
      cmin[icoord] = milieu;
      return (N->droite != NULL) ? integrator_adt2d_node_delete_data(N->droite, ++iniv,
							cmin, cmax, pt, d ) : 0;
    }
}

void * ns_mesh_localization_kill(void * localization_)
{
  if (localization_)
    {
      integrator_adt2d_kill(localization_);
    }
  return NULL;
}


void* ns_mesh_localization(const ns_mesh*mesh_)
{
  R xmin =mesh_->coo[0],xmax=mesh_->coo[0];
  R ymin =mesh_->coo[1],ymax=mesh_->coo[1];
  { I i;
    for (i=1;i<mesh_->nvertex;++i)
      {
	xmin = MIN(xmin,mesh_->coo[2*i]);
	xmax = MAX(xmax,mesh_->coo[2*i]);
	ymin = MIN(ymin,mesh_->coo[2*i+1]);
	ymax = MAX(ymax,mesh_->coo[2*i+1]);
      } }
  R CooCoins[8];
  CooCoins[0]=xmin-0.5;
  CooCoins[1]=ymin-0.5;
  CooCoins[2]=xmax+0.5;
  CooCoins[3]=ymin-0.5;
  CooCoins[4]=xmax+0.5;
  CooCoins[5]=ymax+0.5;
  CooCoins[6]=xmin-0.5;
  CooCoins[7]=ymax+0.5;
  void *localization_data = integrator_adt2d_new(CooCoins,CooCoins+4);
  { I i;
    for (i=0;i<mesh_->nelm;++i)
      {
	R p0[3] = {0.0,0.0,0.0};
	R p1[3] = {0.0,0.0,0.0};
	p0[0] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)];
	p0[1] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)+1];
	p1[0] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)];
	p1[1] = mesh_->coo[2*(mesh_->cnc[6*i+0]-1)+1];
	p0[0] = MIN(p0[0],mesh_->coo[2*(mesh_->cnc[6*i+1]-1)]);
	p0[1] = MIN(p0[1],mesh_->coo[2*(mesh_->cnc[6*i+1]-1)+1]);
	p0[0] = MIN(p0[0],mesh_->coo[2*(mesh_->cnc[6*i+2]-1)]);
	p0[1] = MIN(p0[1],mesh_->coo[2*(mesh_->cnc[6*i+2]-1)+1]);
	p1[0] = MAX(p1[0],mesh_->coo[2*(mesh_->cnc[6*i+1]-1)]);
	p1[1] = MAX(p1[1],mesh_->coo[2*(mesh_->cnc[6*i+1]-1)+1]);
	p1[0] = MAX(p1[0],mesh_->coo[2*(mesh_->cnc[6*i+2]-1)]);
	p1[1] = MAX(p1[1],mesh_->coo[2*(mesh_->cnc[6*i+2]-1)+1]);
	integrator_adt2d_insert(localization_data,p0,p1,i);
      } }
  return localization_data;
}


void  ns_mesh_localize(void   * 	localize,
		       const ns_mesh * 	mesh_,
		       pR 		vert,
		       pI 		select_,
		       pI		select_n_)
{
  nsFLAG ok = (nsFLAG)NO;
  I i,j;
  select_n_[0]=(I)0;  
  R zoom = (R)1.0;
  while (select_n_[0]==0)
    {
      R cc1[2] = {vert[0]-0.005*zoom,vert[1]-0.005*zoom};
      R cc2[2] = {vert[0]+0.005*zoom,vert[1]+0.005*zoom};
      integrator_adt2d_select(localize,&cc1[0],&cc2[0],select_,select_n_);
#ifndef NDEBUG
#if 0
      if (NOT select_n_[0])
	{
	  ns_warn("__rmkZ_findelm_and_project_ifneeded:unzoom localization");
	}
#endif
#if 0
      else
	{
	  printf("mkZ_localize:select_="ifmt"\n",select_[0]);
	}
#endif

#endif
      zoom *=(R)2.0;
    }
  if (select_n_[0]>0)
    {
      for (i=0;i<select_n_[0];++i)
	{
	  R cooelm[6];
	  cooelm[0]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+0]-1)];
	  cooelm[1]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+1]-1)];
	  cooelm[2]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+2]-1)];
	  cooelm[3]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+0]-1)+1];
	  cooelm[4]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+1]-1)+1];
	  cooelm[5]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+2]-1)+1];
#if 0
	  printf("################# search  adt "ifmt" "ifmt"\n%e %e\n",select_n_[0],select_[i],cooelm[0],cooelm[1]);
	  printf("%e %e\n",cooelm[2],cooelm[3]);
	  printf("%e %e\n",cooelm[4],cooelm[5]);
#endif
	  if (elm_contains(cooelm,vert))
	    {

	      ok = (nsFLAG)YES;
	      break;
	    }
	}
      if (i<select_n_[0])
	{
	  j 		= select_[i];
	  select_n_[0] 	= (I)1;
	  select_[0]   	= j;
	}
      else 
	{
#if 0	  
#ifndef NDEBUG
	  ns_warn("mkmshzone_findelm_and_project_ifneeded:perform boundary projection for "rfmt","rfmt"",vert[0],vert[1]);
#endif
#endif
	  R q[3] = {0.0,0.0,0.0};
	  I ind=0;
	  for (i=0;i<select_n_[0];++i)
	    {

#if 0
	  printf("selection  "ifmt"\n",select_n_[0]);
	      printf("elment "ifmt"\n",select_[i]);
#endif
	      R t,p[3],dist2,dist=(R)1.0e+10;
	      R c1[2];
	      R c2[2];
	      R vc1[2];
	      R vc2[2];
#if 0
	      R cooelm[6];
	      cooelm[0]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+0]-1)];
	      cooelm[1]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+1]-1)];
	      cooelm[2]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+2]-1)];
	      cooelm[3]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+0]-1)+1];
	      cooelm[4]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+1]-1)+1];
	      cooelm[5]=mesh_->coo[2*(mesh_->cnc[6*(select_[i])+2]-1)+1];
#endif
	      for (j=0;j<3;++j)
		{
#if 0
		  printf("elment "ifmt" adj "ifmt"\n",select_[i],mesh_->adj[(select_[i])*3+j]);
#endif
		  if (NOT mesh_->adj[(select_[i])*3+j])
		    {
		      c1[0] = mesh_->coo[2*(mesh_->cnc[(select_[i])*6+j]-1)];
		      c1[1] = mesh_->coo[2*(mesh_->cnc[(select_[i])*6+j]-1)+1];
		      c2[0] = mesh_->coo[2*(mesh_->cnc[(select_[i])*6+(j+1)%3]-1)];
		      c2[1] = mesh_->coo[2*(mesh_->cnc[(select_[i])*6+(j+1)%3]-1)+1];
		      /* 
			 on calcule la projection de vert sur l'arete 
		      */
		      vc1[0] = c2[0] - c1[0];
		      vc1[1] = c2[1] - c1[1];
		      vc2[0] = vert[0] - c1[0];
		      vc2[1] = vert[1] - c1[1];      
		      t = (vc1[0]*vc2[0] + vc1[1]*vc2[1]) / (vc1[0]*vc1[0] + vc1[1]*vc1[1]);
		      if (t<((R)0.0))
			t=(R)0.0;      
		      if (t>((R)1.0))
			t=(R)1.0;      
		      p[0] = vc2[0] - t * vc1[0];
		      p[1] = vc2[1] - t * vc1[1];      
		      dist2 = nsSQRT(p[0]*p[0]+p[1]*p[1]);      	  
		      p[0] = c1[0]  + t * vc1[0];
		      p[1] = c1[1]  + t * vc1[1];
		      if (dist2<dist)
			{
			  dist = dist2;
			  q[0] = p[0];
			  q[1] = p[1];
			  ok = YES;
			  ind=i;
			}	  	  		      
		    }
		}	      
	    }


	  if (ok)
	    {
	      vert[0] = q[0];
	      vert[1] = q[1];
	      j = select_[ind];
	      select_n_[0] = (I)1;
	      select_[0]   = j;
	    }
	  else
	    {
	      fprintf(stderr,"no element found "rfmt","rfmt"\n",vert[0],vert[1]);
	      select_n_[0] = (I)0;
	      exit(1);
	    }

	}
    }
  if (NOT ok)
    fprintf(stderr,"__rmpszone_findelm_and_project_ifneeded do not find element\n");
  return;
}/* end function */




void  integrator_adt_localize_recover(void   * 	localize,
			       cst_pR		cooelm_,
			       pI 		select_,
			       pI		select_n_)
{
  R xmin=cooelm_[0];
  R xmax=cooelm_[0];
  R ymin=cooelm_[3];
  R ymax=cooelm_[3];

  ymin = MIN(cooelm_[4],ymin);
  ymin = MIN(cooelm_[5],ymin);
  ymax = MAX(cooelm_[4],ymax);
  ymax = MAX(cooelm_[5],ymax);

  xmin = MIN(cooelm_[1],xmin);
  xmin = MIN(cooelm_[2],xmin);
  xmax = MAX(cooelm_[1],xmax);
  xmax = MAX(cooelm_[2],xmax);

  const R xmil = (xmax-xmin)*((R)0.5);
  const R ymil = (ymax-ymin)*((R)0.5);
  
  select_n_[0] = (I)0;
  R zoom  = (R)1.0;
  while (select_n_[0]==0)
    {
      R cc1[2] = {xmil-(xmil-xmin)*zoom,ymil-(ymil-ymin)*zoom};
      R cc2[2] = {xmil+(xmax-xmil)*zoom,ymil+(ymax-ymil)*zoom};
      integrator_adt2d_select(localize,&cc1[0],&cc2[0],select_,select_n_);
      zoom *=(R)2.0;
    }
  return;
}/* end function */

void ns_mesh_localize_recover(void   * 		localize,
			       cst_pR		cooelm_,
			       pI 		select_,
			       pI		select_n_)
{
  integrator_adt_localize_recover(localize,cooelm_,select_,select_n_);
}
