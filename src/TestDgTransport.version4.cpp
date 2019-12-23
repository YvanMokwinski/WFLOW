
#include <signal.h>
#include <pthread.h>
#include <stdio.h>
#include <math.h>

#include "Cmdline.h"
#include "Monitor.h"
#include "Blas.h"
#include "Type.h"
#include "ns_sys.h"
#include "ns_mesh.h"
#include "mkS.h"
#include "ns_config_lapack.h"
#include "ns_constantes.h"
#include "ExternPardiso.h"
#include "ensBASIS.h"
#include "SmoothedHeaviside.h"

#if 0
typedef struct
{
  I n;
} ns_capturing;

#define __nvp2tria 	6  

enum __ensc { __ensc_maxlevel = 0
		,__ensc_nelms
		,__ensc_ptslvel
		,__ensc_nshapes
		,__ensc_id
		,__ensc_start_hvisu
		,__ensc_start_visutol
		,__ensc_start_eval
		,__ensc_start_vij
		,__ensc_start_fils
		,__ensc_start_idx
		,__ensc_start_reach
		,__ensc_start_visible
		,__ensc_start_iso_visible
		,__ensc_start_level
		,__ensc_len_eval
		,__ensc_len_vij
		,__ensc_len_fils
		,__ensc_len_idx
		,__ensc_len_reach
		,__ensc_len_visible
		,__ensc_len_iso_visible
		,__ensc_len_level
		,__ensc_start_map, 
		__ensc_len_map,
	      __ensc_n};

I  	nsc_map[__ensc_n];
I  	nsc_all4[32];
R 	nsc_cooelm[12];
pR 	nsc_eval 	= NULL;
pR 	nsc_vij		= NULL;
R	nsc_hvisu	= ((R)0.0);
R 	nsc_visutol 	= ((R)0.0);
pI 	nsc_fils	= NULL;
pI 	nsc_idx 	= NULL;
pI 	nsc_iso_visible	= NULL;
pI 	nsc_iso_sign 	= NULL;
pI 	nsc_level	= NULL;
void* 	nsc_ann_data	= NULL;
I  	nsc_countedge	= (I)0;
R 	nsc_iso_value	= ((R)0.0);

#define nsc_set_cooelm(_c) nsc_cooelm[0]=(_c)[0];nsc_cooelm[1]=(_c)[1];nsc_cooelm[2]=(_c)[2];nsc_cooelm[3]=(_c)[3];nsc_cooelm[4]=(_c)[4];nsc_cooelm[5]=(_c)[5]

#define nsc_coo_get(_id,_refcooelm)					\
  _refcooelm[0] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*0+1]-1))*nsc_hvisu; \
  _refcooelm[3] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*0]-1))*nsc_hvisu; \
  									\
  _refcooelm[1] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*1+1]-1))*nsc_hvisu; \
  _refcooelm[4] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*1]-1))*nsc_hvisu; \
  									\
  _refcooelm[2] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*2+1]-1))*nsc_hvisu; \
  _refcooelm[5] 	= ((R)(nsc_idx[__nvp2tria*(_id)+2*2]-1))*nsc_hvisu


#define nsc_sign(_id) (  (nsc_iso_visible[_id]) ? 0 : ( (nsc_iso_sign[_id]==9) ? ((I)1) : ((I)-1) ) )
  



R nsc_calc_jacobian(pR p_)
{
  pR pt1 	= &p_[0];
  pR pt2 	= &p_[1];
  pR pt3 	= &p_[2];
  R d1 	= sqrt( (p_[0]-p_[1])*(p_[0]-p_[1]) + (p_[3]-p_[4])*(p_[3]-p_[4]) );
  R d2 	= sqrt( (p_[1]-p_[2])*(p_[1]-p_[2]) + (p_[4]-p_[5])*(p_[4]-p_[5]) );
  R d3 	= sqrt( (p_[2]-p_[0])*(p_[2]-p_[0]) + (p_[5]-p_[3])*(p_[5]-p_[3]) );
  if (d2<d3)
    {
      R tmp = d2;d2=d3;d3=tmp;
    }
  if (d1<d2)
    {
      R tmp=d1;d1=d2;d2=tmp;
    }
  if (d2<d3)
    {
      R tmp=d2;d2=d3;d3=tmp;
    } 
  return sqrt((d1+d2+d3)*(d3-(d1-d2))*(d3+(d1-d2) )*(d1+(d2-d3) ))*regal1quart;
}


static void nsc_makeij(pI ind_,
				pI res)
{        
  I Pos1[2],Pos2[2],Pos3[2];
  Pos1[0] = (ind_[0] + ind_[2])/2;
  Pos1[1] = (ind_[1] + ind_[3])/2;
  Pos2[0] = (ind_[2] + ind_[4])/2;
  Pos2[1] = (ind_[3] + ind_[5])/2;
  Pos3[0] = (ind_[0] + ind_[4])/2;
  Pos3[1] = (ind_[1] + ind_[5])/2;

  res[0] = ind_[0];
  res[1] = ind_[1];
  res[2] = Pos1[0];
  res[3] = Pos1[1];
  res[4] = Pos3[0];
  res[5] = Pos3[1];
  
  res[6]  = Pos1[0];
  res[7]  = Pos1[1];
  res[8]  = ind_[2];
  res[9]  = ind_[3];
  res[10] = Pos2[0];
  res[11] = Pos2[1];
  
  res[12] = Pos1[0];
  res[13] = Pos1[1];
  res[14] = Pos2[0];
  res[15] = Pos2[1];
  res[16] = Pos3[0];
  res[17] = Pos3[1];

  res[18] = Pos2[0];
  res[19] = Pos2[1];
  res[20] = ind_[4];
  res[21] = ind_[5];
  res[22] = Pos3[0];
  res[23] = Pos3[1];  
}







void nsc_clear()
{
  I stack[64];
  I stack_n 	= 0;
  stack[stack_n++] 	= 0;
 unstack_clear:
  {
    const I id 	= stack[--stack_n];
    nsc_iso_visible[id]	= (I)0;
    nsc_iso_sign   [id]	= (I)0;
    if (nsc_fils[id*4+0])
      {
	stack[stack_n++] = nsc_fils[4*id+3];
	stack[stack_n++] = nsc_fils[4*id+2];
	stack[stack_n++] = nsc_fils[4*id+1];
	stack[stack_n++] = nsc_fils[4*id+0];
      }
  }
  if (stack_n>0)
    goto unstack_clear;    
}


#define nsc_set(_wei)  nsblas_dgemv(transN,&nsc_map[__ensc_len_vij],&nsc_map[__ensc_nshapes],&regal1,&nsc_map[__ensc_nshapes],&nsc_map[__ensc_len_vij],(_wei),&negal1,&regal0,nsc_vij,&negal1)


bool nsc_codelm(const I id)
{
  I cc=0;
  { I j;
    for (j=0;j<3;j++)
      {
	const I myj = (nsc_idx[__nvp2tria*id+2*j+1]-1)*nsc_map[__ensc_ptslvel] - ((nsc_idx[__nvp2tria*id+2*j+1]-1)*(nsc_idx[__nvp2tria*id+2*j+1]-2))/2;
	if (nsc_iso_value>=nsc_vij[myj+nsc_idx[__nvp2tria*id+2*j]-1])
	  ++cc;
      } }
  if (cc<=1)
    {
      return 9;
    }
  else
    {
      return 8;
    }    
}

void nsc_collect_isovisible( const I 	id_,
			     pI c,
			     pI nc,
			     const bool flag)
{
  I stack[64];
  I stack_n = 0;
  
  if (flag==0)
    {
      stack[stack_n++]=0;
    unstack_visible:
      {
	I id = stack[--stack_n];    
	if (nsc_iso_visible[id]!=flag)
	  {
	    if (nsc_fils[4*id])    
	      {
		stack[stack_n++] = nsc_fils[4*id+3];
		stack[stack_n++] = nsc_fils[4*id+2];
		stack[stack_n++] = nsc_fils[4*id+1];
		stack[stack_n++] = nsc_fils[4*id+0];
	      }
	  }
	else c[nc[0]++] = id;
      }
      if (stack_n>0)
	goto unstack_visible;
    }
  else
    {
      stack[stack_n++]=0;
    unstack_visible2:
      {
	I id = stack[--stack_n];    
	if (nsc_iso_visible[id]==flag)
	  {
	    if (nsc_fils[4*id])    
	      {
		stack[stack_n++] = nsc_fils[4*id+3];
		stack[stack_n++] = nsc_fils[4*id+2];
		stack[stack_n++] = nsc_fils[4*id+1];
		stack[stack_n++] = nsc_fils[4*id+0];
	      }
	    else c[nc[0]++] = id;
	  }
      }
      if (stack_n>0)
	goto unstack_visible2;
    }
}


bool nsc_cutelm(R * rtmp,
			   R * vertex)
{  
  I j;
  R r[3],newvertex[4],t,pos[2],pos2[2];
  for (j=0;j<3;++j)
    r[j] = nsc_iso_value-rtmp[j];
    /* cas ou on ne separe pas */
  if ( 
      (rtmp[0]<nsc_iso_value)
      AND
      (rtmp[1]<nsc_iso_value)
      AND
      (rtmp[2]<nsc_iso_value) )
    {
      return 0;
    }
  if ( 
      (rtmp[0]>nsc_iso_value)
      AND
      (rtmp[1]>nsc_iso_value)
      AND
      (rtmp[2]>nsc_iso_value) )
    {
      return 0;
    }
#if 0
  for (j=0;j<3;++j)
    if (rtmp[j]>nsc_iso_value)
      break;
  if (j>=3)
    {
      /* on est completement a l exterieur du fluide */
      return 0;
    }
  for (j=0;j<3;++j)
    if (rtmp[j]<nsc_iso_value)
      break;
  if (j>=3)
    {
      /* on est completement a l interieur du fluide */
      return 0;
    }
#endif

  /* cas ou on separe le triangle en 1, 
     une arete partage l'interface */
  for (j=0;j<3;++j)
    {
      if ( (nsFABS(r[j])<1.0e-14) AND (nsFABS(r[(j+1)%3])<1.0e-14) )
	{
	  if (rtmp[(j+2)%3]>=nsc_iso_value)
	    {
	      /* on est a l interieur */
	      pos[0] 	= vertex[j];
	      pos[1] 	= vertex[3+j];	      
	      pos2[0] 	= vertex[(j+1)%3];
	      pos2[1] 	= vertex[3+(j+1)%3];
	      if (nsc_ann_data) 
		mok_lib_ann_addedge	(nsc_ann_data,pos,pos2,0);	 
	      ++nsc_countedge;
	    }
	  break;
	}
    }
  
  if (j<3)
    {
      return 1;
    }

  /* cas ou on separe le triangle en 2 */
  for (j=0;j<3;++j)
    {

      if (
	  ( 
	   ( ( nsc_iso_value<rtmp[j] )AND (nsc_iso_value>rtmp[(j+1)%3]) )
	   OR
	   ( ( nsc_iso_value>rtmp[j] )AND (nsc_iso_value<rtmp[(j+1)%3]) )
	    )
	  AND
	  (nsFABS(r[(j+2)%3])<1.0e-14) 
	  )
	{
	  /* on subdivise avec le sommet (j+2)%3 */	  
	  /* il faut trouver l intersection dans l arete r[j],r[(j+1)%3] */ 
	  if (nsc_ann_data) 
	    {
	      t = r[j]/(r[j]-r[(j+1)%3]);
	      newvertex[0] = (1.0-t)*vertex[j]   + t*vertex[(j+1)%3];
	      newvertex[1] = (1.0-t)*vertex[3+j] + t*vertex[3+(j+1)%3];		
	      pos[0] = vertex[(j+2)%3];
	      pos[1] = vertex[3+(j+2)%3];	      
	      /* 
		 attention on oriente l arete 
		 on veut que le fluide > 0.5 soit a gauche
	      */	      
	      if (nsc_iso_value>rtmp[j])
		{
		  mok_lib_ann_addedge(nsc_ann_data,pos,newvertex,0);	
		}
	      else
		{
		  mok_lib_ann_addedge(nsc_ann_data,newvertex,pos,0);	  
		}
	    }
	  ++nsc_countedge;	  
	  break;
	}
    }

  if (j<3)
    {
      return 1;
    }
  
  /* cas ou on separe le triangle en 3 */
  for (j=0;j<3;++j)
    {
      if (
	  ( 
	   ( ( nsc_iso_value<rtmp[j] )AND (nsc_iso_value>rtmp[(j+1)%3]) )
	   OR
	   ( ( nsc_iso_value>rtmp[j] )AND (nsc_iso_value<rtmp[(j+1)%3]) )
	    )
	  AND
	  ( 
	   ( ( nsc_iso_value<rtmp[(j+1)%3] )AND (nsc_iso_value>rtmp[(j+2)%3]) )
	   OR
	   ( ( nsc_iso_value>rtmp[(j+1)%3] )AND (nsc_iso_value<rtmp[(j+2)%3]) )
	    )
	  )
	{
	  if (nsc_ann_data) 
	    {
	      t = r[j]/(r[j]-r[(j+1)%3]);
	      newvertex[0] = (((R)1.0)-t)*vertex[j] + t*vertex[(j+1)%3];
	      newvertex[1] = (((R)1.0)-t)*vertex[3+j] + t*vertex[3+(j+1)%3];
	      
	      t = r[(j+1)%3]/(r[(j+1)%3]-r[(j+2)%3]);
	      newvertex[2] 	= (((R)1.0)-t)*vertex[(j+1)%3] + t*vertex[(j+2)%3];
	      newvertex[3] 	= (((R)1.0)-t)*vertex[3+(j+1)%3] + t*vertex[3+(j+2)%3];
	      pos[0] 		= newvertex[0];
	      pos[1] 		= newvertex[1];
	      
	      if (nsc_iso_value<=rtmp[j])
		mok_lib_ann_addedge(nsc_ann_data,pos,&newvertex[2],0);	  
	      else
		mok_lib_ann_addedge(nsc_ann_data,&newvertex[2],pos,0);	  		  
	    }
	  ++nsc_countedge;	  
	  break;
	}
    }
  if (j<3)
    {
      return 1;
    }
  return 0;
}
  

void ns_detect_interface(const I id_,I detect_[1])
{
  I ind[__nvp2tria];
  R rtmp[__nvp2tria];
  I stack[64];
  I stack_n = 0;
  detect_[0] = 0;
  stack[stack_n++]=id_;
 unstack_detect:
  {
    I id = stack[--stack_n];    
    ind[0]	= (nsc_idx[__nvp2tria*id+0]+nsc_idx[__nvp2tria*id+2])/2;
    ind[1]	= (nsc_idx[__nvp2tria*id+1]+nsc_idx[__nvp2tria*id+3])/2;
    ind[2]	= (nsc_idx[__nvp2tria*id+2]+nsc_idx[__nvp2tria*id+4])/2;
    ind[3]	= (nsc_idx[__nvp2tria*id+3]+nsc_idx[__nvp2tria*id+5])/2;
    ind[4]	= (nsc_idx[__nvp2tria*id+4]+nsc_idx[__nvp2tria*id+0])/2;
    ind[5]	= (nsc_idx[__nvp2tria*id+5]+nsc_idx[__nvp2tria*id+1])/2;    
    rtmp[0] 	= nsc_vij[(nsc_idx[__nvp2tria*id+2*0+1]-1)*nsc_map[__ensc_ptslvel] - ( (nsc_idx[__nvp2tria*id+2*0+1]-1)*(nsc_idx[__nvp2tria*id+2*0+1]-2) )/2+nsc_idx[__nvp2tria*id+2*0]-1];
    rtmp[1] 	= nsc_vij[(nsc_idx[__nvp2tria*id+2*1+1]-1)*nsc_map[__ensc_ptslvel] - ( (nsc_idx[__nvp2tria*id+2*1+1]-1)*(nsc_idx[__nvp2tria*id+2*1+1]-2) )/2+nsc_idx[__nvp2tria*id+2*1]-1];
    rtmp[2] 	= nsc_vij[(nsc_idx[__nvp2tria*id+2*2+1]-1)*nsc_map[__ensc_ptslvel] - ( (nsc_idx[__nvp2tria*id+2*2+1]-1)*(nsc_idx[__nvp2tria*id+2*2+1]-2) )/2+nsc_idx[__nvp2tria*id+2*2]-1];    
    rtmp[3] 	= nsc_vij[(ind[2*0+1]-1)*nsc_map[__ensc_ptslvel] - ( (ind[2*0+1]-1)*(ind[2*0+1]-2) ) /2+ind[2*0]-1];
    rtmp[4] 	= nsc_vij[(ind[2*1+1]-1)*nsc_map[__ensc_ptslvel] - ( (ind[2*1+1]-1)*(ind[2*1+1]-2) ) /2+ind[2*1]-1];
    rtmp[5] 	= nsc_vij[(ind[2*2+1]-1)*nsc_map[__ensc_ptslvel] - ( (ind[2*2+1]-1)*(ind[2*2+1]-2) ) /2+ind[2*2]-1];
    if ( ( (rtmp[0]<=nsc_iso_value)
	   OR
	   (rtmp[1]<=nsc_iso_value)
	   OR
	   (rtmp[2]<=nsc_iso_value)
	   OR
	   (rtmp[3]<=nsc_iso_value)
	   OR
	   (rtmp[4]<=nsc_iso_value)
	   OR
	   (rtmp[5]<=nsc_iso_value) )
	 AND
	 ( (rtmp[0]>=nsc_iso_value)
	   OR
	   (rtmp[1]>=nsc_iso_value)
	   OR
	   (rtmp[2]>=nsc_iso_value)
	   OR
	   (rtmp[3]>=nsc_iso_value)
	   OR
	   (rtmp[4]>=nsc_iso_value)
	   OR
	   (rtmp[5]>=nsc_iso_value) ) )
      {
	detect_[0] = 1;
      }
    else
      {
	if (nsc_fils[4*id+0])
	  {
	    stack[stack_n++] = nsc_fils[4*id+3];
	    stack[stack_n++] = nsc_fils[4*id+2];
	    stack[stack_n++] = nsc_fils[4*id+1];
	    stack[stack_n++] = nsc_fils[4*id+0];
	  }
      }
  }    
  if (stack_n>0)
    goto unstack_detect;  
}

  

bool nsc_compute_interface(const I id)
{
  R vertex[6],localcoos[2];
  R rtmp[__nvp2tria];
  /* les nsc_idx */
  { I j;
    for (j=0;j<3;j++)
      {
	const I myj = (nsc_idx[__nvp2tria*id+2*j+1]-1)*nsc_map[__ensc_ptslvel] - ( (nsc_idx[__nvp2tria*id+2*j+1]-1)*(nsc_idx[__nvp2tria*id+2*j+1]-2) )/2;
	rtmp[j] = nsc_vij[myj+nsc_idx[__nvp2tria*id+2*j]-1];
      } }
  { I j;
    for (j=0;j<3;j++)
      {
	localcoos[0] 	= ((R)(nsc_idx[6*id+2*j+1]-1))*nsc_hvisu;
	localcoos[1] 	= ((R)(nsc_idx[6*id+2*j]-1))*nsc_hvisu;
	vertex[j]	= nsc_cooelm[0]*(regal1-localcoos[0]-localcoos[1])+nsc_cooelm[1]*localcoos[0]+nsc_cooelm[2]*localcoos[1];
	vertex[3+j]	= nsc_cooelm[3]*(regal1-localcoos[0]-localcoos[1])+nsc_cooelm[4]*localcoos[0]+nsc_cooelm[5]*localcoos[1];
      } }
  return nsc_cutelm(rtmp,vertex);
}


void nsc_cut( const I 	id,
	      pR 		cutelm,
	      pR 		val,
	      pI 		ncutelm,
	      bool*		inside)
{  
  R r[3],newvertex[4],t;
  I j,myj;
  R vertex[6];
  R rtmp[__nvp2tria];
  R localcoos[2];
  for (j=0;j<3;j++)
    {
      myj 	= (nsc_idx[__nvp2tria*id+2*j+1]-1)*nsc_map[__ensc_ptslvel] - ( (nsc_idx[__nvp2tria*id+2*j+1]-1)*(nsc_idx[__nvp2tria*id+2*j+1]-2) )/2;
      rtmp[j] 	= nsc_vij[myj+nsc_idx[__nvp2tria*id+2*j]-1];
    }

  for (j=0;j<3;j++)
    {
      localcoos[0] 	= ((R)(nsc_idx[6*id+2*j+1]-1))*nsc_hvisu;
      localcoos[1] 	= ((R)(nsc_idx[6*id+2*j]-1))*nsc_hvisu;
      vertex[j]		= cooelm[0]*(regal1-localcoos[0]-localcoos[1])+cooelm[1]*localcoos[0]+cooelm[2]*localcoos[1];
      vertex[3+j]	= cooelm[3]*(regal1-localcoos[0]-localcoos[1])+cooelm[4]*localcoos[0]+cooelm[5]*localcoos[1];
      vertex[j]		= localcoos[0];
      vertex[3+j]	= localcoos[1];
    }      

  r[0] = nsc_iso_value-rtmp[0];
  r[1] = nsc_iso_value-rtmp[1];
  r[2] = nsc_iso_value-rtmp[2];
  for (j=0;j<3;++j)
    if (nsc_iso_value>rtmp[j])
      break;
  if (j>=3)
    {
      ncutelm[0] 	= 1;
      cutelm[0] 	= vertex[0];
      cutelm[1] 	= vertex[1];
      cutelm[2] 	= vertex[2];
      cutelm[3] 	= vertex[3];
      cutelm[4] 	= vertex[4];
      cutelm[5] 	= vertex[5];
      inside[0] 	= 1;
      return;
    }
  for (j=0;j<3;++j)
    if (nsc_iso_value<=rtmp[j])
      break;
  if (j>=3)
    {
#if 0
      mk_err("nmon2");
      vertex[j]		= cooelm[0]*(regal1-localcoos[0]-localcoos[1])+cooelm[1]*localcoos[0]+cooelm[2]*localcoos[1];
      vertex[3+j]	= cooelm[3]*(regal1-localcoos[0]-localcoos[1])+cooelm[4]*localcoos[0]+cooelm[5]*localcoos[1];
#endif
      ncutelm[0] = 1;
      cutelm[0] 	= vertex[0];
      cutelm[1] 	= vertex[1];
      cutelm[2] 	= vertex[2];
      cutelm[3] 	= vertex[3];
      cutelm[4] 	= vertex[4];
      cutelm[5] 	= vertex[5];
      inside[0] = 0;	      
      return;
    }
  
  /* cas ou on separe le triangle en 1, une arete partage l'interface */
#if 1
  for (j=0;j<3;++j)
    {
      if ( (nsFABS(r[j])<1.0e-15) AND (nsFABS(r[(j+1)%3])<1.0e-15) )
	{
	  ncutelm[0] = 1;
	  cutelm[0] 	= vertex[0];
	  cutelm[1] 	= vertex[1];
	  cutelm[2] 	= vertex[2];
	  cutelm[3] 	= vertex[3];
	  cutelm[4] 	= vertex[4];
	  cutelm[5] 	= vertex[5];
	  if (r[(j+2)%3]<=regal0)
	    {
	      inside[0] = 1;
	      /* on est a l interieur */
	    }  
	  else
	    {
	      inside[0] = 0;
	    }
	  val[0] = rtmp[0];
	  val[1] = rtmp[1];	
	  val[2] = rtmp[2];
	  break;
	}
    }
  
  if (j<3)
    {
      return;
    }
#endif

  /* cas ou on separe le triangle en 2 */
  for (j=0;j<3;++j)
    {
      if ( ((r[j]*r[(j+1)%3])<regal0) AND (nsFABS(r[(j+2)%3])<1.0e-15) )
	{
	  /* on subdivise avec le sommet (j+2)%3 */
	  
	  /* il faut trouver l intersection dans l arete r[j],r[(j+1)%3] */ 
	  t = r[j]/(r[j]-r[(j+1)%3]);
	  newvertex[0] = (1.0-t)*vertex[j] + t*vertex[(j+1)%3];
	  newvertex[1] = (1.0-t)*vertex[3+j] + t*vertex[3+(j+1)%3];
	  
	  cutelm[0]    = vertex[j];
	  cutelm[1]    = newvertex[0];
	  cutelm[2]    = vertex[(j+2)%3];
	  cutelm[3]    = vertex[3+j];
	  cutelm[4]    = newvertex[1];
	  cutelm[5]    = vertex[3+(j+2)%3];
	  
	  cutelm[6+0]    = newvertex[0];
	  cutelm[6+1]    = vertex[(j+1)%3];
	  cutelm[6+2]    = vertex[(j+2)%3];
	  cutelm[6+3]    = newvertex[1];
	  cutelm[6+4]    = vertex[3+(j+1)%3];
	  cutelm[6+5]    = vertex[3+(j+2)%3];
	  
	  val[0] = rtmp[j];
	  val[1] = nsc_iso_value;	
	  val[2] = rtmp[(j+2)%3];
	  
	  val[3] = nsc_iso_value;
	  val[4] = rtmp[(j+1)%3];	
	  val[5] = rtmp[(j+2)%3];
	  
	  if (r[j]<=regal0)
	    {
	      inside[0] = 1;
	      inside[1] = 0;
	      /* le triangle j new (j+2%3) est a l interieur */	      
	      
	    }
	  else
	    {
	      inside[0] = 0;
	      inside[1] = 1;
	      /* le triangle new (j+1)%3 (j+2)%3 est a l interieur */
	    }
	  ncutelm[0] = 2;
	  break;
	}
    }
  if (j<3)
    {
      return;
    }
  /* cas ou on separe le triangle en 3 */
  for (j=0;j<3;++j)
    {
      if ( ((r[j]*r[(j+1)%3])<regal0) AND ((r[(j+1)%3]*r[(j+2)%3])<regal0) )
	{
	  t = r[j]/(r[j]-r[(j+1)%3]);
	  newvertex[0] = (regal1-t)*vertex[j] + t*vertex[(j+1)%3];
	  newvertex[1] = (regal1-t)*vertex[3+j] + t*vertex[3+(j+1)%3];
	  	  
	  t 			= r[(j+1)%3]/(r[(j+1)%3]-r[(j+2)%3]);
	  newvertex[2] 	= (regal1-t)*vertex[(j+1)%3] + t*vertex[(j+2)%3];
	  newvertex[3] 	= (regal1-t)*vertex[3+(j+1)%3] + t*vertex[3+(j+2)%3];
	  
	  cutelm[0]    	= newvertex[0];
	  cutelm[1]    	= vertex[(j+1)%3];
	  cutelm[2]    	= newvertex[2];
	  cutelm[3]    	= newvertex[1];
	  cutelm[4]    	= vertex[3+(j+1)%3];
	  cutelm[5]    	= newvertex[3];
	  
	  cutelm[6]    	= vertex[j];
	  cutelm[7]    	= newvertex[0];
	  cutelm[8]    	= vertex[(j+2)%3];
	  cutelm[9]    	= vertex[3+j];
	  cutelm[10]    	= newvertex[1];
	  cutelm[11]    	= vertex[3+(j+2)%3];
	  	  
	  cutelm[12]    	= newvertex[0];
	  cutelm[13]    	= newvertex[2];
	  cutelm[14]    	= vertex[(j+2)%3];
	  cutelm[15]    	= newvertex[1];
	  cutelm[16]    	= newvertex[3];
	  cutelm[17]    	= vertex[3+(j+2)%3];
	    
	  val[0] = nsc_iso_value;
	  val[1] = rtmp[(j+1)%3];
	  val[2] = nsc_iso_value;	
	  
	  val[3] = rtmp[j];
	  val[4] = nsc_iso_value;	
	  val[5] = rtmp[(j+2)%3];

	  val[6] = nsc_iso_value;
	  val[7] = nsc_iso_value;	
	  val[8] = rtmp[(j+2)%3];

	  if (r[(j+1)%3]<=regal0)
	    {
	      /* le triangle est a l interieur, le quadrangle est a l exterieur */
	      inside[0] = 1;
	      inside[1] = 0;
	      inside[2] = 0;
	    }
	  else
	    {
	      inside[0] = 0;
	      inside[1] = 1;
	      inside[2] = 1;
	      /* le triangle est a l exterieur, le quadrangle est a l interieur */
	    }
	  ncutelm[0] = 3;	  
	  break;
	}
    }
}



void nsc_buildtree()
{    
  const I nbv = (I)6;
  I stack[64];
  I stack_n 	= 0;
  stack[stack_n++] 	= 0;
 unstack_build:
  {
    const I id = stack[--stack_n];
    if (nsc_level[id]>=nsc_map[__ensc_maxlevel])
      {
	nsc_fils[id*4+0] = (I)0;
	nsc_fils[id*4+1] = (I)0;
	nsc_fils[id*4+2] = (I)0;
	nsc_fils[id*4+3] = (I)0;
      }  
    else
      {
	nsc_fils[id*4+0] = ++nsc_map[__ensc_id];	    
	nsc_fils[id*4+1] = ++nsc_map[__ensc_id];	    
	nsc_fils[id*4+2] = ++nsc_map[__ensc_id];	    
	nsc_fils[id*4+3] = ++nsc_map[__ensc_id];	    
	nsc_makeij(&nsc_idx[6*id],nsc_all4);	    
	nsc_level[nsc_fils[id*4+0]] = nsc_level[id]+1;      
	nsc_level[nsc_fils[id*4+1]] = nsc_level[id]+1;      
	nsc_level[nsc_fils[id*4+2]] = nsc_level[id]+1;      
	nsc_level[nsc_fils[id*4+3]] = nsc_level[id]+1;      

	nsc_idx[nsc_fils[4*id+0]*nbv+0] = nsc_all4[0*nbv+0];
	nsc_idx[nsc_fils[4*id+0]*nbv+1] = nsc_all4[0*nbv+1];
	nsc_idx[nsc_fils[4*id+0]*nbv+2] = nsc_all4[0*nbv+2];
	nsc_idx[nsc_fils[4*id+0]*nbv+3] = nsc_all4[0*nbv+3];
	nsc_idx[nsc_fils[4*id+0]*nbv+4] = nsc_all4[0*nbv+4];
	nsc_idx[nsc_fils[4*id+0]*nbv+5] = nsc_all4[0*nbv+5];
	    
	nsc_idx[nsc_fils[4*id+1]*nbv+0] = nsc_all4[1*nbv+0];
	nsc_idx[nsc_fils[4*id+1]*nbv+1] = nsc_all4[1*nbv+1];
	nsc_idx[nsc_fils[4*id+1]*nbv+2] = nsc_all4[1*nbv+2];
	nsc_idx[nsc_fils[4*id+1]*nbv+3] = nsc_all4[1*nbv+3];
	nsc_idx[nsc_fils[4*id+1]*nbv+4] = nsc_all4[1*nbv+4];
	nsc_idx[nsc_fils[4*id+1]*nbv+5] = nsc_all4[1*nbv+5];
	    
	nsc_idx[nsc_fils[4*id+2]*nbv+0] = nsc_all4[2*nbv+0];
	nsc_idx[nsc_fils[4*id+2]*nbv+1] = nsc_all4[2*nbv+1];
	nsc_idx[nsc_fils[4*id+2]*nbv+2] = nsc_all4[2*nbv+2];
	nsc_idx[nsc_fils[4*id+2]*nbv+3] = nsc_all4[2*nbv+3];
	nsc_idx[nsc_fils[4*id+2]*nbv+4] = nsc_all4[2*nbv+4];
	nsc_idx[nsc_fils[4*id+2]*nbv+5] = nsc_all4[2*nbv+5];
	    
	nsc_idx[nsc_fils[4*id+3]*nbv+0] = nsc_all4[3*nbv+0];
	nsc_idx[nsc_fils[4*id+3]*nbv+1] = nsc_all4[3*nbv+1];
	nsc_idx[nsc_fils[4*id+3]*nbv+2] = nsc_all4[3*nbv+2];
	nsc_idx[nsc_fils[4*id+3]*nbv+3] = nsc_all4[3*nbv+3];
	nsc_idx[nsc_fils[4*id+3]*nbv+4] = nsc_all4[3*nbv+4];
	nsc_idx[nsc_fils[4*id+3]*nbv+5] = nsc_all4[3*nbv+5];

	stack[stack_n++] 	= nsc_fils[4*id+3];
	stack[stack_n++] 	= nsc_fils[4*id+2];
	stack[stack_n++] 	= nsc_fils[4*id+1];
	stack[stack_n++] 	= nsc_fils[4*id+0];
      }
  }
  if (stack_n>0)
    goto unstack_build;
  return;
}





void nsc_analyze()
{
  I stack[64];
  I stack_n 	= 0;
  stack[stack_n++] 	= 0;
 unstack:
    {
      const I id = stack[--stack_n];
      if (nsc_fils[4*id+0])
	{
	  ns_detect_interface(id,&nsc_iso_visible[id]);
	  if (nsc_iso_visible[id])
	    {
	      stack[stack_n++] 	= nsc_fils[4*id+3];
	      stack[stack_n++] 	= nsc_fils[4*id+2];
	      stack[stack_n++] 	= nsc_fils[4*id+1];
	      stack[stack_n++] 	= nsc_fils[4*id+0];
	    }
	  else
	    {	
	      nsc_iso_sign[id] 			= nsc_codelm(id);	  
	      nsc_iso_sign[nsc_fils[4*id+0]] 	= nsc_codelm(nsc_fils[4*id+0]);
	      nsc_iso_sign[nsc_fils[4*id+1]] 	= nsc_codelm(nsc_fils[4*id+1]);
	      nsc_iso_sign[nsc_fils[4*id+2]] 	= nsc_codelm(nsc_fils[4*id+2]);
	      nsc_iso_sign[nsc_fils[4*id+3]] 	= nsc_codelm(nsc_fils[4*id+3]);
	    }
	}
      else
	{
	  nsc_iso_visible[id] 	= nsc_compute_interface(id);
	  nsc_iso_sign[id] 	= nsc_codelm(id);
	}
    }
    if  (stack_n>0)
      goto unstack;    
}
  
  





void nsc_quadrature(ns_capturing 	hd_,
		    cst_pI   		q_n_,
		    cst_pR  		q_lc2gl_,
		    cst_pI   		q_lc2gloff_,
		    cst_pR  		q_w_,
		    cst_pI   		q_mxn_,
		    cst_pR  		cooelm_,
		    cst_pR  		F_,				  
		    pI               	qi_n_,
		    pR               	qi_w_,
		    pR               	qi_p_,
		    pI               	qo_n_,
		    pR               	qo_w_,
		    pR               	qo_p_,
		    I*		err_)
{

  nsc_set_cooelm(cooelm_);    

  nsc_clear();

  nsc_set(F_);

  nsc_analyze();

  pR q_p0_ = qi_p_;
  pR q_p1_ = qo_p_;
  pR q_w0_ = qi_w_;
  pR q_w1_ = qo_w_;
  
  bool 	inside[64];
  I 	i,j,k,negal3=(I)3,negal1=(I)1,negal2=(I)2;
  I 	collect1[2048],collect2[2048];
  I 	i_ncutelm = (I)0,o_ncutelm=(I)0,ncollect1 = 0,ncollect2 = 0;
  R  	cutelm[2048],val[2048];
  R x;

  nsc_collect_isovisible((I)0,collect1,&ncollect1,0);
  nsc_collect_isovisible((I)0,collect2,&ncollect2,1);  

  for (j=0;j<ncollect1;++j)
    {
      if (nsc_sign(collect1[j])<0)
	o_ncutelm++;
      else
	i_ncutelm++;
    }

  
  for (i=0;i<ncollect2;++i)
    {
      k=0;
      nsc_cut(collect2[i],cutelm,val,&k,inside);
      for (j=0;j<k;++j)
	{
	  if (inside[j]>0)
	    {
	      i_ncutelm++;
	    }
	  else
	    {
	      o_ncutelm++;
	    }
	}
    }
  
  qi_n_[0] = (i_ncutelm>0) ? q_n_[0]*i_ncutelm : 0;
  qo_n_[0] = (o_ncutelm>0) ? q_n_[0]*o_ncutelm : 0;
  
  if (q_mxn_[0]<qi_n_[0])
    {
      err_[0] = 1;
      fprintf(stderr,"integrator_capturing:qmxn_ = " ifmt " < qi_n_ = " ifmt "\n",q_mxn_[0],qi_n_[0]);
      return;
    }
  if (q_mxn_[0]<qo_n_[0])
    {
      err_[0] = 1;
      fprintf(stderr,"integrator_capturing:qmxn_ = " ifmt " < qo_n_ = " ifmt "\n",q_mxn_[0],qo_n_[0]);
      return;
    }
  
  const I o_off_coo 	= o_ncutelm*q_n_[0];
  const I i_off_coo 	= i_ncutelm*q_n_[0];
  I  i_ind 	= 0;
  I  o_ind 	= 0;
  for (j=0;j<ncollect1;++j)
    {
      nsc_coo_get(collect1[j],cutelm);
      if (nsc_sign(collect1[j])<0)
	{	
	  nsblas_dgemm(__vmps_blas_transT,
		       __vmps_blas_transN,
		       q_n_,
		       &negal2,
		       &negal3,
		       &__vmps_blas_regal1,
		       q_lc2gl_,
		       q_lc2gloff_,
		       cutelm,
		       &negal3,
		       &__vmps_blas_regal0,
		       &q_p1_[o_ind*q_n_[0]],
		       &o_off_coo);
	  nsblas_dcopy(q_n_,q_w_,&negal1,&q_w1_[o_ind*q_n_[0]],&negal1);
	  x = nsPOW(2.0,-2.0*((R)nsc_level[collect1[j]]));
	  nsblas_dscal(q_n_,&x,&q_w1_[o_ind*q_n_[0]],&negal1);	      
	  ++o_ind;
	}
      else
	{
	  nsblas_dgemm(__vmps_blas_transT,
		       __vmps_blas_transN,
		       q_n_,
		       &negal2,
		       &negal3,
		       &__vmps_blas_regal1,
		       q_lc2gl_,
		       q_lc2gloff_,
		       cutelm,
		       &negal3,
		       &__vmps_blas_regal0,
		       &q_p0_[i_ind*q_n_[0]],
		       &i_off_coo);
	  nsblas_dcopy(q_n_,
		       q_w_,
		       &negal1,
		       &q_w0_[i_ind*q_n_[0]],
		       &negal1);
	  x = nsPOW(((R)2.0),-((R)2.0)*((R)nsc_level[collect1[j]]));
	  nsblas_dscal(q_n_,
		       &x,
		       &q_w0_[i_ind*q_n_[0]],
		       &negal1);	      
	  ++i_ind;
	}
    }
  for (j=0;j<ncollect2;++j)
    {
      k=0;
      nsc_cut(collect2[j],cutelm,val,&k,inside);
      for (i=0;i<k;++i)
	{
	  /* on calcule le jacobien du triangle
	   */
	  if (inside[i]>0)
	    {
	      if (i_ncutelm>0)
		{
		  nsblas_dgemm(__vmps_blas_transT,
			       __vmps_blas_transN,
			       q_n_,
			       &negal2,
			       &negal3,
			       &__vmps_blas_regal1,
			       q_lc2gl_,
			       q_lc2gloff_,
			       &cutelm[6*i],
			       &negal3,
			       &__vmps_blas_regal0,
			       &q_p0_[i_ind*q_n_[0]],
			       &i_off_coo);
		  x = nsc_calc_jacobian(&cutelm[6*i]);
		  x*=2.0;
		  nsblas_dcopy(q_n_,q_w_,&negal1,&q_w0_[i_ind*q_n_[0]],&negal1);	      
		  nsblas_dscal(q_n_,&x,&q_w0_[i_ind*q_n_[0]],&negal1);
		  ++i_ind;
		}
	    }
	  else
	    {
	      if (o_ncutelm)
		{
		  nsblas_dgemm(__vmps_blas_transT,
			       __vmps_blas_transN,
			       q_n_,&negal2,&negal3,&__vmps_blas_regal1,
			       q_lc2gl_,
			       q_lc2gloff_,
			       &cutelm[6*i],
			       &negal3,
			       &__vmps_blas_regal0,
			       &q_p1_[o_ind*q_n_[0]],
			       &o_off_coo);
		  nsblas_dcopy(q_n_,q_w_,&negal1,&q_w1_[o_ind*q_n_[0]],&negal1);
		  x = nsc_calc_jacobian(&cutelm[6*i]);
		  x*=((R)2.0);
		  nsblas_dscal(q_n_,&x,&q_w1_[o_ind*q_n_[0]],&negal1);
		  ++o_ind;
		}
	    }
	}
    }

}


void nsc_def(cst_pI nf_,
	     cst_pR iso_,
	     cst_pR tol_,
	     cst_pI lev_,
	     I*err_)
{
  nsc_countedge=0;

  I i,j,n,N,nvisuptslevel;
  err_[0] 			= 0;

  nsc_map[__ensc_id] 			= negal0;     
  nsc_map[__ensc_nshapes] 		= nf_[0];
  nsc_map[__ensc_maxlevel] 		= lev_[0];
  nsc_map[__ensc_nelms]    		= ((I)nsPOW(regal4,(R)lev_[0]+1)-1)/negal3;
  nsc_map[__ensc_ptslvel]  		= (I)nsPOW(regal2,(R)lev_[0])+1;
  
  nsc_map[__ensc_len_idx] 		= __nvp2tria * nsc_map[__ensc_nelms];
  nsc_map[__ensc_len_fils] 		= 4 * nsc_map[__ensc_nelms];
  nsc_map[__ensc_len_visible] 		= nsc_map[__ensc_nelms];
  nsc_map[__ensc_len_level] 		= nsc_map[__ensc_nelms];
  nsc_map[__ensc_len_vij] 		= (nsc_map[__ensc_ptslvel]*(nsc_map[__ensc_ptslvel]+1))/2;
  nsc_map[__ensc_len_eval] 		= nf_[0]*nsc_map[__ensc_len_vij];
  nsc_map[__ensc_len_iso_visible] 	= nsc_map[__ensc_nelms];
  
  nsc_idx 				= imalloc(nsc_map[__ensc_len_idx]);      
  nsc_iso_visible 			= imalloc(nsc_map[__ensc_len_iso_visible]);      
  nsc_fils 				= imalloc(nsc_map[__ensc_len_fils]);                     
  nsc_iso_sign 				= imalloc(nsc_map[__ensc_len_iso_visible]);      
  nsc_level				= imalloc(nsc_map[__ensc_len_level]);      
  nsc_vij      				= rmalloc(nsc_map[__ensc_len_vij]);
  nsc_eval     				= rmalloc(nsc_map[__ensc_len_eval]);      

  pR localcoos 			= rmalloc(nsc_map[__ensc_ptslvel]*dimension);

  clr_ivect(nsc_map[__ensc_len_level],nsc_level);

  nsc_idx[0] = (I)1;
  nsc_idx[1] = (I)1;
  nsc_idx[2] = nsc_map[__ensc_ptslvel];
  nsc_idx[3] = (I)1;
  nsc_idx[4] = (I)1;
  nsc_idx[5] = nsc_map[__ensc_ptslvel];
  nsc_hvisu  = nsPOW(regal2,mregal1*((R)nsc_map[__ensc_maxlevel]));            

  nvisuptslevel = nsc_map[__ensc_ptslvel];      
  for (j=0;j<nvisuptslevel;j++)
    localcoos[j]  = regal0;      
  N = negal0;
  for (i=negal1;i<=nvisuptslevel;i++)
    {
      n = nvisuptslevel+1-i;
      for (j=negal0;j<n;j++)	
	localcoos[nvisuptslevel+j] = ((R)j)*nsc_hvisu;

      ns_basis_triangle_l	(shape_f,
				 &n,
				 localcoos,
				 &nvisuptslevel,
				 &nsc_eval[N],
				 &nsc_map[__ensc_len_vij]);

      N += n;	      
      for (j=0;j<n-1;j++)	
	localcoos[j] = ((R)i)*nsc_hvisu;
    }
  nsc_buildtree();
  nsc_iso_value		= iso_[0];
  nsc_visutol 		= tol_[0];
  nsc_ann_data   	= NULL;

  rfree(localcoos);
}


void nsc_free()
{
  ifree(nsc_idx);
  ifree(nsc_iso_visible);   
  ifree(nsc_fils);
  ifree(nsc_iso_sign);
  ifree(nsc_level);
  rfree(nsc_vij);
  rfree(nsc_eval);
}



void nsc_info()
{
  printf("*** Info high definition algorithm\n");
  printf("*** nsc_iso_value " rfmt "\n",nsc_iso_value);
  printf("*** tol       " rfmt "\n",nsc_visutol);
  printf("*** max level " ifmt "\n",nsc_map[__ensc_maxlevel]);
  printf("*** nelm      " ifmt "\n",nsc_map[__ensc_nelms]);
  printf("*** nedge     " ifmt "\n",nsc_countedge);
  printf("*** len       " ifmt "\n",nsc_map[__ensc_len_vij]);
}
#endif

#if 0
struct DG_ZONE
{
  ns_mesh * m_mesh;
};

struct DG_HADAPT
{
  DG_HADAPT * m_previous;
  ns_mesh *   m_mesh;
  DG_HADAPT * m_next;
};

struct DG_VAR
{
  DG_ZONE * 	m_zone;
  mkS_st	m_shape;
  pR 		m_values;
  matrix_handle m_matrix;
  DG_VAR(cst_mkS 	shape,
	 I 		nelm)
  {
    this->m_values = (pR)malloc(sizeof(I)*nelm*mkS_n(shape));
    matrix_handle_def(&this->m_matrix,
		      mkS_n(shape),
		      nelm,
		      m_values,
		      mkS_n(shape));
  };

};

struct DG_VAR_SPECTRAL
{
  mkS_st	m_shape;
  DG_VAR *      m_spectral;

  void up(I ielm)
  {
    ++m_qshape_k[ielm];
  }

  bool remove(I k,I ielm)
  {
    
  }
  
  bool add(I k,I ielm)
  {
    if (m_hash_ielm)
      {
	
      }
  }

  
  I find(I level,I ielm)
  {
    if (level ==0)
      {
	return ielm;
      }
    else
      {
	for (I i=this->m_hash_link[level][ielm%79];i!=-777;i=this->m_hash_ielm[level][2*i+1])
	  {
	    I k = this->m_hash_ielm[level][2*(-i-1)+0];
	    if (k==ielm)
	      {
		return i;
	      }
	  }	
      }
    return -1;
  };

  
  DG_VAR_PADAPT(cst_mkS 	shape,
		I 		nelm)
  {
    this->m_qshape_n[0] = 1;
    this->m_qshape_values[0] = (pR)malloc(sizeof(R)*nelm);
    matrix_handle_def(&this->m_matrix,
		      mkS_n(shape),
		      nelm,
		      m_values,
		      mkS_n(shape));
    
    I degree = mkS_k(shape);
    for (I i=1;i<=degree;++i)
      {
	this->m_qshape_n[i] = i+1;
	this->m_qshape_values[i] = (pR)malloc(sizeof(R)*(i+1)*nelm);	
      }
  };

  I ndofs(const I ielm)
  {
    return 1;
  };

  void dofs(const I 		ielm_,
	    vector_handle& 	dofs_)
  {
    
    I at = 0;
    dofs_.x[at] = m_qshape_values[0][ielm_];
    for (I i=1;i<=m_qshape_k[ielm_];++i)
      {
	I jelm = find(i,ielm_);
	if (jelm + 1 > 0)
	  {
	    for (I j=0;j<this->m_qshape_n[i];++j)
	      {
		dofs_.x[at++] = this->m_qshape_values[i][this->m_qshape_n[i]*jelm+j];
	      }
	  }
      }    
    return 1;
  };
  

};
#endif




ns_mesh * mesh;

struct vector_handle;
struct matrix_handle;

struct temp_gemm
{
  const R a;
  const char transA;
  const char transB;
  const matrix_handle&A;
  const matrix_handle&B;
};

struct temp_gemv
{
  const R a;
  const char trans;
  const matrix_handle&A;
  const vector_handle&b;
};
struct temp_scal;
struct temp_transpose
{
  const char trans;
  const matrix_handle&A;
  inline temp_scal operator * (const R&v_) const;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
};

struct temp_scal
{
  const R a;
  const char trans;
  const matrix_handle&A;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
  inline temp_gemm operator * (const temp_transpose&v_) const;
};

inline temp_gemm temp_scal::operator * (const temp_transpose&v_) const
{
  return {a,trans,v_.trans,A,v_.A};
}

inline temp_gemv temp_scal::operator * (const vector_handle&v_) const
{
  return {a,trans,A, v_};
};
inline temp_gemm temp_scal::operator * (const matrix_handle&v_) const
{
  return {a, trans,'N',A, v_};
};

inline temp_scal temp_transpose::operator * (const R&v_) const
{
  return {v_,trans,A};
};
inline temp_gemv temp_transpose::operator * (const vector_handle&v_) const
{
  return {1.0,trans,A, v_};
};
  
inline temp_gemm temp_transpose::operator * (const matrix_handle&v_) const
{
  return {1.0, trans,'N', A, v_};
};

struct matrix_handle
{
  pR  x{};
  I   n{};
  I   m{};
  I   ld{};
  matrix_handle();
  matrix_handle(I n_, I m_,pR  x_,I ld_);
  inline temp_gemm operator * (const matrix_handle&v_)const;  
  inline temp_gemv operator * (const vector_handle&v_)const;  
  inline temp_scal operator * (const R&v_)const;  
  inline temp_transpose transpose() const
  {
    return {'T',*this};
  };


  inline matrix_handle& operator=(const temp_gemm&temp)
  {
    const R r0 = 0.0;
    I k = (temp.transA == 'N') ? temp.A.m : temp.A.n;
    Blas_dgemm(&temp.transA,
	       &temp.transB,
	       &n,
	       &m,
	       &k,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.B.x,
	       &temp.B.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };

};


temp_scal operator * (const R&v_,const matrix_handle &m) 
{
  return {v_,'N',m};
};

matrix_handle::matrix_handle(){};
matrix_handle::matrix_handle(I n_, I m_,pR  x_,I ld_) : x(x_), n(n_), m(m_), ld(ld_)
  {
  };

  inline temp_gemv matrix_handle::operator * (const vector_handle&v_) const
  {
    return {1.0,'N',*this, v_};
  };

inline temp_gemm matrix_handle::operator * (const matrix_handle&v_) const
  {
    return {1.0,'N','N',*this, v_};
  };

inline temp_scal matrix_handle::operator * (const R&v_) const
  {
    return {v_,'N',*this};
  };
  

struct vector_handle
{
  pR  x{};
  I   n{};
  I   ld{};
  vector_handle(){};
  vector_handle(I n_, pR  x_,I ld_)
    : x(x_),n(n_),ld(ld_)
  {};

  template <typename F>
  inline void apply(F f)
  {
    for (I j=0;j<n;++j)
      {
	pR e = x + j*ld;
	*e = f(*e);
      }    
  };
  
  inline vector_handle& operator=(const temp_gemv&temp)
  {
    const R r0 = 0.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };
  
  inline vector_handle& operator+=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };

  inline vector_handle& operator-=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    R a = -temp.a;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };


};

void vector_handle_def(struct vector_handle * h,I n_, pR  x_,I ld_)
{
  h->n  = n_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_def(struct matrix_handle * h,I n_, I m_,pR  x_,I ld_)
{
  h->n  = n_;
  h->m  = m_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_gemv_low(const struct matrix_handle * h,const char * trans,cst_pR a,cst_pR x,cst_pI xoff,cst_pR ry,pR y,cst_pI yoff)
{

  Blas_dgemv(trans,
	     &h->n,
	     &h->m,
	     a,
	     h->x,
	     &h->ld,
	     x,
	     xoff,
	     ry,
	     y,
	     yoff);

}

void matrix_handle_gemv(const struct matrix_handle * h,const char * trans,cst_pR a,const struct vector_handle * x,cst_pR ry,struct vector_handle * y)
{
  matrix_handle_gemv_low(h,trans,a,x->x,&x->ld,ry,y->x,&y->ld);
}


void matrix_handle_gemm(const struct matrix_handle * h,const char * transA,const char * transB,cst_pR a,const struct matrix_handle * x,cst_pR ry,struct matrix_handle * y)
{
  const I s = transA[0]=='N' ? h->m : h->n;
  Blas_dgemm(transA,transB,&y->n,&y->m,&s,a,h->x,&h->ld,x->x,&x->ld,ry,y->x,&y->ld);

}

void vector_handle_print(const struct vector_handle * h,FILE * f)
{
  for (I i=0;i<h->n;++i)
    {
      fprintf(f," " rfmt,h->x[h->ld*i]);	  
      fprintf(f,"\n");
    }

}

void matrix_handle_print(const struct matrix_handle * h,FILE * f)
{
  for (I i=0;i<h->n;++i)
    {
      for (I j=0;j<h->m;++j)
	{
	  fprintf(f," " rfmt,h->x[h->ld*j+i]);	  
	}
      fprintf(f,"\n");
    }

}

void matrix_handle_gesv(const struct matrix_handle * h,struct vector_handle * rhs,pI lcperm)
{
  I n1=1;
  I info_lapack;
  dgesv(&h->n,
	&n1,
	h->x,
	&h->ld,
	lcperm,
	rhs->x,
	&h->n,
	&info_lapack);
}




#define DGERR_MEMORY  2
#define DGERR_USER    3

#define DG_r_lc        			0
#define DG_r_n         			1

#define DG_ires_err 			0 
#define DG_ires_convergence 		1
#define DG_ires_iter_gauss_seidel 	2
#define DG_ires_required_iw_n 		3
#define DG_ires_required_rw_n 		4
#define DG_ires_n 			5

#define DG_rres_max 			0 
#define DG_rres_nrmL2 			1 
#define DG_rres_nrmLInf			2 
#define DG_rres_areaL1 			3 
#define DG_rres_jumpL2 			4 
#define DG_rres_johnson			5 
#define DG_rres_n 			6


  template <typename impl_t,
	    eTopology 	_faceShape,
	    unsigned int 	_degree>
  class treilli2d_t;
  
  template <eTopology _faceShape,unsigned int _degree> struct treilli2d_traits_t
  {
  };
  template <unsigned int _degree> struct treilli2d_traits_t<__eTopology_TRIANGLE,_degree>
  {
    static constexpr const unsigned int nnodes 		= ( (_degree+1)*(_degree+2))/2;
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim	= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 3;
  };

    
  template <eTopology _faceShape,unsigned int _degree> struct treilli2d_utils
  {	
  };
  
  template <unsigned int _degree> struct treilli2d_utils<__eTopology_TRIANGLE,_degree>
  {
    using traits_t = treilli2d_traits_t<__eTopology_TRIANGLE,_degree>;
    
    using faces_to_nodes_t 		= std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
    using nodes_int_t 			= std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
    using nodes_real_t 	= std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
    
    template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& 	w_,														nodes_real_t& 		p_) noexcept
    {
      
      static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t zero(0.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  unsigned int startedge0 = 3;
	  unsigned int startedge1 = 3 + numPointsOnEdge + numPointsOnEdge-1;  
	  unsigned int startedge2 = 3 + numPointsOnEdge*2 + numPointsOnEdge-1;
	  unsigned int startInterior = 3*numPointsOnEdge+3;
	  
	  p_[0][0] 	= zero;
	  p_[0][1] 	= zero;
	  p_[1][0] 	= one;
	  p_[1][1] 	= zero;
	  p_[2][0] 	= zero;
	  p_[2][1] 	= one;

	  //
	  // Third edge
	  //
	  for (unsigned int j=1;j<_degree;++j)
	    {
	      const auto wj 	= w_[j];
	      const auto wk 	= w_[_degree-j];
	      p_[startedge2][0] 	= zero;
	      p_[startedge2][1] 	= (three + two * wj - wk) / six;
	      //	  std::cout << p_[startedge2][1] << std::endl;
	      //	  std::cout << (three + two * wj - wk) / six << std::endl;
	      --startedge2;
	    }
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi 	= w_[i];
	      const auto wk 	= w_[_degree - i];
	      p_[startedge0][0] 	= (three + two * wi - wk) / six;
	      p_[startedge0][1] 	= zero;
	      ++startedge0;
	    }
	  
	  startedge0=3;
	  // 
	  // Second edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[startedge1][0] 	= p_[startedge0++][0]; 
	      p_[startedge1][1] 	= p_[++startedge2][1];
	      --startedge1;
	    }
	  //
	  // Interior
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi = w_[i];
	      for (unsigned int j=1;j<_degree-i;++j)
		{
		  const auto wj 	= w_[j];
		  const auto wk 	= w_[_degree-i-j];
		  p_[startInterior][0] 	= ( two* (one + wi) - (wj + wk) ) / six;
		  p_[startInterior][1] 	= ( two* (one + wj) - (wi + wk) ) / six;
		  startInterior++;
		}
	    }
	  
	};

	
	static inline void ComputeSubcnc(faces_to_nodes_t & subcnc_)
	{
	  unsigned int subCellIndex = 0;
#define _dec(_i,_j) (( (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	  for (unsigned int j=0;j<_degree;j++)
	    {
	    
	      for (unsigned int i=0;i<j;i++)
		{
		 
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-i);
		  subcnc_[subCellIndex][2] = _dec(i,j+1-i); 
		  ++subCellIndex;
		
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-1-i);
		  subcnc_[subCellIndex][2] = _dec(i+1,j-i);
		  ++subCellIndex;

		}
	    
	      subcnc_[subCellIndex][0] = _dec(j,0);
	      subcnc_[subCellIndex][1] = _dec(j+1,0);
	      subcnc_[subCellIndex][2] = _dec(j,1);
	      ++subCellIndex;
	    }
#undef _dec  
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  // COMPUTE GRID     
	  //	6 
	  //	3 7
	  //	1 4 8 
	  //	0 2 5 9
	  unsigned int nodeIndex = 0;
	  for (unsigned int i=0;i<_degree+1;i++)
	    {
	      for (unsigned int j=0;j<=i;j++)
		{
		  icoo_[nodeIndex][0] = j;
		  icoo_[nodeIndex][1] = i-j;
		  ++nodeIndex;
		}
	    }
	};

	
	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	  // VERTEX 0 
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;
	
	  // VERTEX 1 
	  icoo_[nodeIndex][0] = _degree;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;

	  // VERTEX 2
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = _degree;
	  ++nodeIndex;
	    
	  // EDGE 0 
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      icoo_[nodeIndex][0] = i+1;
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	  // EDGE 1 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree-(i+1); 
	      icoo_[nodeIndex][1] = i+1;
	      ++nodeIndex;
	    }
	  // EDGE 2 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree-(i+1);
	      ++nodeIndex;
	    }
	  // INTERIOR
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      for (unsigned int j=0;j<_degree-i-2;++j)
		{
		  icoo_[nodeIndex][0] = i + 1;
		  icoo_[nodeIndex][1] = j + 1;
		  ++nodeIndex;
		}
	    }
	};
	  
      };


#if 0
  template <unsigned int _degree> struct treilli2d_utils<FaceType::Quadrilateral,_degree>
      {
	using traits_t = treilli2d_traits_t<FaceType::Quadrilateral,_degree>;		
	using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

	template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& w_,													    nodes_real_t& 	p_) noexcept
	{
	  static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t mone(-1.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  p_[0][0] 	= mone;
	  p_[0][1] 	= mone;
	  p_[1][0] 	= one;
	  p_[1][1] 	= mone;
	  p_[2][0] 	= one;
	  p_[2][1] 	= one;
	  p_[3][0] 	= mone;
	  p_[3][1] 	= one;

	  unsigned int pointIndex = 4;
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[i];
	      p_[pointIndex][1] 	= mone;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= one;
	      p_[pointIndex][1] 	= w_[i];
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[_degree - i];
	      p_[pointIndex][1] 	= one;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= mone;
	      p_[pointIndex][1] 	= w_[_degree - i];
	      ++pointIndex;
	    }
	  
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  p_[pointIndex][0] = w_[i];
		  p_[pointIndex][1] = w_[j];	
		  ++pointIndex;
		} 
	    }
	};

	static inline void ComputeSubcnc(faces_to_nodes_t& subcnc_)
	{
	  unsigned int subCellIndex = 0;
	  // COMPUTE CNC
#define _dec(_i,_j)   (_degree+1) * (_i) + (_j) 
	  { 
	    for (unsigned int i=0;i<_degree;i++)
	      {
		for (unsigned int j=0;j<_degree;j++)
		  {		    
		    subcnc_[subCellIndex][0] = _dec( (i+1), (j+1) );
		    subcnc_[subCellIndex][1] = _dec( (i), (j+1) );
		    subcnc_[subCellIndex][2] = _dec( (i), (j) );
		    subcnc_[subCellIndex][3] = _dec( (i+1), (j) );
		    ++subCellIndex;
		  } 
	      }
	  }
#undef _dec
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  for (unsigned int i=0;i<_degree+1;++i)
	    {
	      for (unsigned int j=0;j<_degree+1;j++)
		{
		  icoo_[i*(_degree+1)+j][0] = i;	  
		  icoo_[i*(_degree+1)+j][1] = j;	
		}
	    }
	};


	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
	
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  /* EDGE 0 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = i;	  
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	
	  /* EDGE 1 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = _degree;	  
	      icoo_[nodeIndex][1] = i;
	      ++nodeIndex;
	    } 
	
	  /* EDGE 2 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree - i - 1;	  
	      icoo_[nodeIndex][1] = _degree;	
	      ++nodeIndex;
	    }

	  /* EDGE 3 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree - i - 1;	
	      ++nodeIndex;
	    }
	
	  /* INSIDE */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  icoo_[nodeIndex][0] = i;
		  icoo_[nodeIndex][1] = j;	
		  ++nodeIndex;
		} 
	    }
	  
	};

  };
#endif
      template <typename impl_t, eTopology _faceShape,unsigned int _degree>
      class treilli2d_t // : public CRTP< treilli2d_t<impl_t,_faceShape,_degree> > 
      {
      private: using traits_t 	= treilli2d_traits_t<_faceShape,_degree>;
	
      private: using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	
      protected: static constexpr const unsigned int s_numNodes 	= traits_t::nnodes;
      protected: static constexpr const unsigned int s_numElements 	= traits_t::nsubcells;
      protected: static constexpr const unsigned int s_dimension 	= traits_t::dim;
      protected: static constexpr const unsigned int s_degree 		= traits_t::degree;
      protected: static constexpr const unsigned int s_numNodesInFace	= traits_t::nnodesInFace;
	
      public: inline unsigned int 	degree() 		const noexcept { return s_degree; };
      public: inline unsigned int 	dimension() 		const noexcept { return s_dimension; };
      public: inline unsigned int 	nnodes() 		const noexcept { return s_numNodes; };    
      public: inline unsigned int 	nsubcells() 		const noexcept { return s_numElements; };
      public: inline eTopology	celltype()		const noexcept { return _faceShape; };
      public: inline unsigned int 	nnodesincell() 		const noexcept { return s_numNodesInFace; };    
	
      protected: faces_to_nodes_t m_subcnc;

      public: template <typename _float_type>
      inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				       const unsigned int dimensionIndex_) const noexcept
	{
	  return static_cast<const impl_t&>(*this).GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
	};
	
      public: inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
					       const unsigned int&localNodeIndex_) const noexcept
	{
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,subElementIndex_,(unsigned int)0,this->s_numElements-1);
//	  Debug::IsInRange(__TRACE__,localNodeIndex_,(unsigned int)0,this->s_numNodesInFace-1);
//#endif
	  return this->m_subcnc[subElementIndex_][localNodeIndex_];
	};
	
      protected: inline treilli2d_t(nodes_int_t& icoo) noexcept
	{
	  //
	  // Compute integer coordinates.
	  //
	  treilli2d_utils<_faceShape,_degree>::ComputeCoordinates(icoo);
	  std::array<unsigned int, (_degree+1)*(_degree+1)> perm;	  
	  nodes_int_t  ilagr;
	  {
	  //
	  // Compute integer coordinates.
	  //
	    treilli2d_utils<_faceShape,_degree>::ComputeCoordinates2(ilagr);	  
	    for (unsigned int i=0;i<s_numNodes;++i)
	      {	
		perm[ icoo[i][0] * (_degree+1) + icoo[i][1] ] = i+1;
	      } 
	  
	    treilli2d_utils<_faceShape,_degree>::ComputeSubcnc(this->m_subcnc);
	  
	    for (unsigned int subElementIndex=0;subElementIndex<s_numElements;++subElementIndex)
	      {
		for (unsigned int iv=0;iv<s_numNodesInFace;++iv)
		  {	  
		    const unsigned int l 				= m_subcnc[subElementIndex][iv];
		    const unsigned int i 				= ilagr[l][0];
		    const unsigned int j 				= ilagr[l][1];
		    m_subcnc[subElementIndex][iv] = perm[(_degree+1)*i+j]-1;
		  }
	      }
	  }
	
	};
    
	inline ~treilli2d_t() noexcept
	{
	};
    
      };
    

      template <eTopology _faceShape,unsigned int _degree>
      class Uniform : public treilli2d_t<Uniform<_faceShape,_degree> ,_faceShape,_degree>
      {
      
      private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
      private: nodes_int_t m_icoo;

      public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									       const unsigned int dimensionIndex_) const noexcept
	{
	  static constexpr const _float_type idegree = _float_type(1.0) / _float_type(_degree);
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
//	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
//#endif      
	  return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
	};
            
      public: inline Uniform() noexcept : treilli2d_t<Uniform<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
	{  	
	};
	
	inline ~Uniform() noexcept
	{
	};
    
      };


  template <eTopology _faceShape,unsigned int _degree>
  class Generator : public treilli2d_t<Generator<_faceShape,_degree> ,_faceShape,_degree>
  {
    
  private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
  protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
  protected: using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

  private: nodes_int_t m_icoo;
  private: nodes_real_t m_rcoo;
    
    //!
    //! @brief Get the coordinates.
    //!
  public: template <typename _float_type>
  inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				   const unsigned int dimensionIndex_) const noexcept
    {
      //#ifndef NDEBUG
      //	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
      //	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
      //#endif      
      return this->m_rcoo[nodeIndex_][dimensionIndex_];
    };
    
    
  public: inline Generator(const double * p) noexcept : treilli2d_t<Generator<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
    {
      //      Quadrature::Edge::Legendre<double,_degree-1> l;
      std::array<double,_degree+1> w;
      w[0]=-1.0;
      for (int i=0;i<_degree-1;++i)
	{
	  w[1+i] = p[i]; // l.GetPosition(i,0);
	}	  
      w[_degree]=1.0;
      treilli2d_utils<_faceShape,_degree>::LobattoPoints(w, 
						       this->m_rcoo);
      
    };
    
    inline ~Generator() noexcept
    {
    };
  };

  





cst_mkS 		mkS_derivative	(mkS shape_,const I idim)
{
  if (idim==0)
    {
      return mkS_dx(shape_);
    }
  else if (idim==1)
    {
      return mkS_dy(shape_);
    }
  else if (idim==2)
    {
      return mkS_dz(shape_);
    }
  else
    {
      fprintf(stderr,"mkS_derivative error idim = " ifmt "\n",idim);
      exit(1);
    }
  return NULL;
}


void dg_print_sol(mkS shape,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }

#if 1
    mkS_st shapeP;
    I degree = mkS_k(shape);
    Err err;
    mkS_definit	(&shapeP,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 degree,
		 __emk_discontinuous,
		 &err);

    double rst[64];
    double beval[64];
    
    R rwork[1024*2];
    I rwork_n = 1024*2;

    I np = mkS_n(&shapeP);
    I ee;
    mkS_lagrange_localspl_tria(&degree,
			       rst,
			       &n);
    
    mkS_basis(mkS_b(shape),
	      &np,	      
	      beval,
	      &n,
	      rst,
	      &n,
	      rwork,
	      &rwork_n,
	      &ee);
    
#if 0
    for (I i=0;i<n;++i)
      {
	for (I j=0;j<n;++j)
	  {
	    std::cout << " " <<beval[j*n+i] ;
	  }
	std::cout    << std::endl;

      }
#endif
#endif
    
    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*n);
    { I j;
      for (j=0;j<N;++j)
	{
	      double r1=1.0;
	      double r0=0.0;
	      I n1=1;
	      Blas_dgemv("T",
			 &n,
			 &n,
			 &r1,
			 beval,
			 &n,
			 &sol[n*j],
			 &n1,
			 &r0,
			 tmp,
			 &n1);
	  for (I i=0;i<n;++i)
	    {
	      fprintf(fil,"" rfmt "\n",tmp[i]);
	    }
	} } 
    fclose(fil); }
}

void dg_print_mod(mkS shape,I kmod,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }

#if 1
    mkS_st shapeP;
    I degree = mkS_k(shape);
    Err err;
    mkS_definit	(&shapeP,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 degree,
		 __emk_discontinuous,
		 &err);

    double rst[64];
    double beval[64];
    
    R rwork[1024*2];
    I rwork_n = 1024*2;

    I np = mkS_n(&shapeP);
    I ee;
    mkS_lagrange_localspl_tria(&degree,
			       rst,
			       &n);
    double tt[64];
    mkS_basis(mkS_b(shape),
	      &np,	      
	      beval,
	      &n,
	      rst,
	      &n,
	      rwork,
	      &rwork_n,
	      &ee);
    
#if 0
    for (I i=0;i<n;++i)
      {
	for (I j=0;j<n;++j)
	  {
	    std::cout << " " <<beval[j*n+i] ;
	  }
	std::cout    << std::endl;

      }
#endif
#endif
    
    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*n);
    { I j;
      for (j=0;j<N;++j)
	{
	      double r1=1.0;
	      double r0=0.0;
	      I n1=1;
	      for (I i=0;i<n;++i)
		{
		  tt[i] = sol[n*j+i];
		}

	      // 0 [0 0[ [1,n[
	      // 1 [0 1[ [3,n[
	      // 2 [0 3[ [6,n[
	      // 3 [0 6[ [10,n[
	      I low = kmod==0 ? 0 : kmod-1;	      
	      I start = 0;
	      I bound = (kmod==0) ? 0 : ((low+1)*(low+2))/2;
	      for (I i=start;i<bound;++i)
		{
		  tt[i] = 0.0;
		}
	      start = ((kmod+1)*(kmod+2))/2;

	      for (I i=start;i<n;++i)
		{
		  tt[i] = 0.0;
		}


	      
	      Blas_dgemv("T",
			 &n,
			 &n,
			 &r1,
			 beval,
			 &n,
			 tt,
			 &n1,
			 &r0,
			 tmp,
			 &n1);
		   
	  for (I i=0;i<n;++i)
	    {
	      fprintf(fil,"" rfmt "\n",tmp[i]);
	    }
	} } 
    fclose(fil); }
}

void dg_print_mesh(const ns_mesh*s_,const char * name_,...)
{
  //  const I numNodes = ns_mesh_get_numNodes(s_);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	} } 

    fprintf(fil,"Triangles\n" ifmt "\n",s_->nelm); 
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,((I)0));
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void dg_print_mesh(const ns_mesh*s_,cst_pI codelm,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	} } 

    fprintf(fil,"Triangles\n" ifmt "\n",s_->nelm); 
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void dg_print_mesh(I numNodes,I nelm,pI cnc,I cncoff,pR coo,pI codnodes,I codnodesoff,pI codelm,I codelmoff,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",nelm*3);
    { I i;
      for (i=0;i<nelm;++i)
	{
	  for (I k=0;k<3;++k)
	    {
	      cncelm[k] = cnc[cncoff*i+k];
	    }

	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",coo[cncelm[j]*2+0],coo[cncelm[j]*2+1],codnodes[codnodesoff * cncelm[j]]); } }
	} } 
    
    fprintf(fil,"Triangles\n" ifmt "\n",nelm); 
    { I i;
      for (i=0;i<nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i*codelmoff]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}



struct DG_VIEW
{
  static constexpr I dim = 2;
  //  Generator<__eTopology_TRIANGLE,_degree> *generator;
  Uniform<__eTopology_TRIANGLE,5> *generator;


  I 		m_nsubcells;
  R 		m_mem[4096];

  matrix_handle m_rst;
  matrix_handle m_beval_xyz;
  matrix_handle m_beval_f;
  matrix_handle m_xyz;
  vector_handle m_f;

  R p[32];

  DG_VIEW(mkS shape_x,mkS shape_f)
  {
    //Uniform<__eTopology_TRIANGLE,4> treilli;(p);
    p[0]=-2.0/3.0;
    p[1]=-1.0/3.0;
    p[2]=0;
    p[3]=1.0/3.0;
    p[4]=2.0/3.0;

    p[0]=-0.9061798459386639927976268;
    p[1]=    -0.5384693101056830910363;
    p[2]=0;
    p[3]=0.5384693101056830910363144207;
    p[4]=0.9061798459386639927976268;



p[0] =     -0.973906528517171720077964;
p[1] =     -0.8650633666889845;
p[2] =     -0.679409568299024;
p[3] =     -0.4333953941292471907992659;
p[4] =     -0.14887433898163121088482600112971;
p[5] =     0.1488743389816312108848260011;
p[6] =     0.433395394129247190799;
p[7] =     0.67940956829902440;
p[8] =     0.865063366688984;
p[9] =     0.9739065285171717200779640;



#if 0
    p[0] = -0.774596669241483377035853079956;
    p[1] = 0;
    p[2] = 0.774596669241483377035853079956;
    p[0] = -0.5;
    p[1] = 0;
    p[2] = 0.6;
#endif
    
    // generator = new Generator<__eTopology_TRIANGLE,_degree> (p);
    
    generator 	= new Uniform<__eTopology_TRIANGLE,5>();
    I nx 	= mkS_n(shape_x);
    I n 	= mkS_n(shape_f);
    I nsubcells = generator->nsubcells();    
    I at     	= 0;
    I nnodes 	= generator->nnodes();
    m_nsubcells = nsubcells;
    
    matrix_handle_def(&m_rst,nnodes,dim,&m_mem[at],nnodes);
    at += nnodes*dim;
    matrix_handle_def(&m_beval_xyz,nx,nnodes,&m_mem[at],nx);
    at += nnodes*nx;
    matrix_handle_def(&m_beval_f,n,nnodes,&m_mem[at],n);
    at += n*nnodes;
    matrix_handle_def(&m_xyz,nnodes,dim,&m_mem[at],nnodes);
    at += nnodes*dim;
    vector_handle_def(&m_f,nnodes,&m_mem[at],1);
    at += nnodes;

    for (I i=0;i<nnodes;++i)
      {
	m_rst.x[i] = generator->GetCoordinate<double>(i,0);
      }
    for (I i=0;i<nnodes;++i)
      {
	m_rst.x[m_rst.ld+i] = generator->GetCoordinate<double>(i,1);
      }

    I err;
    R rwork[4096];
    I rwork_n = 4096;
    
    mkS_basis(mkS_b(shape_x),
	      &nnodes,	      
	      this->m_beval_xyz.x,
	      &this->m_beval_xyz.ld,
	      this->m_rst.x,
	      &this->m_rst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);

    mkS_basis(mkS_b(shape_f),
	      &nnodes,	      
	      this->m_beval_f.x,
	      &this->m_beval_f.ld,
	      this->m_rst.x,
	      &this->m_rst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
  };

  void xyz(I subcell,
	   matrix_handle&xyz)
  {
    for (I j=0;j<dim;++j)
      {
 	for (I i=0;i<3;++i)
	  {
	    xyz.x[xyz.ld*j+i] = m_xyz.x[ m_xyz.ld * j + generator->GetNodeIndex(subcell,i)];
	  }
      }
  };
  void f(I subcell,
	 vector_handle&f)
  {
    for (I i=0;i<3;++i)
      {
	f.x[f.ld*i] = this->m_f.x[ this->m_f.ld * generator->GetNodeIndex(subcell,i)];
      }
  };
  
  void xyz(const matrix_handle&xyz)
  {
    this->m_xyz = this->m_beval_xyz.transpose() * xyz;
  };
  
  void f(const vector_handle&f)
  {
    this->m_f = this->m_beval_f.transpose() * f;
    //    vector_handle_print(&f,stdout);
  };
  
};


void dg_print_sol(DG_VIEW& 		dg_view_,
		  mkS shape_x,mkS shape,I N,cst_pR sol,const char * name_,...)
{
  I n = mkS_n(shape);
  I nx = mkS_n(shape_x);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }
    //    R tmp[64];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N*dg_view_.m_nsubcells*3);
    vector_handle f;
    R fx[128];
    vector_handle_def(&f,n,fx,1);
    vector_handle f2;
    R fx2[128];
    vector_handle_def(&f2,nx,fx2,1);
    { I j;
      for (j=0;j<N;++j)
	{
	  for (I i=0;i<n;++i)
	    {
	      fx[i] = sol[n*j+i];
	    }
	  
	  dg_view_.f(f);
	  for (I k=0;k<dg_view_.m_nsubcells;++k)
	    {
	      dg_view_.f(k,f2);
	      for (I i=0;i<nx;++i)
		{
		  fprintf(fil,"" rfmt "\n",f2.x[i]);
		}
	    }
	} } 
    fclose(fil); }
}


void dg_print_mesh(const ns_mesh*	s_,
		   DG_VIEW& 		dg_view_,
		   const char * 	name_,...)
{
  
  //  const I numNodes = ns_mesh_get_numNodes(s_);

  
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }

    I nsubcells = dg_view_.m_nsubcells;
    I nelm = s_->nelm*nsubcells;
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*nsubcells*3);
    
    R cooelm[6];
    matrix_handle xyz;
    matrix_handle_def(&xyz,3,2,cooelm,3);
    R cooelm2[6];
    matrix_handle xyz2;
    matrix_handle_def(&xyz2,3,2,cooelm2,3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);	  
	  ns_mesh_cooelm(s_,
			 &i,
			 cooelm);	  
	  dg_view_.xyz(xyz);
	  
	  for (I k=0;k<nsubcells;++k)
	    {	      
	      dg_view_.xyz(k,xyz2);
	      {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",xyz2.x[xyz2.ld*0+j],xyz2.x[xyz2.ld*1+j],((I)0)); }}
	    }	  
	} } 

    
    fprintf(fil,"Triangles\n" ifmt "\n",nelm); 
    { I i;
      for (i=0;i<nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,((I)0));
	} }

    
    fprintf(fil,"End\n");						
    fclose(fil); }  
}

static const double tria_L1_L1[9] = {
((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)

,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)

};

static const double tria_L2_L2[36] = {
  ((double)1.6666666666666666E-2)
  ,((double)-2.7777777777777779E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)0.)
  ,((double)-1.1111111111111112E-2)
  ,((double)0.)
,((double)-2.7777777777777779E-3)
,((double)1.6666666666666666E-2)
,((double)-2.7777777777777779E-3)
,((double)0.)
,((double)0.)
,((double)-1.1111111111111112E-2)

,((double)-2.7777777777777779E-3)
,((double)-2.7777777777777779E-3)
,((double)1.6666666666666666E-2)
,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)0.)

,((double)0.)
,((double)0.)
,((double)-1.1111111111111112E-2)
,((double)8.8888888888888892E-2)
,((double)4.4444444444444446E-2)
,((double)4.4444444444444446E-2)

,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)0.)
,((double)4.4444444444444446E-2)
,((double)8.8888888888888892E-2)
,((double)4.4444444444444446E-2)

,((double)0.)
,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)4.4444444444444446E-2)
,((double)4.4444444444444446E-2)
,((double)8.8888888888888892E-2)

};

static const double tria_L2p_L2p[49] = {
((double)0.2388888888888889)
,((double)0.10833333333333334)
,((double)0.10833333333333334)
,((double)-0.12222222222222222)
,((double)0.12222222222222222)
,((double)-0.12222222222222222)
,((double)0.16071428571428573)
,((double)0.10833333333333334)
,((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)
,((double)-6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)7.4999999999999997E-2)

,((double)0.10833333333333334)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)
,((double)-6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)7.4999999999999997E-2)

,((double)-0.12222222222222222)
,((double)-6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.5714285714285715E-2)

,((double)0.12222222222222222)
,((double)6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.5714285714285715E-2)

,((double)-0.12222222222222222)
,((double)-6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.5714285714285715E-2)

,((double)0.16071428571428573)
,((double)7.4999999999999997E-2)
,((double)7.4999999999999997E-2)
,((double)-8.5714285714285715E-2)
,((double)8.5714285714285715E-2)
,((double)-8.5714285714285715E-2)
,((double)0.14464285714285716)

};

static const double tria_L3_L3[100] = {
((double)5.6547619047619046E-3)
,((double)8.1845238095238097E-4)
,((double)8.1845238095238097E-4)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)8.1845238095238097E-4)
,((double)5.6547619047619046E-3)
,((double)8.1845238095238097E-4)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)2.6785714285714286E-3)

,((double)8.1845238095238097E-4)
,((double)8.1845238095238097E-4)
,((double)5.6547619047619046E-3)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.6785714285714286E-3)

,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)1.2053571428571429E-2)

,((double)0.)
,((double)1.3392857142857143E-3)
,((double)2.0089285714285712E-3)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)1.2053571428571429E-2)

,((double)2.0089285714285712E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)1.2053571428571429E-2)

,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)1.2053571428571429E-2)

,((double)0.)
,((double)2.0089285714285712E-3)
,((double)1.3392857142857143E-3)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)1.2053571428571429E-2)

,((double)1.3392857142857143E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)1.2053571428571429E-2)

,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.14464285714285716)

};

static const double tria_L1_L1_L1[27] = {
((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)

,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)

};

static const double tria_L2_L2_L2[216] = {
((double)7.1428571428571426E-3)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)4.7619047619047623E-3)
,((double)1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)7.1428571428571426E-3)
,((double)-7.9365079365079365E-4)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)

,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)7.1428571428571426E-3)
,((double)1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)

,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-6.3492063492063492E-3)
,((double)5.7142857142857141E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)

,((double)1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)-6.3492063492063492E-3)
,((double)0.)
,((double)0.)
,((double)1.9047619047619049E-2)
,((double)5.7142857142857141E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)

,((double)4.7619047619047623E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)0.)
,((double)-6.3492063492063492E-3)
,((double)0.)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)5.7142857142857141E-2)

};

static const double tria_L2p_L2p_L2p[343] = {
((double)0.17857142857142858)
,((double)7.6984126984126988E-2)
,((double)7.6984126984126988E-2)
,((double)-9.3650793650793651E-2)
,((double)9.3650793650793651E-2)
,((double)-9.3650793650793651E-2)
,((double)0.11785714285714285)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)9.3650793650793651E-2)
,((double)4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.4285714285714279E-2)
,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)0.11785714285714285)
,((double)5.3571428571428568E-2)
,((double)5.3571428571428568E-2)
,((double)-6.4285714285714279E-2)
,((double)6.4285714285714279E-2)
,((double)-6.4285714285714279E-2)
,((double)0.10607142857142857)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)5.0000000000000003E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)

,((double)9.3650793650793651E-2)
,((double)4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.4285714285714279E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.4285714285714279E-2)
,((double)3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)5.7857142857142857E-2)

,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)

,((double)0.11785714285714285)
,((double)5.3571428571428568E-2)
,((double)5.3571428571428568E-2)
,((double)-6.4285714285714279E-2)
,((double)6.4285714285714279E-2)
,((double)-6.4285714285714279E-2)
,((double)0.10607142857142857)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)
,((double)6.4285714285714279E-2)
,((double)3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)5.7857142857142857E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)
,((double)0.10607142857142857)
,((double)4.8214285714285716E-2)
,((double)4.8214285714285716E-2)
,((double)-5.7857142857142857E-2)
,((double)5.7857142857142857E-2)
,((double)-5.7857142857142857E-2)
,((double)0.10650974025974026)

};

static const double tria_L3_L3_L3[1000] = {
((double)2.5162337662337662E-3)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)7.3051948051948055E-4)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)-1.2175324675324675E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)9.8620129870129864E-4)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)9.8620129870129864E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)6.5746753246753243E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)-1.0957792207792207E-3)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)8.7662337662337668E-4)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)2.5162337662337662E-3)
,((double)1.8939393939393939E-4)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)7.3051948051948055E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.2175324675324675E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)6.5746753246753243E-4)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.0957792207792207E-3)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)6.5746753246753243E-4)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)9.8620129870129864E-4)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)9.8620129870129864E-4)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)8.7662337662337668E-4)

,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)-1.2175324675324675E-5)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)2.5162337662337662E-3)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)7.3051948051948055E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)6.5746753246753243E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)8.7662337662337668E-4)

,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)9.8620129870129864E-4)
,((double)2.1915584415584414E-3)
,((double)1.3149350649350649E-3)
,((double)1.497564935064935E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)1.2820616883116883E-2)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.9310064935064934E-3)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.3011363636363637E-3)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.2735389610389614E-3)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)8.5470779220779228E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)5.9172077922077923E-3)

,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)9.8620129870129864E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.9310064935064934E-3)
,((double)1.3149350649350649E-3)
,((double)2.1915584415584414E-3)
,((double)1.497564935064935E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)1.2820616883116883E-2)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)8.5470779220779228E-3)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)-2.3011363636363637E-3)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-4.2735389610389614E-3)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)

,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)9.8620129870129864E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)6.5746753246753243E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)-4.2735389610389614E-3)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)8.5470779220779228E-3)
,((double)1.497564935064935E-3)
,((double)2.1915584415584414E-3)
,((double)1.3149350649350649E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)1.2820616883116883E-2)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-2.3011363636363637E-3)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)5.9172077922077923E-3)

,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-1.0957792207792207E-3)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.3011363636363637E-3)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)1.497564935064935E-3)
,((double)1.3149350649350649E-3)
,((double)2.1915584415584414E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)1.2820616883116883E-2)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)8.5470779220779228E-3)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-4.2735389610389614E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)

,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)6.5746753246753243E-4)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)9.8620129870129864E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.2735389610389614E-3)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)-2.3011363636363637E-3)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)-4.2735389610389614E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)8.5470779220779228E-3)
,((double)1.3149350649350649E-3)
,((double)1.497564935064935E-3)
,((double)2.1915584415584414E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.2820616883116883E-2)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)5.9172077922077923E-3)

,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)-1.0957792207792207E-3)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)9.8620129870129864E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)6.5746753246753243E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)8.5470779220779228E-3)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-4.2735389610389614E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)2.1915584415584414E-3)
,((double)1.497564935064935E-3)
,((double)1.3149350649350649E-3)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)1.2820616883116883E-2)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)5.9172077922077923E-3)

,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)8.7662337662337668E-4)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)8.7662337662337668E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)8.7662337662337668E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)5.9172077922077923E-3)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)5.9172077922077923E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)5.9172077922077923E-3)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)5.9172077922077923E-3)
,((double)8.7662337662337668E-4)
,((double)8.7662337662337668E-4)
,((double)8.7662337662337668E-4)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)0.10650974025974026)

};

static const double tria_L1_L2_L2[108] = {
((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)2.3809523809523812E-3)
,((double)3.9682539682539683E-4)
,((double)-3.1746031746031746E-3)
,((double)-1.5873015873015873E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)2.3809523809523812E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)3.8095238095238099E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)4.7619047619047623E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)3.8095238095238099E-2)
,((double)2.3809523809523812E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)2.3809523809523812E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)-4.7619047619047623E-3)
,((double)3.8095238095238099E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)-4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)3.8095238095238099E-2)
,((double)1.2698412698412698E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)

,((double)2.3809523809523812E-3)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)3.9682539682539683E-4)
,((double)2.3809523809523812E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)1.2698412698412698E-2)
,((double)3.8095238095238099E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)3.8095238095238099E-2)

};

static const double tria_L1_L2p_L2p[147] = {
((double)8.4920634920634924E-2)
,((double)2.7777777777777776E-2)
,((double)2.7777777777777776E-2)
,((double)-2.3809523809523808E-2)
,((double)2.3809523809523808E-2)
,((double)-2.3809523809523808E-2)
,((double)5.3571428571428568E-2)
,((double)2.7777777777777776E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)-1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)2.1428571428571429E-2)
,((double)2.7777777777777776E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)-1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)2.1428571428571429E-2)
,((double)-2.3809523809523808E-2)
,((double)-1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-2.1428571428571429E-2)
,((double)2.3809523809523808E-2)
,((double)1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)2.1428571428571429E-2)
,((double)-2.3809523809523808E-2)
,((double)-1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-2.1428571428571429E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)2.1428571428571429E-2)
,((double)-2.1428571428571429E-2)
,((double)2.1428571428571429E-2)
,((double)-2.1428571428571429E-2)
,((double)4.8214285714285716E-2)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)5.0000000000000003E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

};

static const double tria_L1_L3_L3[300] = {
((double)4.464285714285714E-3)
,((double)3.720238095238095E-4)
,((double)3.720238095238095E-4)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)0.)
,((double)3.720238095238095E-4)
,((double)5.9523809523809529E-4)
,((double)7.4404761904761911E-5)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)7.4404761904761911E-5)
,((double)5.9523809523809529E-4)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)-1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)2.2321428571428571E-4)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)-8.9285714285714283E-4)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)2.6785714285714286E-3)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)4.8214285714285716E-2)
,((double)5.9523809523809529E-4)
,((double)3.720238095238095E-4)
,((double)7.4404761904761911E-5)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)4.464285714285714E-3)
,((double)3.720238095238095E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)0.)
,((double)7.4404761904761911E-5)
,((double)3.720238095238095E-4)
,((double)5.9523809523809529E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)1.1160714285714285E-3)
,((double)2.6785714285714286E-3)
,((double)1.5625000000000001E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)1.2053571428571429E-2)
,((double)1.5625000000000001E-3)
,((double)2.6785714285714286E-3)
,((double)1.1160714285714285E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)1.2053571428571429E-2)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)4.8214285714285716E-2)

,((double)5.9523809523809529E-4)
,((double)7.4404761904761911E-5)
,((double)3.720238095238095E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)1.3392857142857143E-3)
,((double)7.4404761904761911E-5)
,((double)5.9523809523809529E-4)
,((double)3.720238095238095E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)3.720238095238095E-4)
,((double)4.464285714285714E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)0.)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-8.9285714285714283E-4)
,((double)-1.3392857142857143E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)0.)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)2.6785714285714286E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.6785714285714286E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)-8.9285714285714283E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)4.8214285714285716E-2)

};

static const double tria_L2_L3_L3[600] = {
  ((double)3.2738095238095239E-3)
  ,((double)2.2321428571428571E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)0.)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.9841269841269841E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-1.9841269841269841E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-1.9841269841269841E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.2738095238095239E-3)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)3.4722222222222222E-5)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.9841269841269841E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.3392857142857144E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)0.)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)2.6785714285714287E-4)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)-1.9841269841269841E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)2.2321428571428571E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-1.9841269841269841E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.2738095238095239E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-2.2321428571428571E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)1.1904761904761906E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.5793650793650796E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.142857142857143E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)5.7857142857142857E-2)

  ,((double)3.9682539682539683E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-5)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)2.142857142857143E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)1.1904761904761906E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)2.5793650793650796E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714283E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.7857142857142857E-2)

  ,((double)1.1904761904761906E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)2.5793650793650796E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)2.142857142857143E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)1.1904761904761906E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)5.7857142857142857E-2)

};

static const double tria_L2_L2p_L2p[294] = {
  ((double)7.1428571428571426E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-8.7301587301587304E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)1.6666666666666666E-2)
  ,((double)0.)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)-5.9523809523809521E-3)
  ,((double)0.)
  ,((double)-5.5555555555555558E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)0.)
  ,((double)-5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)-8.7301587301587304E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)0.)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)1.6666666666666666E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)1.5873015873015873E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)7.7777777777777779E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)5.7857142857142857E-2)

  ,((double)9.3650793650793651E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)6.9841269841269843E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)2.2222222222222223E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.214285714285714E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)2.2222222222222223E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)6.9841269841269843E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)5.7857142857142857E-2)

  ,((double)7.7777777777777779E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)5.7857142857142857E-2)

};

static const double tria_L3_L2p_L2p[490] = {
  ((double)1.1111111111111112E-2)
  ,((double)2.7777777777777779E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)2.142857142857143E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)2.142857142857143E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)2.142857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)2.142857142857143E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)8.7662337662337668E-4)
  ,((double)3.1746031746031746E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)1.984126984126984E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)3.9682539682539683E-4)
  ,((double)-3.9682539682539683E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)1.984126984126984E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)2.3809523809523812E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.7662337662337668E-4)

  ,((double)3.1746031746031746E-3)
  ,((double)1.984126984126984E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)3.9682539682539683E-4)
  ,((double)-3.9682539682539683E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.984126984126984E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)7.1428571428571426E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.7662337662337668E-4)

  ,((double)2.3214285714285715E-2)
  ,((double)1.7857142857142857E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)5.3571428571428572E-3)
  ,((double)1.7857142857142856E-2)
  ,((double)0.)
  ,((double)-8.9285714285714281E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.7857142857142856E-2)
  ,((double)2.6785714285714284E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.0714285714285714E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)0.)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)-1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)1.0714285714285714E-2)
  ,((double)1.7857142857142857E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)-1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)2.5000000000000001E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)2.3214285714285715E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)2.5000000000000001E-2)
  ,((double)2.6785714285714284E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)-2.3214285714285715E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.3214285714285715E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)2.5000000000000001E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)2.3214285714285715E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)2.5000000000000001E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.3214285714285715E-2)
  ,((double)5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)1.7857142857142856E-2)
  ,((double)-8.9285714285714281E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)0.)
  ,((double)1.7857142857142856E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.0714285714285714E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)1.0714285714285714E-2)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)0.)
  ,((double)1.7857142857142857E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)0.)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)0.11785714285714285)
  ,((double)5.3571428571428568E-2)
  ,((double)5.3571428571428568E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)0.10607142857142857)
  ,((double)5.3571428571428568E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)5.3571428571428568E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)5.7857142857142857E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)0.10607142857142857)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)5.7857142857142857E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)0.10650974025974026)

};

static const double tria_L1dx_L1[9] = {
  ((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

};

static const double tria_L2dx_L2[36] = {
  ((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.10000000000000001)
  ,((double)0.10000000000000001)
  ,((double)-3.3333333333333333E-2)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)0.10000000000000001)
  ,((double)-0.10000000000000001)
  ,((double)0.)
  ,((double)0.)
  ,((double)-0.13333333333333333)
  ,((double)0.13333333333333333)

  ,((double)-3.3333333333333333E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)0.13333333333333333)
  ,((double)0.26666666666666666)
  ,((double)0.26666666666666666)

  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-6.6666666666666666E-2)
  ,((double)-0.13333333333333333)
  ,((double)-0.26666666666666666)
  ,((double)-0.26666666666666666)

};

static const double tria_L3dx_L3[100] = {
  ((double)-3.8095238095238099E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)1.607142857142857E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)6.160714285714286E-2)
  ,((double)6.160714285714286E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-1.607142857142857E-2)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-9.6428571428571433E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)-0.14464285714285716)

  ,((double)-3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)0.)
  ,((double)9.6428571428571433E-2)
  ,((double)0.)
  ,((double)-9.6428571428571433E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.14464285714285716)

  ,((double)1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)0.19285714285714287)
  ,((double)2.4107142857142858E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.24107142857142858)

  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)5.8928571428571427E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)0.19285714285714287)
  ,((double)0.19285714285714287)
  ,((double)-7.2321428571428578E-2)
  ,((double)9.6428571428571433E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-5.8928571428571427E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.19285714285714287)
  ,((double)7.2321428571428578E-2)
  ,((double)-9.6428571428571433E-2)

  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.24107142857142858)

  ,((double)-1.607142857142857E-2)
  ,((double)1.607142857142857E-2)
  ,((double)0.)
  ,((double)0.14464285714285716)
  ,((double)-0.14464285714285716)
  ,((double)-0.24107142857142858)
  ,((double)-9.6428571428571433E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)0.24107142857142858)
  ,((double)0.)

};

static const double tria_L1_L1dx_L1[27] = {
  ((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

};

static const double tria_L2_L2dx_L2[216] = {
  ((double)-3.0952380952380953E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-7.1428571428571426E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)4.3650793650793652E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.8095238095238099E-2)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)3.1746031746031744E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)9.5238095238095247E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-9.5238095238095247E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.0952380952380953E-2)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.8095238095238099E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.5873015873015873E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.5873015873015873E-3)
  ,((double)-9.5238095238095247E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)

  ,((double)3.5714285714285713E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)4.3650793650793652E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047616E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-4.7619047619047616E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)6.3492063492063489E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-2.5396825396825397E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)

  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063489E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.2698412698412698E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-2.5396825396825397E-2)
  ,((double)-0.10158730158730159)
  ,((double)0.)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)0.15238095238095239)
  ,((double)7.6190476190476197E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-0.15238095238095239)
  ,((double)-7.6190476190476197E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)2.5396825396825397E-2)
  ,((double)0.)
  ,((double)0.10158730158730159)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)7.6190476190476197E-2)
  ,((double)0.15238095238095239)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-7.6190476190476197E-2)
  ,((double)-0.15238095238095239)

};

static const double tria_L3_L3dx_L3[1000] = {
  ((double)-1.7708333333333333E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)2.976190476190476E-3)
  ,((double)1.488095238095238E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)3.3482142857142857E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6785714285714284E-2)
  ,((double)2.2767857142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.0089285714285716E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)2.2098214285714287E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.2767857142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.0133928571428573E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-6.6964285714285718E-5)
  ,((double)-1.4062499999999999E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.4375000000000006E-3)
  ,((double)0.)
  ,((double)4.0178571428571428E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-6.6964285714285718E-5)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.2098214285714286E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-2.976190476190476E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)1.488095238095238E-3)
  ,((double)1.7708333333333333E-2)
  ,((double)1.488095238095238E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6785714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.2767857142857143E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-6.4285714285714285E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-2.6785714285714284E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)-2.2098214285714287E-2)
  ,((double)1.1049107142857144E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)5.3571428571428572E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.4062499999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-5.3571428571428572E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.4464285714285714E-2)

  ,((double)-1.488095238095238E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.976190476190476E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)1.488095238095238E-3)
  ,((double)2.976190476190476E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.2187500000000003E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)0.)
  ,((double)4.2187500000000003E-3)
  ,((double)0.)
  ,((double)-3.8169642857142855E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)3.8169642857142855E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)3.214285714285714E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.214285714285714E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-6.6964285714285718E-5)
  ,((double)1.4062499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-3.8169642857142855E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-4.0178571428571432E-2)
  ,((double)1.40625E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)3.3482142857142857E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.40625E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)6.8303571428571432E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.0089285714285716E-2)
  ,((double)1.6071428571428571E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)4.2187500000000003E-3)
  ,((double)4.3392857142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-9.0401785714285723E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-2.3504464285714285E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)6.5089285714285711E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.40625E-2)
  ,((double)3.2142857142857142E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-6.8303571428571432E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.40625E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.0312500000000002E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-4.2187500000000003E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.8080357142857145E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.3392857142857143E-2)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.3504464285714285E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-6.5089285714285711E-2)

  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.4866071428571428E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)2.6116071428571429E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)-1.40625E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.2053571428571428E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.700892857142857E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)-2.2098214285714287E-2)
  ,((double)-3.8169642857142855E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.8080357142857145E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-2.3504464285714285E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)0.10848214285714286)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.0089285714285712E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-0.13017857142857142)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-8.6785714285714285E-2)

  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.4866071428571428E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)6.8303571428571432E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)-1.40625E-2)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.2053571428571428E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)3.0133928571428573E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.424107142857143E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)0.10848214285714286)
  ,((double)5.424107142857143E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-0.10848214285714286)
  ,((double)-5.424107142857143E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.424107142857143E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.10848214285714286)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.0044642857142858E-2)
  ,((double)-6.8303571428571432E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)3.2142857142857142E-3)
  ,((double)1.40625E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.4866071428571428E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.1049107142857144E-2)
  ,((double)-3.0133928571428573E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)6.4285714285714285E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.2098214285714286E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-5.424107142857143E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)0.10848214285714286)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-0.10848214285714286)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-5.424107142857143E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)0.10848214285714286)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)1.0044642857142858E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.40625E-2)
  ,((double)-4.0178571428571432E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4866071428571428E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.2098214285714287E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-1.4464285714285714E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-4.700892857142857E-2)
  ,((double)0.)
  ,((double)-2.0089285714285712E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)2.3504464285714285E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-0.10848214285714286)
  ,((double)-6.5089285714285711E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)0.)
  ,((double)3.2544642857142855E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.13017857142857142)
  ,((double)8.6785714285714285E-2)

  ,((double)-2.6785714285714286E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.8928571428571428E-2)
  ,((double)1.4732142857142858E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.8928571428571428E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)0.13017857142857142)
  ,((double)8.4375000000000006E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.17357142857142857)
  ,((double)4.8214285714285711E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.424107142857143E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-4.3392857142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)-0.17357142857142857)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)0.)
  ,((double)6.5089285714285711E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)-8.6785714285714285E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)8.6785714285714285E-2)
  ,((double)0.)

};

static const double tria_L1dy_L1[9] = {
  ((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)

};

static const double tria_L2dy_L2[36] = {
  ((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-3.3333333333333333E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.10000000000000001)
  ,((double)0.10000000000000001)

  ,((double)3.3333333333333333E-2)
  ,((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.26666666666666666)
  ,((double)-0.26666666666666666)
  ,((double)-0.13333333333333333)

  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.26666666666666666)
  ,((double)0.26666666666666666)
  ,((double)0.13333333333333333)

  ,((double)0.10000000000000001)
  ,((double)0.)
  ,((double)-0.10000000000000001)
  ,((double)0.13333333333333333)
  ,((double)-0.13333333333333333)
  ,((double)0.)

};

static const double tria_L3dy_L3[100] = {
  ((double)-3.8095238095238099E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)1.607142857142857E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)1.1309523809523809E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)6.160714285714286E-2)
  ,((double)6.160714285714286E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)-1.607142857142857E-2)

  ,((double)2.6785714285714286E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-0.19285714285714287)
  ,((double)-2.4107142857142858E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)-0.24107142857142858)

  ,((double)-1.3392857142857142E-2)
  ,((double)-5.8928571428571427E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.19285714285714287)
  ,((double)7.2321428571428578E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-9.6428571428571433E-2)

  ,((double)1.3392857142857142E-2)
  ,((double)5.8928571428571427E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)0.19285714285714287)
  ,((double)0.19285714285714287)
  ,((double)-7.2321428571428578E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)9.6428571428571433E-2)

  ,((double)1.3392857142857142E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)-4.8214285714285716E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)0.19285714285714287)
  ,((double)9.6428571428571433E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.24107142857142858)

  ,((double)-3.4821428571428573E-2)
  ,((double)0.)
  ,((double)-6.160714285714286E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)0.)
  ,((double)9.6428571428571433E-2)
  ,((double)0.14464285714285716)

  ,((double)6.160714285714286E-2)
  ,((double)0.)
  ,((double)3.4821428571428573E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)0.)
  ,((double)-0.14464285714285716)

  ,((double)-1.607142857142857E-2)
  ,((double)0.)
  ,((double)1.607142857142857E-2)
  ,((double)0.24107142857142858)
  ,((double)9.6428571428571433E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)-0.24107142857142858)
  ,((double)-0.14464285714285716)
  ,((double)0.14464285714285716)
  ,((double)0.)

};

static const double tria_L1_L1dy_L1[27] = {
  ((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)

  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)

};

static const double tria_L2_L2dy_L2[216] = {
  ((double)-3.0952380952380953E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.1428571428571426E-3)
  ,((double)4.3650793650793652E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-9.5238095238095247E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.3650793650793652E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-4.7619047619047616E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047616E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)0.)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)

  ,((double)3.5714285714285713E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.0952380952380953E-2)
  ,((double)4.7619047619047623E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.5873015873015873E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-9.5238095238095247E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)-3.1746031746031744E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-0.15238095238095239)
  ,((double)-7.6190476190476197E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)0.15238095238095239)
  ,((double)7.6190476190476197E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)0.10158730158730159)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)

  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063489E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-7.6190476190476197E-2)
  ,((double)-0.15238095238095239)
  ,((double)-5.0793650793650794E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)7.6190476190476197E-2)
  ,((double)0.15238095238095239)
  ,((double)5.0793650793650794E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-2.5396825396825397E-2)
  ,((double)0.)
  ,((double)-0.10158730158730159)
  ,((double)-2.5396825396825397E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063489E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)0.)
  ,((double)-3.1746031746031744E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)0.)

};

static const double tria_L3_L3dy_L3[1000] = {
  ((double)-1.7708333333333333E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.976190476190476E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)1.488095238095238E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)1.4732142857142858E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)6.6964285714285718E-5)
  ,((double)-4.0178571428571425E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)5.0223214285714289E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-5.0223214285714289E-3)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.4062499999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)8.4375000000000006E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)6.4285714285714285E-3)
  ,((double)3.0133928571428573E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)2.4107142857142856E-3)
  ,((double)2.2767857142857143E-3)
  ,((double)2.2098214285714287E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)2.0089285714285716E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-2.976190476190476E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.8541666666666662E-4)
  ,((double)2.976190476190476E-3)
  ,((double)1.488095238095238E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)1.2053571428571429E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-3.8169642857142855E-3)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.214285714285714E-2)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)3.214285714285714E-2)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.6964285714285718E-5)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)0.)
  ,((double)-1.8080357142857143E-3)
  ,((double)0.)
  ,((double)-2.4107142857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)0.)
  ,((double)4.2187500000000003E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)0.)
  ,((double)1.8080357142857143E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-4.2187500000000003E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)0.)

  ,((double)-1.488095238095238E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.976190476190476E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.488095238095238E-3)
  ,((double)1.488095238095238E-3)
  ,((double)1.7708333333333333E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)1.4062499999999999E-3)
  ,((double)0.)
  ,((double)-7.4330357142857141E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)0.)
  ,((double)2.8124999999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.8214285714285711E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-2.8124999999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)6.6964285714285718E-5)
  ,((double)5.3571428571428572E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.2098214285714286E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6785714285714284E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-2.2098214285714287E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)7.0312500000000002E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.2767857142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)7.0312500000000002E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.4732142857142858E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)1.4464285714285714E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-4.0178571428571432E-2)
  ,((double)1.40625E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)0.)
  ,((double)1.4866071428571428E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-0.10848214285714286)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)2.3504464285714285E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-2.0089285714285712E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.2053571428571428E-3)
  ,((double)-4.700892857142857E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)0.)
  ,((double)2.2098214285714287E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)7.2321428571428578E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.8080357142857145E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)7.2321428571428571E-3)
  ,((double)0.13017857142857142)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)8.6785714285714285E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.40625E-2)
  ,((double)3.2142857142857142E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-6.8303571428571432E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4866071428571428E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-0.10848214285714286)
  ,((double)-5.424107142857143E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)7.633928571428571E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)0.10848214285714286)
  ,((double)5.424107142857143E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.424107142857143E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.4285714285714285E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.10848214285714286)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)

  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.4866071428571428E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.40625E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)6.8303571428571432E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)5.424107142857143E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-0.10848214285714286)
  ,((double)2.7120535714285715E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)7.633928571428571E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)0.10848214285714286)
  ,((double)-2.7120535714285715E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)5.424107142857143E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)5.424107142857143E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)3.0133928571428573E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)0.)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)4.8214285714285711E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-0.10848214285714286)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)

  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.4866071428571428E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.40625E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.10848214285714286)
  ,((double)3.6160714285714289E-2)
  ,((double)-2.3504464285714285E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-2.2098214285714287E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.4464285714285714E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)-1.8080357142857145E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.1049107142857144E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.700892857142857E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)0.)
  ,((double)-7.2321428571428571E-3)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-3.2544642857142855E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-8.6785714285714285E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.0044642857142858E-2)
  ,((double)-6.8303571428571432E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)3.2142857142857142E-3)
  ,((double)1.40625E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)-1.40625E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)2.3504464285714285E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-1.6071428571428571E-3)
  ,((double)0.)
  ,((double)-2.0089285714285716E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857145E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)-4.2187500000000003E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)9.0401785714285723E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-4.3392857142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-6.5089285714285711E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)1.0044642857142858E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.40625E-2)
  ,((double)-4.0178571428571432E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.3482142857142857E-4)
  ,((double)-1.8080357142857143E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)6.8303571428571432E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)-1.40625E-2)
  ,((double)-3.2142857142857142E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.3504464285714285E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)4.2187500000000003E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)0.)
  ,((double)1.6071428571428571E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)6.5089285714285711E-2)

  ,((double)-2.6785714285714286E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.8928571428571428E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)-2.8928571428571428E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-8.4375000000000006E-3)
  ,((double)-6.5089285714285711E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)-0.17357142857142857)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)4.8214285714285711E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.424107142857143E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.17357142857142857)
  ,((double)-4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)0.13017857142857142)
  ,((double)6.0267857142857146E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-3.2544642857142855E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)1.4464285714285714E-2)
  ,((double)8.6785714285714285E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-8.6785714285714285E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)0.)

};

#define dim 		((I)2)
#define nfaceinelm 	((I)3)

int comp(const void * a ,const void * b)
{
  cst_pI a_ = (cst_pI)a;
  cst_pI b_ = (cst_pI)b;
  if (a_[0] < b_[0] )
    {
      return -1;
    }
  else if (a_[0] > b_[0] )
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

struct DG_JACOBIAN
{
  I  m_nelm;
  I  m_nfaceinelm;
  I  m_size_block;
  I  m_size_blockXsize_block;
  pR m_x;
  
  I  m_n;
  I  m_nc;
  pI m_begin;
  pI m_index;
  pR m_values;

  pR m_idiagonal {};
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


  void inverse_diagonal()
  {
    const I n = m_size_block;
    if ( NULL == m_idiagonal)
      {
	this->m_idiagonal = (pR)malloc(sizeof(R)*m_nelm * n*n);
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

#if 0
  void inverse_diagonal_gemv(cst_pR a,cst_pR x,cst_pI xoff,cst_pR ry,pR y,cst_pI yoff)
  {
    for (I ielm=0;ielm<m_nelm;++ielm)
      {
	Blas_dgemv("N",
		   &this->m_size_block,
		   &this->m_size_block,
		   a,
		   &this->m_idiagonal[ielm*this->m_size_blockXsize_block],
		   &this->m_size_block,
		   x,
		   xoff,
		   ry,
		   y,
		   yoff);
      }
  };
#endif

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
    this->m_x 				= (pR)malloc(sizeof(R) * this->m_size_blockXsize_block * (this->m_nfaceinelm + 1) * this->m_nelm);
    
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
  };



  DG_JACOBIAN(I 	nelm_,
	      I 	nfaceinelm_,
	      I 	size_block_,
	      cst_pI 	adj_,
	      cst_pI 	adjoff_,
	      cst_pI    perm_)
  {
    this->m_size_block 			= size_block_;
    this->m_size_blockXsize_block 	= size_block_*size_block_;
    this->m_nelm			= nelm_;
    this->m_nfaceinelm			= nfaceinelm_;
    this->m_x 				= (pR)malloc(sizeof(R) * this->m_size_blockXsize_block * (this->m_nfaceinelm + 1) * this->m_nelm);
    
    this->m_n 				= this->m_size_block * this->m_nelm;
    this->m_nc 				= this->m_size_block * this->m_size_block * (this->m_nfaceinelm+1) * this->m_nelm;
    this->m_begin 			= (pI)calloc(this->m_n+1,sizeof(I));
    this->m_index 			= (pI)malloc(this->m_nc*sizeof(I));
    this->m_values 			= (pR)malloc(this->m_nc*sizeof(R));
    
    for (I ielm=0;ielm<nelm_;++ielm)
      {
	I jelm = perm_[ielm+1]-1;

	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_size_block;++k)
	  {
	    m_begin[jelm * m_size_block + k + 1] += m_size_block;
	  }
#if 1
	//
	// Extra diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*jelm+localFaceIndex];
	    if (nei)
	      {
		for (I k = 0;k <  this->m_size_block;++k)
		  {
		    m_begin[jelm * m_size_block + k + 1] += m_size_block;
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
	I jelm = perm_[ielm+1]-1;
#if 1
	//
	// Before diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*jelm+localFaceIndex];
	    if (nei)
	      {
		nei = perm_[nei];
		if (nei-1 < jelm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[jelm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[jelm * m_size_block + k]+=1;
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
		m_index[m_begin[jelm * m_size_block + k]] = jelm * m_size_block + j;
		m_begin[jelm * m_size_block + k]+=1;
	      }
	  }
#if 1
	//
	// After diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*jelm+localFaceIndex];
	    if (nei)
	      {
		nei = perm_[nei];
		if (nei - 1 > jelm)
		  {
		    for (I k = 0;k <  this->m_size_block;++k)
		      {
			for (I j = 0;j <  this->m_size_block;++j)
			  {					    
			    m_index[m_begin[jelm * m_size_block + k]] = (nei-1) * m_size_block + j;
			    m_begin[jelm * m_size_block + k]+=1;
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
    fprintf(stdout,"hello\n");
    //    this->spy("dgperm.txt");
    //    exit(1);
  };

  
  void addelm(I 	ielm_,
	      I 	jelm_,
	      R         s_,
	      pR 	block_x_,
	      I  	block_off_)
  {
#if 0
    printf("addelm\n");
    for (I i = 0;i<m_size_block;++i)
      {
	for (I k=0;k<m_size_block;++k)
	  {
	    printf("    " rfmt "",block_x_[block_off_*k+i]);
	  }
	printf("\n");
      }
    printf("addelm done\n");

#endif
    for (I i = 0;i<m_size_block;++i)
      {
	I bound = m_begin[ielm_*m_size_block + i + 1];
	for (I at = m_begin[ielm_*m_size_block + i];at<bound;++at)
	  {
	    if (m_index[at] == jelm_ * m_size_block)
	      {
		for (I k=0;k<m_size_block;++k)
		  {
		    m_values[at+k] = block_x_[block_off_*k+i] + s_ * m_values[at+k];
		  }
		break;
	      }
	  }
      }
  };


#if 0  

  void jacobian(const matrix_handle&	xelm,
		const matrix_handle& 	aelm,
		const matrix_handle& 	uelm,
		const matrix_handle& 	felm,
		const matrix_handle& 	felmnei[],
		matrix_handle& 		jacobianelm)
  {
    auto mass = build_mass * aelm;
    jacobianelm += mass * felm;
    
    auto convection = build_convection * uelm;
    jacobianelm += convection * felm;

    auto flux = build_flux * uelm;
    residuelm += flux * felm;
    
  };

  void residual(const matrix_handle&	xelm,
		const matrix_handle& 	aelm,
		const matrix_handle& 	uelm,
		const matrix_handle& 	felm,
		matrix_handle& 		residuelm)
  {    
    auto mass = build_mass * aelm;
    residuelm += mass * felm;
    
    auto convection = build_convection * uelm;
    residuelm += convection * felm;

    auto flux = build_flux * uelm;
    residuelm += flux * felm;
  };
  
  void residual(I nelm_,
		cst_pI adj_,
		cst_pI adjoff_)
  {
    for (I ielm=0;ielm<nelm_;++ielm)
      {
	diag_residual(ielm);
      }    
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_[0]*ielm+localFaceIndex];
	    if (nei)
	      {

	      }
	  }
      }
	
  };
  
  void extra_diag_residual(I 			ielm_,
			   const double *	xelm,
			   const double * 	uelm,
			   const double * 	felm,
			   const double * 	felmnei,
			   double * 		residu)
  {    
  };
  
  void diag_residual(I 			ielm_,
		     const double *	xelm,
		     const double * 	aelm,
		     const double * 	uelm,
		     const double * 	felm,
		     double * 		residu)
  {
  };

  
  void diag(I 				ielm_,
	    const double *		xelm,
	    const double * 		aelm,
	    const double * 		uelm,
	    double * 			diag)
  {  
  };
  
  void extra_diag(I 			ielm_,
		  I 			jelm_,
		  const double *	xelm,
		  const double * 	uelm,
		  double * 		extra_diag)
  {
  };
#endif  
};





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

  R nrmelm[dim * nfaceinelm];

  R lcrhs[128];
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



void solve_gauss_seidel(DG_JACOBIAN&jacobian,I mxiter,I degree,pR x_,pR b_)
{

#if 0
  
  I niter = 40;
  for (I iter=0;iter < niter;++iter)
    {

      
      solve_gauss_seidel(&jacobian[degree],7,degree,corr_level[degree],rhs_level[degree]);

      
      {
	
	//
	// Project the residu on degree k.
	//
	I n = mkS_n(&shape_EE[k]);
	I N = mkS_n(&shape_EE[k+1]);
	for (I jelm=0;jelm<mesh->nelm;++jelm)
	  {
	    for (I j=0;j<n;++j)
	      {
		corr_level[degree_solve-1][jelm*n + j] = corr_level[degree_solve+1][jelm* N + j];
	      }
	  }	  	  
      }
      
      {
	
	//
	// Project the residu on degree k.
	//
	I n = mkS_n(&shape_EE[k]);
	I N = mkS_n(&shape_EE[k+1]);
	for (I jelm=0;jelm<mesh->nelm;++jelm)
	  {
	    for (I j=0;j<n;++j)
	      {
		rhs_level[degree_solve-1][jelm*n + j] = rhs_level[degree_solve+1][jelm* N + j];
	      }
	  }	  	  
      }

      solve_gauss_seidel(&jacobian[degree-1],700,degree-1,corr_level[degree-1],rhs_level[degree-1]);

    }

#endif
  jacobian.inverse_diagonal();
  jacobian.inverse_lower_gemv(b_,&jacobian.m_size_block);
  
  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);

  
  for (I i=0;i<mxiter;++i)
    {
      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);

      
      //
      // Compute the correction.
      //
      jacobian.extra_gemv_upper(tmpx,tmp2);
      jacobian.inverse_lower_gemv(tmpx,&jacobian.m_size_block);      
      daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);
      daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);

      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<mesh->nelm;++ielm)
	      {
		double h1 = 0.0;
		I start = ((ideg) * (ideg+1) )/2;
		I bound = ((ideg+1) * (ideg+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*jacobian.m_size_block+j];
		    h1 += x*x;
		  }
		h1 *= mesh->jacelm[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }
#if 0      
      double mr1=-1.0;
      daxpy(&jacobian.m_n,&mr1,corr_,&n1,tmp2,&n1);
      double h = dnrm2(&jacobian.m_n,tmp2,&n1);
      double h0 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block];
	  h0 += mesh->jacelm[ielm]*x*x;
	}
      double h1 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+1];
	  double y = tmp2[ielm*jacobian.m_size_block+2];
	  h1 += mesh->jacelm[ielm]*(x*x + y*y);
	}
      double h2 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+3];
	  double y = tmp2[ielm*jacobian.m_size_block+4];
	  double z = tmp2[ielm*jacobian.m_size_block+5];
	  h2 += mesh->jacelm[ielm]*(x*x + y*y + z*z);
	}
      double h3 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+6];
	  double y = tmp2[ielm*jacobian.m_size_block+7];
	  double z = tmp2[ielm*jacobian.m_size_block+8];
	  double t = tmp2[ielm*jacobian.m_size_block+9];
	  h3 += mesh->jacelm[ielm]*(x*x + y*y + z*z + t*t);
	}
      printf(" " ifmt " %e %e %e %e %e\n",i,h,sqrt(h0),sqrt(h1),sqrt(h2),sqrt(h3));
#endif

      if (h < 1.0e-12)
	{
	  break;
	}
      
#if 0
	    dg_print_mesh( mesh,
			   view,
			   "dgjacobi." ifmt,i);  
	    dg_print_sol(view,
			 shape_X,
			 shape_E,
			 jacobian.m_nelm,
			 corr_,
			 "dgjacobi." ifmt,i);
#endif
	}
      //
      //
      // corr = extra * corr;
      // corr = diag^{-1} * corr
      // corr += rhs;
      //
      //
    }



void solve_gauss_seidel_select(DG_JACOBIAN&jacobian,I q,I degree,pR x_,pR b_)
{
  jacobian.inverse_diagonal();
  jacobian.inverse_lower_gemv(b_,&jacobian.m_size_block);

  I n1=1;
  //  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpx = (pR)calloc(jacobian.m_n,sizeof(R));
  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  
  for (I i=0;i<10;++i)
    {
      for (I j=0;j<jacobian.m_n;++j) tmpx[i]=0.0;
      
      //      dcopy(&jacobian.m_n,x_,&n1,tmpx,&n1);
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*jacobian.m_size_block+j] = x_[ielm*jacobian.m_size_block + j];
	    }
	}
      //
      // Compute the correction.
      //
      jacobian.extra_gemv_upper(tmpx,tmp2);
      
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  
	  for (I j = 0;j<jacobian.m_size_block;++j)
	    {
	      if (j<start || j>=bound)
		{
		  tmpx[ielm*jacobian.m_size_block+j] = 0.0;
		}
	    }
	}

      jacobian.inverse_lower_gemv(tmpx,&jacobian.m_size_block);      
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  
	  for (I j = 0;j<jacobian.m_size_block;++j)
	    {
	      if (j<start || j>=bound)
		{
		  tmpx[ielm*jacobian.m_size_block+j] = 0.0;
		}
	    }
	}

      // daxpy(&jacobian.m_n,&mr1,x_,&n1,tmpx,&n1);
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*jacobian.m_size_block+j] -= x_[ielm*jacobian.m_size_block + j];
	    }
	}

      
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      tmpx[ielm*jacobian.m_size_block+j] += b_[ielm*jacobian.m_size_block + j];
	    }
	}
      //      daxpy(&jacobian.m_n,&r1,b_,&n1,tmpx,&n1);

      //
      //
      //
      // daxpy(&jacobian.m_n,&r1,tmpx,&n1,x_,&n1);      
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  //	  double h1 = 0.0;
	  I start = ((q) * (q+1) )/2;
	  I bound = ((q+1) * (q+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      x_[ielm*jacobian.m_size_block + j] += tmpx[ielm*jacobian.m_size_block+j];
	    }
	}

      
      double h = dnrm2(&jacobian.m_n,tmpx,&n1);
      printf(" " ifmt " %e",i,h);
      {
	// for (I ideg=0;ideg<=degree;++ideg)
	  {
	    double nrm=0.0;
	    for (I ielm=0;ielm<mesh->nelm;++ielm)
	      {
		double h1 = 0.0;
		I start = ((q) * (q+1) )/2;
		I bound = ((q+1) * (q+2) )/2;
		for (I j = start;j<bound;++j)
		  {
		    double x = tmpx[ielm*jacobian.m_size_block+j];
		    h1 += x*x;
		  }
		h1 *= mesh->jacelm[ielm];
		nrm+=h1;
	      }
	    printf(" %e",sqrt(nrm));
	  }
	printf("\n");
      }
#if 0      
      double mr1=-1.0;
      daxpy(&jacobian.m_n,&mr1,corr_,&n1,tmp2,&n1);
      double h = dnrm2(&jacobian.m_n,tmp2,&n1);
      double h0 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block];
	  h0 += mesh->jacelm[ielm]*x*x;
	}
      double h1 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+1];
	  double y = tmp2[ielm*jacobian.m_size_block+2];
	  h1 += mesh->jacelm[ielm]*(x*x + y*y);
	}
      double h2 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+3];
	  double y = tmp2[ielm*jacobian.m_size_block+4];
	  double z = tmp2[ielm*jacobian.m_size_block+5];
	  h2 += mesh->jacelm[ielm]*(x*x + y*y + z*z);
	}
      double h3 =0.0;
      for (I ielm=0;ielm<mesh->nelm;++ielm)
	{
	  double x = tmp2[ielm*jacobian.m_size_block+6];
	  double y = tmp2[ielm*jacobian.m_size_block+7];
	  double z = tmp2[ielm*jacobian.m_size_block+8];
	  double t = tmp2[ielm*jacobian.m_size_block+9];
	  h3 += mesh->jacelm[ielm]*(x*x + y*y + z*z + t*t);
	}
      printf(" " ifmt " %e %e %e %e %e\n",i,h,sqrt(h0),sqrt(h1),sqrt(h2),sqrt(h3));
#endif

      if (h < 1.0e-12)
	{
	  break;
	}
      
#if 0
	    dg_print_mesh( mesh,
			   view,
			   "dgjacobi." ifmt,i);  
	    dg_print_sol(view,
			 shape_X,
			 shape_E,
			 jacobian.m_nelm,
			 corr_,
			 "dgjacobi." ifmt,i);
#endif
	}
      //
      //
      // corr = extra * corr;
      // corr = diag^{-1} * corr
      // corr += rhs;
      //
      //
    }






void solve_gauss_seidel_mg(DG_JACOBIAN&jacobian,I degree,pR x_,pR b_)
{
  jacobian.inverse_diagonal();
  jacobian.inverse_lower_gemv(b_,&jacobian.m_size_block);
  
  //  R r1=1.0;I n1=1;R mr1=-1.0;
  pR tmpb = (pR)malloc(sizeof(R)*jacobian.m_n);
  //  pR tmpx = (pR)malloc(sizeof(R)*jacobian.m_n);
  //  pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
  I n1=1;
  dcopy(&jacobian.m_n,b_,&n1,tmpb,&n1);
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      for (I j = 1;j<jacobian.m_size_block;++j)
	{
	  tmpb[ielm*jacobian.m_size_block+j]=0.0;
	}
    }
  
  solve_gauss_seidel_select(jacobian, 0, degree, x_, tmpb);
#if 0
  dcopy(&jacobian.m_n,b_,&n1,tmpb,&n1);
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      for (I j = 3;j<jacobian.m_size_block;++j)
	{
	  tmpb[ielm*jacobian.m_size_block+j]=0.0;
	}
    }
  solve_gauss_seidel_select(jacobian,1,degree,x_,tmpb);
#endif
  //  solve_gauss_seidel(jacobian,degree,x_,b_);
}

struct DG
{  

  
public:  
  typedef enum enum_DG
    {
    
      IA_lc=0,
      IA_lc_elm=1,
      IA_lc_face=2,
      IA_lc_elm_a=3,
      IA_lc_u=4,
      IA_lc_elm_u0=5,
      IA_lc_elm_u1=6,
      I_npoints=11,
      I_npoints_element=12,
      I_npoints_boundary=13,
      RA_lc=14,
      I_lc_len=16,
      RA_bmat=17,
      RA_bmatx=18,
      I_bmat_n=19,
      I_bmat_m=20,
      I_bmat_len=21,
      ERR=22,
      RA_EVAL_TETA_U=23,

      I_QFACE_N=28,
      I_QELM_N=29,

      I_TRIAL_FAMILY=30,
      I_TRIAL_DEGREE=31,
      I_TRIAL_NBASIS=32,

    
      I_TEST_FAMILY=33,
      I_TEST_DEGREE=34,
      I_TEST_NBASIS=35,
    
      I_TETA_FAMILY=36,
      I_TETA_DEGREE=37,
      I_TETA_NBASIS=38,
    
      I_TETA_U_FAMILY=39,
      I_TETA_U_DEGREE=40,
      I_TETA_U_NBASIS=41,
    
      I_TETA_A_FAMILY=42,
      I_TETA_A_DEGREE=43,
      I_TETA_A_NBASIS=44,

      I_n=64

    } info_t;

struct INFO
{
  I 		dg_iinfo[DG::I_n];
  R		dg_rres	[DG_rres_n];
  I		dg_ires	[DG_ires_n];
  I 		dg_rinfo_n;
  I 		dg_iinfo_n;
  pR		dg_rinfo;
  I 		dg_rwork_n;
  pR		dg_rwork;
  I 		dg_iwork_n;
  pI		dg_iwork;
  INFO()
  {
    dg_iinfo[DG::I_QELM_N]	= 10;
    dg_iinfo[DG::I_QFACE_N]  	= 10;  
    dg_rinfo_n			= (I)640000;
    dg_rinfo 			= (pR)malloc(sizeof(R)*dg_rinfo_n);
    dg_rwork_n			= (I)1280000;
    dg_rwork 			= (pR)malloc(sizeof(R)*dg_rwork_n);
    dg_iwork_n			= (I)1280000;
    dg_iwork 			= (pI)malloc(sizeof(I)*dg_iwork_n);
  };
  
  static void define(INFO*self)
  {
    printf("I_QFACE_N %d\n",DG::I_QFACE_N);
    self->dg_iinfo_n = DG::I_n;
    self->dg_iinfo[DG::I_QELM_N]	= 10;
    self->dg_iinfo[DG::I_QFACE_N]  	= 10;  
    self->dg_rinfo_n			= (I)640000;
    self->dg_rinfo 			= (pR)malloc(sizeof(R)*self->dg_rinfo_n);
    self->dg_rwork_n			= (I)1280000;
    self->dg_rwork 			= (pR)malloc(sizeof(R)*self->dg_rwork_n);
    self->dg_iwork_n			= (I)1280000;
    self->dg_iwork 			= (pI)malloc(sizeof(I)*self->dg_iwork_n);
  };
  
};
  
public: static DG_DATA * create_data(cst_mkS 		s_teta_a_,
				     cst_mkS 		s_teta_u_,
				     cst_mkS 		s_teta_,
				     cst_mkS 		s_test_,
				     cst_mkS 		s_trial_)
  {
    return new DG_DATA(mkS_n(s_teta_),
		       mkS_n(s_trial_),
		       mkS_n(s_test_),
		       mkS_n(s_teta_u_));
  };
  
public: static void dgadvection_init_with_quadrature(cst_mkS 		s_teta_a_,
						     cst_mkS 		s_teta_u_,
						     cst_mkS 		s_teta_,
						     cst_mkS 		s_test_,
						     cst_mkS 		s_trial_,
						     
						     cst_pI 		qelm_n_,
						     cst_pR 		qelm_p_,
						     cst_pR 		qelm_w_,
						     
						     cst_pI 		qface_n_,
						     cst_pR 		qface_p_,
						     cst_pR 		qface_w_,
						     
						     cst_pI           	iinfo_n_,
						     pI			iinfo_,
						     cst_pI		rinfo_n_,
						     pR 		rinfo_,
						     cst_pI 		rwork_n_,
						     pR			rwork_)
  {
  
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_n_[0] < DG::I_n)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_N]<1\n");
	exit(1);
      }
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }
  
    iinfo_[DG::I_TETA_A_NBASIS]  = (s_teta_a_) ? mkS_n(s_teta_a_):((I)0);
    iinfo_[DG::I_TETA_U_NBASIS]	 = mkS_n(s_teta_u_);
    iinfo_[DG::I_TETA_NBASIS]	 = mkS_n(s_teta_);
    iinfo_[DG::I_TRIAL_NBASIS]	 = mkS_n(s_trial_);
    iinfo_[DG::I_TEST_NBASIS]	 = mkS_n(s_test_);
    iinfo_[DG::I_QELM_N] 	 = qelm_n_[0];  

    if (ferr_[0]) 
      {
	return;
      }
  
    if (iinfo_[DG::I_TRIAL_NBASIS]!=iinfo_[DG::I_TEST_NBASIS])   
      { 
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_TEST_NBASIS]!=IINFO[DG::I_TRIAL_NBASIS]\n");    
      }

    mkS_st s_teta_a;
    mkS_st s_teta_u;
    mkS_st s_teta;
    mkS_st s_test;
    mkS_st s_trial;

    mkS_copy(&s_teta_a,s_teta_a_);
    mkS_copy(&s_teta_u,s_teta_u_);
    mkS_copy(&s_teta,s_teta_);
    mkS_copy(&s_test,s_test_);
    mkS_copy(&s_trial,s_trial_);

    const I s_teta_a_degree = mkS_k(s_teta_a_);
    const I s_teta_u_degree = mkS_k(s_teta_u_);

    const I trial_n	= iinfo_[DG::I_TRIAL_NBASIS];
    const I test_n	= iinfo_[DG::I_TEST_NBASIS];
    const I teta_n	= iinfo_[DG::I_TETA_NBASIS];
    const I teta_a_n 	= iinfo_[DG::I_TETA_A_NBASIS];
    const I teta_u_n 	= iinfo_[DG::I_TETA_U_NBASIS];

    fprintf(stdout,"  Discontinuous Galerkin Definition\n");
    fprintf(stdout,"  trial_n   " ifmt "\n",trial_n);
    fprintf(stdout,"  test_n    " ifmt "\n",test_n);
    fprintf(stdout,"  teta_n    " ifmt "\n",teta_n);
    fprintf(stdout,"  teta_a_n  " ifmt "\n",teta_a_n);
    fprintf(stdout,"  teta_u_n  " ifmt "\n",teta_u_n);
  
    const I trial_nXtrial_n 	= trial_n * trial_n;

    
    const I npoints_element	= teta_a_n + dim * teta_u_n;
    const I npoints_boundary	= nfaceinelm * qface_n_[0];
    const I npoints 	 	= npoints_element + npoints_boundary;  
    const I bmat_n 		= trial_n * test_n + teta_n * test_n;
    //    const I bmat_m 		= npoints;
    
    const I bmatx_size 		= bmat_n*(nfaceinelm * qface_n_[0] * nfaceinelm);
    const I evalu_size 		= teta_u_n * (dim*teta_u_n + nfaceinelm * qface_n_[0]);  
    if (rinfo_n_[0] < 2 * npoints + bmat_n*npoints + bmatx_size + evalu_size)
      {
	fprintf(stderr,"*** DGI too small rinfo_ array (" ifmt "<" ifmt ")\n",
		rinfo_n_[0],
		2*npoints+bmat_n*npoints+bmatx_size+evalu_size);
	ferr_[0] = DGERR_MEMORY;
	return;
      }
  
    /* on recopie la quadrature pour les faces */  
    /* on recopie les interpolants */
    
    iinfo_[DG::I_npoints]     		= npoints;
    iinfo_[DG::I_npoints_element]  	= npoints_element;
    iinfo_[DG::I_npoints_boundary] 	= npoints_boundary;
    
    iinfo_[DG::I_bmat_n] 		= bmat_n;
    iinfo_[DG::I_bmat_m] 		= iinfo_[DG::I_npoints];
    iinfo_[DG::I_bmat_len] 		= iinfo_[DG::I_bmat_n]*iinfo_[DG::I_npoints];

    iinfo_[DG::I_lc_len] 		= 2*iinfo_[DG::I_npoints];
  
    /* matrice de flux */  
    /* les points d evaluation en coordonnees locales */
    iinfo_[DG::IA_lc]        	= 0;
    iinfo_[DG::IA_lc_elm]    	= 0;
    iinfo_[DG::IA_lc_elm_a]  	= 0;
    iinfo_[DG::IA_lc_elm_u0] 	= iinfo_[DG::I_TETA_A_NBASIS];
    iinfo_[DG::IA_lc_elm_u1] 	= iinfo_[DG::I_TETA_A_NBASIS]+iinfo_[DG::I_TETA_U_NBASIS];
    iinfo_[DG::IA_lc_face]   	= npoints_element;
    iinfo_[DG::IA_lc_u]      	= iinfo_[DG::I_TETA_A_NBASIS];
  
    iinfo_[DG::RA_bmat] 		= 2*npoints;
    iinfo_[DG::RA_bmatx] 		= iinfo_[DG::RA_bmat]  + iinfo_[DG::I_bmat_len];
    iinfo_[DG::RA_EVAL_TETA_U]    = iinfo_[DG::RA_bmatx] + bmatx_size;

    

    const I npoints_velocity_involved 	= iinfo_[DG::I_npoints]-iinfo_[DG::I_TETA_A_NBASIS];  
    
    pR bmat_x = &rinfo_[iinfo_[DG::RA_bmat]];  
    const I bmat_ld = bmat_n;
    pR teta_a_trial_test 	= &bmat_x[0];

    pR teta_u_nabla_trial_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_trial_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld];
      }
    
    //    pR teta_u_dxtrial_test 	= &bmat_x[teta_a_n*bmat_ld];
    //    pR teta_u_dytrial_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_ld];
    pR teta_u_nabla_teta_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_teta_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld+trial_nXtrial_n];
      }

    pR teta_a_teta_test  	= &bmat_x[trial_nXtrial_n];
    //    pR teta_u_dxteta_test 	= &bmat_x[teta_a_n*bmat_n+trial_nXtrial_n];
    //    pR teta_u_dyteta_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_n+trial_nXtrial_n];
    
    if (teta_a_n>0)
      {
	
	//
	// Form the matrix.
	//
	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_trial),
		mkS_b(&s_test),
		teta_a_trial_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);

#if 0
	I i;for (i=0;i<36;++i) { printf("hh %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }      
#endif

	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_teta),
		mkS_b(&s_test),
		teta_a_teta_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
#if 0
	for (i=0;i<36;++i) { printf("gg %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }
#endif
      
      }

  
    if (teta_u_n>0)
      {
	
	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_trial,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_trial_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }

	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_teta,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_teta_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }
      }  

    /* -------------------------------------------------------------------------------------- */
    /* LOCAL COORDINATES -------------------------------------------------------------------  */
    /* -------------------------------------------------------------------------------------- */            
    { pR lc = &rinfo_[iinfo_[DG::IA_lc]];

      if (iinfo_[DG::I_TETA_A_NBASIS]>0)
	{
	  mkS_lagrange_localspl_tria(&s_teta_a_degree,
				     &lc[iinfo_[DG::IA_lc_elm_a]],
				     &npoints);
	}
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u0]],
					 &npoints);
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u1]],
					 &npoints);
    
      mkS_bmapping		(qface_n_[0],
				 &lc[iinfo_[DG::IA_lc_face]],
				 &npoints,
				 qface_p_);
    
      /* -------------------------------------------------------------------------------------- */
      /* CALCUL DE L EVALUATION DE TETA_U POUR TOUS LES POINTS D INTEGRATION ELM + 3*NFACES */
      /* -------------------------------------------------------------------------------------- */            
      mkS_basis(mkS_b(&s_teta_u),
		&npoints_velocity_involved,
		&rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		&iinfo_[DG::I_TETA_U_NBASIS],
		&rinfo_[iinfo_[DG::IA_lc_u]],
		&iinfo_[DG::I_npoints],// DG::I_npoints			   			   
		rwork_,
		rwork_n_,
		ferr_);  }
  
#if 0
    {
      I i,j;pR lc = &rinfo_[iinfo_[DG::IA_lc]];
      for (i=0;i<npoints_velocity_involved;++i)
	{
	  printf("allo %e %e " ifmt "\n",lc[iinfo_[DG::IA_lc_face]+i],lc[iinfo_[DG::IA_lc_face]+npoints+i],iinfo_[DG::I_TETA_U_NBASIS]);
	  for (j=0;j<iinfo_[DG::I_TETA_U_NBASIS];++j)
	    {
	      printf("%e\n",rinfo_[iinfo_[DG::RA_EVAL_TETA_U] + iinfo_[DG::I_TETA_U_NBASIS]*i+j]);
	    }
	}
      exit(1);}
#endif


    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD */
    /* -------------------------------------------------------------------------------------- */        

    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_trial),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n],
	     bmat_n,
	     ferr_);
  
    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_teta),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n+trial_nXtrial_n],
	     bmat_n,
	     ferr_);

  
    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD CONTRAPOSEE AVEC LE BORD */
    /* -------------------------------------------------------------------------------------- */            

    { pR bmatflux = &rinfo_[iinfo_[DG::RA_bmatx]];

      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_trial),
		   mkS_b(&s_test),
		   bmatflux,
		   bmat_n,
		   ferr_);    
    
      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_teta),
		   mkS_b(&s_test),
		   &bmatflux[trial_nXtrial_n],
		   bmat_n,
		   ferr_); }
    return;
#if 0
    I s_teta_a[4]      = {iinfo_[DG::I_TETA_A_FAMILY],0,iinfo_[DG::I_TETA_A_DEGREE],iinfo_[DG::I_TETA_A_NBASIS]};
    I s_teta_u[4]      = {iinfo_[DG::I_TETA_U_FAMILY],0,iinfo_[DG::I_TETA_U_DEGREE],iinfo_[DG::I_TETA_U_NBASIS]};
    I s_teta[4]        = {iinfo_[DG::I_TETA_FAMILY],0,iinfo_[DG::I_TETA_DEGREE],iinfo_[DG::I_TETA_NBASIS]};
    I s_test[4]        = {iinfo_[DG::I_TEST_FAMILY],0,iinfo_[DG::I_TEST_DEGREE],iinfo_[DG::I_TEST_NBASIS]};
    I s_trial[4]       = {iinfo_[DG::I_TRIAL_FAMILY],0,iinfo_[DG::I_TRIAL_DEGREE],iinfo_[DG::I_TRIAL_NBASIS]};
#endif
  };

  static void define(DG_HANDLE* 	handle_,
		     cst_mkS 		s_teta_a_,
		     cst_mkS 		s_teta_u_,
		     cst_mkS 		s_teta_,
		     cst_mkS 		s_test_,
		     cst_mkS 		s_trial_,
		     cst_pI         	iinfo_n_,
		     pI			iinfo_,		  
		     cst_pI		rinfo_n_,
		     pR 		rinfo_,
		     cst_pI 		rwork_n_,
		     pR			rwork_)
  {
    
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }

    const I qelm_n_	[1] 	= {iinfo_[DG::I_QELM_N]*iinfo_[DG::I_QELM_N]};
    const I qface_n_	[1] 	= {iinfo_[DG::I_QFACE_N]};
    
    R qelm_p_	[64*64*2];
    R qelm_w_	[64*64];
    R qface_p_	[64];
    R qface_w_	[64];
  
    fprintf(stdout,"generate quadrature " ifmt " " ifmt "\n",
	    qelm_n_[0],
	    qface_n_[0]);
  
    {
      const I rwork_n = (rwork_n_[0]>3*iinfo_[DG::I_QELM_N])?rwork_n_[0]-3*iinfo_[DG::I_QELM_N]:(I)0;
      if (rwork_n<1)
	{
	  fprintf(stderr,"*** MOK_DGI_INIT ERROR TOO SMALL RWORK_N\n");
	  ferr_[0] = DGERR_MEMORY;
	  return;
	}
      mkQ_legendre_interval_(qface_n_,
			     qface_p_,
			     &negal1,
			     qface_w_,
			     &negal1,
			     rwork_,
			     rwork_n_,
			     ferr_);
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:first mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      mkQ_legendre_interval_(&iinfo_[DG::I_QELM_N],
			     rwork_,
			     &negal1,
			     &rwork_[2*iinfo_[DG::I_QELM_N]],
			     &negal1,
			     &rwork_[3*iinfo_[DG::I_QELM_N]],
			     rwork_n_,
			     ferr_);    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:second mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      I qelm_n=0;
      mkQ_collapse_		(&iinfo_[DG::I_QELM_N],
				 rwork_,
				 &negal1,
				 &rwork_[2*iinfo_[DG::I_QELM_N]],
				 &negal1,
				 &qelm_n,
				 qelm_p_,
				 qelm_n_,
				 qelm_w_,
				 &negal1,
				 ferr_);
    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:mkQ_legendre_interval_ failed\n");
	  return;
	}
    }
    fprintf(stdout,"generate quadrature done\n");



    dgadvection_init_with_quadrature(s_teta_a_,
				     s_teta_u_,
				     s_teta_,
				     s_test_,
				     s_trial_,
		  
				     qelm_n_,
				     qelm_p_,
				     qelm_w_,
		  
				     qface_n_,
				     qface_p_,
				     qface_w_,

				     iinfo_n_,
				     iinfo_,

				     rinfo_n_,
				     rinfo_,

				     rwork_n_,
				     rwork_);


    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	na      = &iinfo_[DG::I_TETA_A_NBASIS];
    const I bmat_n 		= iinfo_[DG::I_bmat_n];  
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  


    
    pR brhs      		= &rwork_[0];
    pR brhs_uelm 		= &rwork_[1];
    I  npts_involving_u 	= bmat_m-na[0];/*nTot-teta_a_n*/    
    pR tmpbrhs        		= &rwork_[bmat_m]; 
    pR tmpbrhs_uelm   		= &tmpbrhs[na[0]]; 
    pR tmpbrhs_u 		= tmpbrhs_uelm;
    const I tmpbrhs_ufaceoff 	= bmat_m;
    
    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    //    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    matrix_handle_def(&handle_->m_BMAT,
		      bmat_n,
		      bmat_m,(pR)bmat,
		      bmat_n);
    

    matrix_handle_def(&handle_->m_mat_tmpbrhs_uelm,
		      nu[0],
		      dim,
		      tmpbrhs_uelm,
		      tmpbrhs_ufaceoff);

    matrix_handle_def(&handle_->m_UVWDOFS,
		      nu[0],
		      dim,
		      brhs_uelm,
		      nu[0]);

    matrix_handle_def(&handle_->m_brhs_uvw,
		      npts_involving_u,
		      dim,
		      tmpbrhs_u,
		      tmpbrhs_ufaceoff);
    
    vector_handle_def(&handle_->m_BRHS,
		      bmat_m,
		      brhs,
		      1);

    matrix_handle_def(&handle_->m_EVALU,
		      nu[0],
		      npts_involving_u,
		      &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		      iinfo_[DG::I_TETA_U_NBASIS]);
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {    
	matrix_handle_def(&handle_->m_mat_tmpbrhs_uface[localFaceIndex],
			  qface_n_[0],
			  dim,
			  &tmpbrhs[na[0]+dim*nu[0]+localFaceIndex*qface_n_[0]],
			  tmpbrhs_ufaceoff);
      }
    
    
  };


  
  static void solve(DG_HANDLE*  handle_,
		    DG_DATA*    data_,
		    cst_pR 	xa_,
		    cst_pR 	xu_,
		    cst_pR	rhs_,
		    cst_pI 	rhsoff_,
		    
		    cst_pI  	cnc_u_,
		    cst_pI  	cncoff_u_,
		    cst_pR 	data_u_,		       
		    cst_pR	data_v_,
		    
		    cst_pR	sol_,
		    cst_pI 	soloff_,
		    
		    pR		corr_,
		    cst_pI 	corroff_,
		    
		    cst_pR 	t_,
		    cst_pI 	nelm_,
		    cst_pR 	coo_,
		    cst_pI 	cooff_,
		    cst_pI 	cnc_,
		    cst_pI 	cncoff_,
		    cst_pI 	adj_,
		    cst_pI 	adjoff_,
		    cst_pI 	vcod_,
		    cst_pI 	noboundary_cod_,
		    
		    cst_pI 	rwork_n_,
		    pR  	rwork_,
		    cst_pI 	iwork_n_,
		    pI	iwork_,
		    
		    cst_pR 	rinfo_,
		    const I  	iinfo_[DG::I_n],
		    R 		rres_[DG_rres_n],
		    I 		ires_[DG_ires_n])
  {
    
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    cst_pR bmat 	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1,
      n2 	= (I)2;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/

    I marker_size	= 0;
  
    pI
      lcperm 	= NULL,
      perm 	= NULL,
      graph 	= NULL,
      marker 	= NULL,
      marker_begin= NULL,
      marker_flag = NULL;

    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  
    R mx;  
    I
      id,
      ielm,
      iter_gauss_seidel,
      b1,
      b2,
      pp,
      n,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R
      jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1],
      lcsol[21];
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*nu[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*nu[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/



    //
    //
    //

#if 0
    R 	matflux[21*(21+21)];/*trial_n*(trial_n+teta_n)*/
    pR 		matflux_corr 	= &matflux[0];
    pR 		matflux_residu 	= &matflux[trial_n[0]*trial_n[0]];

    matrix_handle_def(&handle_->m_FC,trial_n[0],trial_n[0],matflux_corr,trial_n[0]);
    matrix_handle_def(&handle_->m_FR,trial_n[0],teta_n[0],matflux_residu,trial_n[0]);
#endif

    // R lcrhs[21*2];
    const I renum_iwork_n = 3*nelm_[0] + 4*(nelm_[0]+1) + trial_n[0];
    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
  
    const I tmpbrhs_ufaceoff 	= bmat_m;
    pR tmpbrhs        	= &rwork_[bmat_m]; 
    pR tmpbrhs_ufaces 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface0 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface1 	= &tmpbrhs[1+2*nu[0]+qface_n[0]];
    pR tmpbrhs_uface2 	= &tmpbrhs[1+2*nu[0]+2*qface_n[0]];
    
    //  pR tmpbrhs_uelm   	= &tmpbrhs[1]; 
    //  pR tmpbrhs_u	= tmpbrhs_uelm;
    /*_____________________________________________________*/
    if (renum_iwork_n>iwork_n_[0])
      {
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_iw_n] = renum_iwork_n;
	fprintf(stderr,"too small iwork_n_ " ifmt " " ifmt "\n",renum_iwork_n,iwork_n_[0]);
	return;
      }

    lcperm  	= &iwork_[0];
    perm  	= &iwork_[trial_n[0]];
    graph 	= &iwork_[(nelm_[0]+1) + trial_n[0] ];
    marker 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]];
    marker_begin 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ (nelm_[0]+1)];
    marker_flag 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ 2*(nelm_[0]+1)];
    /*_____________________________________________________*/

    goto __state_renum;

    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  
  __state_renum:
    b1=1;
    b2=0;    
    { I j;for (j=trial_n[0];j<renum_iwork_n;++j) iwork_[j]=(I)0;}  
    /* */
    /* on donne un sens a l adjacence */
    { I jelm;
      for (jelm=0;jelm<nelm_[0];++jelm)
	{        
	  /* compute jaface and nrmelm */
	  neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];
	  cnc[0]     = cnc_[jelm*cncoff_[0]+0];cnc[1]     = cnc_[jelm*cncoff_[0]+1];cnc[2]     = cnc_[jelm*cncoff_[0]+2];
	  cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	  cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];

	  jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	  jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	  jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  

	  longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];  
	  data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];				
	  x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;						  
	  data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	  x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;  
	  data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	  x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	  jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	  /*____ INFO ELEMENT DONE */ 	
	  /*___________________________________________________________________________________________________________________*/  
	  /* 
	     eval u faces
	  */      
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
  
	  /*___________________________________________________________________________________________________________________*/
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);


	  
	  Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	
	  { I k;
	    for (k=0;k<qface_n[0];++k) if (uface0[k]<((R)0.0)) break;
	    graph[jelm*3+0] = (neids[0]) ? ( (k<qface_n[0])?-neids[0]:neids[0] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	  { I k;
	    for (k=0;k<qface_n[0];++k)if (uface1[k]<((R)0.0))break;
	    graph[jelm*3+1] = (neids[1]) ? ( (k<qface_n[0])?-neids[1]:neids[1] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	  { I k;
	    for (k=0;k<qface_n[0];++k)if (uface2[k]<((R)0.0))break;
	    graph[jelm*3+2] = (neids[2]) ? ( (k<qface_n[0])?-neids[2]:neids[2] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }	
	} }
    /* on met tous ceux qui n'ont pas d antecedent en premier*/
    { I jelm;
      for (jelm=0;jelm<nelm_[0];++jelm)
	{ 
	  if (graph[jelm*3]<-nelm_[0]) {perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	  if (graph[jelm*3+1]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	  if (graph[jelm*3+2]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	} }
    /*________________________________________________________________________________________________________________*/  
    pp	= b2;
    n	= (I)1;
  __state_renum_do1:
    marker_size=0;
    /* INJECTION */
    n=0;
    { I i;	
      for (i=b1;i<=b2;++i)
	{
	  { I j;
	    for (j=0;j<3;++j)
	      {
		const I s = graph[(perm[i]-1)*3+j];
		if ( (s>0) AND ( (NOT marker[s]) AND (NOT marker_flag[s]) ))
		  {
		    marker_size=marker_size+1; marker_begin[marker_size] = s;
		    marker[s] = marker_size;n=n+1;
		  }
	      } }
	} }
    
    /* FIN INJECTION */
    { I i;	    
      for (i=1;i<=marker_size;++i)
	{
	  marker_flag[marker_begin[i]]=1;perm[b2+i]=marker_begin[i];
	  marker[marker_begin[i]]=0;marker_begin[i]=0;
	} }
  
    pp+=n;b1 = b2+1;b2 = b2+n;
    if (b1>b2)
      {
	if (pp==nelm_[0])
	  goto __state_renum_do1_done;
	else
	  {
	    /* on cherche un element*/
	    { I jelm;
	      for (jelm=1;jelm<=nelm_[0];++jelm)
		{
		  if (NOT marker_flag[jelm])
		    {
		      marker_flag[jelm]=(I)1;
		      perm[b1] = jelm;
		      ++n;b2 = b1;
		      break;
		    }
		} }
	  }
      }
    if (n>0)
      {
	/* on continue a injecter */
	goto __state_renum_do1;
      }
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_renum_do1_done:
    if (b2!=nelm_[0])
      {      
	fprintf(stderr,"*** DGERR RENUM FAILED " ifmt "!=" ifmt "",b2,nelm_[0]);
	goto __state_error;
      }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    iter_gauss_seidel = 0;
  __state_next_iter_gauss_seidel:

    fprintf(stdout,"next iter gauss seidel \n");
    mx = (R)0.0;
    /*___________________________________________________________________________________________________________________*/  
    ielm=0;      
  __state_next_elm:
    id 		= perm[ielm+1]-1;
    neids[0]	= adj_[id*adjoff_[0]];
    neids[1]	= adj_[id*adjoff_[0]+1];
    neids[2]	= adj_[id*adjoff_[0]+2];

    cnc[0]        = cnc_[id*cncoff_[0]+0];
    cnc[1]        = cnc_[id*cncoff_[0]+1];
    cnc[2]        = cnc_[id*cncoff_[0]+2];

#if 0    
    vcod[0] 	= vcod_[cnc[0]-1];
    vcod[1] 	= vcod_[cnc[1]-1];
    vcod[2] 	= vcod_[cnc[2]-1];
    codface[0] 	= (neids[0])? noboundary_cod_[0] :  MAX(vcod_[0],vcod_[1]);
    codface[1] 	= (neids[1])? noboundary_cod_[0] :  MAX(vcod_[1],vcod_[2]);
    codface[2] 	= (neids[2])? noboundary_cod_[0] :  MAX(vcod_[2],vcod_[0]);
#endif
    
    neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
    /*---*/
    if (neids[0])
      {	if (neids[0]<0) neids[0]=-neids[0];
	if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
    if (neids[1])
      {	if (neids[1]<0) neids[1]=-neids[1];
	if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
    if (neids[2])
      {	if (neids[2]<0) neids[2]=-neids[2];
	if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
    /*---*/
    cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
    cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
    cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
    cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
    cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
    cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
    jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
    jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
    jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
    longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
    jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
    data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
    data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
    data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
    x = data_->mat_belm.x[2];
    data_->mat_belm.x[2]=data_->mat_belm.x[1];
    data_->mat_belm.x[1]=x;
    
    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
    /*____ INFO ELEMENT DONE */ 
    /*____ COPY LOCAL RHS */ 
    {I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	


    if (data_a_)
      {
	for (I k=0;k<teta_a_n[0];++k)
	  {
	    brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	  }
      }
    else
      {
	brhs_a[0] = jacelm[0] * xa_[0];
      }
      
  
    /*___________________________________________________________________________________________________________________*/  
    /* 
       eval u : elm and faces
    */  
    /*___________________________________________________________________________________________________________________*/  

    //
    // Copy the velocity.
    //
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }

    //
    // Wrap the solution.
    //
    vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
  
    //
    // Evaluate u at all the required coordinates.
    //
    handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
    
    /*____ TRANSFORM ON UELM */
    handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

    //
    // On the boundary.
    //
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
      }
    
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex].apply([xu_](const R r)
					       {
						 return (r<0.0) ? -r*xu_[0] : ((R)0.0);
					       });
      }
      
    /*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	if (neids[localFaceIndex])
	  {
	    Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);

	    
	    vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
	    vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);

	    
	    data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
	    data_->m_local_rhs += data_->m_flux_part_sol  * data_->vec_neisol;	 
	  }
      }
  
  
    /*--- COMPUTE MATRICES  */
    data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
  
    /*--- COMPUTE LOCAL RESIDUAL  */
    data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
  
    /*--- COMPUTE LOCAL CORRECTION   */
#if 0
    std::cout << "ielm " << ielm << std::endl;
    matrix_handle_print(&data_->m_local_matrix_part_corr,stdout);
    std::cout << " rhs " << ielm << std::endl;
    vector_handle_print(&data_->m_local_rhs,stdout);
#endif    
    matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
  
    /*--- COMPUTE NRMS OF LOCAL CORRECTION */
    { I j; for (j=0;j<test_n[0];++j) {lcsol[j] = corr_[corroff_[0]*id+j];}}
    { I j;for (j=0;j<test_n[0];++j) corr_[corroff_[0]*id+j] = data_->m_local_rhs.x[j];}
    { I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j]-=lcsol[j];}
    { R xx=(R)0.0; {I j;for (j=0;j<test_n[0];++j) xx+=data_->m_local_rhs.x[j]*data_->m_local_rhs.x[j];}
      if (mx<xx) mx = xx; }		    
  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    ++ielm;
    if (ielm<nelm_[0])
      {
	goto  __state_next_elm;
      }  
    else
      {
	goto  __state_next_elm_done;
      }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_next_elm_done:
    {
      fprintf(stdout,"*** DGINFO time %.5f iter " ifmt "/" ifmt " tol = %8.15e\n",t_[0],iter_gauss_seidel,nelm_[0],mx);
      if ( (mx>((R)1.0e-24)) AND (iter_gauss_seidel<nelm_[0]) )
	{
	  ++iter_gauss_seidel;
	  goto __state_next_iter_gauss_seidel;
	}
      else
	{	
	  goto __state_computation_done;
	}
    }
    
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
  __state_computation_done:
    {    
      ires_[DG_ires_iter_gauss_seidel] 	= iter_gauss_seidel;
      rres_[DG_rres_max] 			= mx;
      if (iter_gauss_seidel<nelm_[0])
	{
	  ires_[DG_ires_convergence] = (I)1;
	}
      else
	{
	  ires_[DG_ires_convergence] = (I)0;
	}

      rres_[DG_rres_areaL1]	= 0.0; 
      rres_[DG_rres_nrmL2]	= 0.0; 
      rres_[DG_rres_nrmLInf]	= 0.0; 
      rres_[DG_rres_jumpL2]	= 0.0; 
      rres_[DG_rres_johnson]	= 0.0;   
      { I jelm;
	for (jelm=0;jelm<nelm_[0];++jelm)
	  { 
	    /**/
	    neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];

	    cnc[0]        = cnc_[jelm*cncoff_[0]+0];
	    cnc[1]        = cnc_[jelm*cncoff_[0]+1];
	    cnc[2]        = cnc_[jelm*cncoff_[0]+2];	  
	    cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	    cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	    cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	    cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	    cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	    cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
	    jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	    jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	    jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	    longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	    jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  
	    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];data_->mat_belm.x[1] = cooelm[0]-cooelm[2];data_->mat_belm.x[2] = cooelm[3]-cooelm[4];data_->mat_belm.x[3] = cooelm[1]-cooelm[0];  	  
	    x = data_->mat_belm.x[2];data_->mat_belm.x[2]=data_->mat_belm.x[1];data_->mat_belm.x[1]=x;	  
	    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	    /**/	  

	    Blas_dgemv("N",trial_n,trial_n,jacelm,bmat,trial_n,&corr_[jelm*corroff_[0]],&n1,&r0,data_->m_local_rhs.x,&n1);

	    rres_[DG_rres_nrmL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);

	    { I j;
	      for (j=0;j<trial_n[0];++j)
		{
		  rres_[DG_rres_areaL1] += data_->m_local_rhs.x[j];
		  if (rres_[DG_rres_nrmLInf]<nsFABS(corr_[corroff_[0]*jelm+j]))
		    {
		      rres_[DG_rres_nrmLInf]=nsFABS(corr_[corroff_[0]*jelm+j]);
		    }
		} }

	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	  
	    /* u au points de gauss de l element */
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);
	  
	    /* u.n sur les faces */
	    Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	  
	    /* mise en a zero */
	    { I j;for (j=0;j<qface_n[0];++j) uface0[j] = (uface0[j]<0.0) ? uface0[j]*(-xu_[0]) : (R)0.0;}       
	    { I j;for (j=0;j<qface_n[0];++j) uface1[j] = (uface1[j]<0.0) ? uface1[j]*(-xu_[0]) : (R)0.0;}	
	    { I j;for (j=0;j<qface_n[0];++j) uface2[j] = (uface2[j]<0.0) ? uface2[j]*(-xu_[0]) : (R)0.0;}
	  	  
	    { I k;
	      for (k=0;k<3;++k)
		{
		  if ( neids[k] )
		    { 
		      /* 
			 on forme la matrice de flux pour l'arete
		      */
		      Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[k][neids_face[k]],&bmat_n,(k==2)?uface2:((k==1)?uface1:uface0),&n1,&r0,data_->m_flux_memory,&n1);
		      
		      /* 
			 on calcule le produit scalaire
		      */
		      Blas_dgemv("N",trial_n,trial_n,&jacface[k],data_->m_flux_part_corr.x,trial_n,&corr_[(neids[k]-1)*corroff_[0]],&n1,&r1,data_->m_local_rhs.x,&n1);	   
		      rres_[DG_rres_jumpL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);/* saut */	
		    } 
		} }
	  
	  } }

      printf("compute nrms\n"
	     "   areaL1  : %e\n"
	     "   nrmL2   : %e\n"
	     "   nrmLInf : %e\n"
	     "   jumpL2  : %e\n"
	     "   jhonson : %e\n",
	     rres_[DG_rres_areaL1],
	     rres_[DG_rres_nrmL2],
	     rres_[DG_rres_nrmLInf],
	     rres_[DG_rres_jumpL2],
	     rres_[DG_rres_johnson]);
    
      return;
    }  
    /*_________________________________________________________________________________________________________________*/    
  __state_error:
    {
      printf("out\n");
      ires_[DG_ires_err] = (I)1;
      return;
    }
  }

// x_{n+1} = (D-L)^{-1} * U * x_n + (D-L)^{-1} * b_;
// x_{n+1} - x_n = ( (D-L)^{-1} * U - I) x_n + (D-L)^{-1} * b_;





  static void solve_global(I jacdeg,
			   DG_JACOBIAN* jacobians[],
			   mkS shape_X,
			   mkS shape_E,
			   DG_VIEW& view,
			   DG_HANDLE*  handle_,
			   DG_DATA*    data_,
			   cst_pR 	xa_,
			   cst_pR 	xu_,
			   pR	rhs_,
			   cst_pI 	rhsoff_,
			   
			   cst_pI  	cnc_u_,
			   cst_pI  	cncoff_u_,
			   cst_pR 	data_u_,		       
			   cst_pR	data_v_,
			   
			   cst_pR	sol_,
			   cst_pI 	soloff_,
		    
			   pR		corr_,
			   cst_pI 	corroff_,
			   
			   cst_pR 	t_,
			   cst_pI 	nelm_,
			   cst_pR 	coo_,
			   cst_pI 	cooff_,
			   cst_pI 	cnc_,
			   cst_pI 	cncoff_,
			   cst_pI 	adj_,
			   cst_pI 	adjoff_,
			   cst_pI 	vcod_,
			   cst_pI 	noboundary_cod_,
			   
			   cst_pI 	rwork_n_,
			   pR  	rwork_,
			   cst_pI 	iwork_n_,
			   pI	iwork_,
			   
			   cst_pR 	rinfo_,
			   const I  	iinfo_[DG::I_n],
			   R 		rres_[DG_rres_n],
			   I 		ires_[DG_ires_n])
  {

    //    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];

#define JJ (*jacobians[jacdeg])
    
#if 0
    DG_JACOBIAN jacobian(nelm_[0],
			 nfaceinelm,
			 trial_n[0],
			 adj_,
			 adjoff_);
#endif    
    
    //    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    //    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];

    cst_pI	teta_n      = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	test_n      = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_u_n      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    //const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    //    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    //    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    //    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      //      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/


    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  

    I
      id,
      ielm,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1];
      
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*teta_u_n[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*teta_u_n[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/

    //
    //
    //
    // R lcrhs[21*2];

    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
      
    /*___________________________________________________________________________________________________________________*/  
    for (ielm = 0; ielm < nelm_[0];++ielm)
      {
	id = ielm;
	{I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	
	
	
	//	std::cout << "ielm " << ielm << std::endl;
	neids[0]	= adj_[id*adjoff_[0]];
	neids[1]	= adj_[id*adjoff_[0]+1];
	neids[2]	= adj_[id*adjoff_[0]+2];
	cnc[0]        = cnc_[id*cncoff_[0]+0];
	cnc[1]        = cnc_[id*cncoff_[0]+1];
	cnc[2]        = cnc_[id*cncoff_[0]+2];

	
	neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
	/*---*/
	if (neids[0])
	  {	if (neids[0]<0) neids[0]=-neids[0];
	    if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
	if (neids[1])
	  {	if (neids[1]<0) neids[1]=-neids[1];
	    if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
	if (neids[2])
	  {	if (neids[2]<0) neids[2]=-neids[2];
	    if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
	/*---*/
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
	jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

	data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
	data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
	data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
	data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
	x = data_->mat_belm.x[2];
	data_->mat_belm.x[2]=data_->mat_belm.x[1];
	data_->mat_belm.x[1]=x;
    
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  

	//
	//
	//
	if (data_a_)
	  {
	    for (I k=0;k<teta_a_n[0];++k)
	      {
		brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	      }
	  }
	else
	  {
	    brhs_a[0] = jacelm[0] * xa_[0];
	  }
  
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u : elm and faces
	*/  
	/*___________________________________________________________________________________________________________________*/  

	//
	// Copy the velocity.
	//
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	
	//
	// Evaluate u at all the required coordinates.
	//
	handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
	
	/*____ TRANSFORM ON UELM */
	handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

	//
	// On the boundary.
	//
	for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
	  }
    
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex].apply([xu_](const R r)
						   {
						     return (r<0.0) ? -r*xu_[0] : ((R)0.0);
						   });
	  }
	
	/*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    if (neids[localFaceIndex])
	      {
		//		printf("  flux edge " ifmt " \n",localFaceIndex);
		double mr1= -1.0;
		Blas_dgemv("N",&bmat_n,qface_n,&mr1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);
		// vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
		// vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);		
		JJ.addelm(id,neids[localFaceIndex]-1,0.0,data_->m_flux_part_corr.x,data_->m_flux_part_corr.ld);
		
		//
		//
		//
		//		data_->m_flux_part_corr;
		
		//
		//
		//
		//		data_->m_flux_part_sol;

		
		// data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
		vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);
		data_->m_local_rhs -= data_->m_flux_part_sol  * data_->vec_neisol;	 
	      }
	  }
	
	/*--- COMPUTE MATRICES  */
	{
	  // matrix_handle msub;
	  // vector_handle rsub;
	  
	  // matrix_handle_def(&msub,handle_->m_BMAT.n,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BMAT.x,handle_->m_BMAT.ld);
	  // vector_handle_def(&rsub,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BRHS.x,1);
	  
	  data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
	  // data_->m_local_matrices = msub * rsub;
	  JJ.addelm(id,id,0.0,data_->m_local_matrix_part_corr.x,data_->m_local_matrix_part_corr.ld);

	  vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
	  data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	  //	  printf("addelmcorr\n");
	}

	
	{I j;for (j=0;j<test_n[0];++j) rhs_[rhsoff_[0]*id+j] = data_->m_local_rhs.x[j];}	
	//	{I j;for (j=0;j<test_n[0];++j) printf(" sisis" ifmt " %e\n",j, rhs_[rhsoff_[0]*id+j]);}
	//
	// aE + u.nabla(E) + h(u,E) = -aF - u.nabla(F) - h(u,F)
	// JACOBIAN(trial,test) E = -R(F,test)
	// 
	// corr P3 F P2  Operator P2-P3 Operator P3-P3
	//
	// P0
	//   P1\P0
	//        P2\P0
	//
	//
	//
	//
	//
	//
	//
	//
	//

	//	data_->m_local_matrix_part_corr;
	/*--- COMPUTE LOCAL RESIDUAL  */
	// data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	/*--- COMPUTE LOCAL CORRECTION   */
	//	matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
	// matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
      }
    
    //    jacobian.spy("before.txt");
#if 0
    pSparse matrix = Sparse_build(JJ.m_n,
				  JJ.m_n,
				  JJ.m_nc,
				  JJ.m_begin,
				  JJ.m_index,
				  JJ.m_values);
#endif
    
    solve_gauss_seidel(JJ,
		       700,
		       mkS_k(shape_E),
		       corr_,
		       rhs_);

#if 0
#if 1
#if 1
    {
      printf("inverse diagonal\n");
      jacobian.inverse_diagonal();
      
      // Jacobi
      // corr_{n+1} = (diag + L){-1} * U * corr_{n} + (diag + L){-1} rhs
      // (I - invD * extra) corr = invD * rhs
      // step-1
      // rhs   = diag^{-1} * rhs

      printf("update rhs\n");
      jacobian.inverse_lower_gemv(rhs_,&jacobian.m_size_block);
      
      R r1=1.0;I n1=1;
      pR tmp = (pR)malloc(sizeof(R)*jacobian.m_n);
      pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
      printf("extra_gemv corr\n");
      
      for (I i=0;i<10;++i)
	{
	  dcopy(&jacobian.m_n,corr_,&n1,tmp2,&n1);
	  
	  jacobian.extra_gemv_upper(corr_,tmp);
	  //jacobian.extra_diagonal_gemv(corr,&jacobian.m_size_block);
	  //	  printf("diag_gemv corr\n");
	  jacobian.inverse_lower_gemv(corr_,&jacobian.m_size_block);
	  //	  printf("daxpy\n");
	  daxpy(&jacobian.m_n,&r1,rhs_,&n1,corr_,&n1);
	  double mr1=-1.0;
	  daxpy(&jacobian.m_n,&mr1,corr_,&n1,tmp2,&n1);
	  double h = dnrm2(&jacobian.m_n,tmp2,&n1);
	  double h0 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block];
	      h0 += mesh->jacelm[ielm]*x*x;
	    }
	  double h1 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+1];
	      double y = tmp2[ielm*jacobian.m_size_block+2];
	      h1 += mesh->jacelm[ielm]*(x*x + y*y);
	    }
	  double h2 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+3];
	      double y = tmp2[ielm*jacobian.m_size_block+4];
	      double z = tmp2[ielm*jacobian.m_size_block+5];
	      h2 += mesh->jacelm[ielm]*(x*x + y*y + z*z);
	    }
	  double h3 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+6];
	      double y = tmp2[ielm*jacobian.m_size_block+7];
	      double z = tmp2[ielm*jacobian.m_size_block+8];
	      double t = tmp2[ielm*jacobian.m_size_block+9];
	      h3 += mesh->jacelm[ielm]*(x*x + y*y + z*z + t*t);
	    }
	  
	  printf(" " ifmt " %e %e %e %e %e\n",i,h,sqrt(h0),sqrt(h1),sqrt(h2),sqrt(h3));
	  if (h < 1.0e-12)
	    {
	      break;
	    }
#if 0
	    dg_print_mesh( mesh,
			   view,
			   "dgjacobi." ifmt,i);  
	    dg_print_sol(view,
			 shape_X,
			 shape_E,
			 jacobian.m_nelm,
			 corr_,
			 "dgjacobi." ifmt,i);
#endif
	}
      //
      //
      // corr = extra * corr;
      // corr = diag^{-1} * corr
      // corr += rhs;
      //
      //
    }

#else
    {
      printf("inverse diagonal\n");
      jacobian.inverse_diagonal();
      
      // Jacobi
      // corr_{n+1} = diag^{-1} * extra * corr_{n} + diag^{-1} * rhs
      // (I - invD * extra) corr = invD * rhs
      // step-1
      // rhs   = diag^{-1} * rhs

      printf("update rhs\n");
      jacobian.inverse_diagonal_gemv(rhs_,&jacobian.m_size_block);
      
      R r1=1.0;I n1=1;
      pR tmp = (pR)malloc(sizeof(R)*jacobian.m_n);
      pR tmp2 = (pR)malloc(sizeof(R)*jacobian.m_n);
      printf("extra_gemv corr\n");
      
      for (I i=0;i<600;++i)
	{
	  dcopy(&jacobian.m_n,corr_,&n1,tmp2,&n1);
	  
	  jacobian.extra_gemv(corr_,tmp);
	  //jacobian.extra_diagonal_gemv(corr,&jacobian.m_size_block);
	  //	  printf("diag_gemv corr\n");
	  jacobian.inverse_diagonal_gemv(corr_,&jacobian.m_size_block);
	  //	  printf("daxpy\n");
	  daxpy(&jacobian.m_n,&r1,rhs_,&n1,corr_,&n1);
	  double mr1=-1.0;
	  daxpy(&jacobian.m_n,&mr1,corr_,&n1,tmp2,&n1);
	  double h = dnrm2(&jacobian.m_n,tmp2,&n1);
	  double h0 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block];
	      h0 += mesh->jacelm[ielm]*x*x;
	    }
	  double h1 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+1];
	      double y = tmp2[ielm*jacobian.m_size_block+2];
	      h1 += mesh->jacelm[ielm]*(x*x + y*y);
	    }
	  double h2 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+3];
	      double y = tmp2[ielm*jacobian.m_size_block+4];
	      double z = tmp2[ielm*jacobian.m_size_block+5];
	      h2 += mesh->jacelm[ielm]*(x*x + y*y + z*z);
	    }
	  double h3 =0.0;
	  for (I ielm=0;ielm<mesh->nelm;++ielm)
	    {
	      double x = tmp2[ielm*jacobian.m_size_block+6];
	      double y = tmp2[ielm*jacobian.m_size_block+7];
	      double z = tmp2[ielm*jacobian.m_size_block+8];
	      double t = tmp2[ielm*jacobian.m_size_block+9];
	      h3 += mesh->jacelm[ielm]*(x*x + y*y + z*z + t*t);
	    }
	  
	  printf(" " ifmt " %e %e %e %e %e\n",i,h,sqrt(h0),sqrt(h1),sqrt(h2),sqrt(h3));
	  if (h < 1.0e-12)
	    {
	      break;
	    }
#if 0	    
	    dg_print_mesh( mesh,
			   view,
			   "dgjacobi." ifmt,i);  
	    dg_print_sol(view,
			 shape_X,
			 shape_E,
			 jacobian.m_nelm,
			 corr_,
			 "dgjacobi." ifmt,i);
#endif
	}
      //
      //
      // corr = extra * corr;
      // corr = diag^{-1} * corr
      // corr += rhs;
      //
      //
    }
#endif
#else
    for (I i=0;i<=jacobian.m_n;++i) jacobian.m_begin[i]+=1;
    for (I i=0;i<jacobian.m_nc;++i) jacobian.m_index[i]+=1;    
    matrix->format = 1;
    {
      pExternPardiso pardiso = ExternPardiso_new();    

      std::cout << "precompute: symbolic phase ... " << std::endl;
      std::cout << "            n  = " << jacobian.m_n << std::endl;
      std::cout << "            nc = " << jacobian.m_nc << std::endl;
      ExternPardiso_precompute(pardiso,
			       matrix);
      
      std::cout << "precompute done. " << std::endl;
      
      
      std::cout << "compute: numerical phase ..." << std::endl;
      ExternPardiso_compute(pardiso,
			    "N",
			    corr_,
			    rhs_);
      std::cout << "compute done. " << std::endl;
      ExternPardiso_kill(pardiso);
    }
    
    for (I i=0;i<=jacobian.m_n;++i) jacobian.m_begin[i]-=1;
    for (I i=0;i<jacobian.m_nc;++i) jacobian.m_index[i]-=1;
#endif
#endif
  }









  
  static DG_JACOBIAN* create_jacobian(DG_HANDLE* 	handle_,				      
				      DG::INFO * 	dgi,
				      DG_DATA*   	data_,
				      cst_pR 		xa_,
				      cst_pR 		xu_,
				      cst_pI  		cnc_u_,
				      cst_pI  		cncoff_u_,
				      cst_pR 		data_u_,		       
				      cst_pR		data_v_,
				      cst_pI 		nelm_,
				      cst_pR 		coo_,
				      cst_pI 		cooff_,
				      cst_pI 		cnc_,
				      cst_pI 		cncoff_,
				      cst_pI 		adj_,
				      cst_pI 		adjoff_,
				      cst_pI 		vcod_,
				      cst_pI 		noboundary_cod_)
  {
    cst_pI 	rwork_n_ 	= &dgi->dg_rwork_n;      
    pR  	rwork_ 		= &dgi->dg_rwork[0];
    //    cst_pI 	iwork_n_	= &dgi->dg_iwork_n;
    //    pI		iwork_		= &dgi->dg_iwork[0];			      
    cst_pR 	rinfo_ 		= &dgi->dg_rinfo[0];
    pI  	iinfo_ 		= &dgi->dg_iinfo[0];
    //    pR 		rres_ 		= &dgi->dg_rres[0];
    pI 		ires_		= &dgi->dg_ires[0];


    
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    
    DG_JACOBIAN* jacobian = new DG_JACOBIAN(nelm_[0],
					    nfaceinelm,
					    trial_n[0],
					    adj_,
					    adjoff_);
    
    //    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    //    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];

    //    cst_pI	teta_n      = &iinfo_[DG::I_TETA_NBASIS];
    //    cst_pI	test_n      = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_u_n      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    //const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    //    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    //    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    //    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      //      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/


    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  

    I
      id,
      ielm,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1];
      
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*teta_u_n[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*teta_u_n[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/

    //
    //
    //
    // R lcrhs[21*2];

    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return NULL;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
      
    /*___________________________________________________________________________________________________________________*/  
    for (ielm = 0; ielm < nelm_[0];++ielm)
      {
	id = ielm;

	//	std::cout << "ielm " << ielm << std::endl;
	neids[0]	= adj_[id*adjoff_[0]];
	neids[1]	= adj_[id*adjoff_[0]+1];
	neids[2]	= adj_[id*adjoff_[0]+2];
	cnc[0]        = cnc_[id*cncoff_[0]+0];
	cnc[1]        = cnc_[id*cncoff_[0]+1];
	cnc[2]        = cnc_[id*cncoff_[0]+2];

	
	neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
	/*---*/
	if (neids[0])
	  {	if (neids[0]<0) neids[0]=-neids[0];
	    if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
	if (neids[1])
	  {	if (neids[1]<0) neids[1]=-neids[1];
	    if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
	if (neids[2])
	  {	if (neids[2]<0) neids[2]=-neids[2];
	    if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
	/*---*/
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
	jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

	data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
	data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
	data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
	data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
	x = data_->mat_belm.x[2];
	data_->mat_belm.x[2]=data_->mat_belm.x[1];
	data_->mat_belm.x[1]=x;
    
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  

	//
	//
	//
	if (data_a_)
	  {
	    for (I k=0;k<teta_a_n[0];++k)
	      {
		brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	      }
	  }
	else
	  {
	    brhs_a[0] = jacelm[0] * xa_[0];
	  }
  
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u : elm and faces
	*/  
	/*___________________________________________________________________________________________________________________*/  

	//
	// Copy the velocity.
	//
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	
	//
	// Evaluate u at all the required coordinates.
	//
	handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
	
	/*____ TRANSFORM ON UELM */
	handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

	//
	// On the boundary.
	//
	for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
	  }
    
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex].apply([xu_](const R r)
						   {
						     return (r<0.0) ? -r*xu_[0] : ((R)0.0);
						   });
	  }
	
	/*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    if (neids[localFaceIndex])
	      {
		//		printf("  flux edge " ifmt " \n",localFaceIndex);
		double mr1= -1.0;
		Blas_dgemv("N",&bmat_n,qface_n,&mr1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);
		// vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
		// vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);		
		jacobian->addelm(id,neids[localFaceIndex]-1,0.0,data_->m_flux_part_corr.x,data_->m_flux_part_corr.ld);
		
		//
		//
		//
		//		data_->m_flux_part_corr;
		//
		//		data_->m_flux_part_sol;

		// data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
		//	vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);
		//		data_->m_local_rhs -= data_->m_flux_part_sol  * data_->vec_neisol;	 
	      }
	  }
	
	/*--- COMPUTE MATRICES  */
	{
	  // matrix_handle msub;
	  // vector_handle rsub;
	  
	  // matrix_handle_def(&msub,handle_->m_BMAT.n,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BMAT.x,handle_->m_BMAT.ld);
	  // vector_handle_def(&rsub,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BRHS.x,1);
	  
	  data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
	  // data_->m_local_matrices = msub * rsub;
	  jacobian->addelm(id,id,0.0,data_->m_local_matrix_part_corr.x,data_->m_local_matrix_part_corr.ld);

	  //	  vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
	  //	  data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	  //	  printf("addelmcorr\n");
	}
	

	//	{I j;for (j=0;j<test_n[0];++j) printf(" sisis" ifmt " %e\n",j, rhs_[rhsoff_[0]*id+j]);}
	//
	// aE + u.nabla(E) + h(u,E) = -aF - u.nabla(F) - h(u,F)
	// JACOBIAN(trial,test) E = -R(F,test)
	// 
	// corr P3 F P2  Operator P2-P3 Operator P3-P3
	//
	// P0
	//   P1\P0
	//        P2\P0
	//
	//
	//
	//
	//
	//
	//
	//
	//

	//	data_->m_local_matrix_part_corr;
	/*--- COMPUTE LOCAL RESIDUAL  */
	// data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	/*--- COMPUTE LOCAL CORRECTION   */
	//	matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
	// matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
      }
    return jacobian;
  }

  static void compute_residual(DG_HANDLE* 	handle_,				      
			       DG::INFO * 	dgi,
			       DG_DATA*   	data_,
			       cst_pR 		xa_,
			       cst_pR 		xu_,
			       cst_pI  		cnc_u_,
			       cst_pI  		cncoff_u_,
			       cst_pR 		data_u_,		       
			       cst_pR		data_v_,
			       cst_pI 		nelm_,
			       cst_pR 		coo_,
			       cst_pI 		cooff_,
			       cst_pI 		cnc_,
			       cst_pI 		cncoff_,
			       cst_pI 		adj_,
			       cst_pI 		adjoff_,
			       cst_pI 		vcod_,
			       cst_pI 		noboundary_cod_)
#if 0
			       mkS 	shape_X,
			       mkS 	shape_E,
			       DG_HANDLE* handle_,
			       DG_DATA*   data_,
			       cst_pR 	xa_,
			       cst_pR 	xu_,
			       pR		rhs_,
			       cst_pI 	rhsoff_,
			       
			       cst_pI  	cnc_u_,
			       cst_pI  	cncoff_u_,
			       cst_pR 	data_u_,		       
			       cst_pR	data_v_,
				     
			       cst_pR	sol_,
			       cst_pI 	soloff_,
			      
			       pR		corr_,
			       cst_pI 	corroff_,
			      
			       cst_pR 	t_,
			       cst_pI 	nelm_,
			       cst_pR 	coo_,
			       cst_pI 	cooff_,
			       cst_pI 	cnc_,
			       cst_pI 	cncoff_,
			       cst_pI 	adj_,
			       cst_pI 	adjoff_,
			       cst_pI 	vcod_,
			       cst_pI 	noboundary_cod_,
			      
			       cst_pI 	rwork_n_,
			       pR  	rwork_,
			       cst_pI 	iwork_n_,
			       pI	iwork_,			      
			       cst_pR 	rinfo_,
			       const I  	iinfo_[DG::I_n],
			       R 		rres_[DG_rres_n],
			       I 		ires_[DG_ires_n]
#endif
			       
  {
    
    //    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    
    //    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    //    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];

    cst_pI	teta_n      = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	test_n      = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_u_n      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    //const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    //    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    //    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    //    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      //      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1;
  
    /*_____________________________________________________*/
    const R mxu=xu_[0]*((R)1.0);

    /*_____________________________________________________*/


    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  

    I
      id,
      ielm,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1];
      
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*teta_u_n[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*teta_u_n[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/

    //
    //
    //
    // R lcrhs[21*2];

    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
      
    /*___________________________________________________________________________________________________________________*/  
    for (ielm = 0; ielm < nelm_[0];++ielm)
      {
	id = ielm;
	{I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	
	
	
	//	std::cout << "ielm " << ielm << std::endl;
	neids[0]	= adj_[id*adjoff_[0]];
	neids[1]	= adj_[id*adjoff_[0]+1];
	neids[2]	= adj_[id*adjoff_[0]+2];
	cnc[0]        = cnc_[id*cncoff_[0]+0];
	cnc[1]        = cnc_[id*cncoff_[0]+1];
	cnc[2]        = cnc_[id*cncoff_[0]+2];

	
	neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
	/*---*/
	if (neids[0])
	  {	if (neids[0]<0) neids[0]=-neids[0];
	    if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
	    else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
	if (neids[1])
	  {	if (neids[1]<0) neids[1]=-neids[1];
	    if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
	    else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
	if (neids[2])
	  {	if (neids[2]<0) neids[2]=-neids[2];
	    if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
	    else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
	/*---*/
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
	cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
	cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
	cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
	cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
    
	jacface[0] = sqrt( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = sqrt( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = sqrt( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
	data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	jacelm[0]  = sqrt((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

	data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
	data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
	data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
	data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
	x = data_->mat_belm.x[2];
	data_->mat_belm.x[2]=data_->mat_belm.x[1];
	data_->mat_belm.x[1]=x;
    
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  

	//
	//
	//
	if (data_a_)
	  {
	    for (I k=0;k<teta_a_n[0];++k)
	      {
		brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	      }
	  }
	else
	  {
	    brhs_a[0] = jacelm[0] * xa_[0];
	  }
  
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u : elm and faces
	*/  
	/*___________________________________________________________________________________________________________________*/  

	//
	// Copy the velocity.
	//
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	{ I k;
	  for (k=0;k<teta_u_n[0];++k)
	    {
	      data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	    } }
	
	//
	// Evaluate u at all the required coordinates.
	//
	handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
	
	/*____ TRANSFORM ON UELM */
	handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

	//
	// On the boundary.
	//
	for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
	  }
    
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    data_->vec_uface[localFaceIndex].apply([xu_](const R r)
						   {
						     return (r<0.0) ? -r*xu_[0] : ((R)0.0);
						   });
	  }
	
	/*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
	for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	  {
	    if (neids[localFaceIndex])
	      {
		//		printf("  flux edge " ifmt " \n",localFaceIndex);
		double mr1= -1.0;
		Blas_dgemv("N",&bmat_n,qface_n,&mr1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);
		// vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
		// vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);		

		//		jacobian->addelm(id,neids[localFaceIndex]-1,0.0,data_->m_flux_part_corr.x,data_->m_flux_part_corr.ld);

		//
		//
		//
		//		data_->m_flux_part_corr;
		//
		//		data_->m_flux_part_sol;

		// data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
		vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);
		data_->m_local_rhs -= data_->m_flux_part_sol  * data_->vec_neisol;	 
	      }
	  }
	
	/*--- COMPUTE MATRICES  */
	{
	  // matrix_handle msub;
	  // vector_handle rsub;
	  
	  // matrix_handle_def(&msub,handle_->m_BMAT.n,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BMAT.x,handle_->m_BMAT.ld);
	  // vector_handle_def(&rsub,teta_a_n[0] + dim * teta_u_n[0],handle_->m_BRHS.x,1);
	  
	  data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;

	  // data_->m_local_matrices = msub * rsub;
	  // jacobian->addelm(id,id,0.0,data_->m_local_matrix_part_corr.x,data_->m_local_matrix_part_corr.ld);

	  vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
	  data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	  //	  printf("addelmcorr\n");
	}
	
	{I j;for (j=0;j<test_n[0];++j) rhs_[rhsoff_[0]*id+j] = data_->m_local_rhs.x[j];}	
	//	{I j;for (j=0;j<test_n[0];++j) printf(" sisis" ifmt " %e\n",j, rhs_[rhsoff_[0]*id+j]);}
	//
	// aE + u.nabla(E) + h(u,E) = -aF - u.nabla(F) - h(u,F)
	// JACOBIAN(trial,test) E = -R(F,test)
	// 
	// corr P3 F P2  Operator P2-P3 Operator P3-P3
	//
	// P0
	//   P1\P0
	//        P2\P0
	//
	//
	//
	//
	//
	//
	//
	//
	//

	//	data_->m_local_matrix_part_corr;
	/*--- COMPUTE LOCAL RESIDUAL  */
	// data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
	/*--- COMPUTE LOCAL CORRECTION   */
	//	matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
	// matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
      }

  }

  
};



static void coordinates2d(long long int				N_, 
			  double				box_[],
			  double* 				vertex_,
			  long long int                         vertexoff_,
			  unsigned long long int *  		cod2_,
			  long long int				codoff_)
{
  static const unsigned long long m = 1LL<<62;
  static const int BitTab[2] 	= {1,2};
  static const int GeoCod[4]	= {1,2,0,3};
  static const int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
  static const double len 	= 4.611686018427387904e+18;
  double box[4],dbl;
  unsigned long long int IntCrd[2],cod;
  int rot[4];
  box[0] 	= box_[0];
  box[1] 	= box_[1];
  box[2] 	= len / (box_[2] - box_[0]);
  box[3] 	= len / (box_[3] - box_[1]);
  double loc[2];
  { int i;
    for(i=0; i<N_; i++)
      {
	loc[0] = vertex_[vertexoff_*i+0];
	loc[1] = vertex_[vertexoff_*i+1];
	/* Convert double precision coordinates to integers */
	dbl 	= (loc[0] - box[0]) * box[0+2];
	IntCrd[0] = (unsigned long long int)dbl;
	dbl 	= (loc[1] - box[1]) * box[1+2];
	IntCrd[1] = (unsigned long long int)dbl;
	/* Binary hilbert renumbering loop */
	cod = 0;
	rot[0] = GeoCod[0];
	rot[1] = GeoCod[1];
	rot[2] = GeoCod[2];
	rot[3] = GeoCod[3];
	{ int b;
	  for(b=0;b<31;b++)
	    {
	      int GeoWrd = 0;

	      if(IntCrd[0] & m)
		GeoWrd |= BitTab[0];
	      IntCrd[0] = IntCrd[0]<<1;

	      if(IntCrd[1] & m)
		GeoWrd |= BitTab[1];
	      IntCrd[1] = IntCrd[1]<<1;

	      const int NewWrd = rot[ GeoWrd ];
	      cod = cod<<2 | NewWrd;
	      rot[0] = HilCod[ NewWrd ][ rot[0] ];
	      rot[1] = HilCod[ NewWrd ][ rot[1] ];
	      rot[2] = HilCod[ NewWrd ][ rot[2] ];
	      rot[3] = HilCod[ NewWrd ][ rot[3] ];
	    } }
	cod2_[codoff_*i+0] 	= cod;
	cod2_[codoff_*i+1] 	= i;
      } }

}

struct hilbert_reordering_t
{
  unsigned long long int cod;
  unsigned long long int id;
};


int hilbert_comp_t(const void * a_,const void * b_)
{
  const I * a = (const I*)a_;
  const I * b = (const I*)b_;
  if (*a < *b)
    {
      return -1;
    }
  else if (*a > *b)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

void hilbert_reordering(ns_mesh * 	mesh_,
			pI 		perm_)
{
  double * vertices = (double *)malloc(sizeof(double) * 2*mesh_->nelm);
  I * tmp = (I*)malloc(sizeof(I) * 2*mesh_->nelm);
  for (I ielm=0;ielm<mesh_->nelm;++ielm)
    {
      //      double x = 0.0;
      //      double y = 0.0;

      I i0 = mesh_->cnc[6*ielm+0]-1;
      I i1 = mesh_->cnc[6*ielm+1]-1;
      I i2 = mesh_->cnc[6*ielm+2]-1;

      cst_pR p0 = &mesh_->coo[2*i0];
      cst_pR p1 = &mesh_->coo[2*i1];
      cst_pR p2 = &mesh_->coo[2*i2];
      
      vertices[2*ielm+0]  = (p0[0]+p1[0]+p2[0])/3.0;
      vertices[2*ielm+1]  = (p0[1]+p1[1]+p2[1])/3.0;
    }
  
  //
  // Compute the centroids of the elements.
  //
  
  double bb[4];
  bb[0] = vertices[0];
  bb[1] = vertices[1];
  bb[2] = vertices[0];
  bb[3] = vertices[1];
  for (I ielm=1;ielm<mesh_->nelm;++ielm)
    {
  
  bb[0] = std::min(bb[0],vertices[2*ielm+0]);
  bb[1] = std::min(bb[1],vertices[2*ielm+1]);
  bb[2] = std::max(bb[2],vertices[2*ielm+0]);
  bb[3] = std::max(bb[3],vertices[2*ielm+1]);
 
}


  bb[0] = 0;
  bb[1] = 0;
  bb[2] = 1;
  bb[3] = 2;
#if 1
  coordinates2d(mesh_->nelm, 
    bb,
    vertices,
    2,
    (unsigned long long int*)tmp,
    2);

  qsort(tmp,mesh_->nelm,2*sizeof(I),hilbert_comp_t);
#endif
  for (I ielm=0;ielm<mesh_->nelm;++ielm)
    {
  // perm_[tmp[2*ielm+1]] = ielm;
  perm_[ielm] = tmp[2*ielm+1];
}


  FILE * f = fopen("curve.txt","w");
  for (I i=0;i<mesh_->nelm;++i)
    {
  fprintf(f," %e %e\n",vertices[2*perm_[i]+0],vertices[2*perm_[i]+1]);
}
  fclose(f);
  
  for (I ielm=0;ielm<mesh_->nelm;++ielm)
    {
  // perm_[tmp[2*ielm+1]] = ielm;
  perm_[tmp[2*ielm+1]] = ielm;
}
  
}


struct topology_elem
{
  I m_partid[1];
  I m_cnc[6];
  I m_adj[3];
  I m_gid[1];
  I m_h[1];
};


int comp_part(const void * a ,const void * b)
{
  cst_pI a_ = (cst_pI)a;
  cst_pI b_ = (cst_pI)b;
  if (a_[0] < b_[0] )
    {
      return -1;
    }
  else if (a_[0] > b_[0] )
    {
      return 1;
    }
  else
    {
      //
      // inside the partitions.
      //
      
      
      

      
      return 0;
    }
}

struct partition
{
  I  m_nelm;
  I  m_nelm_interior_ghost;
  I  m_nelm_exterior_ghost;

  pI m_cnc;
  I  m_cncoff;
  pI m_adj;
  I  m_adjoff;
  
  pI m_coo;
  pI m_gidelm;
  pI m_gidnode;
};

 void partition(ns_mesh * mesh,I nparts,pI partid,I partid_off)
 {
  //
  // Space filling curve reordering of the elements.
  //
  pI perm = (pI)malloc(sizeof(I)*mesh->nelm);
  hilbert_reordering(mesh,perm);

  I nelm_per_partitions = mesh->nelm / nparts;
  pI part_begin = (pI)calloc(nparts+1,sizeof(I));
  part_begin[0] = 0;
  for (I ipart=1;ipart<nparts;++ipart)
    {
      part_begin[ipart] = part_begin[ipart-1]  + nelm_per_partitions;
    }
  part_begin[nparts] = part_begin[nparts-1]  + nelm_per_partitions + (mesh->nelm % nparts);

  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      I jelm = perm[ielm];
      I k = 1;
      for (;k<=nparts;++k)
	{
	  if (jelm < part_begin[k])
	    {
	      partid[ielm] = k;
	      break;
	    }
	}
    }


  //
  //
  //
  topology_elem * elm = (topology_elem *)malloc(sizeof(topology_elem)*mesh->nelm);
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      for (I i=0;i<6;++i)
	elm[ielm].m_cnc[i] = mesh->cnc[6*ielm+i];
      for (I i=0;i<3;++i)
	elm[ielm].m_adj[i] = mesh->adj[3*ielm+i];
      elm[ielm].m_gid[0] = ielm;
      elm[ielm].m_partid[0] = partid[ielm];
      elm[ielm].m_h[0] = perm[ielm];
    }
  
  qsort(elm,mesh->nelm,sizeof(topology_elem),
	[](const void * a,const void * b)
	{
	  cst_pI a_ = (cst_pI)a;
	  cst_pI b_ = (cst_pI)b;
	  if (a_[0] < b_[0] )
	    {
	      return -1;
	    }
	  else if (a_[0] > b_[0] )
	    {
	      return 1;
	    }
	  else
	    {
	      
	      if (a_[11] < b_[11] )
		{
		  return -1;
		}
	      else if (a_[11] > b_[11] )
		{
		  return 1;
		}
	      
#if 0
	      //
	      // inside the partitions.
	      //
	      bool a_is = false;
	      for (I i=0;i<3;++i)
		{
		  if (a_[7+i] && partid[ a_[7+i] ] != a_[0] )
		    {
		      a_is = true;
		      break;
		    }
		}

	      bool b_is = false;
	      for (I i=0;i<3;++i)
		{
		  if (b_[7+i] && partid[ b_[7+i] ] != b_[0] )
		    {
		      b_is = true;
		      break;
		    }
		}

	      if (a_is && !b_is)
		{
		  return 1;
		}
	      else if (!a_is && b_is)
		{
		  return -1;
		}
	      else
		{
		  return 0;
		}
#endif
	      return 0;
	    }
	  
	});

  perm = (pI)malloc(sizeof(I)*mesh->nelm);
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      perm[elm[ielm].m_gid[0]] = ielm;
    }
  
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      for (I i=0;i<6;++i)
	mesh->own_cnc[6*ielm+i] = elm[ielm].m_cnc[i];
      for (I i=0;i<3;++i)	
	mesh->own_adj[3*ielm+i] = (elm[ielm].m_adj[i]==0) ? 0 : perm[elm[ielm].m_adj[i]-1] + 1;
      partid[ielm] = elm[ielm].m_partid[0];
    }

  free(perm);
  perm = NULL;

  
  pI marker = (pI)calloc(mesh->nelm,sizeof(I));
  pI select = (pI)calloc(mesh->nelm,sizeof(I));
  I select_n=0;
  pI part_adj    = (pI)calloc(nparts*nparts,sizeof(I));
  pI part_nghost = (pI)calloc(nparts,sizeof(I));


  //  pI part_cncs[32];
  for (I ipart = 0;ipart<nparts;++ipart)
    {
      
      I start = part_begin[ipart];      
      I bound = part_begin[ipart+1];
      
      for (I k=0;k<select_n;++k)
	{
	  marker[select[k]]=0;	  
	}
      select_n=0;
      
      //
      // Flag the ghost.
      //
      for (;start<bound;++start)
	{
	  I ielm = start;
	  for (I k=0;k<3;++k)
	    {
	      I nei = elm[ielm].m_adj[k];
	      if (nei > 0)
		{
		  I neipart = partid[nei-1];
		  if (neipart != ipart)
		    {
		      if (0==marker[nei-1])
			{
			  part_adj[nparts * ipart + neipart] += 1;
			  
			  select[select_n] = nei-1;
			  marker[nei-1] = ++select_n;
			  
			  part_nghost[ipart] +=1;				  
			}		  
		    }
		}
	    }
	}

      
      I nelm = part_nghost[ipart] + part_begin[ipart+1] - part_begin[ipart];
      topology_elem * pelm = (topology_elem *)malloc(sizeof(topology_elem)*nelm);
      memcpy(pelm,&elm[part_begin[ipart]],sizeof(topology_elem)*(part_begin[ipart+1] - part_begin[ipart]));
      start = part_begin[ipart];

      
      for (I k=0;k<select_n;++k)
	{
	  marker[select[k]]=0;	  
	}
      select_n=0;

      
      I iii = part_begin[ipart+1] - part_begin[ipart];
      for (;start<bound;++start)
	{
	  I ielm = start;
	  for (I k=0;k<3;++k)
	    {
	      I nei = elm[ielm].m_adj[k];
	      if (nei > 0)
		{
		  I neipart = partid[nei-1];
		  if (neipart != ipart)
		    {
		      part_adj[nparts * ipart + neipart] = 1;
		      part_adj[nparts * neipart + ipart] = 1;

		      if (0==marker[nei-1])
			{
			  memcpy(&pelm[iii],&elm[nei-1],sizeof(topology_elem));
			  elm[iii].m_partid[0] = neipart;
			  //			  +iii;
			  
			  select[select_n] = nei-1;
			  marker[nei-1] = ++select_n;			  
			}		  
		    }
		}
	    }
	}

      
    }

  

  
  
  {
    void dg_print_mesh(I numNodes,I nelm,pI cnc,I cncoff,pR coo,pI codnodes,I codnodesoff,pI codelm,I codelmoff,const char * name_,...);
  }
  
#if 0  
  
  partition * p = (partition *)malloc(sizeof(partition)*nparts); 
  for (I k = 0;k<nparts;++k)
    {
      I  part_nelm = part_begin[k+1]-part_begin[k];
      pI part_cnc  = malloc(6*sizeof(I)*part_nelm);
      I  part_nelm_exterior_ghost = 0;
      I  part_nelm_interior_ghost = 0;
    }
#endif  

  
#if 0  
struct partition
{
  I  m_nelm;
  I  m_nelm_interior_ghost;
  I  m_nelm_exterior_ghost;
  pI m_cnc;
  pI m_coo;
  pI m_adj;
  pI m_gidelm;
  pI m_gidnode;
};
#endif


  free(part_begin);
 }


struct DG_BOUNDARY_CONDITION
{
  double m_mem[16000];
  I m_qn;
  I m_nu;
  I m_nx;
  I m_ntest;
  //
  // Matrix to build the evaluation of u.
  //
  matrix_handle m_beval_uvw;
  matrix_handle m_beval_uvw_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  matrix_handle m_beval_xyz;
  matrix_handle m_beval_xyz_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  matrix_handle m_beval_test;
  matrix_handle m_beval_test_part[nfaceinelm];

  
  matrix_handle m_qrst;
  matrix_handle m_qrst_part[nfaceinelm];
  
#if 0  
  matrix_handle eval_uvw;
  matrix_handle nrmelm;
  matrix_handle cooelm;
  matrix_handle uvwelm;
  vector_handle udotn;
  matrix_handle qxyz;

  vector_handle f;
#endif
  DG_BOUNDARY_CONDITION(){};
  DG_BOUNDARY_CONDITION(mkS s_test_,
			mkS s_teta_u_,
			mkS s_teta_x_,
			I 	qn,
			cst_pR 	qw,
			cst_pR 	qp)
  {
    const I qnXnfaceinelm = qn * nfaceinelm;
    const I nu = mkS_n(s_teta_u_);
    const I nx = mkS_n(s_teta_x_);
    const I ntest = mkS_n(s_test_);
    m_qn = qn;
    m_nu = nu;
    m_nx = nx;
    m_ntest = ntest;
    
    I at = 0;

    //
    // Reference quadrature points 
    //
    {
      I qrst_size = qnXnfaceinelm*dim;
      matrix_handle_def(&this->m_qrst,
			qnXnfaceinelm,
			dim,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    this->m_qrst.x + localFaceIndex * qn,
			    this->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   this->m_qrst.x,
		   &this->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      matrix_handle_def(&this->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&this->m_mem[at],
			nx);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    this->m_beval_xyz.x + (localFaceIndex * qn) * this->m_beval_xyz.ld,
			    this->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      matrix_handle_def(&this->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&this->m_mem[at],
			ntest);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    this->m_beval_test.x + (localFaceIndex * qn) * this->m_beval_test.ld,
			    this->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      matrix_handle_def(&this->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    this->m_beval_uvw.x + localFaceIndex * qn * this->m_beval_uvw.ld,
			    this->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    I err;
    R rwork[4096];
    I rwork_n = 4096;

    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      this->m_beval_xyz.x,
	      &this->m_beval_xyz.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      this->m_beval_uvw.x,
	      &this->m_beval_uvw.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      this->m_beval_test.x,
	      &this->m_beval_test.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    for (I l=0;l<nfaceinelm;++l)
      {
	for (I j=0;j<qn;++j)
	  {
	    for (I i=0;i<ntest;++i)
	      {
		m_beval_test_part[l].x[m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }



  static void define(DG_BOUNDARY_CONDITION*self,
		     mkS s_test_,
		     mkS s_teta_u_,
		     mkS s_teta_x_,
		     I 	qn,
		     cst_pR 	qw,
		     cst_pR 	qp)
  {
    const I qnXnfaceinelm = qn * nfaceinelm;
    const I nu = mkS_n(s_teta_u_);
    const I nx = mkS_n(s_teta_x_);
    const I ntest = mkS_n(s_test_);
    self->m_qn = qn;
    self->m_nu = nu;
    self->m_nx = nx;
    self->m_ntest = ntest;
    
    I at = 0;

    //
    // Reference quadrature points 
    //
    {
      I qrst_size = qnXnfaceinelm*dim;
      matrix_handle_def(&self->m_qrst,
			qnXnfaceinelm,
			dim,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    self->m_qrst.x + localFaceIndex * qn,
			    self->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   self->m_qrst.x,
		   &self->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      matrix_handle_def(&self->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&self->m_mem[at],
			nx);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    self->m_beval_xyz.x + (localFaceIndex * qn) * self->m_beval_xyz.ld,
			    self->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      matrix_handle_def(&self->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&self->m_mem[at],
			ntest);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    self->m_beval_test.x + (localFaceIndex * qn) * self->m_beval_test.ld,
			    self->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      matrix_handle_def(&self->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    self->m_beval_uvw.x + localFaceIndex * qn * self->m_beval_uvw.ld,
			    self->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    I err;
    R rwork[4096];
    I rwork_n = 4096;

    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      self->m_beval_xyz.x,
	      &self->m_beval_xyz.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      self->m_beval_uvw.x,
	      &self->m_beval_uvw.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      self->m_beval_test.x,
	      &self->m_beval_test.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    for (I l=0;l<nfaceinelm;++l)
      {
	for (I j=0;j<qn;++j)
	  {
	    for (I i=0;i<ntest;++i)
	      {
		self->m_beval_test_part[l].x[self->m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }

  

  struct DATA
  {
    double m_mem[512];
    matrix_handle eval_xyz;
    matrix_handle eval_uvw;    
    vector_handle eval_udotn;
    vector_handle eval_f;
    vector_handle lrhs;

    DATA(I qn,I ntest)
    {
      // qn * (2 * dim + 3)
      I at = 0;
      matrix_handle_def(&eval_xyz,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      matrix_handle_def(&eval_uvw,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      vector_handle_def(&eval_udotn,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&eval_f,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&lrhs,ntest,&m_mem[at],1);
      at+=ntest;
    };
  };

  DATA * CreateData()
  {
    return new DATA(this->m_qn,this->m_ntest);
  };
  R m_eps;
  
  SmoothedHeaviside m_hea;
  template <typename userfct_t>
  void boundary_condition(const I 		localFaceIndex,
			  cst_pR xu_,
			  const vector_handle&	normal,
			  const matrix_handle&	xyz,
			  const matrix_handle&	uvw,
			  DATA*			data,
			  userfct_t             userfct)
  {
    m_eps = 0.25;
    {
      Err e;
      SmoothedHeaviside_def	( &this->m_hea,
				  __eHeaviside_m4p3,
				  &m_eps,
				  &e);
    }
    
    
    //
    // Evaluate xyz.
    //
    data->eval_xyz = this->m_beval_xyz_part[localFaceIndex].transpose() * xyz;
#if 0
    
    fprintf(stdout,"xyz " ifmt " " ifmt " " ifmt " " ifmt "\n",data->eval_xyz.n,data->eval_xyz.m,data->eval_xyz.ld,data->eval_f.ld);
    matrix_handle_print(&data->eval_xyz,stdout);
#endif

    //
    // User function.
    //


    //
    // eval_f pouvait ne pas etre initialise dans userfct
    //
    userfct(data->eval_xyz.n,
	    data->eval_xyz.x,
	    data->eval_xyz.ld,
	    data->eval_f.x,
	    data->eval_f.ld);
    
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif
    //
    // Evaluate uvw.
    //
#if 0
    fprintf(stdout,"uvw\n");
    matrix_handle_print(&uvw,stdout);
#endif
    
    data->eval_uvw = this->m_beval_uvw_part[localFaceIndex].transpose() * uvw;
    
#if 0
    fprintf(stdout,"eval_uvw\n");
    matrix_handle_print(&data->eval_uvw,stdout);
#endif    
    //
    // Compute the dot product with the normal
    //
    data->eval_udotn = data->eval_uvw * normal;
#if 0
    fprintf(stdout,"eval_udotn\n");
    vector_handle_print(&data->eval_udotn,stdout);
#endif

    //
    //
    //
    
    //
    // Form the flux functions.
    //
    
    for (I i=0;i<data->eval_f.n;++i)
      {	data->eval_f.x[i*data->eval_f.ld]=1.0;
	data->eval_f.x[i*data->eval_f.ld]
	  *= ( data->eval_udotn.x[i*data->eval_udotn.ld] < ((R)0.0) ) ? -data->eval_udotn.x[i*data->eval_udotn.ld] * xu_[0] : ((R)0.0) ;
      }
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif	      
    //
    // Form the residu.
    //
    data->lrhs = this->m_beval_test_part[localFaceIndex] * data->eval_f;
#if 0
    fprintf(stdout,"lrhs\n");
    vector_handle_print(&data->lrhs,stdout);
#endif
  }

  
  void boundary_condition(ns_mesh * 	mesh,
			  cst_pI  	cnc_u_,
			  cst_pI  	cncoff_u_,
			  cst_pR 	data_u_,		       
			  cst_pR	data_v_,
			  pR 		rhs_,
			  I 		rhsoff_)
  {    

    vector_handle normal;
    matrix_handle xyz;
    matrix_handle uvw;
    double normal_values[6];
    double cooelm[32];
    double uvw_values[32];
    matrix_handle_def(&xyz,m_nx,dim,cooelm,m_nx);
    matrix_handle_def(&uvw,m_nu,dim,uvw_values,m_nu);
    vector_handle_def(&normal,2,normal_values,1);
    
    DATA* data = this->CreateData();
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {

		//
		//
		//
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const I jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		normal_values[0]	= mesh->normaledge[2*jedge+0];
		normal_values[1] 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];
		//

		//
		//
		//
		
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];;
		  }
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[m_nu+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		  }
		//
		
		boundary_condition(jadj,
				   &xjac,
				   normal,
				   xyz,
				   uvw,
				   data,
				   [this](I n_,cst_pR pos,I posoff,pR f,I foff)
				   {
				     for (I k=0;k<n_;++k)
				       {
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];

					 
					 if ( (x < 1.0e-13) && (y>0.5 && y<=1.0) )
					   {
					 f[foff*k] = 0.0;
				       }
					 else if ( (x < 1.0e-13) && (y<=0.5 && y >0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x>0.5 && x<=1.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x<=0.5 && x >=0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (x > 1.0-1.0e-13) )
					   {
					 double z  = sin(acos(-1.0)*y);
					 f[foff*k] = z*z*z;
				       }
					 else
					   {
					 f[foff*k] = 0.0;
				       }
				       }
#if 0
				     
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = y - 0.5; // exp(-sin(32.0*y)*y);
				       }
				     Err e;
				     SmoothedHeaviside_eval	(&this->m_hea,
								 n_,
								 f,
								 &e);
#endif
#if 0
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = exp(-sin(32.0*y)*y);
				       }
#endif
				   });
		
		for (I k=0;k<m_ntest;++k)
		  {
		    rhs_[jelm * rhsoff_ + k] += data->lrhs.x[k];
		  }		
	      }
	  }
      }
    
#if 0
    printf("ddddddddddddddddddd\n");
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (I k=0;k<m_ntest;++k)
	  {
	    std::cout << rhs_[jelm * rhsoff_ + k] << std::endl;
	  }
	printf("\n");
      }
    printf("ddddddddddddddddddd\n");
#endif    

    
#if 0
    double mem[2048];
    matrix_handle beval_uvw[nfaceinelm];
    matrix_handle eval_uvw;
    matrix_handle beval_xyz[nfaceinelm];
    matrix_handle nrmelm;
    matrix_handle cooelm;
    matrix_handle uvwelm;
    vector_handle udotn;
    matrix_handle qxyz;
    vector_handle f;
  
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const I jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		const R nx	= mesh->normaledge[2*jedge+0];
		const R ny 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];	    
		
	      //
	      // Interpolation of u * nx + v*ny + w * nz 
	      //

	      
	      //
	      // Evaluation of u on the edge.
	      //
	      eval_uvw = beval_uvw[jadj] * dofs_uvw;

	      //
	      // Evaluation of u on the edge.
	      //
	      udotn = eval_uvw * nrmedge;

	      //
	      //
	      // Cartesian coordinates of the quadrature points.
	      //
	      //
	      qxyz = build_qxyz * cooelm;

	      //
	      // Evaluation of the user function.
	      //
	      

	      //
	      // Form the flux functions.
	      //
	      for (I i=0;i<nq;++i)
		{
		  f.x[i] *= ( udotn.x[i] < ((R)0.0) ) ? -udotn.x[i] * xu_[0] : ((R)0.0) ;
		}
	      
	      //
	      // Form the residu.
	      //
	      lrhs = qbasis * f;

	    }
	}
    }
#endif
};
};



struct DG_OPERATOR
{
  DG_JACOBIAN * 		m_jacobian;
  DG_BOUNDARY_CONDITION * 	m_bc;
  DG_DATA * 			m_dgdata;
  DG_HANDLE * 			m_dghandle;
  DG::INFO * 			m_dginfo;
  
  DG_OPERATOR(mkS shape_X,mkS shape_A,mkS shape_U,mkS shape_F,mkS shape_E,mkS shape_Test,
	      cst_pR 		xa_,
				      cst_pR 		xu_,
				      cst_pI  		cnc_u_,
				      cst_pI  		cncoff_u_,
				      cst_pR 		data_u_,		       
				      cst_pR		data_v_,
				      cst_pI 		nelm_,
				      cst_pR 		coo_,
				      cst_pI 		cooff_,
				      cst_pI 		cnc_,
				      cst_pI 		cncoff_,
				      cst_pI 		adj_,
				      cst_pI 		adjoff_,
				      cst_pI 		vcod_,
				      cst_pI 		noboundary_cod_)
  {
    
    double iw[6];
    double ip[6];
    iw[0] = 1.713244923791704e-1;
    iw[1] = 3.607615730481386e-1;
    iw[2] = 4.679139345726911e-1;
    iw[3] = 4.679139345726911e-1;
    iw[4] = 3.607615730481386e-1;
    iw[5] = 1.713244923791704e-1;
  
    ip[0] = -9.32469514203152e-1;
    ip[1] = -6.612093864662645e-1;
    ip[2] = -2.386191860831969e-1;
    ip[3] = 2.386191860831969e-1;
    ip[4] = 6.612093864662645e-1;
    ip[5] = 9.32469514203152e-1;
    //    I rhsoff 	= mkS_n(shape_E);
    //    I soloff 	= mkS_n(shape_F);
    //    I corroff 	= mkS_n(shape_E);
    //    I adjoff=3;
    //    I cncoff=6;

    this->m_bc = (DG_BOUNDARY_CONDITION *)calloc(1,sizeof(DG_BOUNDARY_CONDITION));
    this->m_dghandle = (DG_HANDLE *)calloc(1,sizeof(DG_HANDLE));
    this->m_dginfo = (DG::INFO *)calloc(1,sizeof(DG::INFO));
    printf("info define\n");
    DG::INFO::define(this->m_dginfo);

    printf("bc define\n");
    DG_BOUNDARY_CONDITION::define(this->m_bc,
				  shape_E,
				  shape_U,
				  shape_X,
				  6,iw,ip);
    printf("solver define\n");
    DG::define(this->m_dghandle,
	       shape_A,
	       shape_U,
	       shape_F,
	       shape_E,
	       shape_E,
	       &this->m_dginfo->dg_iinfo_n,
	       this->m_dginfo->dg_iinfo,
	       &this->m_dginfo->dg_rinfo_n,
	       this->m_dginfo->dg_rinfo,
	       &this->m_dginfo->dg_rwork_n,
	       this->m_dginfo->dg_rwork);
    
    this->m_dgdata 	= DG::create_data(shape_A,
					  shape_U,
					  shape_F,
					  shape_E,
					  shape_Test);



    this->m_jacobian 	= DG::create_jacobian(this->m_dghandle,
					      this->m_dginfo,
					      this->m_dgdata,
					      xa_,
					      xu_,					   
					      cnc_u_,
					      cncoff_u_,
					      data_u_,
					      data_v_,
					      nelm_,
					      coo_,
					      cooff_,
					      cnc_,
					      cncoff_,
					      adj_,
					      adjoff_,
					      vcod_,
					      noboundary_cod_);
    
  };

  void compute_residual()
  {
    
  }
  
  DG::compute_residual(mkS 	shape_X,
			       mkS 	shape_E,
			       DG_HANDLE* handle_,
			       DG_DATA*   data_,
			       cst_pR 	xa_,
			       cst_pR 	xu_,
			       pR		rhs_,
			       cst_pI 	rhsoff_,
			       
			       cst_pI  	cnc_u_,
			       cst_pI  	cncoff_u_,
			       cst_pR 	data_u_,		       
			       cst_pR	data_v_,
				     
			       cst_pR	sol_,
			       cst_pI 	soloff_,
			      
			       pR		corr_,
			       cst_pI 	corroff_,
			      
			       cst_pR 	t_,
			       cst_pI 	nelm_,
			       cst_pR 	coo_,
			       cst_pI 	cooff_,
			       cst_pI 	cnc_,
			       cst_pI 	cncoff_,
			       cst_pI 	adj_,
			       cst_pI 	adjoff_,
			       cst_pI 	vcod_,
			       cst_pI 	noboundary_cod_,
			      
			       cst_pI 	rwork_n_,
			       pR  	rwork_,
			       cst_pI 	iwork_n_,
			       pI	iwork_,			      
			       cst_pR 	rinfo_,
			       const I  	iinfo_[DG::I_n],
			       R 		rres_[DG_rres_n],
			       I 		ires_[DG_ires_n])

  
};





int main(int 		argc,
	 const char**	argv)
{



  

  //
  // Set up the monitor
  //
#if 0
  Monitor_def(0,
	      argv[0],
	      MonitorMode_STD);
#endif  
  mesh = (ns_mesh*)calloc(1,sizeof(ns_mesh));
  STR errmsg;
  Err err;
  ns_mesh_read(mesh,
	       errmsg,
	       &err,
	       argv[1]);

#if 0
  I npartitions = 5;

  pI codelm = (pI)malloc(sizeof(I)*mesh->nelm);
  partition(mesh,npartitions,codelm,1);

  dg_print_mesh(mesh,
		codelm,
		"dg");



  



  
  
  printf("hilbert done.");
  exit(1);
#endif
  

  
#if 0
  ns_mesh_write_medit(mesh,
		      out);
#endif
  
  

  I degree = 3;
  
  mkS_st shape_A;
  mkS_st shape_F;
  mkS_st shape_U;
  mkS_st shape_E;
  mkS_st shape_X;
  mkS_st shape_EE[32];
  mkS_st shape_FF[32];
    
mkS_definit	(&shape_A,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 0,
		 __emk_discontinuous,
		 &err);
  mkS_definit	(&shape_F,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_l2orthonormal,
		 degree,
		 __emk_discontinuous,
		 &err);
  
  mkS_definit	(&shape_E,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_l2orthonormal,
		 degree,
		 __emk_discontinuous,
		 &err);

  for (I i=0;i<=degree;++i)
    {
      mkS_definit	(&shape_EE[i],
			 __eTopology_TRIANGLE,
			 __emkS_FAMILY_l2orthonormal,
			 i,
			 __emk_discontinuous,
			 &err);
      mkS_definit	(&shape_FF[i],
			 __eTopology_TRIANGLE,
			 __emkS_FAMILY_l2orthonormal,
			 i,
			 __emk_discontinuous,
			 &err);
    }

  mkS_definit	(&shape_U,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 2,
		 __emk_discontinuous,
		 &err);		    

  mkS_definit	(&shape_X,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 1,
		 __emk_discontinuous,
		 &err);

  I nddlu 		= mesh->nddlP6;
  I nddlv 		= mesh->nddlP6;
  
  double * u 		= (double*)calloc(nddlu,sizeof(R)); 
  double * v 		= (double*)calloc(nddlv,sizeof(R));

  //
  // On interpo
  //
  for (int i=0;i<nddlu;++i)
    {
      double x= mesh->coo[2*i+0];
      double y= mesh->coo[2*i+1];
      u[i] = 10.0*y*y-12.0*x+1.0;
    }
  for (int i=0;i<nddlv;++i)
    {
      double y= mesh->coo[2*i+1];
      v[i] = 1.0 + y;
    }

    

  DG_OPERATOR*dg[32];

  const I noboundary_cod  	= 100;
  
  const I  cooff = 2;
  const I cncoff = 6;
  const I adjoff = 3;
  
  cst_pI cnc_u = mesh->cnc;
  double xa=0.0,xu=1.0;
  for (I ideg=0;ideg<=degree;++ideg)
    {
      printf("############################\n");
      dg[ideg] = new DG_OPERATOR(&shape_X,
				 &shape_A,
				 &shape_U,
				 &shape_FF[ideg],
				 &shape_EE[ideg],
				 &shape_EE[ideg],
				 &xa,
				 &xu,		  
				 cnc_u,
				 &cncoff,
				 u,
				 v,
				 &mesh->nelm,
				 mesh->coo,
				 &cooff,
				 mesh->cnc,
				 &cncoff,
				 mesh->adj,
				 &adjoff,
				 mesh->cod,
				 &noboundary_cod);
      if (dg[ideg])
	{

	}
      printf("############################\n");
    }

  I rhs_level_n[32];
  pR rhs_level[32];
  

  for (I i=0;i<=degree;++i)
    {
      rhs_level_n[i] = mesh->nelm * mkS_n(&shape_EE[i]);
      rhs_level[i] = (double*)calloc(rhs_level_n[i],sizeof(R));
    }

  

  
  dg[degree]->m_bc->boundary_condition(mesh,
				       mesh->cnc,
				       &cncoff,
				       u,
				       v,			   
				       rhs_level[degree],
				       mkS_n(&shape_EE[degree]));
  
#if 1
  I sol_level_n[32];
  I corr_level_n[32];
  pR sol_level[32];
  for (I i=0;i<=degree;++i)
    {
      sol_level_n[i] = mesh->nelm * mkS_n(&shape_FF[i]);
      sol_level[i] = (pR)calloc(sol_level_n[i],sizeof(R));
    }
  pR corr_level[32];
  for (I i=0;i<=degree;++i)
    {
      corr_level_n[i] = mesh->nelm * mkS_n(&shape_EE[i]);
      corr_level[i] = (pR)calloc(corr_level_n[i],sizeof(R));
    }

#endif

  solve_gauss_seidel(*dg[degree]->m_jacobian,
		     700,
		     degree,
		     corr_level[degree],
		     rhs_level[degree]);

  DG_VIEW view(&shape_X,&shape_EE[degree]);

  dg_print_mesh( mesh,
		 view,
		 "dgjacobi." ifmt,degree);  

  dg_print_sol(view,
	       &shape_X,
	       &shape_E,
	       mesh->nelm,
	       corr_level[degree],
	       "dgjacobi." ifmt,degree);

  
  
#define NONO 0
#if NONO
  
#if 0  
  //
  // Evaluate RHS
  //
  dg_bc[degree].boundary_condition(mesh,
					 cnc_u,
					 &cncoff_u,
					 u,
					 v,			   
					 rhs_level[degree_solve],
					 mkS_n(&shape_EE[degree_solve]));

  
  

  DG_DATA *data  = DG::create_data(&shape_A,
				       &shape_U,
				       &shape_FF[degree_solve],
				       &shape_EE[degree_solve],
				       &shape_EE[degree_solve]);
      DG_VIEW view(&shape_X,&shape_EE[degree_solve]);
      
      I rhsoff=mkS_n(&shape_EE[degree_solve]);
      I soloff=mkS_n(&shape_FF[degree_solve]);
      I corroff=mkS_n(&shape_EE[degree_solve]);
#endif

      





  
  double iw[6];
  double ip[6];
  iw[0] = 1.713244923791704e-1;
  iw[1] = 3.607615730481386e-1;
  iw[2] = 4.679139345726911e-1;
  iw[3] = 4.679139345726911e-1;
  iw[4] = 3.607615730481386e-1;
  iw[5] = 1.713244923791704e-1;
  
  ip[0] = -9.32469514203152e-1;
  ip[1] = -6.612093864662645e-1;
  ip[2] = -2.386191860831969e-1;
  ip[3] = 2.386191860831969e-1;
  ip[4] = 6.612093864662645e-1;
  ip[5] = 9.32469514203152e-1;


  DG_HANDLE h[32];
  
  DG_BOUNDARY_CONDITION dg_bc[32];
  DG::INFO dg_info[32];

  for (I ideg=0;ideg<=degree;++ideg)
    {
      printf("degree " ifmt " \n",ideg);
      DG::define(&h[ideg],
		 &shape_A,
		 &shape_U,
		 &shape_FF[ideg],
		 &shape_EE[ideg],
		 &shape_EE[ideg],
		 &dg_info[ideg].dg_iinfo_n,
		 dg_info[ideg].dg_iinfo,
		 &dg_info[ideg].dg_rinfo_n,
		 dg_info[ideg].dg_rinfo,
		 &dg_info[ideg].dg_rwork_n,
		 dg_info[ideg].dg_rwork);
      
      DG_BOUNDARY_CONDITION::define(&dg_bc[ideg],
				    &shape_EE[ideg],
				    &shape_U,
				    &shape_X,
				    6,iw,ip);
    }

  
  const double xu 		= 1.0;
  const double t 		= 1.0;
  const I noboundary_cod  	= 100;

  const I  cooff = 2;
  const I cncoff = 6;
  const I adjoff = 3;

  cst_pI cnc_u = mesh->cnc;
  const I cncoff_u = 6;
  //    exit(1);


  
  
  I degree_solve = 0;
  const I nelm 		= ns_mesh_nelm(mesh);
  const I numNodes 	= ns_mesh_get_numNodes(mesh);
#if 0
  const I fn 		= mkS_n(&shape_FF[degree_solve]);
  const I en 		= mkS_n(&shape_EE[degree_solve]);
  const I flen 		= fn * nelm;
  const I elen 		= en * nelm;
  double * sol		= (double*)calloc(flen,sizeof(R)); 
  double * rhs 		= (double*)calloc(elen,sizeof(R)); 
  double * corr 	= (double*)calloc(elen,sizeof(R));
  const I rhsoff 	= fn;
  const I soloff 	= fn;
  const I corroff 	= fn;
#endif
  I sol_level_n[32];
  I rhs_level_n[32];
  I corr_level_n[32];
  double * sol_level[32];
  for (I i=0;i<=degree;++i)
    {
      sol_level_n[i] = nelm * mkS_n(&shape_FF[i]);
      sol_level[i] = (double*)calloc(sol_level_n[i],sizeof(R));
    }
  
  double * rhs_level[32];
  for (I i=0;i<=degree;++i)
    {
      rhs_level_n[i] = nelm * mkS_n(&shape_EE[i]);
      rhs_level[i] = (double*)calloc(rhs_level_n[i],sizeof(R));
    }

  double * corr_level[32];
  for (I i=0;i<=degree;++i)
    {
      corr_level_n[i] = nelm * mkS_n(&shape_EE[i]);
      corr_level[i] = (double*)calloc(corr_level_n[i],sizeof(R));
    }

  DG_JACOBIAN * jacobians[32];
    DG_DATA *data[32]; 
#if 0
#endif
  for (I i=0;i<=degree;++i)
    {
      
      data[i] = DG::create_data(&shape_A,
				&shape_U,
				&shape_FF[i],
				&shape_EE[i],
				&shape_EE[i]);

      double xa = 0.0;
      I rhsoff = mkS_n(&shape_EE[i]);
      I soloff = mkS_n(&shape_FF[i]);
      I corroff = mkS_n(&shape_EE[i]);
#if 0
      jacobians[i] = DG::create_jacobian(&shape_X,
					 &shape_EE[i],
					 &h[i],
					 data[i],
					 &xa,
					 &xu,
					 rhs_level[i],
					 &rhsoff,
					       
					 cnc_u,
					 &cncoff_u,
					 u,
					 v,
		       
					 sol_level[i],
					 &soloff,
					       
					 corr_level[i],
					 &corroff,
					       
					 &t,
					       
					 &nelm,
					 mesh->coo,
					 &cooff,
					 mesh->cnc,
					 &cncoff,
					 mesh->adj,
					 &adjoff,
					 mesh->cod,
					 &noboundary_cod,
					       
					 &dg_info[degree_solve].dg_rwork_n,
					 dg_info[degree_solve].dg_rwork,
					 &dg_info[degree_solve].dg_iwork_n,
					 dg_info[degree_solve].dg_iwork,
					 dg_info[degree_solve].dg_rinfo,
					 dg_info[degree_solve].dg_iinfo,
					 &dg_info[degree_solve].dg_rres[0],
					 &dg_info[degree_solve].dg_ires[0]);

#endif      
      //
      //      new DG_JACOBIAN(nelm_[0],
      //				     nfaceinelm,
      //				     &shape_EE[i],
      //				     adj_,
      //				     adjoff_);
    }

  //
  // Form rhs of the system of the highest degree.
  //
  dg_bc[degree].boundary_condition(mesh,
				   cnc_u,
				   &cncoff_u,
				   u,
				   v,			   
				   rhs_level[degree],
				   mkS_n(&shape_EE[degree]));
  
  //
  // Residual_k(sol_k)
  //
  //
  // Projection on P_{k-1}
  //
  //
  // Projection on P_{k-2}
  //
  //
  // Projection on P_{k-3}
  //

  

  
#if 0

#endif


  
  double xa = 0.0;
  degree_solve=degree;
  for (;degree_solve <=degree;++degree_solve)
    {
  
  dg_bc[degree_solve].boundary_condition(mesh,
    cnc_u,
    &cncoff_u,
    u,
    v,			   
    rhs_level[degree_solve],
    mkS_n(&shape_EE[degree_solve]));
  

  DG_DATA *data  = DG::create_data(&shape_A,
				       &shape_U,
				       &shape_FF[degree_solve],
				       &shape_EE[degree_solve],
				       &shape_EE[degree_solve]);
      DG_VIEW view(&shape_X,&shape_EE[degree_solve]);
      
      I rhsoff=mkS_n(&shape_EE[degree_solve]);
      I soloff=mkS_n(&shape_FF[degree_solve]);
      I corroff=mkS_n(&shape_EE[degree_solve]);


      
//  static void solve_global(I jacdeg,
//			   DG_JACOBIAN* jacobians,
//			   mkS shape_X,
//			   mkS shape_E,
//			   DG_VIEW& view,
//			   DG_HANDLE*  handle_,
//			   DG_DATA*    data_,
//			   cst_pR 	xa_,
//			   cst_pR 	xu_,
//			   pR	rhs_,
//			   cst_pI 	rhsoff_,
//
      DG::solve_global(degree_solve,
		       jacobians,

		       
		       &shape_X,
		       &shape_EE[degree_solve],
		       view,
		       &h[degree_solve],
		       &data[degree_solve],
		       &xa,
		       &xu,
		       rhs_level[degree_solve],
		       &rhsoff,
		       
		       cnc_u,
		       &cncoff_u,
		       u,
		       v,
		       
		       sol_level[degree_solve],
		       &soloff,
		   
		   corr_level[degree_solve],
		   &corroff,
		   
		   &t,
		      
		   &nelm,
		   mesh->coo,
		   &cooff,
		   mesh->cnc,
		   &cncoff,
		   mesh->adj,
		   &adjoff,
		   mesh->cod,
		   &noboundary_cod,

		   &dg_info[degree_solve].dg_rwork_n,
		   dg_info[degree_solve].dg_rwork,
		   &dg_info[degree_solve].dg_iwork_n,
		   dg_info[degree_solve].dg_iwork,
		   dg_info[degree_solve].dg_rinfo,
		   dg_info[degree_solve].dg_iinfo,
		   &dg_info[degree_solve].dg_rres[0],
		   &dg_info[degree_solve].dg_ires[0]);
  
      I n1=1;
      R r1=1.0;
      daxpy(&corr_level_n[degree_solve],&r1,corr_level[degree_solve],&n1,sol_level[degree_solve],&n1);
      
      
      //
      // Copy the solution into the next solution
      //
      if (degree_solve < degree)
	{
	  I n = mkS_n(&shape_FF[degree_solve]);
	  I N = mkS_n(&shape_FF[degree_solve+1]);
	  for (I jelm=0;jelm<mesh->nelm;++jelm)
	    {
	      for (I j=0;j<n;++j)
		{
		  sol_level[degree_solve+1][jelm* N + j] = sol_level[degree_solve][jelm*n + j]    ;
		}
	      
	    }
	}

  

      

  
#if 1
  std::cout << "post-processing ..." << std::endl;
  std::cout << "   export dgview.mesh ..." << std::endl;

  dg_print_mesh( mesh,
		 view,
		 "dgview");  
  std::cout << "   export dgview.bb ..." << std::endl;
  dg_print_sol(	view,
		&shape_X,
		&shape_EE[degree_solve],
		nelm,
		sol_level[degree_solve],
		"dgview");
  std::cout << "post-processing done." << std::endl;
#endif



    }


  
  
#if 0
  dg_print_mesh(mesh,
		view,
		"dgview.0");
  for (I i=0;i<nelm;++i)
    {
      for (I j=0;j<6;++j)
	{
	  corr[i*10+j]=0.0;
	}
#if 0
      for (I j=6;j<10;++j)
	{
	  corr[i*10+j]=0.0;
	}
#endif
    }
  dg_print_sol(	view,
		&shape_X,&shape_E,
		nelm,corr,
		"dgview.0");
#endif


#if 0
  //
  // Make a selection.
  //
  I  nselect = 10;
  pI select  = malloc(sizeof(I)*nselect);
  for (I i=0;i<nselect;++i)
    {
      select[i]=i+1;
    }

  
  //
  // Form a submesh with split.
  //
  I split_nelm 	= 4*nselect;
  pI * cnc 	= malloc(sizeof(I)*6*split_nelm);
  pI * adj 	= malloc(sizeof(I)*3*split_nelm);
  I lcnc[6];
  for (I i=0;i<nselect;++i)
    {
      I jelm = select[i];
      for (I j=0;j<6;++j)
	{
	  cnc[(4*i)*6+j] = mesh->cnc[6*jelm+j];
	}      
      for (I j=0;j<3;++j)
	{
	  adj[(4*i)*3+j] = mesh->adj[3*jelm+j];
	}      
    }

  
  for (I i=0;i<nselect;++i)
    {
      I jelm = select[i];
      for (I j=0;j<6;++j)
	{
	  lcnc[j] = cnc[(4*i)*6+j];
	}
      
      cnc[(4*i)*6+0] = lcnc[3];
      cnc[(4*i)*6+1] = lcnc[4];
      cnc[(4*i)*6+2] = lcnc[5];
      
      cnc[(4*i + 1)*6+0] = lcnc[0];
      cnc[(4*i + 1)*6+1] = lcnc[3];
      cnc[(4*i + 1)*6+2] = lcnc[5];

      cnc[(4*i + 2)*6+0] = lcnc[1];
      cnc[(4*i + 2)*6+1] = lcnc[4];
      cnc[(4*i + 2)*6+2] = lcnc[3];

      cnc[(4*i + 3)*6+0] = lcnc[2];
      cnc[(4*i + 3)*6+1] = lcnc[5];
      cnc[(4*i + 3)*6+2] = lcnc[4];
      
      for (I j=0;j<3;++j)
	{
	  ladj[j] = adj[(4*i)*3+j];
	}

      adj[(4*i)*3+0] 	= mesh->nelm+4*i+1+1;
      adj[(4*i)*3+1] 	= mesh->nelm+4*i+2+1;
      adj[(4*i)*3+2] 	= mesh->nelm+4*i+3+1;

      adj[(4*i+1)*3+0] 	= -ladj[0];
      adj[(4*i+1)*3+1] 	= mesh->nelm + 4*i+0+1;
      adj[(4*i+1)*3+2] 	= -ladj[2];

      adj[(4*i+2)*3+0] 	= -ladj[1];
      adj[(4*i+2)*3+1] 	= mesh->nelm + 4*i+0+1;
      adj[(4*i+2)*3+2] 	= -ladj[0];
      
      adj[(4*i+3)*3+0] 	= -ladj[2];
      adj[(4*i+3)*3+1] 	= mesh->nelm + 4*i+0+1;
      adj[(4*i+3)*3+2] 	= -ladj[1];

    }

  
  
  //
  // Interpolate velocity.
  //

  //
  // Interpolate correction.
  //
  
  //
  // Interpolate solution.
  //

  
  //
  // Calculate boundary conditions from level-1 in the residual.
  //
  
  //
  // Solve
  //

  //
  // Force correction to be zero where it is refined.
  //

  //
  // Calculate boundary conditions from level+1 in the residual.
  //

  //
  // Solve
  //
#endif  
  std::cout << "#####" << std::endl;
#endif  
  return 0;
}



