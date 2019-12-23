#include "Config.h"
#include "Type.h"
#include "Medit.h"
#include "extern_medit.h"
#include "eVolume.h"
#include "eFace.h"

typedef struct
{
  Medit 		self_empty;
  int 			version;
  int 			inm;
  eTopologyDimension 	topologyDimension;
  eDim 			dim;
  I 			nbVolumes[__eVolume_ALL];
  I 			nbFaces[__eFace_ALL];
  I 			nbEdges;
  I 			nbVertices;
} MeditImpl,* RESTRICT pMeditImpl;

typedef const MeditImpl	* RESTRICT cst_pMeditImpl;


static const enum GmfKwdCod __eFace_TO_GmfKwCod[ __eFace_ALL ] 		= {GmfReserved1,
									   GmfTriangles,
									   GmfQuadrilaterals,
									   GmfReserved1};

static const enum GmfKwdCod __eVolume_TO_GmfKwCod[ __eVolume_ALL ] 	= {GmfReserved1,
									   GmfTetrahedra,
									   GmfHexahedra,
									   GmfWedges,
									   GmfPyramids,
									   GmfReserved1};


eDim 			Medit_get_eDim			(cst_pMedit 	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  cst_pMeditImpl self = (cst_pMeditImpl)self_;
  return self->dim;
}

eTopologyDimension 	Medit_get_topologyDimension	(cst_pMedit   	const	self_)
{  
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  pMeditImpl  self = (pMeditImpl)self_;
  return self->topologyDimension;
}


pMedit  		Medit_kill			(pMedit		const	self_)
{
  if (self_)
    {
      pMeditImpl self = (pMeditImpl)self_;
      if (self->inm)
	{
	  self->inm = GmfCloseMesh(self->inm);     
	}
      self->inm 	= 0;
      self->dim 	= __eDim_ERROR;
      self->version 	= 0;
      free(self);
    }
  return NULL;
}


void  			Medit_set_precision		(pMedit		const 	self_,
							 cst_ePrecision		precision_)
{
#if 0
  pMeditImpl self = malloc(sizeof(MeditImpl));
#endif
}

I 			Medit_get_nbFaces		(cst_pMedit 	const 	self_,
							 cst_eFace 		face_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  cst_pMeditImpl self = (cst_pMeditImpl)self_;
  return self->nbFaces[face_];
}

I 			Medit_get_nbEdges		(cst_pMedit 	const	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  cst_pMeditImpl self = (cst_pMeditImpl)self_;
  return self->nbEdges;
}

I 			Medit_get_nbVertices		(cst_pMedit 	const	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  cst_pMeditImpl self = (cst_pMeditImpl)self_;
  return self->nbVertices;
}

I 			Medit_get_nbVolumes		(cst_pMedit 	const 	self_,
							 cst_eVolume 		volume_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  cst_pMeditImpl self = (cst_pMeditImpl)self_;
  return self->nbVolumes[volume_];
}


pMedit  		Medit_new			(cst_pS const filename_,
							 ...)
{
#ifndef NDEBUG
  DebugVerif(filename_);
#endif

  STR filename;

  { va_list args;
    va_start(args,filename_);
    vsprintf(filename,filename_,args);
    va_end(args); }

  pMeditImpl 	self 	= (pMeditImpl)calloc(1,sizeof(MeditImpl));
  self->inm 		= GmfOpenMesh(filename,GmfRead,&self->version,&self->dim);

  if (!self->inm)
    return Medit_kill(&self->self_empty);    
  
  { eVolume volume = __eVolume_ERROR;
    for (++volume;volume<__eVolume_ALL;++volume)
      {
	self->nbVolumes[volume] = GmfStatKwd(self->inm,__eVolume_TO_GmfKwCod[volume]);
      } }
  
  { eFace face = __eFace_ERROR;
    for (++face;face<__eFace_ALL;++face)
      {
	self->nbFaces[face] = GmfStatKwd(self->inm,__eFace_TO_GmfKwCod[face]);
      } }

  self->nbEdges 	= GmfStatKwd(self->inm,GmfEdges);
  self->nbVertices 	= GmfStatKwd(self->inm,GmfVertices);

  self->topologyDimension = __eTopologyDimension_Point;  
  { eVolume volume = __eVolume_ERROR;
    for (++volume;volume<__eVolume_ALL;++volume)
      {
	if (self->nbVolumes[volume]>0)
	  {
	    self->topologyDimension = __eTopologyDimension_Volume;
	    break;	    
	  }
      } }

  if (self->topologyDimension==__eTopologyDimension_Point)
    {
      { eFace face = __eFace_ERROR;
	for (++face;face<__eFace_ALL;++face)
	  {
	    if (self->nbFaces[face]>0)
	      {
		self->topologyDimension = __eTopologyDimension_Face;
		break;	    
	      }
	  } }
    }

  if (self->topologyDimension==__eTopologyDimension_Point)
    {
      if (self->nbEdges > 0)
	{
	  self->topologyDimension = __eTopologyDimension_Edge;
	}
    }
  return &self->self_empty;
}




#define NC     32
#define NC2    (NC*NC)
#define NC3    (NC2*NC)
#define KA     31
#define KB     57
#define KC     79
static int idir[5]  = {0,1,2,0,1};

/* very sioux! (09/2002) auteur : Frey, package medit */
void hash(int nt,int**adja,unsigned char**voy,int * medit_cnc) 
{
  int       k,kk,l,ll,MINs,MAXs,MINs1,MAXs1,hsize;
  int      *hcode,*link,inival,iadr,pp;
  char     *hvoy;
  unsigned char     i,i1,i2,ii;
  unsigned int key;
  /* avoid building again! */
  if ( 4*sizeof(char) != sizeof(int) )  {
    printf("aaa\n");
    exit(1);
  }
  /* memory alloc */
  hcode = (int*)calloc(MAX(1,3*nt/4)+1,sizeof(int));
  link  = (int*)calloc(3*nt+1,sizeof(int));
  hsize = MAX(2,3*nt/4-1);
  hvoy = (char*)hcode;
  /* init */
  inival = 2147483647;
  for (k=0; k<=3*nt/4; k++)
    hcode[k] = -inival;
  /* build hash table */
  for (k=1; k<=nt; k++) {

    for (i=0; i<3; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];

      MINs = MIN(medit_cnc[(k-1)*3+i1],medit_cnc[(k-1)*3+i2]);
      MAXs = MAX(medit_cnc[(k-1)*3+i1],medit_cnc[(k-1)*3+i2]);
     
      /* compute key */
      key = KA*MINs + KB*MAXs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  /* set adjacency */
  for (l=3*nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = idir[i+1];
    i2 = idir[i+2];
    MINs = MIN(medit_cnc[(k-1)*3+i1],medit_cnc[(k-1)*3+i2]);
    MAXs = MAX(medit_cnc[(k-1)*3+i1],medit_cnc[(k-1)*3+i2]);


    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = idir[ii+1];
      i2 = idir[ii+2];

      MINs1 = MIN(medit_cnc[(kk-1)*3+i1],medit_cnc[(kk-1)*3+i2]);
      MAXs1 = MAX(medit_cnc[(kk-1)*3+i1],medit_cnc[(kk-1)*3+i2]);

      /* adjacent found */
      if ( MINs1 == MINs && MAXs1 == MAXs ) {
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = kk;
        hvoy[l] = ii;
        link[ll]= k;
        hvoy[ll]= i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  adja[0] = (int*)link;
  voy[0]  = (unsigned char*)hcode;
#if 1
  for (l=0;l<nt; ++l) 
    {
      int tmp;

#if 1
      unsigned char ctmp;
      ctmp = voy[0][1+3*l+0];
      voy[0][1+3*l+0]=voy[0][1+3*l+2];
      voy[0][1+3*l+2]=ctmp;

      ctmp = voy[0][1+3*l+1];
      voy[0][1+3*l+1]=voy[0][1+3*l+2];
      voy[0][1+3*l+2]=ctmp;
#endif



      tmp = link[1+3*l+0];
      link[1+3*l+0]=link[1+3*l+2];
      link[1+3*l+2]=tmp;

      tmp = link[1+3*l+1];
      link[1+3*l+1]=link[1+3*l+2];
      link[1+3*l+2]=tmp;

    }  
#endif
#undef MIN
#undef MAX
}
#undef NC
#undef NC2
#undef NC3
#undef KA 
#undef KB 
#undef KC 


void Medit_gotoSectionCell(pMedit const self_,const eVolume volume_)
{
  pMeditImpl 	self 	= (pMeditImpl)self_;
  GmfGotoKwd(self->inm,__eVolume_TO_GmfKwCod [volume_]); 
}

void Medit_get_cellToNodes	(pMedit const self_,const eVolume volume_,pI const cellToNodeIndices_)
{
  pMeditImpl 	self 	= (pMeditImpl)self_;
  const enum GmfKwdCod gmfcod 	=  __eVolume_TO_GmfKwCod [volume_];
  int inm = self->inm;
  int localcod;
  int ci[8];
  switch(volume_)
    {
    case __eVolume_TETRAHEDRON:
      {
	GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&localcod);
	cellToNodeIndices_[0] = ci[0]-1;
	cellToNodeIndices_[1] = ci[1]-1;
	cellToNodeIndices_[2] = ci[2]-1;
	cellToNodeIndices_[3] = ci[3]-1;
	break;
      }
    case __eVolume_HEXAHEDRON:
      {
	GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&ci[6],&ci[7],&localcod);
	cellToNodeIndices_[0] = ci[0]-1;
	cellToNodeIndices_[1] = ci[1]-1;
	cellToNodeIndices_[2] = ci[2]-1;
	cellToNodeIndices_[3] = ci[3]-1;
	cellToNodeIndices_[4] = ci[4]-1;
	cellToNodeIndices_[5] = ci[5]-1;
	cellToNodeIndices_[6] = ci[6]-1;
	cellToNodeIndices_[7] = ci[7]-1;
	break;
      }
    case __eVolume_PYRAMID:
      {
	GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&localcod);
	cellToNodeIndices_[0] = ci[0]-1;
	cellToNodeIndices_[1] = ci[1]-1;
	cellToNodeIndices_[2] = ci[2]-1;
	cellToNodeIndices_[3] = ci[3]-1;
	cellToNodeIndices_[4] = ci[4]-1;
	break;
      }
    case __eVolume_WEDGE:
      {
	GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&localcod);
	cellToNodeIndices_[0] = ci[0]-1;
	cellToNodeIndices_[1] = ci[1]-1;
	cellToNodeIndices_[2] = ci[2]-1;
	cellToNodeIndices_[3] = ci[3]-1;
	cellToNodeIndices_[4] = ci[4]-1;
	cellToNodeIndices_[5] = ci[5]-1;
	break;
      }
    case __eVolume_POLYHEDRON:
    case __eVolume_ERROR:
    case __eVolume_ALL:
      {
	break;
      }
    }

}

void Medit_gotoSectionVertices(pMedit 	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  pMeditImpl 	self 	= (pMeditImpl)self_;
  GmfGotoKwd(self->inm,GmfVertices); 
}

void Medit_get_GridVertex(pMedit 	const 	self_,
			  pR 		const 	coo_)
{ 
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(coo_);
#endif
  int icod;
  pMeditImpl 	self 	= (pMeditImpl)self_;
  GmfGetLin(self->inm,GmfVertices,&coo_[0],&coo_[1],&coo_[2],&icod);
  /*  printf("%e %e %e\n",coo_[0],coo_[1],coo_[2]);*/
}


void Medit_get_cncv(pMedit 	const self_,
		    cst_pI	const cootr_,
		    pR 		const coo_,
		    cst_pI	const cooff_,
		    pI		const cod_,
		    cst_pI	const codoff_)
{ 
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(cootr_);
  DebugVerif(coo_);
  DebugVerif(cooff_);
  DebugVerif(cod_);
  DebugVerif(codoff_);
#endif

  pMeditImpl 	self 	= (pMeditImpl)self_;
  cst_eDim 	dim 	= self->dim;
  int 		inm 	= self->inm;
  const I 	nv 	= self->nbVertices;
  const I 	cootr	= cootr_[0];
  const I 	cooff	= cooff_[0];
  const I 	codoff	= codoff_[0];
  pI		cod	= cod_;
  double 	c0[3];
  int 		icod;
  GmfGotoKwd(self->inm,GmfVertices); 
  if (!cootr)
    {
      pR coo = coo_;
      if (dim==__eDim_3)
	{
	  { I i;
	    for (i=0;i<nv;++i)
	      {
		GmfGetLin(inm,GmfVertices,&c0[0],&c0[1],&c0[2],&icod);
		coo[0]=c0[0];
		coo[1]=c0[1];
		coo[2]=c0[2];
		cod[0]=icod;
		coo += cooff;
		cod += codoff;
	      } }
	}
      else if (dim==__eDim_2)
	{
	  { I i;
	    for (i=0;i<nv;++i)
	      {
		GmfGetLin(inm,GmfVertices,&c0[0],&c0[1],&icod);
		coo[0]=c0[0];
		coo[1]=c0[1];
		cod[0]=icod;
		cod+=codoff;
		coo += cooff;
	      } }
	}
    }
  else
    {     
      NotYetImplemented;
    }
}


void Medit_get_Points(pMedit 	self_,
		      pPoints 	xyz_,
		      pI	cod_,
		      cst_pI	codoff_)
{ 
  pMeditImpl 	self 	= (pMeditImpl)self_;
  cst_eDim 	dim 	= self->dim;
  int 		inm 	= self->inm;
  const I 	nbVertices 	= self->nbVertices;
  const I 	codoff	= codoff_[0];
  pI		cod	= cod_;
  double 	c0[3];
  int 		icod;
  GmfGotoKwd(self->inm,GmfVertices); 
  switch(dim)
    {
    case __eDim_3:
      {
	{ I i;
	  for (i=0;i<nbVertices;++i)
	    {
	      GmfGetLin(inm,GmfVertices,&c0[0],&c0[1],&c0[2],&icod);
	      Points_set(xyz_,i,c0);
	      cod[codoff*i+0]=icod;
	    } }
	break;
      }
    case __eDim_2:
      {
	{ I i;
	  for (i=0;i<nbVertices;++i)
	    {
	      GmfGetLin(inm,GmfVertices,&c0[0],&c0[1],&icod);
	      Points_set(xyz_,i,c0);
	      cod[codoff*i+0]=icod;
	    } }
	break;
      }
    case __eDim_1:
      {
	NotYetImplemented;
	break;	
      }
    case __eDim_ERROR:
    case __eDim_ALL:
      {
	break;
      }
    }
}


void Medit_get_cncVolume(pMedit 	const 	self_,
			 pI 		const 	cnc_,
			 const I 		cncoff_,
			 pI 		const 	cod_,
			 const I 		codoff_,
			 cst_eVolume 		volume_)
{ 
  pI cnc 			= cnc_;
  pI cod 			= cod_;
  pMeditImpl 	self 		= (pMeditImpl)self_;
  int inm			= self->inm;
  const I nbVolumes 		= self->nbVolumes[volume_];
  const enum GmfKwdCod gmfcod 	=  __eVolume_TO_GmfKwCod [volume_];
  GmfGotoKwd(inm,gmfcod); 
  int localcod;
  int ci[8];
  switch(volume_)
    {
    case __eVolume_TETRAHEDRON:
      {
	{ I i;
	  for (i=0;i<nbVolumes;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc[3] = ci[3]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod   = localcod;
		  cod   += codoff_;
		}
	    } }
	break;
      }
    case __eVolume_HEXAHEDRON:
      {
	{ I i;
	  for (i=0;i<nbVolumes;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&ci[6],&ci[7],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc[3] = ci[3]-1;
	      cnc[4] = ci[4]-1;
	      cnc[5] = ci[5]-1;
	      cnc[6] = ci[6]-1;
	      cnc[7] = ci[7]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod = localcod;
		  cod += codoff_;
		}

	    } }
	break;
      }
    case __eVolume_PYRAMID:
      {
	{ I i;
	  for (i=0;i<nbVolumes;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc[3] = ci[3]-1;
	      cnc[4] = ci[4]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod = localcod;
		  cod += codoff_;
		}

	    } }
	break;
      }
    case __eVolume_WEDGE:
      {
	{ I i;
	  for (i=0;i<nbVolumes;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc[3] = ci[3]-1;
	      cnc[4] = ci[4]-1;
	      cnc[5] = ci[5]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod = localcod;
		  cod += codoff_;
		}

	    } }
	break;
      }
    case __eVolume_POLYHEDRON:
      {
	NotYetImplemented;
	break;
      }
    case __eVolume_ERROR:
    case __eVolume_ALL:
      {
	break;
      }
    }
}



eVolume Medit_get_volumeKind(cst_pMedit 	const 		self_,
			     const I 				volumeIndex_)
{
    cst_pMeditImpl self = (cst_pMeditImpl)self_;

  I n = 0;
  { eVolume volume = __eVolume_ERROR;
    for (++volume;volume<__eVolume_ALL;++volume)
      {
	n += self->nbVolumes[volume];
	if (volumeIndex_ < n)
	  {
	    return volume;
	  }
      } }
  return __eVolume_ERROR;  
}

#if 0
int MeditVolumeIterator_next	(pMeditVolumeIterator self_);
void MeditVolumeIterator_def	(pMeditVolumeIterator self_,pMedit const medit_);
void MeditVolumeIterator_free	(pMeditVolumeIterator self_);

void MeditVolumeIterator_free	(pMeditVolumeIterator self_)
{
  if (self_)
    {
      self_->m_medit = 0;
    }
}

void MeditVolumeIterator_def(pMeditVolumeIterator self_,pMedit const medit_)
{
  eVolume volume = __eVolume_ERROR;
  ++volume;  
  const enum GmfKwdCod gmfcod = __eVolume_TO_GmfKwCod [volume];
  GmfGotoKwd(self_->m_medit->inm,gmfcod); 
  self_->VolumeIndex 	= -1;
  self_->VolumeKind 	= volume;
  self_->Indices	= malloc(sizeof(I)*32);
  self_->NumIndices 	= 32;
  self_->m_medit 	= medit_;
}

I MeditVolumeIterator_get(cst_pMeditVolumeIterator self_,pI const indices_)
{
  { I k;
    for (k=0;k<self_->NumIndices;++k)
      {
	indices_[k] = self_->Indices[k];
      } }
  return self_->NumIndices;
}


int MeditVolumeIterator_next(pMeditVolumeIterator self_)
{
  ++self_->VolumeIndex;
  if (self_->VolumeIndex < self_->NumVolumes)
    {
      const I volumeKind = Medit_get_volumeKind(self_->m_medit,self_->VolumeIndex);
      if (volumeKind != self_->VolumeKind)
	{
	  const enum GmfKwdCod gmfcod =  __eVolume_TO_GmfKwCod [volume];
	  GmfGotoKwd(inm,gmfcod);
	}
      self_->VolumeKind = volumeKind;
      
      int localcod;
      int ci[8];
      switch(volumeKind)
	{
	case __eVolume_TETRAHEDRON:
	  {
	    self_->NumIndices = 4;
	    GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&localcod);
	    break;
	  }
	case __eVolume_HEXAHEDRON:
	  {
	    self_->NumIndices = 8;
	    GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&ci[6],&ci[7],&localcod);
	    break;
	  }
	case __eVolume_PYRAMID:
	  {
	    self_->NumIndices = 5;	    
	    GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&localcod);
	    break;
	  }
	case __eVolume_WEDGE:
	  {
	    self_->NumIndices  =6;
	    GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&ci[4],&ci[5],&localcod);
	    break;
	  }
	case __eVolume_POLYHEDRON:
	  {
	    int numFaces;
	    GmfGetLin(inm,gmfcod,&numFaces,ci,&localcod);
	    self_->NumIndices  = numFaces;
	    break;
	  }
	case __eVolume_ERROR:
	case __eVolume_ALL:
	  {
	    break;
	  }      
	}
      
      self_->VolumeCod = localcod;
      { I k;
	for (k=0;k<self_->m_numIndicesInVolume;++k)
	  {
	    self_->m_indices[k] = ci[k];
	  } }

      return 1;
    }
  else
    {
      return 0;
    }
}

pMeditVolumeIterator Medit_get_volumeIterator(pMedit 	const	self_)
{
  eVolume volume = __eVolume_ERROR;
  ++volume;
  const enum GmfKwdCod gmfcod 	=  __eVolume_TO_GmfKwCod [volume];
  GmfGotoKwd(inm,gmfcod); 
}
#endif

void Medit_get_cncFace(pMedit 	const 	self_,
		       pI 	const 	cnc_,
		       const I 		cncoff_,
		       pI 	const 	cod_,
		       const I 		codoff_,
		       cst_eFace 	face_)
{ 
  pI cnc 			= cnc_;
  pI cod 			= cod_;
  pMeditImpl 	self 		= (pMeditImpl)self_;
  int inm			= self->inm;
  const I nbFaces 		= self->nbFaces[face_];
  const enum GmfKwdCod gmfcod 	=  __eFace_TO_GmfKwCod [face_];
  GmfGotoKwd(inm,gmfcod); 
  int localcod;
  int ci[32];
  switch(face_)
    {
    case __eFace_TRIANGLE:
      {
	{ I i;
	  for (i=0;i<nbFaces;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod = localcod;
		  cod += codoff_;
		}

	    } }
	break;
      }
    case __eFace_QUADRILATERAL:
      {
	{ I i;
	  for (i=0;i<nbFaces;++i)
	    {
	      GmfGetLin(inm,gmfcod,&ci[0],&ci[1],&ci[2],&ci[3],&localcod);
	      cnc[0] = ci[0]-1;
	      cnc[1] = ci[1]-1;
	      cnc[2] = ci[2]-1;
	      cnc[3] = ci[3]-1;
	      cnc += cncoff_;
	      if (cod)
		{
		  *cod = localcod;
		  cod += codoff_;
		}

	    } }
	break;
      }
    case __eFace_POLYGON:
      {
	NotYetImplemented;
	break;
      }
    case __eFace_ERROR:
    case __eFace_ALL:
      {
	break;
      }
    }
}




void Medit_get_cncEdge(pMedit 		const	self_,
		       pI 		const	cnc_,
		       const I 			cncoff_,
		       pI 		const 	cod_,
		       const I 			codoff_)
{
  pI 		cnc		= cnc_;
  pI 		cod		= cod_;
  pMeditImpl 	self 		= (pMeditImpl)self_;
  int 		inm		= self->inm;
  const I 	nbEdges		= self->nbEdges;
  const enum GmfKwdCod gmfcod 	= GmfEdges;
  GmfGotoKwd(inm,gmfcod); 
  int localcod;
  int ci[32];
  { I i;
    for (i=0;i<nbEdges;++i)
      {
	GmfGetLin(inm,gmfcod,
		  &ci[0],
		  &ci[1],
		  &localcod);
	cnc[0] = ci[0]-1;
	cnc[1] = ci[1]-1;
	cnc += cncoff_;
	if (cod)
	  {
	    *cod = localcod;
	    cod += codoff_;
	  }	
      } }
}



unsigned char * Medit_get_adjcncs		(pMedit 	self_,
						 const L 	CIndexation_,
						 pI 		cnc_,
						 cst_pI		cncoff_,
						 pI 		cod_,
						 cst_pI		codoff_,
						 pI 		adj_,
						 cst_pI 	adjoff_)
{
  const I cncIndexation = (CIndexation_) ? ((I)1) : ((I)0);
  pMeditImpl 	self 	= (pMeditImpl)self_;
  int cod;
  int ci[8];
  cst_eDim dim = self->dim;
  int inm	= self->inm;
  const I ns 	= self->nbFaces[__eFace_TRIANGLE];
  const I cncoff = cncoff_[0];
  const I adjoff = adjoff_[0];
  
  unsigned char * medit_voy = NULL;    
  if (dim==__eDim_2)
    {
      GmfGotoKwd(inm,GmfTriangles); 
      
      int medit_ns 	= (int)ns;
      int * medit_adj 	= NULL;
      int * medit_cnc	= (int*)malloc(sizeof(int)*3*medit_ns);
      if (cod_)
	{
	  const I codoff = codoff_[0];
	  { I i;
	    for (i=0;i<ns;++i)
	      {
		GmfGetLin(inm,GmfTriangles,&medit_cnc[i*3+0],&medit_cnc[i*3+1],&medit_cnc[i*3+2],&cod);
		cnc_[cncoff * i + 0] = medit_cnc[i*3+0]-cncIndexation;
		cnc_[cncoff * i + 1] = medit_cnc[i*3+1]-cncIndexation;
		cnc_[cncoff * i + 2] = medit_cnc[i*3+2]-cncIndexation;
		cod_[codoff * i + 0] = cod;
	      } }	  
	}
      else
	{
	  { I i;
	    for (i=0;i<ns;++i)
	      {
		GmfGetLin(inm,GmfTriangles,&medit_cnc[i*3+0],&medit_cnc[i*3+1],&medit_cnc[i*3+2],&cod);
		cnc_[cncoff * i + 0] = medit_cnc[i*3+0]-cncIndexation;
		cnc_[cncoff * i + 1] = medit_cnc[i*3+1]-cncIndexation;
		cnc_[cncoff * i + 2] = medit_cnc[i*3+2]-cncIndexation;
	      } }	  
	}

      if (ns > 1)
	{
	  hash(medit_ns,&medit_adj,&medit_voy,medit_cnc);
	  free(medit_cnc);
	  medit_cnc = NULL;
	  { I k;
	    for (k=0; k<ns; k++) 
	      {
		adj_[adjoff*k+0]=medit_adj[1+3*k+0];
		adj_[adjoff*k+1]=medit_adj[1+3*k+1];
		adj_[adjoff*k+2]=medit_adj[1+3*k+2]; 
	      } }
	  free(medit_adj);
	  medit_adj = NULL;
	}
      else
	{
	  { I k;
	    for (k=0; k<ns; k++) 
	      {
		adj_[adjoff*k+0]=0;
		adj_[adjoff*k+1]=0;
		adj_[adjoff*k+2]=0; 
	      } }
	}
    }
  else if (dim==__eDim_3)
    {
      GmfGotoKwd(inm,GmfTetrahedra);   
      if (cod_)
	{
	  const I codoff = codoff_[0];
	  { I i;
	    for (i=0;i<ns;++i)
	      {
		GmfGetLin(inm,GmfTetrahedra,&ci[0],&ci[1],&ci[2],&ci[3],&cod);
		cnc_[cncoff * i + 0] = ci[0]-cncIndexation;
		cnc_[cncoff * i + 1] = ci[1]-cncIndexation;
		cnc_[cncoff * i + 2] = ci[2]-cncIndexation;
		cnc_[cncoff * i + 3] = ci[3]-cncIndexation;
		cod_[codoff * i + 0] = cod;
	      } }
	}
      else
	{
	  { I i;
	    for (i=0;i<ns;++i)
	      {
		GmfGetLin(inm,GmfTetrahedra,&ci[0],&ci[1],&ci[2],&ci[3],&cod);
		cnc_[cncoff * i + 0] = ci[0]-cncIndexation;
		cnc_[cncoff * i + 1] = ci[1]-cncIndexation;
		cnc_[cncoff * i + 2] = ci[2]-cncIndexation;
		cnc_[cncoff * i + 3] = ci[3]-cncIndexation;
	      } }
	}
    } 
  return medit_voy;
}

