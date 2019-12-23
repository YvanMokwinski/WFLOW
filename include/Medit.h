#ifndef __header_Medit_h__
#define __header_Medit_h__


#include "ePrecision.h"
#include "Points.h"
#include "eVolume.h"
#include "eTopologyDimension.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
} Medit,* RESTRICT pMedit;

typedef const Medit 	* RESTRICT 		cst_pMedit;

  
#if 0
  typedef struct
{
  I 		NumVolumes;
  I 		VolumeIndex;
  I 		VolumeKind;
  I 		VolumeCod;
  pMedit 	m_medit;
  pI		Indices;
  pI		NumIndices;
} MeditVolumeIterator,*RESTRICT pMeditVolumeIterator;
#endif

unsigned char * Medit_get_adjcncs		(pMedit 	self_,
						 const L 	CIndexation_,
						 pI 		cnc_,
						 cst_pI		cncoff_,
						 pI 		cod_,
						 cst_pI		codoff_,
						 pI 		adj_,
						 cst_pI 	adjoff_);

void Medit_get_cellToNodes	(pMedit const self_,const eVolume volume_,pI const cellToNodeIndices_);
void Medit_gotoSectionCell(pMedit const self_,const eVolume volume_);
  void Medit_gotoSectionVertices(pMedit 	const 	self_);
  void Medit_get_GridVertex(pMedit 	const 	self_,
			    pR 		const 	coo_);
  
  eVolume 		Medit_get_volumeKind(cst_pMedit 	const 		self_,
					     const I 				volumeIndex_);

eTopologyDimension 	Medit_get_topologyDimension	(cst_pMedit 	const	self_);

eDim 			Medit_get_eDim			(cst_pMedit 	const	self_);

I 			Medit_get_nbVertices		(cst_pMedit 	const 	self_);
I 			Medit_get_nbEdges		(cst_pMedit 	const 	self_);

I 			Medit_get_nbFaces		(cst_pMedit 	const 	self_,
							 cst_eFace 		face_);

I 			Medit_get_nbVolumes		(cst_pMedit 	const 	self_,
							 cst_eVolume 		volume_);

void 			Medit_get_cncVolume		(pMedit 	const 	self_,
							 pI 		const 	cnc_,
							 const I 		cncoff_,
							 pI 		const 	cod_,
							 const I 		codoff_,
							 cst_eVolume 		volume_);

void 			Medit_get_cncFace		(pMedit 	const 	self_,
							 pI 		const	cnc_,
							 const I 		cncoff_,
							 pI 		const 	cod_,
							 const I 		codoff_,
							 cst_eFace 		face_);


void 			Medit_get_cncEdge		(pMedit 	const 	self_,
							 pI 		const	cnc_,
							 const I 		cncoff_,
							 pI 		const 	cod_,
							 const I 		codoff_);


pMedit  		Medit_kill			(pMedit		const	self_);

pMedit  		Medit_new			(cst_pS		const 	filename_,
							 ...);

void  			Medit_set_precision		(pMedit 	const 	self_,
							 cst_ePrecision 	precision_);

void 			Medit_get_Points		(pMedit 	const	self_,
							 pPoints 	const 	xyz_,
							 pI		const 	cod_,
							 cst_pI		const 	codoff_);

void 			Medit_get_cncv			(pMedit 	const self_,
							 cst_pI		cootr_,
							 pR 		coo_,
							 cst_pI		cooff_,
							 pI		cod_,
							 cst_pI		codoff_);


  

#ifdef __cplusplus
}
#endif
#endif
