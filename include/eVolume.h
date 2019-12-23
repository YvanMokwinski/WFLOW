#ifndef __HEADER_eVolume_H__
#define __HEADER_eVolume_H__

#include "eFace.h"

GridEnum(eVolume,
	 __eVolume_TETRAHEDRON,		/**< value for a tetrahedra */
	 __eVolume_HEXAHEDRON,		/**< value for a hexahedra */
	 __eVolume_WEDGE,		/**< value for a wedge */
	 __eVolume_PYRAMID,		/**< value for a pyramid */
	 __eVolume_POLYHEDRON		/**< value for a polyhedron */ );

#ifdef __cplusplus
extern "C"
{
#endif

  eTopology	eVolume_get_eTopology		(cst_eVolume 		self_);
  eVolume	eTopology_get_eVolume		(cst_eTopology		self_);

  I 		eVolume_get_maxNbNodesInFace	(cst_eVolume 		self_);
  I 		eVolume_get_nbNodes		(cst_eVolume 		self_);
  I 		eVolume_get_nbFaces		(cst_eVolume 		self_);
  I 		eVolume_get_nbEdges		(cst_eVolume 		self_);
  eFace		eVolume_get_eFace		(cst_eVolume 		self_,
						 const I 		iface_);
  
  I 		eVolume_get_nbKindOfFaces	(cst_eVolume 		self_);
  void 		eVolume_get_faceConnectivity	(cst_eVolume 		self_,
						 const I 		ilocface_,
						 pI 		const 	cnc_);
  void 		eVolume_get_edgeConnectivity	(cst_eVolume 		self_,
						 const I 		ilocedge_,
						 pI 		const 	cnc_);
#ifdef __cplusplus
}
#endif

#endif




