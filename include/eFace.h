#ifndef __HEADER_eFace_H__
#define __HEADER_eFace_H__

#include "eEdge.h"

GridEnum(eFace,
	   __eFace_TRIANGLE,	/**< value for a triangle  */
	   __eFace_QUADRILATERAL,	/**< value for a quad */
	   __eFace_POLYGON	/**< value for a polygon */);

#ifdef __cplusplus
extern "C"
{
#endif
  eTopology	eFace_get_eTopology		(cst_eFace 	self_);
  eFace		eTopology_get_eFace		(cst_eTopology	self_);
  I 		eFace_get_nbNodes		(cst_eFace 	self_);
  I 		eFace_get_nbEdges		(cst_eFace 	self_);
  void 		eFace_get_edgeConnectivity	(cst_eFace 	self_,
						 const I 	iedge_,
						 pI 		locncedge_);

  void 		eFace_orientation		(cst_eFace	face_,
						 cst_pI		icnc_,
						 cst_pI		jcnc_,
						 pS 		irot_);

  void  	eFace_orientationPermutation	(cst_eFace 	face_,
						 const I 	sign_,
						 const I	degree_,
						 pI 		permloc_);

#ifdef __cplusplus
}
#endif

#endif




