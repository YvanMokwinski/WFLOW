#ifndef __HEADER_eEdge_H__
#define __HEADER_eEdge_H__

#include "eTopology.h"

GridEnum(eEdge,__eEdge_SEGMENT,__eEdge_POLYLINE);

eTopology	eEdge_get_eTopology		(cst_eEdge 		self_);
eEdge		eTopology_get_eEdge		(cst_eTopology 		self_);

#ifdef __cplusplus
extern "C"
{
#endif
  
#ifdef __cplusplus
}
#endif
  
#endif




