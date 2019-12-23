#ifndef __HEADER_eTopology_H__
#define __HEADER_eTopology_H__

#include "GridEnum.h"
#include "eTopologyDimension.h"

GridEnum(eTopology,
	 __eTopology_NODE,
	 __eTopology_SEGMENT,
	 __eTopology_POLYLINE,
	 __eTopology_TRIANGLE,
	 __eTopology_QUADRILATERAL,
	 __eTopology_POLYGON,
	 __eTopology_TETRAHEDRON,
	 __eTopology_HEXAHEDRON,
	 __eTopology_WEDGE,
	 __eTopology_PYRAMID,
	 __eTopology_POLYHEDRON);

#ifdef __cplusplus
extern "C"
{
#endif

  eTopologyDimension 	eTopology_get_eTopologyDimension		(cst_eTopology 		self_);

#ifdef __cplusplus
}
#endif

#endif




