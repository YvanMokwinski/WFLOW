#ifndef __HEADER_eTopologyDimension_H__
#define __HEADER_eTopologyDimension_H__

typedef enum __eTopologyDimension{__eTopologyDimension_Point=0,
				  __eTopologyDimension_Edge,
				  __eTopologyDimension_Face,
				  __eTopologyDimension_Volume} eTopologyDimension;

typedef const eTopologyDimension cst_eTopologyDimension;

#endif




