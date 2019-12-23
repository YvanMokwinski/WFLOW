#ifndef __HEADER_EMK_ZELM_H__
#define __HEADER_EMK_ZELM_H__

#include "emkDIM.h"

#if __mok_pure_enum__

enum __emkZELM 	{   __emkZELM_error = 0, 		      
		    __emkZELM_interval,
		    __emkZELM_triangle, 
		    __emkZELM_quadrangle, 
		    __emkZELM_tetrahedra, 		      
		    __emkZELM_n }; 

typedef enum __emkZELM emkZELM;

#else

#define __emkZELM_error      	0
#define __emkZELM_interval   	1
#define __emkZELM_triangle   	2
#define __emkZELM_quadrangle 	3
#define __emkZELM_tetrahedra 	4
#define __emkZELM_n 		5

typedef nsINT emkZELM;

#endif

typedef const emkZELM cst_emkZELM;

#include "emk_template.h"

__emk_template_headers(emkZELM);

#if __mk_debug__

emkDIM      	__emkZELM_mindim	(cst_emkZELM elm_);
nsINT        	__emkZELM_nedge		(cst_emkZELM elm_);
nsINT        	__emkZELM_nface		(cst_emkZELM elm_);
nsINT        	__emkZELM_nvertex	(cst_emkZELM elm_);
emkZELM 	__emkZELM_kind_face	(cst_emkZELM elm_);
cst_mkR 	__emkZELM_refcoo	(cst_emkZELM elm_);
emk_logic 	__emkZELM_isSimplex	(cst_emkZELM elm_);

#else

extern cst_emk_logic 	__vmkZELM_isSimplex	[__emkZELM_n];
extern cst_emkDIM 	__vmkZELM_minimum_dimension_to_be_defined[__emkZELM_n];
extern const nsINT	__vmkZELM_nedge		[__emkZELM_n];
extern const nsINT	__vmkZELM_nface		[__emkZELM_n];
extern const nsINT	__vmkZELM_nvertex	[__emkZELM_n];
extern cst_emkZELM	__vmkZELM_kind_face	[__emkZELM_n];
extern const cst_mkR   	__vmkZELM_refcoo	[__emkZELM_n];

#define __emkZELM_mindim(_e)    __vmkZELM_minimum_dimension_to_be_defined[(_e)]
#define __emkZELM_nedge(_e) 	__vmkZELM_nedge	[(_e)]
#define __emkZELM_nface(_e) 	__vmkZELM_nface	[(_e)]
#define __emkZELM_nvertex(_e) 	__vmkZELM_nvertex[(_e)]
#define __emkZELM_kind_face(_e)	__vmkZELM_kind_face[(_e)]
#define __emkZELM_refcoo(_e) 	__vmkZELM_refcoo[(_e)]
#define __emkZELM_refcoo(_e) 	__vmkZELM_refcoo[(_e)]
#define __emkZELM_isSimplex(_e) __vmkZELM_isSimplex[(_e)]

#endif

#define __mmkZELM_triangle_name  	"triangle"
#define __mmkZELM_quadrangle_name  	"quadrangle"
#define __mmkZELM_interval_name  	"interval"
#define __mmkZELM_tetrahedra_name  	"tetrahedra"


#define __mmkZELM_interval_refcoo_x0 ((nsREAL)-1.0)
#define __mmkZELM_interval_refcoo_x1 ((nsREAL)1.0)

#define __mmkZELM_triangle_refcoo_x0 ((nsREAL)0.0)
#define __mmkZELM_triangle_refcoo_y0 ((nsREAL)0.0)
#define __mmkZELM_triangle_refcoo_x1 ((nsREAL)1.0)
#define __mmkZELM_triangle_refcoo_y1 ((nsREAL)0.0)
#define __mmkZELM_triangle_refcoo_x2 ((nsREAL)0.0)
#define __mmkZELM_triangle_refcoo_y2 ((nsREAL)1.0)

#define __mmkZELM_quadrangle_refcoo_x0 ((nsREAL)1.0)
#define __mmkZELM_quadrangle_refcoo_y0 ((nsREAL)1.0)
#define __mmkZELM_quadrangle_refcoo_x1 ((nsREAL)-1.0)
#define __mmkZELM_quadrangle_refcoo_y1 ((nsREAL)1.0)
#define __mmkZELM_quadrangle_refcoo_x2 ((nsREAL)-1.0)
#define __mmkZELM_quadrangle_refcoo_y2 ((nsREAL)-1.0)
#define __mmkZELM_quadrangle_refcoo_x3 ((nsREAL)1.0)
#define __mmkZELM_quadrangle_refcoo_y3 ((nsREAL)-1.0)

#define __mmkZELM_tetrahedra_refcoo_x0 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_y0 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_z0 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_x1 ((nsREAL)1.0)
#define __mmkZELM_tetrahedra_refcoo_y1 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_z1 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_x2 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_y2 ((nsREAL)1.0)
#define __mmkZELM_tetrahedra_refcoo_z2 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_x3 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_y3 ((nsREAL)0.0)
#define __mmkZELM_tetrahedra_refcoo_z3 ((nsREAL)1.0)


emk_logic  		__emkZELM_contains			(cst_mkR 	cooelm,
								 cst_mkR 	coo_);
void 			__emkZELM_belement_triangle		(cst_mkR 	array_coo,
								 mkR 		Belement);
void 			__emkZELM_bmapping			(cst_emkZELM 	elm_,
								 cst_mkI 	n_,
								 mkR 		r_,
								 cst_mkI 	roff_,
								 cst_mkR 	p_,
								 emk_err*	err);
void 			__emkZELM_cncedge			(cst_emkZELM 	elm_,
								 mkI 		loc_nedge, 
								 nsINT 		loc_edge[16]);


#endif
