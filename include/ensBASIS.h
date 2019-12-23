#ifndef __header_ensBASIS_h__
#define __header_ensBASIS_h__

#include "eFace.h"
#define __NS_SHAPE_NMAX__ 64

#include "ConfigEnum.h"

ConfigEnum(ensBASIS,
	   __ensBASIS_LAGRANGE_0,
	   __ensBASIS_LAGRANGE_1,
	   __ensBASIS_LAGRANGE_2,
	   __ensBASIS_LAGRANGE_3,
	   __ensBASIS_L2ORTHONORMAL_0,
	   __ensBASIS_L2ORTHONORMAL_1,
	   __ensBASIS_L2ORTHONORMAL_2,
	   __ensBASIS_L2ORTHONORMAL_3,
	   __ensBASIS_LAGRANGEBUBBLE_0,
	   __ensBASIS_LAGRANGEBUBBLE_1,
	   __ensBASIS_LAGRANGEBUBBLE_2);

/*
enum __ensBASIS 		{ __ensBASIS_error=0,
				  __ensBASIS_LAGRANGE_0,
				  __ensBASIS_LAGRANGE_1,
				  __ensBASIS_LAGRANGE_2,
				  __ensBASIS_LAGRANGE_3,
				  __ensBASIS_L2ORTHONORMAL_0,
				  __ensBASIS_L2ORTHONORMAL_1,
				  __ensBASIS_L2ORTHONORMAL_2,
				  __ensBASIS_L2ORTHONORMAL_3,
				  __ensBASIS_LAGRANGEBUBBLE_0,
				  __ensBASIS_LAGRANGEBUBBLE_1,
				  __ensBASIS_LAGRANGEBUBBLE_2,			    
				  __ensBASIS_n};
typedef enum 		__ensBASIS ensBASIS;
typedef const enum 	__ensBASIS cst_ensBASIS;
*/

#ifdef __cplusplus
extern "C"
{
#endif
  I	 	ensBASIS_n			(cst_ensBASIS 	b_,
						 cst_eFace 	elm_);
  I 		ensBASIS_n_onface		(cst_ensBASIS 	basis_,
						 cst_eFace 	elm_);

  cst_pI 	ensBASIS_cnctreilli		(cst_ensBASIS 	shape_,
						 cst_eFace 	elm_);

  I  		ensBASIS_nelmtreilli		(cst_ensBASIS 	shape_,
						 cst_eFace 	elm_);

  cst_pR 	ensBASIS_cootreilli		(cst_ensBASIS 	shape_,
						 cst_eFace 	elm_);

  cst_pR 	ensBASIS_lc2gltreilli_ptr	(cst_ensBASIS 	shape_,
						 cst_eFace 	elm_);

  void 		ensBASIS_lc2gltreilli		(cst_ensBASIS 	shape_,
						 cst_eFace 	elm_,
						 cst_pR 	cooelm_,
						 cst_pI 	cooelmoff_,
						 pR 		gl_,
						 cst_pI 	gloff_);

  I	 	ensBASIS_degree			(cst_ensBASIS b_);

  ensBASIS 	ensBASIS_down			(cst_ensBASIS b_);
  ensBASIS 	ensBASIS_up			(cst_ensBASIS b_);


#define COMMON_HEADER    (cst_ensBASIS basis_,cst_eFace elm_,cst_pS ptr_,cst_pS btr_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_)
  void 	 	ensBASIS_b			COMMON_HEADER;
  void 	 	ensBASIS_dr			COMMON_HEADER;
  void 	 	ensBASIS_ds			COMMON_HEADER;
  void 	 	ensBASIS_dt			COMMON_HEADER;
#undef COMMON_HEADER

#ifdef __cplusplus
};
std::ostream& operator<<(std::ostream&s,cst_ensBASIS & b);
#endif

#endif
