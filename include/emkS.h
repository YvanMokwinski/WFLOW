#ifndef __HEADER_EMKS_H__
#define __HEADER_EMKS_H__

#include "emkZELM.h"
#include "emkCONTINUITY.h"

#include "emk_template.h"

#if __mok_pure_enum__

enum  __emkS_FAMILY { __emkS_FAMILY_error=0,
		       __emkS_FAMILY_lagrange,
		       __emkS_FAMILY_canonic,
		       __emkS_FAMILY_l2orthonormal,
		       __emkS_FAMILY_lagrangebubble,
		       __emkS_FAMILY_n };

typedef enum __emkS_FAMILY emkS_FAMILY;

#else

#define   __emkS_FAMILY_error 		0
#define   __emkS_FAMILY_lagrange 	1
#define   __emkS_FAMILY_canonic 	2
#define   __emkS_FAMILY_l2orthonormal 	3
#define   __emkS_FAMILY_lagrangebubble 4
#define   __emkS_FAMILY_n 		5
typedef nsINT emkS_FAMILY;

#endif

typedef const emkS_FAMILY cst_emkS_FAMILY;

__emk_template_headers(emkS_FAMILY);

#define __mmk_s_family_lagrange_name 		"Lagrange"
#define __mmk_s_family_canonic_name 		"Canonic"
#define __mmk_s_family_l2orthonormal_name 	"L2Orthonormal"
#define __mmk_s_family_lagrangebubble_name 	"LagrangeBubble"

#if __mok_pure_enum__

enum  __emk_s_diff { 
  __emk_s_diff_no=0,
  __emk_s_diff_dx,
  __emk_s_diff_dy,
  __emk_s_diff_dz,
  __emk_s_diff_n };

typedef enum __emk_s_diff emk_s_diff;

#else

#define   __emk_s_diff_no 0
#define   __emk_s_diff_dx 1
#define   __emk_s_diff_dy 2
#define   __emk_s_diff_dz 3
#define   __emk_s_diff_n  4

typedef nsINT emk_s_diff;

#endif

typedef const emk_s_diff cst_emk_s_diff;

__emk_template_headers(emk_s_diff);

enum __emk_s { __emk_s_error = 0, 					    
	       __emk_s_continuous_lagrange_interval_02,
	       __emk_s_continuous_lagrange_interval_03,
	       __emk_s_continuous_lagrange_interval_04,
	       __emk_s_continuous_lagrange_interval_05,
	       __emk_s_continuous_lagrange_interval_06,
	       __emk_s_continuous_lagrange_interval_07,
	       __emk_s_continuous_lagrange_interval_08,
	       __emk_s_continuous_lagrange_triangle_03,
	       __emk_s_continuous_lagrangebubble_triangle_04,
	       __emk_s_continuous_lagrange_triangle_06,
	       __emk_s_continuous_lagrangebubble_triangle_07,
	       __emk_s_continuous_lagrange_triangle_10,
	       __emk_s_continuous_lagrange_triangle_15,
	       __emk_s_continuous_lagrange_triangle_21,
	       __emk_s_continuous_lagrange_triangle_28,
	       __emk_s_continuous_lagrange_triangle_35,
	       __emk_s_continuous_lagrange_quadrangle_04,
	       __emk_s_continuous_lagrange_quadrangle_09,
	       __emk_s_continuous_lagrange_quadrangle_16,
	       __emk_s_continuous_lagrange_quadrangle_25,
	       __emk_s_continuous_lagrange_quadrangle_36,
	       __emk_s_continuous_lagrange_quadrangle_49,
	       __emk_s_continuous_lagrange_quadrangle_64,
	       __emk_s_continuous_lagrange_tetrahedra_04,
	       __emk_s_continuous_lagrange_tetrahedra_10,

	       __emk_s_discontinuous_lagrange_interval_01,
	       __emk_s_discontinuous_lagrange_interval_02,
	       __emk_s_discontinuous_lagrange_interval_03,
	       __emk_s_discontinuous_lagrange_interval_04,
	       __emk_s_discontinuous_lagrange_interval_05,
	       __emk_s_discontinuous_lagrange_interval_06,
	       __emk_s_discontinuous_lagrange_interval_07,
	       __emk_s_discontinuous_lagrange_interval_08,

	       __emk_s_discontinuous_lagrange_triangle_01,
	       __emk_s_discontinuous_lagrange_triangle_03,
	       __emk_s_discontinuous_lagrange_triangle_06,
	       __emk_s_discontinuous_lagrange_triangle_10,
	       __emk_s_discontinuous_lagrange_triangle_15,
	       __emk_s_discontinuous_lagrange_triangle_21,
	       __emk_s_discontinuous_lagrange_triangle_28,
	       __emk_s_discontinuous_lagrange_triangle_35,

	       __emk_s_discontinuous_lagrange_quadrangle_01,
	       __emk_s_discontinuous_lagrange_quadrangle_04,
	       __emk_s_discontinuous_lagrange_quadrangle_09,
	       __emk_s_discontinuous_lagrange_quadrangle_16,
	       __emk_s_discontinuous_lagrange_quadrangle_25,
	       __emk_s_discontinuous_lagrange_quadrangle_36,
	       __emk_s_discontinuous_lagrange_quadrangle_49,
	       __emk_s_discontinuous_lagrange_quadrangle_64,

	       __emk_s_discontinuous_lagrange_tetrahedra_01,
 	       __emk_s_discontinuous_lagrange_tetrahedra_04,
	       __emk_s_discontinuous_lagrange_tetrahedra_10,

	       __emk_s_discontinuous_canonic_interval_01,
	       __emk_s_discontinuous_canonic_interval_02,
	       __emk_s_discontinuous_canonic_interval_03,
	       __emk_s_discontinuous_canonic_interval_04,
	       __emk_s_discontinuous_canonic_interval_05,
	       __emk_s_discontinuous_canonic_interval_06,
	       __emk_s_discontinuous_canonic_interval_07,
	       __emk_s_discontinuous_canonic_interval_08,

	       __emk_s_discontinuous_canonic_triangle_01,
	       __emk_s_discontinuous_canonic_triangle_03,
	       __emk_s_discontinuous_canonic_triangle_06,
	       __emk_s_discontinuous_canonic_triangle_10,
	       __emk_s_discontinuous_canonic_triangle_15,
	       __emk_s_discontinuous_canonic_triangle_21,
	       __emk_s_discontinuous_canonic_triangle_28,
	       __emk_s_discontinuous_canonic_triangle_35,

	       __emk_s_discontinuous_canonic_quadrangle_01,
	       __emk_s_discontinuous_canonic_quadrangle_04,
	       __emk_s_discontinuous_canonic_quadrangle_09,
	       __emk_s_discontinuous_canonic_quadrangle_16,
	       __emk_s_discontinuous_canonic_quadrangle_25,
	       __emk_s_discontinuous_canonic_quadrangle_36,
	       __emk_s_discontinuous_canonic_quadrangle_49,
	       __emk_s_discontinuous_canonic_quadrangle_64,

	       __emk_s_discontinuous_canonic_tetrahedra_01,
 	       __emk_s_discontinuous_canonic_tetrahedra_04,
	       __emk_s_discontinuous_canonic_tetrahedra_10,	       	       

	       __emk_s_discontinuous_l2orthonormal_interval_01,
	       __emk_s_discontinuous_l2orthonormal_interval_02,
	       __emk_s_discontinuous_l2orthonormal_interval_03,
	       __emk_s_discontinuous_l2orthonormal_interval_04,
	       __emk_s_discontinuous_l2orthonormal_interval_05,
	       __emk_s_discontinuous_l2orthonormal_interval_06,
	       __emk_s_discontinuous_l2orthonormal_interval_07,
	       __emk_s_discontinuous_l2orthonormal_interval_08,

	       __emk_s_discontinuous_l2orthonormal_triangle_01,
	       __emk_s_discontinuous_l2orthonormal_triangle_03,
	       __emk_s_discontinuous_l2orthonormal_triangle_06,
	       __emk_s_discontinuous_l2orthonormal_triangle_10,
	       __emk_s_discontinuous_l2orthonormal_triangle_15,
	       __emk_s_discontinuous_l2orthonormal_triangle_21,
	       __emk_s_discontinuous_l2orthonormal_triangle_28,
	       __emk_s_discontinuous_l2orthonormal_triangle_35,

	       __emk_s_discontinuous_l2orthonormal_quadrangle_01,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_04,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_09,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_16,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_25,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_36,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_49,
	       __emk_s_discontinuous_l2orthonormal_quadrangle_64,

	       __emk_s_discontinuous_l2orthonormal_tetrahedra_01,
 	       __emk_s_discontinuous_l2orthonormal_tetrahedra_04,
	       __emk_s_discontinuous_l2orthonormal_tetrahedra_10,	       	       

	       __emk_s_n }; 

typedef enum __emk_s emk_s;
typedef const enum __emk_s cst_emk_s;

__emk_template_headers(emk_s);


cst_emk_s mkS_emk_s_enum(cst_emkCONTINUITY 	continuity_,
			 cst_emkS_FAMILY 	family_,
			 cst_emkZELM 		element_,
			 const nsINT 		nshape_);




emkS_FAMILY 	__emk_s_emkS_FAMILY	(cst_emk_s kind_);
emkCONTINUITY 	__emk_s_emkCONTINUITY	(cst_emk_s kind_);
emkZELM 	__emk_s_emkZELM		(cst_emk_s kind_);
nsINT 		__emk_s_degree		(cst_emk_s kind_);

nsINT 		__emkS_FAMILY_nshape	(cst_emkS_FAMILY family_,
					 cst_emkZELM 	  elm_,
					 const nsINT 	  degree_);

#endif


#if 0

emk_s 			__emk_s_enum		(cst_emkCONTINUITY 	continuity_,
						 cst_emkS_FAMILY 	family_,
						 cst_emkZELM 		elm_,
						 const nsINT 		degree_);
const char * 		__emk_s_string		(cst_emk_s 		flag_);
#endif
#if 0
const char * 		__emk_s_diff_string	(cst_emk_s_diff diff_);
emk_s_diff 		__emk_s_diff_enum	(const char * name_);
const char * 		__emkS_FAMILY_string 		(cst_emkS_FAMILY family_);
emkS_FAMILY		__emkS_FAMILY_enum		(const char * name_);
void 			__emkS_FAMILY_strings_copy	(const char * names_[__emkS_FAMILY_n]);
#endif
