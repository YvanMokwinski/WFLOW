#include "mkS.h"


I 		__emk_s_degree(cst_emk_s kind_)
{
  I degree=(I)-1;
  switch(kind_)
    {
    case __emk_s_discontinuous_lagrange_interval_01:
    case __emk_s_discontinuous_lagrange_triangle_01:
    case __emk_s_discontinuous_lagrange_quadrangle_01:
    case __emk_s_discontinuous_lagrange_tetrahedra_01:
    case __emk_s_discontinuous_canonic_interval_01:
    case __emk_s_discontinuous_canonic_triangle_01:
    case __emk_s_discontinuous_canonic_quadrangle_01:
    case __emk_s_discontinuous_canonic_tetrahedra_01:
    case __emk_s_discontinuous_l2orthonormal_interval_01:
    case __emk_s_discontinuous_l2orthonormal_triangle_01:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_01:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_01:
      {
	degree=0;
	break;
      }

    case __emk_s_discontinuous_l2orthonormal_tetrahedra_04:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_04:
    case __emk_s_discontinuous_l2orthonormal_triangle_03:
    case __emk_s_discontinuous_l2orthonormal_interval_02:
    case __emk_s_discontinuous_canonic_quadrangle_04:
    case __emk_s_discontinuous_lagrange_interval_02:
    case __emk_s_discontinuous_lagrange_triangle_03:
    case __emk_s_discontinuous_lagrange_quadrangle_04:
    case __emk_s_discontinuous_lagrange_tetrahedra_04:
    case __emk_s_discontinuous_canonic_interval_02:
    case __emk_s_discontinuous_canonic_triangle_03:
    case __emk_s_discontinuous_canonic_tetrahedra_04:
    case __emk_s_continuous_lagrange_tetrahedra_04:
    case __emk_s_continuous_lagrange_quadrangle_04:
    case __emk_s_continuous_lagrange_triangle_03:
    case __emk_s_continuous_lagrange_interval_02:
    case __emk_s_continuous_lagrangebubble_triangle_04:
      {
	degree=1;
	break;
      }
    case __emk_s_continuous_lagrange_tetrahedra_10:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_10:	       	       
    case __emk_s_continuous_lagrange_quadrangle_09:
    case __emk_s_continuous_lagrange_interval_03:
    case __emk_s_continuous_lagrange_triangle_06:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_09:
    case __emk_s_discontinuous_lagrange_triangle_06:
    case __emk_s_discontinuous_canonic_tetrahedra_10:	       	       
    case __emk_s_discontinuous_lagrange_quadrangle_09:
    case __emk_s_discontinuous_l2orthonormal_triangle_06:
    case __emk_s_discontinuous_canonic_interval_03:
    case __emk_s_discontinuous_canonic_triangle_06:
    case __emk_s_discontinuous_canonic_quadrangle_09:
    case __emk_s_discontinuous_lagrange_tetrahedra_10:
    case __emk_s_discontinuous_l2orthonormal_interval_03:
    case __emk_s_discontinuous_lagrange_interval_03:
    case __emk_s_continuous_lagrangebubble_triangle_07:

      {
	degree=2;
	break;
      }
    case __emk_s_continuous_lagrange_interval_04:
    case __emk_s_discontinuous_l2orthonormal_interval_04:
    case __emk_s_discontinuous_canonic_interval_04:
    case __emk_s_discontinuous_lagrange_triangle_10:
    case __emk_s_discontinuous_canonic_quadrangle_16:
    case __emk_s_continuous_lagrange_triangle_10:
    case __emk_s_discontinuous_lagrange_quadrangle_16:
    case __emk_s_discontinuous_canonic_triangle_10:
    case __emk_s_continuous_lagrange_quadrangle_16:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_16:
    case __emk_s_discontinuous_lagrange_interval_04:
    case __emk_s_discontinuous_l2orthonormal_triangle_10:
      {
	degree=3;
	break;
      }
      
    case __emk_s_discontinuous_lagrange_quadrangle_25:
    case __emk_s_continuous_lagrange_quadrangle_25:
    case __emk_s_discontinuous_canonic_quadrangle_25:
    case __emk_s_discontinuous_l2orthonormal_interval_05:
    case __emk_s_continuous_lagrange_triangle_15:
    case __emk_s_continuous_lagrange_interval_05:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_25:
    case __emk_s_discontinuous_lagrange_interval_05:
    case __emk_s_discontinuous_canonic_triangle_15:
    case __emk_s_discontinuous_lagrange_triangle_15:
    case __emk_s_discontinuous_l2orthonormal_triangle_15:
    case __emk_s_discontinuous_canonic_interval_05:
      {
	degree=4;
	break;
      }
    case __emk_s_continuous_lagrange_quadrangle_36:
    case __emk_s_discontinuous_canonic_quadrangle_36:
    case __emk_s_discontinuous_lagrange_quadrangle_36:
    case __emk_s_continuous_lagrange_triangle_21:
    case __emk_s_discontinuous_l2orthonormal_interval_06:
    case __emk_s_continuous_lagrange_interval_06:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_36:
    case __emk_s_discontinuous_lagrange_triangle_21:
    case __emk_s_discontinuous_l2orthonormal_triangle_21:
    case __emk_s_discontinuous_canonic_triangle_21:
    case __emk_s_discontinuous_lagrange_interval_06:
    case __emk_s_discontinuous_canonic_interval_06:
      {
	degree=5;
	break;
      }
    case __emk_s_continuous_lagrange_quadrangle_49:
    case __emk_s_continuous_lagrange_triangle_28:
    case __emk_s_discontinuous_l2orthonormal_interval_07:
    case __emk_s_continuous_lagrange_interval_07:
    case __emk_s_discontinuous_canonic_quadrangle_49:
    case __emk_s_discontinuous_lagrange_quadrangle_49:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_49:
    case __emk_s_discontinuous_lagrange_interval_07:
    case __emk_s_discontinuous_canonic_interval_07:
    case __emk_s_discontinuous_canonic_triangle_28:
    case __emk_s_discontinuous_lagrange_triangle_28:
    case __emk_s_discontinuous_l2orthonormal_triangle_28:
      {
	degree=6;
	break;
      }
    case __emk_s_continuous_lagrange_quadrangle_64:
    case __emk_s_discontinuous_lagrange_quadrangle_64:
    case __emk_s_discontinuous_canonic_quadrangle_64:
    case __emk_s_discontinuous_lagrange_triangle_35:
    case __emk_s_discontinuous_canonic_triangle_35:
    case __emk_s_continuous_lagrange_triangle_35:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_64:
    case __emk_s_discontinuous_lagrange_interval_08:
    case __emk_s_discontinuous_l2orthonormal_triangle_35:
    case __emk_s_discontinuous_l2orthonormal_interval_08:
    case __emk_s_discontinuous_canonic_interval_08:
    case __emk_s_continuous_lagrange_interval_08:
      {
	degree=7;
	break;
      }
    
    case __emk_s_error:
    case __emk_s_n :
      {
	fprintf(stderr,"__emkS_degree:switch failed on __emk_s\n");
	exit(1);
	break;
      }
    }
  return degree;
}

#define __mmk_s_family_lagrange_name 		"Lagrange"
#define __mmk_s_family_canonic_name 		"Canonic"
#define __mmk_s_family_l2orthonormal_name 	"L2Orthonormal"
#define __mmk_s_family_lagrangebubble_name 	"LagrangeBubble"

static const char * __emkS_FAMILY_strings[__emkS_FAMILY_n]={"__emkS_FAMILY_error",
							    __mmk_s_family_lagrange_name,
							    __mmk_s_family_canonic_name,
							    __mmk_s_family_l2orthonormal_name,
							    __mmk_s_family_lagrangebubble_name};

const char * __emkS_FAMILY_string(const emkS_FAMILY family_)
{
  return __emkS_FAMILY_strings[family_];
}


#define __mmk_continuity_continuous_name 	"continuous"
#define __mmk_continuity_discontinuous_name 	"discontinuous"


const char*__emkCONTINUITY_string(const emkCONTINUITY continuity_)
{
  const char * f = NULL;
  switch(continuity_)
    {
    case __emk_continuous:
      {
	f = __mmk_continuity_continuous_name;
	break;
      }
    case __emk_discontinuous:
      {
	f = __mmk_continuity_discontinuous_name;
	break;
      }
    default:
      {
	fprintf(stderr,"emkCONTINUITY_string:switch failed on __emkCONTINUITY_string\n");
	exit(1);
	break;
      }

    }
  return f;
}


#define __mmkZELM_triangle_name  	"triangle"
#define __mmkZELM_quadrangle_name  	"quadrangle"
#define __mmkZELM_interval_name  	"interval"
#define __mmkZELM_tetrahedra_name  	"tetrahedra"

const char * __emkZELM_strings [] = {"emkZELM_error",__mmkZELM_interval_name,__mmkZELM_triangle_name,__mmkZELM_quadrangle_name,__mmkZELM_tetrahedra_name };

const char * emkZELM_string(const eTopology element)
{
  switch(element)
    {
    case __eTopology_TRIANGLE:
      {
	return __emkZELM_strings[2];
      }
    default:
      {
	fprintf(stderr,"ofkofkofkdokeokoekoe\n");
	exit(1);
	break;
      }
    }
  return NULL;
}


emkS_FAMILY __emk_s_emkS_FAMILY(cst_emk_s kind_)
{
  emkS_FAMILY family = __emkS_FAMILY_error;
  switch(kind_)
    {
    case __emk_s_continuous_lagrange_triangle_03:
    case __emk_s_continuous_lagrange_triangle_06:
    case __emk_s_continuous_lagrange_triangle_10:
    case __emk_s_continuous_lagrange_triangle_15:
    case __emk_s_continuous_lagrange_triangle_21:
    case __emk_s_continuous_lagrange_triangle_28:
    case __emk_s_continuous_lagrange_triangle_35:
    case __emk_s_discontinuous_lagrange_triangle_01:
    case __emk_s_discontinuous_lagrange_triangle_03:
    case __emk_s_discontinuous_lagrange_triangle_06:
    case __emk_s_discontinuous_lagrange_triangle_10:
    case __emk_s_discontinuous_lagrange_triangle_15:
    case __emk_s_discontinuous_lagrange_triangle_21:
    case __emk_s_discontinuous_lagrange_triangle_28:
    case __emk_s_discontinuous_lagrange_triangle_35:
    case __emk_s_continuous_lagrange_interval_02:
    case __emk_s_continuous_lagrange_interval_03:
    case __emk_s_continuous_lagrange_interval_04:
    case __emk_s_continuous_lagrange_interval_05:
    case __emk_s_continuous_lagrange_interval_06:
    case __emk_s_continuous_lagrange_interval_07:
    case __emk_s_continuous_lagrange_interval_08:
    case __emk_s_discontinuous_lagrange_interval_01:
    case __emk_s_discontinuous_lagrange_interval_02:
    case __emk_s_discontinuous_lagrange_interval_03:
    case __emk_s_discontinuous_lagrange_interval_04:
    case __emk_s_discontinuous_lagrange_interval_05:
    case __emk_s_discontinuous_lagrange_interval_06:
    case __emk_s_discontinuous_lagrange_interval_07:
    case __emk_s_discontinuous_lagrange_interval_08:
    case __emk_s_continuous_lagrange_quadrangle_04:
    case __emk_s_continuous_lagrange_quadrangle_09:
    case __emk_s_continuous_lagrange_quadrangle_16:
    case __emk_s_continuous_lagrange_quadrangle_25:
    case __emk_s_continuous_lagrange_quadrangle_36:
    case __emk_s_continuous_lagrange_quadrangle_49:
    case __emk_s_continuous_lagrange_quadrangle_64:
    case __emk_s_discontinuous_lagrange_quadrangle_01:
    case __emk_s_discontinuous_lagrange_quadrangle_04:
    case __emk_s_discontinuous_lagrange_quadrangle_09:
    case __emk_s_discontinuous_lagrange_quadrangle_16:
    case __emk_s_discontinuous_lagrange_quadrangle_25:
    case __emk_s_discontinuous_lagrange_quadrangle_36:
    case __emk_s_discontinuous_lagrange_quadrangle_49:
    case __emk_s_discontinuous_lagrange_quadrangle_64:
    case __emk_s_continuous_lagrange_tetrahedra_04:
    case __emk_s_continuous_lagrange_tetrahedra_10:
    case __emk_s_discontinuous_lagrange_tetrahedra_01:
    case __emk_s_discontinuous_lagrange_tetrahedra_04:
    case __emk_s_discontinuous_lagrange_tetrahedra_10:
      {
	family = __emkS_FAMILY_lagrange;
	break;
      }
    case __emk_s_continuous_lagrangebubble_triangle_04:
    case __emk_s_continuous_lagrangebubble_triangle_07:
      {
	family = __emkS_FAMILY_lagrangebubble;
	break;
      }

    case __emk_s_discontinuous_canonic_interval_01:
    case __emk_s_discontinuous_canonic_interval_02:
    case __emk_s_discontinuous_canonic_interval_03:
    case __emk_s_discontinuous_canonic_interval_04:
    case __emk_s_discontinuous_canonic_interval_05:
    case __emk_s_discontinuous_canonic_interval_06:
    case __emk_s_discontinuous_canonic_interval_07:
    case __emk_s_discontinuous_canonic_interval_08:
    case __emk_s_discontinuous_canonic_triangle_01:
    case __emk_s_discontinuous_canonic_triangle_03:
    case __emk_s_discontinuous_canonic_triangle_06:
    case __emk_s_discontinuous_canonic_triangle_10:
    case __emk_s_discontinuous_canonic_triangle_15:
    case __emk_s_discontinuous_canonic_triangle_21:
    case __emk_s_discontinuous_canonic_triangle_28:
    case __emk_s_discontinuous_canonic_triangle_35:
    case __emk_s_discontinuous_canonic_quadrangle_01:
    case __emk_s_discontinuous_canonic_quadrangle_04:
    case __emk_s_discontinuous_canonic_quadrangle_09:
    case __emk_s_discontinuous_canonic_quadrangle_16:
    case __emk_s_discontinuous_canonic_quadrangle_25:
    case __emk_s_discontinuous_canonic_quadrangle_36:
    case __emk_s_discontinuous_canonic_quadrangle_49:
    case __emk_s_discontinuous_canonic_quadrangle_64:
    case __emk_s_discontinuous_canonic_tetrahedra_01:
    case __emk_s_discontinuous_canonic_tetrahedra_04:
    case __emk_s_discontinuous_canonic_tetrahedra_10:  
      {
	family = __emkS_FAMILY_canonic;
	break;
      }
    case __emk_s_discontinuous_l2orthonormal_interval_01:
    case __emk_s_discontinuous_l2orthonormal_interval_02:
    case __emk_s_discontinuous_l2orthonormal_interval_03:
    case __emk_s_discontinuous_l2orthonormal_interval_04:
    case __emk_s_discontinuous_l2orthonormal_interval_05:
    case __emk_s_discontinuous_l2orthonormal_interval_06:
    case __emk_s_discontinuous_l2orthonormal_interval_07:
    case __emk_s_discontinuous_l2orthonormal_interval_08:
    case __emk_s_discontinuous_l2orthonormal_triangle_01:
    case __emk_s_discontinuous_l2orthonormal_triangle_03:
    case __emk_s_discontinuous_l2orthonormal_triangle_06:
    case __emk_s_discontinuous_l2orthonormal_triangle_10:
    case __emk_s_discontinuous_l2orthonormal_triangle_15:
    case __emk_s_discontinuous_l2orthonormal_triangle_21:
    case __emk_s_discontinuous_l2orthonormal_triangle_28:
    case __emk_s_discontinuous_l2orthonormal_triangle_35:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_01:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_04:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_09:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_16:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_25:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_36:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_49:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_64:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_01:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_04:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_10:	       	       
      {
	family = __emkS_FAMILY_l2orthonormal;
	break;
      }
    case __emk_s_error:
    case __emk_s_n :
      {

	fprintf(stderr,"emk_s_emkS_FAMILY:switch failed on __emk_s\n");
	exit(1);
	break;
      }
    }
  return family;
}


eTopology __emk_s_emkZELM(cst_emk_s kind_)
{
  eTopology element = __eTopology_ERROR;
  switch(kind_)
    {

    case __emk_s_discontinuous_canonic_interval_01:
    case __emk_s_discontinuous_canonic_interval_02:
    case __emk_s_discontinuous_canonic_interval_03:
    case __emk_s_discontinuous_canonic_interval_04:
    case __emk_s_discontinuous_canonic_interval_05:
    case __emk_s_discontinuous_canonic_interval_06:
    case __emk_s_discontinuous_canonic_interval_07:
    case __emk_s_discontinuous_canonic_interval_08:
    case __emk_s_continuous_lagrange_interval_02:
    case __emk_s_continuous_lagrange_interval_03:
    case __emk_s_continuous_lagrange_interval_04:
    case __emk_s_continuous_lagrange_interval_05:
    case __emk_s_continuous_lagrange_interval_06:
    case __emk_s_continuous_lagrange_interval_07:
    case __emk_s_continuous_lagrange_interval_08:
    case __emk_s_discontinuous_lagrange_interval_01:
    case __emk_s_discontinuous_lagrange_interval_02:
    case __emk_s_discontinuous_lagrange_interval_03:
    case __emk_s_discontinuous_lagrange_interval_04:
    case __emk_s_discontinuous_lagrange_interval_05:
    case __emk_s_discontinuous_lagrange_interval_06:
    case __emk_s_discontinuous_lagrange_interval_07:
    case __emk_s_discontinuous_lagrange_interval_08:
    case __emk_s_discontinuous_l2orthonormal_interval_01:
    case __emk_s_discontinuous_l2orthonormal_interval_02:
    case __emk_s_discontinuous_l2orthonormal_interval_03:
    case __emk_s_discontinuous_l2orthonormal_interval_04:
    case __emk_s_discontinuous_l2orthonormal_interval_05:
    case __emk_s_discontinuous_l2orthonormal_interval_06:
    case __emk_s_discontinuous_l2orthonormal_interval_07:
    case __emk_s_discontinuous_l2orthonormal_interval_08:
      {
	element = __eTopology_SEGMENT;
	break;
      }

    case __emk_s_discontinuous_lagrange_triangle_01:
    case __emk_s_discontinuous_lagrange_triangle_03:
    case __emk_s_discontinuous_lagrange_triangle_06:
    case __emk_s_discontinuous_lagrange_triangle_10:
    case __emk_s_discontinuous_lagrange_triangle_15:
    case __emk_s_discontinuous_lagrange_triangle_21:
    case __emk_s_discontinuous_lagrange_triangle_28:
    case __emk_s_discontinuous_lagrange_triangle_35:
    case __emk_s_continuous_lagrange_triangle_03:
    case __emk_s_continuous_lagrangebubble_triangle_04:
    case __emk_s_continuous_lagrange_triangle_06:
    case __emk_s_continuous_lagrangebubble_triangle_07:
    case __emk_s_continuous_lagrange_triangle_10:
    case __emk_s_continuous_lagrange_triangle_15:
    case __emk_s_continuous_lagrange_triangle_21:
    case __emk_s_continuous_lagrange_triangle_28:
    case __emk_s_continuous_lagrange_triangle_35:
    case __emk_s_discontinuous_canonic_triangle_01:
    case __emk_s_discontinuous_canonic_triangle_03:
    case __emk_s_discontinuous_canonic_triangle_06:
    case __emk_s_discontinuous_canonic_triangle_10:
    case __emk_s_discontinuous_canonic_triangle_15:
    case __emk_s_discontinuous_canonic_triangle_21:
    case __emk_s_discontinuous_canonic_triangle_28:
    case __emk_s_discontinuous_canonic_triangle_35:
    case __emk_s_discontinuous_l2orthonormal_triangle_01:
    case __emk_s_discontinuous_l2orthonormal_triangle_03:
    case __emk_s_discontinuous_l2orthonormal_triangle_06:
    case __emk_s_discontinuous_l2orthonormal_triangle_10:
    case __emk_s_discontinuous_l2orthonormal_triangle_15:
    case __emk_s_discontinuous_l2orthonormal_triangle_21:
    case __emk_s_discontinuous_l2orthonormal_triangle_28:
    case __emk_s_discontinuous_l2orthonormal_triangle_35:
      {
	element = __eTopology_TRIANGLE;
	break;
      }
    case __emk_s_continuous_lagrange_quadrangle_04:
    case __emk_s_continuous_lagrange_quadrangle_09:
    case __emk_s_continuous_lagrange_quadrangle_16:
    case __emk_s_continuous_lagrange_quadrangle_25:
    case __emk_s_continuous_lagrange_quadrangle_36:
    case __emk_s_continuous_lagrange_quadrangle_49:
    case __emk_s_continuous_lagrange_quadrangle_64:
    case __emk_s_discontinuous_canonic_quadrangle_01:
    case __emk_s_discontinuous_canonic_quadrangle_04:
    case __emk_s_discontinuous_canonic_quadrangle_09:
    case __emk_s_discontinuous_canonic_quadrangle_16:
    case __emk_s_discontinuous_canonic_quadrangle_25:
    case __emk_s_discontinuous_canonic_quadrangle_36:
    case __emk_s_discontinuous_canonic_quadrangle_49:
    case __emk_s_discontinuous_canonic_quadrangle_64:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_01:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_04:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_09:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_16:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_25:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_36:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_49:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_64:
    case __emk_s_discontinuous_lagrange_quadrangle_01:
    case __emk_s_discontinuous_lagrange_quadrangle_04:
    case __emk_s_discontinuous_lagrange_quadrangle_09:
    case __emk_s_discontinuous_lagrange_quadrangle_16:
    case __emk_s_discontinuous_lagrange_quadrangle_25:
    case __emk_s_discontinuous_lagrange_quadrangle_36:
    case __emk_s_discontinuous_lagrange_quadrangle_49:
    case __emk_s_discontinuous_lagrange_quadrangle_64:
      {
	element = __eTopology_QUADRILATERAL;
	break;
      }
    case __emk_s_continuous_lagrange_tetrahedra_04:
    case __emk_s_continuous_lagrange_tetrahedra_10:
    case __emk_s_discontinuous_lagrange_tetrahedra_01:
    case __emk_s_discontinuous_lagrange_tetrahedra_04:
    case __emk_s_discontinuous_lagrange_tetrahedra_10:
    case __emk_s_discontinuous_canonic_tetrahedra_01:
    case __emk_s_discontinuous_canonic_tetrahedra_04:
    case __emk_s_discontinuous_canonic_tetrahedra_10:  
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_01:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_04:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_10:	       	       
      {
	element = __eTopology_TETRAHEDRON;
	break;
      }
    case __emk_s_error:
    case __emk_s_n :
      {
	fprintf(stderr,"emk_s_emkZELM:switch failed on __emk_s\n");
	exit(1);
	break;
      }
    }
  return element;
}


emkCONTINUITY 		__emk_s_emkCONTINUITY(cst_emk_s kind_)
{
  emkCONTINUITY continuity = __emk_continuous;
  switch(kind_)
    {
    case __emk_s_continuous_lagrange_interval_02:
    case __emk_s_continuous_lagrange_interval_03:
    case __emk_s_continuous_lagrange_interval_04:
    case __emk_s_continuous_lagrange_interval_05:
    case __emk_s_continuous_lagrange_interval_06:
    case __emk_s_continuous_lagrange_interval_07:
    case __emk_s_continuous_lagrange_interval_08:
    case __emk_s_continuous_lagrange_triangle_03:
    case __emk_s_continuous_lagrangebubble_triangle_04:
    case __emk_s_continuous_lagrange_triangle_06:
    case __emk_s_continuous_lagrangebubble_triangle_07:
    case __emk_s_continuous_lagrange_triangle_10:
    case __emk_s_continuous_lagrange_triangle_15:
    case __emk_s_continuous_lagrange_triangle_21:
    case __emk_s_continuous_lagrange_triangle_28:
    case __emk_s_continuous_lagrange_triangle_35:
    case __emk_s_continuous_lagrange_quadrangle_04:
    case __emk_s_continuous_lagrange_quadrangle_09:
    case __emk_s_continuous_lagrange_quadrangle_16:
    case __emk_s_continuous_lagrange_quadrangle_25:
    case __emk_s_continuous_lagrange_quadrangle_36:
    case __emk_s_continuous_lagrange_quadrangle_49:
    case __emk_s_continuous_lagrange_quadrangle_64:
    case __emk_s_continuous_lagrange_tetrahedra_04:
    case __emk_s_continuous_lagrange_tetrahedra_10:
      {
	continuity = __emk_continuous;
	break;
      }
    case __emk_s_discontinuous_lagrange_interval_01:
    case __emk_s_discontinuous_lagrange_interval_02:
    case __emk_s_discontinuous_lagrange_interval_03:
    case __emk_s_discontinuous_lagrange_interval_04:
    case __emk_s_discontinuous_lagrange_interval_05:
    case __emk_s_discontinuous_lagrange_interval_06:
    case __emk_s_discontinuous_lagrange_interval_07:
    case __emk_s_discontinuous_lagrange_interval_08:

    case __emk_s_discontinuous_lagrange_triangle_01:
    case __emk_s_discontinuous_lagrange_triangle_03:
    case __emk_s_discontinuous_lagrange_triangle_06:
    case __emk_s_discontinuous_lagrange_triangle_10:
    case __emk_s_discontinuous_lagrange_triangle_15:
    case __emk_s_discontinuous_lagrange_triangle_21:
    case __emk_s_discontinuous_lagrange_triangle_28:
    case __emk_s_discontinuous_lagrange_triangle_35:

    case __emk_s_discontinuous_lagrange_quadrangle_01:
    case __emk_s_discontinuous_lagrange_quadrangle_04:
    case __emk_s_discontinuous_lagrange_quadrangle_09:
    case __emk_s_discontinuous_lagrange_quadrangle_16:
    case __emk_s_discontinuous_lagrange_quadrangle_25:
    case __emk_s_discontinuous_lagrange_quadrangle_36:
    case __emk_s_discontinuous_lagrange_quadrangle_49:
    case __emk_s_discontinuous_lagrange_quadrangle_64:

    case __emk_s_discontinuous_lagrange_tetrahedra_01:
    case __emk_s_discontinuous_lagrange_tetrahedra_04:
    case __emk_s_discontinuous_lagrange_tetrahedra_10:

    case __emk_s_discontinuous_canonic_interval_01:
    case __emk_s_discontinuous_canonic_interval_02:
    case __emk_s_discontinuous_canonic_interval_03:
    case __emk_s_discontinuous_canonic_interval_04:
    case __emk_s_discontinuous_canonic_interval_05:
    case __emk_s_discontinuous_canonic_interval_06:
    case __emk_s_discontinuous_canonic_interval_07:
    case __emk_s_discontinuous_canonic_interval_08:

    case __emk_s_discontinuous_canonic_triangle_01:
    case __emk_s_discontinuous_canonic_triangle_03:
    case __emk_s_discontinuous_canonic_triangle_06:
    case __emk_s_discontinuous_canonic_triangle_10:
    case __emk_s_discontinuous_canonic_triangle_15:
    case __emk_s_discontinuous_canonic_triangle_21:
    case __emk_s_discontinuous_canonic_triangle_28:
    case __emk_s_discontinuous_canonic_triangle_35:

    case __emk_s_discontinuous_canonic_quadrangle_01:
    case __emk_s_discontinuous_canonic_quadrangle_04:
    case __emk_s_discontinuous_canonic_quadrangle_09:
    case __emk_s_discontinuous_canonic_quadrangle_16:
    case __emk_s_discontinuous_canonic_quadrangle_25:
    case __emk_s_discontinuous_canonic_quadrangle_36:
    case __emk_s_discontinuous_canonic_quadrangle_49:
    case __emk_s_discontinuous_canonic_quadrangle_64:

    case __emk_s_discontinuous_canonic_tetrahedra_01:
    case __emk_s_discontinuous_canonic_tetrahedra_04:
    case __emk_s_discontinuous_canonic_tetrahedra_10:	       	       

    case __emk_s_discontinuous_l2orthonormal_interval_01:
    case __emk_s_discontinuous_l2orthonormal_interval_02:
    case __emk_s_discontinuous_l2orthonormal_interval_03:
    case __emk_s_discontinuous_l2orthonormal_interval_04:
    case __emk_s_discontinuous_l2orthonormal_interval_05:
    case __emk_s_discontinuous_l2orthonormal_interval_06:
    case __emk_s_discontinuous_l2orthonormal_interval_07:
    case __emk_s_discontinuous_l2orthonormal_interval_08:

    case __emk_s_discontinuous_l2orthonormal_triangle_01:
    case __emk_s_discontinuous_l2orthonormal_triangle_03:
    case __emk_s_discontinuous_l2orthonormal_triangle_06:
    case __emk_s_discontinuous_l2orthonormal_triangle_10:
    case __emk_s_discontinuous_l2orthonormal_triangle_15:
    case __emk_s_discontinuous_l2orthonormal_triangle_21:
    case __emk_s_discontinuous_l2orthonormal_triangle_28:
    case __emk_s_discontinuous_l2orthonormal_triangle_35:

    case __emk_s_discontinuous_l2orthonormal_quadrangle_01:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_04:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_09:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_16:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_25:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_36:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_49:
    case __emk_s_discontinuous_l2orthonormal_quadrangle_64:

    case __emk_s_discontinuous_l2orthonormal_tetrahedra_01:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_04:
    case __emk_s_discontinuous_l2orthonormal_tetrahedra_10:	       	       
      {
	continuity = __emk_discontinuous;
	break;
      }

    case __emk_s_error:
    case __emk_s_n :
      {
	fprintf(stderr,"emk_s_emkCONTINUITY:switch failed on __emk_s\n");
	exit(1);
	break;
      }
    }
  return continuity;
}


I __emkS_FAMILY_nshape(const emkS_FAMILY 	family_,
		       const eTopology 	elm_,
		       const I 	degree_);

#if 0
#include "mkSYS.h"
#include "mkS_def_interval.h"
#include "mkS_def_tria.h"
#include "mkS_def_quad.h"
#include "mkS_def_tetra.h"
#endif

void mkS_def_transformation	(mkS 			shape_,
				 eTopology		elm_,
				 Err * 		err_)
{
  err_[0]=__eErr_no;
  mkS_def(shape_,__emk_s_continuous_lagrange_triangle_03,err_);
#if 0
	switch(elm_)
    {
    case ____eTopology_TRIANGLE:
      {
	mkS_def(shape_,__emk_s_continuous_lagrange_triangle_03,err_);
	break;
      }
    case ____eTopology_SEGMENT:
      {
	mkS_def(shape_,__emk_s_continuous_lagrange_interval_02,err_);
	break;
      }
    case ____eTopology_TETRAHEDRA:
      {
	mkS_def(shape_,__emk_s_continuous_lagrange_tetrahedra_04,err_);
	break;
      }
    case ____eTopology_QUADRILATERAL:
      {
	mkS_def(shape_,__emk_s_continuous_lagrange_quadrangle_04,err_);
	break;
      }
    case ____eTopology_ALL:
    case ____eTopology_ERROR:
      {
	break;
      }
    }
#endif
	
}

mkS 			mkS_new_transformation		(eTopology	elm_,
							 Err * 		err_)
{
  err_[0] = __eErr_no;
  mkS S = malloc(sizeof(mkS_st));
  if (NOT S)
    {
      err_[0] = __eErr_memory;
      fprintf(stderr,"mkS_new_transformation:malloc failed\n");
      return NULL;
    }
  mkS_def_transformation(S,elm_,err_);
  if (err_[0])
    {
      fprintf(stderr,"mkS_new_transformation:mkS_def_transformation failed\n");
      free(S);
      return NULL;
    }
  return S;
}

mkS mkS_newinit(eTopology    		element_,
		emkS_FAMILY    	family_,
		const I 		degree_,
		emkCONTINUITY   	continuity_,
		Err * 		err_)
{
  err_[0] = __eErr_no;
  mkS S = malloc(sizeof(mkS_st));
  if (NOT S)
    {
      err_[0] = __eErr_memory;
      fprintf(stderr,"mkS_newinit:malloc failed\n");
      return NULL;
    }
  mkS_definit(S,element_,family_,degree_,continuity_,err_);
  if (err_[0])
    {
      fprintf(stderr,"mkS_newinit:mkS_definit failed\n");
      free(S);
      return NULL;
    }
  return S;
}

#define __mmkZELM_triangle_name  	"triangle"
#define __mmkZELM_quadrangle_name  	"quadrangle"
#define __mmkZELM_interval_name  	"interval"
#define __mmkZELM_tetrahedra_name  	"tetrahedra"

const char * __emk_s_strings[__emk_s_n] =  { "__emk_s_error",  
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_02",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_03",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_05",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_06",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_07",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_08",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_03",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrangebubble_name"_04",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_06",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrangebubble_name"_07",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_10",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_15",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_21",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_28",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_35",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_09",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_16",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_25",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_36",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_49",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_64",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_continuous_name"_"__mmk_s_family_lagrange_name"_10",
	       
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_01",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_02",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_03",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_05",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_06",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_07",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_08",

					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_01",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_03",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_06",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_10",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_15",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_21",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_28",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_35",

					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_01",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_09",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_16",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_25",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_36",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_49",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_64",
	       
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_01",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_04",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_lagrange_name"_10",
	       
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_01",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_02",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_03",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_04",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_05",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_06",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_07",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_08",
	       
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_01",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_03",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_06",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_10",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_15",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_21",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_28",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_35",
	       
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_01",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_04",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_09",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_16",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_25",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_36",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_49",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_64",
	       
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_01",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_04",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_canonic_name"_10",	       	       
	       
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_01",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_02",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_03",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_04",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_05",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_06",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_07",
					     ""__mmkZELM_interval_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_08",
	       
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_01",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_03",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_06",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_10",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_15",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_21",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_28",
					     ""__mmkZELM_triangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_35",
	       
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_01",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_04",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_09",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_16",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_25",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_36",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_49",
					     ""__mmkZELM_quadrangle_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_64",
	       
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_01",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_04",
					     ""__mmkZELM_tetrahedra_name"_"__mmk_continuity_discontinuous_name"_"__mmk_s_family_l2orthonormal_name"_10"};	       	       

cst_emk_s __emk_s_enum(const char * ctmp)
{
  int i=0;
  for (i=1;i<__emk_s_n;++i)
    {
      if (!strcmp(ctmp,__emk_s_strings[i]))
	{
	  return i;
	}
    }
  return __emk_s_error;
}


cst_emk_s mkS_emk_s_enum(const emkCONTINUITY 	continuity_,
			 const emkS_FAMILY 	family_,
			 const eTopology	element_,
			 const I 		nshape_)
{

  
  //  printf("warning: come here and correct the format %%Ld\n");
  char ctmp[256];
  sprintf(ctmp,"%s_%s_%s_%.2Ld",
	  emkZELM_string(element_),
	  __emkCONTINUITY_string(continuity_),
	  __emkS_FAMILY_string(family_),
	  nshape_);
  return __emk_s_enum(ctmp);
}


void mkS_definit(mkS 			shape_,
		 eTopology    		element_,
		 emkS_FAMILY    	family_,
		 const I 		degree_,
		 emkCONTINUITY   	continuity_,
		 Err * 		err_)
{
  err_[0] = __eErr_no;
  memset(shape_,0,sizeof(mkS_st));
  shape_->element 	= element_;
  shape_->family	= family_;
  //  shape_->continu 	= continuity_; 
  shape_->degree        = degree_;
  shape_->nshape    	= __emkS_FAMILY_nshape(family_,shape_->element,shape_->degree);
  shape_->shape 	= mkS_emk_s_enum(shape_->continu,
					 shape_->family,
					 shape_->element,
					 shape_->nshape);      
  if (shape_->shape==__emk_s_error)
    {
      err_[0] = __eErr_intern;
      fprintf(stderr,"mkS_definit:mkS_emk_s_enum failed\n");
      return;
    }
}



I __emkS_FAMILY_nshape(const emkS_FAMILY 	family_,
		       const eTopology 	elm_,
		       const I 	degree_)
{

  I n=0;
  switch(elm_)
    {
    case __eTopology_SEGMENT:
      {
	n=degree_+1;
	break;
      }
    case __eTopology_QUADRILATERAL:
      {
	switch(family_)
	  {
	  case __emkS_FAMILY_lagrangebubble:
	  case __emkS_FAMILY_canonic:
	  case __emkS_FAMILY_l2orthonormal:
	  case __emkS_FAMILY_lagrange:
	    {
	      n= (degree_+1)*(degree_+1);
	      break;
	    }
	  case __emkS_FAMILY_n:
	  case __emkS_FAMILY_error:
	    {
	      fprintf(stderr,"mkS_nshapelm:switch failed on __emkZELM\n");
	      break;
	    }
	  }	
	break;
      }
    case __eTopology_TRIANGLE:
      {
	switch(family_)
	  {
	  case __emkS_FAMILY_lagrangebubble:
	    {
	      if (degree_==1)
		n=4;
	      else if (degree_==2)
		n=7;
	      else
		{
		  fprintf(stderr,"__emkS_FAMILY_nshape:lagrangebubble degree="ifmt"\n",degree_);
		  n=0;
		  break;
		}
	      break;
	    }
	  case __emkS_FAMILY_canonic:
	  case __emkS_FAMILY_l2orthonormal:
	  case __emkS_FAMILY_lagrange:
	    {
	      n= ((degree_+1)*(degree_+2))/2;
	      break;
	    }
	  case __emkS_FAMILY_n:
	  case __emkS_FAMILY_error:
	    {
	      fprintf(stderr,"mkS_nshapelm:switch failed on __emkZELM\n");
	      exit(1);
	      break;
	    }
	  }	
	break;
      }
    case __eTopology_TETRAHEDRON:
      {

	switch(family_)
	  {
	  case __emkS_FAMILY_lagrangebubble:
	    {
	      if (degree_==1)
		n=5;
	      else if (degree_==2)
		n=11;
	      else
		{
		  fprintf(stderr,"__emkS_FAMILY_nshape:lagrangebubble degree="ifmt"\n",degree_);
		  n=0;
		  break;
		}
	      break;
	    }
	  case __emkS_FAMILY_canonic:
	  case __emkS_FAMILY_l2orthonormal:
	  case __emkS_FAMILY_lagrange:
	    {
	      if (degree_==1)
		n=4;
	      else if (degree_==2)
		n=10;
	      else
		{
		  fprintf(stderr,"__emkS_FAMILY_nshape:tetrahedra degree="ifmt"\n",degree_);
		  n=0;
		  break;
		}
	      break;
	    }
	  case __emkS_FAMILY_n:
	  case __emkS_FAMILY_error:
	    {
	      fprintf(stderr,"mkS_nshapelm:switch failed on __emkZELM\n");
	      exit(1);
	      break;
	    }
	  }	
	break;

      }
    case __eTopology_ERROR:
    case __eTopology_ALL:
    default:
      {
	fprintf(stderr,"mkS_nshapelm:switch failed on __emkZELM\n");
	exit(1);
	break;
      }
    }
  return n;
}

void mkS_def	(mkS 			shape_,
		 cst_emk_s 		kind_,
		 Err * 		err_)
{
#if __mk_debug__
  __mkS_debug(shape_,mkS_def);
  __mkS_debug(err_,mkS_def);
#endif
  err_[0] = __eErr_no;
  memset(shape_,0,sizeof(mkS_st));
  shape_->shape		= kind_;
  shape_->element 	= __emk_s_emkZELM(kind_);
  shape_->family	= __emk_s_emkS_FAMILY(kind_);
  shape_->continu 	= __emk_s_emkCONTINUITY(kind_);  
  shape_->degree        = __emk_s_degree(kind_);
  shape_->nshape    	= __emkS_FAMILY_nshape(shape_->family,shape_->element,shape_->degree);
}

