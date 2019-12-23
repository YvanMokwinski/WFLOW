
#include "mkS.h"
#if 0

#define __mmkS_lagrange_localspl_interval(_degree,_coo,_off)		\
  {									\
    switch((_degree)[0])						\
      {									\
      case 0:								\
	{								\
	  (_coo)[0] = (__mmkZELM_interval_refcoo_x0+__mmkZELM_interval_refcoo_x1)*( (nsREAL)0.5 ); \
	  break;							\
	}								\
      case 1:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  break;							\
	}								\
      case 2:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] =  __mmkZELM_interval_refcoo_x0*((nsREAL)0.5)+__mmkZELM_interval_refcoo_x1*((nsREAL)0.5); \
	  break;							\
	}								\
      case 3:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((nsREAL)2.0)/((nsREAL)3.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)1.0)/((nsREAL)3.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((nsREAL)1.0)/((nsREAL)3.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)2.0)/((nsREAL)3.0)); \
	  break;							\
	}								\
      case 4:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((nsREAL)3.0)/((nsREAL)4.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)1.0)/((nsREAL)4.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*((nsREAL)0.5)+__mmkZELM_interval_refcoo_x1*((nsREAL)0.5); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*(((nsREAL)1.0)/((nsREAL)4.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)3.0)/((nsREAL)4.0)); \
	  break;							\
	}								\
      case 5:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((nsREAL)4.0)/((nsREAL)5.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)1.0)/((nsREAL)5.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((nsREAL)3.0)/((nsREAL)5.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)2.0)/((nsREAL)5.0)); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*(((nsREAL)2.0)/((nsREAL)5.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)3.0)/((nsREAL)5.0)); \
	  (_coo)[5] = __mmkZELM_interval_refcoo_x0*(((nsREAL)1.0)/((nsREAL)5.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)4.0)/((nsREAL)5.0)); \
	  break;							\
	}								\
      case 6:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((nsREAL)5.0)/((nsREAL)6.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)1.0)/((nsREAL)6.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((nsREAL)4.0)/((nsREAL)6.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)2.0)/((nsREAL)6.0)); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*((nsREAL)0.5)+__mmkZELM_interval_refcoo_x1*((nsREAL)0.5); \
	  (_coo)[5] = __mmkZELM_interval_refcoo_x0*(((nsREAL)2.0)/((nsREAL)6.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)4.0)/((nsREAL)6.0)); \
	  (_coo)[6] = __mmkZELM_interval_refcoo_x0*(((nsREAL)1.0)/((nsREAL)6.0))+__mmkZELM_interval_refcoo_x1*(((nsREAL)5.0)/((nsREAL)6.0)); \
	  break;							\
	}								\
      default:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  nsINT _i;							\
	  const nsREAL _dt = (__mmkZELM_interval_refcoo_x1-__mmkZELM_interval_refcoo_x0)/((nsREAL)((_degree)[0])); \
	  for (_i=1;_i<(_degree)[0];++_i)				\
	    (_coo)[_i+1] = __mmkZELM_interval_refcoo_x0 + _dt*((nsREAL)_i); \
	  break;							\
	}								\
      }									\
  }



#define __mmkS_lagrange_localspl_tetra(_degree,_coo,_off)	\
  {								\
    switch((_degree)[0])					\
      {								\
      case ((nsINT)0):						\
	{							\
	  (_coo)[0] = (__mmkZELM_tetrahedra_refcoo_x0+__mmkZELM_tetrahedra_refcoo_x1+__mmkZELM_tetrahedra_refcoo_x2+__mmkZELM_tetrahedra_refcoo_x3)*((nsREAL)0.25); \
	  (_coo)[1] = (__mmkZELM_tetrahedra_refcoo_y0+__mmkZELM_tetrahedra_refcoo_y1+__mmkZELM_tetrahedra_refcoo_y2+__mmkZELM_tetrahedra_refcoo_y3)*((nsREAL)0.25); \
	  (_coo)[2] = (__mmkZELM_tetrahedra_refcoo_z0+__mmkZELM_tetrahedra_refcoo_z1+__mmkZELM_tetrahedra_refcoo_z2+__mmkZELM_tetrahedra_refcoo_z3)*((nsREAL)0.25); \
	  break;							\
	}								\
      case ((nsINT)1):							\
	{								\
	  (_coo)[0] 		= __mmkZELM_tetrahedra_refcoo_x0;	\
	  (_coo)[(*(_off))+0] 	= __mmkZELM_tetrahedra_refcoo_y0;	\
	  (_coo)[2*(*(_off))+0] = __mmkZELM_tetrahedra_refcoo_z0;	\
	  								\
	  (_coo)[1] 		= __mmkZELM_tetrahedra_refcoo_x1;	\
	  (_coo)[(*(_off))+1] 	= __mmkZELM_tetrahedra_refcoo_y1;	\
	  (_coo)[2*(*(_off))+1] = __mmkZELM_tetrahedra_refcoo_z1;	\
	  								\
	  (_coo)[2] 		= __mmkZELM_tetrahedra_refcoo_x2;	\
	  (_coo)[(*(_off))+2] 	= __mmkZELM_tetrahedra_refcoo_y2;	\
	  (_coo)[2*(*(_off))+2] = __mmkZELM_tetrahedra_refcoo_z2;	\
	  								\
	  (_coo)[3] 		= __mmkZELM_tetrahedra_refcoo_x3;	\
	  (_coo)[(*(_off))+3] 	= __mmkZELM_tetrahedra_refcoo_y3;	\
	  (_coo)[2*(*(_off))+3] = __mmkZELM_tetrahedra_refcoo_z3;	\
	  break;							\
	}								\
									\
      case ((nsINT)2):							\
	{								\
	  (_coo)[0] 		= __mmkZELM_tetrahedra_refcoo_x0;	\
	  (_coo)[(*(_off))+0] 	= __mmkZELM_tetrahedra_refcoo_y0;	\
	  (_coo)[2*(*(_off))+0] = __mmkZELM_tetrahedra_refcoo_z0;	\
	  								\
	  (_coo)[1] 		= __mmkZELM_tetrahedra_refcoo_x1;	\
	  (_coo)[(*(_off))+1] 	= __mmkZELM_tetrahedra_refcoo_y1;	\
	  (_coo)[2*(*(_off))+1] = __mmkZELM_tetrahedra_refcoo_z1;	\
	  								\
	  (_coo)[2] 		= __mmkZELM_tetrahedra_refcoo_x2;	\
	  (_coo)[(*(_off))+2] 	= __mmkZELM_tetrahedra_refcoo_y2;	\
	  (_coo)[2*(*(_off))+2] = __mmkZELM_tetrahedra_refcoo_z2;	\
	  								\
	  (_coo)[3] 		= __mmkZELM_tetrahedra_refcoo_x3;	\
	  (_coo)[(*(_off))+3] 	= __mmkZELM_tetrahedra_refcoo_y3;	\
	  (_coo)[2*(*(_off))+3] = __mmkZELM_tetrahedra_refcoo_z3;	\
	  (_coo)[4] 		= (__mmkZELM_tetrahedra_refcoo_x0+__mmkZELM_tetrahedra_refcoo_x1)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+4] 	= (__mmkZELM_tetrahedra_refcoo_y0+__mmkZELM_tetrahedra_refcoo_y1)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+4] = (__mmkZELM_tetrahedra_refcoo_z0+__mmkZELM_tetrahedra_refcoo_z1)*((nsREAL)0.5);			\
	  (_coo)[5] 		= (__mmkZELM_tetrahedra_refcoo_x1+__mmkZELM_tetrahedra_refcoo_x2)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+5] 	= (__mmkZELM_tetrahedra_refcoo_y1+__mmkZELM_tetrahedra_refcoo_y2)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+5] = (__mmkZELM_tetrahedra_refcoo_z1+__mmkZELM_tetrahedra_refcoo_z2)*((nsREAL)0.5);			\
	  (_coo)[6] 		= (__mmkZELM_tetrahedra_refcoo_x2+__mmkZELM_tetrahedra_refcoo_x0)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+6] 	= (__mmkZELM_tetrahedra_refcoo_y2+__mmkZELM_tetrahedra_refcoo_y0)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+6] = (__mmkZELM_tetrahedra_refcoo_z2+__mmkZELM_tetrahedra_refcoo_z0)*((nsREAL)0.5);			\
	  (_coo)[7] 		= (__mmkZELM_tetrahedra_refcoo_x3+__mmkZELM_tetrahedra_refcoo_x0)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+7] 	= (__mmkZELM_tetrahedra_refcoo_y3+__mmkZELM_tetrahedra_refcoo_y0)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+7] = (__mmkZELM_tetrahedra_refcoo_z3+__mmkZELM_tetrahedra_refcoo_z0)*((nsREAL)0.5);			\
	  (_coo)[8] 		= (__mmkZELM_tetrahedra_refcoo_x3+__mmkZELM_tetrahedra_refcoo_x1)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+8] 	= (__mmkZELM_tetrahedra_refcoo_y3+__mmkZELM_tetrahedra_refcoo_y1)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+8] = (__mmkZELM_tetrahedra_refcoo_z3+__mmkZELM_tetrahedra_refcoo_z1)*((nsREAL)0.5);			\
	  (_coo)[9] 		= (__mmkZELM_tetrahedra_refcoo_x3+__mmkZELM_tetrahedra_refcoo_x2)*((nsREAL)0.5);			\
	  (_coo)[(*(_off))+9] 	= (__mmkZELM_tetrahedra_refcoo_y3+__mmkZELM_tetrahedra_refcoo_y2)*((nsREAL)0.5);			\
	  (_coo)[2*(*(_off))+9] = (__mmkZELM_tetrahedra_refcoo_z3+__mmkZELM_tetrahedra_refcoo_z2)*((nsREAL)0.5);			\
	  break;						\
	}							\
      default:							\
	{							\
	  mk_err("__mmkS_lagrange_localspl_tetra degree>2");	\
	  break;						\
	}							\
      }								\
  }



#define __mmkS_lagrange_localspl_quad(_degree,_coo,_off)		\
  {									\
    switch((_degree)[0])						\
      {									\
      case 0:								\
	{								\
	  (_coo)[0] = (((__mmkZELM_quadrangle_refcoo_x0+__mmkZELM_quadrangle_refcoo_x1)*((nsREAL)0.5))+((__mmkZELM_quadrangle_refcoo_x2+__mmkZELM_quadrangle_refcoo_x3)*((nsREAL)0.5)))*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+0]	= ( ((__mmkZELM_quadrangle_refcoo_y3+__mmkZELM_quadrangle_refcoo_y0)*((nsREAL)0.5)) + ((__mmkZELM_quadrangle_refcoo_y1+__mmkZELM_quadrangle_refcoo_y2)*((nsREAL)0.5)) )*( (nsREAL)0.5 ); \
	  break;							\
	}								\
      case 1:								\
	{								\
	  (_coo)[0] =__mmkZELM_quadrangle_refcoo_x0;			\
	  (_coo)[1] =__mmkZELM_quadrangle_refcoo_x1;			\
	  (_coo)[2] =__mmkZELM_quadrangle_refcoo_x2;			\
	  (_coo)[3] =__mmkZELM_quadrangle_refcoo_x3;			\
	  (_coo)[(_off)[0]+0] =__mmkZELM_quadrangle_refcoo_y0;		\
	  (_coo)[(_off)[0]+1] =__mmkZELM_quadrangle_refcoo_y1;		\
	  (_coo)[(_off)[0]+2] =__mmkZELM_quadrangle_refcoo_y2;		\
	  (_coo)[(_off)[0]+3] =__mmkZELM_quadrangle_refcoo_y3;		\
	  break;							\
	}								\
      case 2:								\
	{								\
	  (_coo)[0] 		= __mmkZELM_quadrangle_refcoo_x0;	\
	  (_coo)[1] 		= __mmkZELM_quadrangle_refcoo_x1;	\
	  (_coo)[2] 		= __mmkZELM_quadrangle_refcoo_x2;	\
	  (_coo)[3] 		= __mmkZELM_quadrangle_refcoo_x3;	\
	  (_coo)[(_off)[0]+0] 	= __mmkZELM_quadrangle_refcoo_y0;	\
	  (_coo)[(_off)[0]+1] 	= __mmkZELM_quadrangle_refcoo_y1;	\
	  (_coo)[(_off)[0]+2] 	= __mmkZELM_quadrangle_refcoo_y2;	\
	  (_coo)[(_off)[0]+3] 	= __mmkZELM_quadrangle_refcoo_y3;	\
	  (_coo)[4] 		= (__mmkZELM_quadrangle_refcoo_x0+__mmkZELM_quadrangle_refcoo_x1)*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+4]	= (__mmkZELM_quadrangle_refcoo_y0+__mmkZELM_quadrangle_refcoo_y1)*((nsREAL)0.5); \
	  (_coo)[5] 		= (__mmkZELM_quadrangle_refcoo_x1+__mmkZELM_quadrangle_refcoo_x2)*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+5]	= (__mmkZELM_quadrangle_refcoo_y1+__mmkZELM_quadrangle_refcoo_y2)*((nsREAL)0.5); \
	  (_coo)[6] 		= (__mmkZELM_quadrangle_refcoo_x2+__mmkZELM_quadrangle_refcoo_x3)*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+6]	= (__mmkZELM_quadrangle_refcoo_y2+__mmkZELM_quadrangle_refcoo_y3)*((nsREAL)0.5); \
	  (_coo)[7] 		= (__mmkZELM_quadrangle_refcoo_x3+__mmkZELM_quadrangle_refcoo_x0)*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+7]	= (__mmkZELM_quadrangle_refcoo_y3+__mmkZELM_quadrangle_refcoo_y0)*((nsREAL)0.5); \
	  (_coo)[(_off)[0]+8] = (((__mmkZELM_quadrangle_refcoo_x0+__mmkZELM_quadrangle_refcoo_x1)*((nsREAL)0.5))+((__mmkZELM_quadrangle_refcoo_x2+__mmkZELM_quadrangle_refcoo_x3)*((nsREAL)0.5)))*((nsREAL)0.5); \
	  (_coo)[8]	= ( ((__mmkZELM_quadrangle_refcoo_y3+__mmkZELM_quadrangle_refcoo_y0)*((nsREAL)0.5)) + ((__mmkZELM_quadrangle_refcoo_y1+__mmkZELM_quadrangle_refcoo_y2)*((nsREAL)0.5)) )*( (nsREAL)0.5 ); \
	  break;							\
	}								\
      default:								\
	{								\
	  const nsINT _n   	= (_degree)[0]-1;			\
	  const nsREAL _dt 	= (__mmkZELM_quadrangle_refcoo_x0-__mmkZELM_quadrangle_refcoo_x1)/((nsREAL)((_degree)[0])); \
	  (_coo)[0] 		= __mmkZELM_quadrangle_refcoo_x0;	\
	  (_coo)[1] 		= __mmkZELM_quadrangle_refcoo_x1;	\
	  (_coo)[2] 		= __mmkZELM_quadrangle_refcoo_x2;	\
	  (_coo)[3] 		= __mmkZELM_quadrangle_refcoo_x3;	\
	  (_coo)[(_off)[0]+0] 	= __mmkZELM_quadrangle_refcoo_y0;	\
	  (_coo)[(_off)[0]+1] 	= __mmkZELM_quadrangle_refcoo_y1;	\
	  (_coo)[(_off)[0]+2] 	= __mmkZELM_quadrangle_refcoo_y2;	\
	  (_coo)[(_off)[0]+3] 	= __mmkZELM_quadrangle_refcoo_y3;	\
	  { nsINT _i; for (_i=0;_i<_n;++_i){ (_coo)[4+_i] = ((nsREAL)1.0) - _dt*((nsREAL)_i+1); (_coo)[(_off)[0]+4+_i] = (nsREAL)1.0; }} \
	  { nsINT _i; for (_i=0;_i<_n;++_i){ (_coo)[(_off)[0]+4+_n+_i] = ((nsREAL)1.0) - _dt*((nsREAL)_i+1); (_coo)[4+_i+_n]	= (nsREAL)-1.0; }} \
	  { nsINT _i; for (_i=0;_i<_n;++_i){ (_coo)[4+_n*2+_i] = ((nsREAL)-1.0) + _dt*((nsREAL)_i+1); (_coo)[(_off)[0]+4+_i+_n*2]	= (nsREAL)-1.0; }} \
	  { nsINT _i; for (_i=0;_i<_n;++_i) {(_coo)[(_off)[0]+4+_n*3+_i] = ((nsREAL)-1.0) + _dt*((nsREAL)_i+1); (_coo)[4+_i+_n*3]	= (nsREAL)1.0; }} \
	  { nsINT _i,_j; for (_i=0;_i<_n;++_i) for (_j=0;_j<_n;++_j) {  (_coo)[(_off)[0]+4+_n*4+_i*_n+_j] = ((nsREAL)-1.0) + _dt*((nsREAL)_j+1); (_coo)[4+_i*_n+_j+_n*4] = ((nsREAL)-1.0) + _dt*((nsREAL)_i+1); } } \
	  								\
	  break;							\
	}								\
      }									\
  }



#if 0
#define __mmkS_lagrange_localspl_tria(_degree,_coo,_off)		\
  {									\
    const nsREAL __mmkS_lagrange_localspl_triax = (nsREAL)1.0/((nsREAL)(_degree)[0]); \
    switch((_degree)[0])						\
      {									\
      case 0:								\
	{								\
	  (_coo)[0] = (nsREAL)1.0/(nsREAL)3.0;				\
	  (_coo)[(_off)[0]] = (nsREAL)1.0/(nsREAL)3.0;			\
	  break;							\
	}								\
      case 1:								\
	{								\
	  (_coo)[(_off)[0]+0] = (nsREAL)0.0;				\
	  (_coo)[1] = (nsREAL)1.0;					\
	  (_coo)[(_off)[0]+1] = (nsREAL)0.0;				\
	  (_coo)[2] = (nsREAL)1.0;					\
	  (_coo)[(_off)[0]+2] = (nsREAL)1.0;				\
	  break;							\
	}								\
      case 2:								\
	{								\
	  (_coo)[0] = (nsREAL)1.0/(nsREAL)3.0;				\
	  (_coo)[(_off)[0]] = (nsREAL)1.0/(nsREAL)3.0;			\
	  break;							\
	}								\
      default:								\
	{								\
	  nsINT __mmkS_lagrange_localspl_triak=(nsINT)0;		\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= (nsREAL)0.0;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= (nsREAL)0.0; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= (nsREAL)1.0;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= (nsREAL)0.0; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= (nsREAL)0.0;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= (nsREAL)1.0; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  { nsINT __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak] = ((nsREAL)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((nsREAL)0.0); \
	      }}							\
	  								\
	  { nsINT __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak] = ((nsREAL)((_degree)[0]-__mmkS_lagrange_localspl_triaj))*__mmkS_lagrange_localspl_triax; \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((nsREAL)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
	      }}							\
	  								\
	  { nsINT __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak]      	= ((nsREAL)0.0); \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((nsREAL)((_degree)[0]-__mmkS_lagrange_localspl_triaj))*__mmkS_lagrange_localspl_triax; \
	      } }							\
									\
	  { nsINT __mmkS_lagrange_localspl_triai;			\
	    for (__mmkS_lagrange_localspl_triai=1;			\
		 __mmkS_lagrange_localspl_triai<(_degree)[0]-1;		\
		 ++__mmkS_lagrange_localspl_triai)			\
	      { nsINT __mmkS_lagrange_localspl_triaj;			\
		for (__mmkS_lagrange_localspl_triaj=1;			\
		     __mmkS_lagrange_localspl_triaj<(_degree)[0]-__mmkS_lagrange_localspl_triai; \
		     ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
		  {							\
		    (_coo)[__mmkS_lagrange_localspl_triak]      = ((nsREAL)__mmkS_lagrange_localspl_triai)*__mmkS_lagrange_localspl_triax; \
		    (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] = ((nsREAL)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
		  } }							\
	  } }								\
	break;								\
      }									\
  }									\

#endif





#define __mmkS_lagrange_localspl( _elm,_degree,_coo,_off)	\
  {								\
    switch(_elm)						\
      {								\
      case __emkZELM_triangle:					\
	{							\
	  __mmkS_lagrange_localspl_tria( _degree,		\
					 _coo,			\
					 _off);			\
	  break;						\
	}							\
      case __emkZELM_quadrangle:				\
	{							\
	  __mmkS_lagrange_localspl_quad( _degree,		\
					 _coo,			\
					 _off);			\
	  break;						\
	}							\
      case __emkZELM_tetrahedra:				\
	{							\
	  __mmkS_lagrange_localspl_tetra( _degree,		\
					  _coo,			\
					  _off);		\
	  break;						\
	}							\
      case __emkZELM_interval:					\
	{							\
	  __mmkS_lagrange_localspl_interval( _degree,		\
					     _coo,		\
					     _off);		\
	  break;						\
	}							\
      case __emkZELM_n:						\
      case __emkZELM_error:					\
	{							\
	  break;						\
	}							\
      }								\
  }





#define __mmkS_lagrangebubble_localspl( _elm,_degree,_coo,_off)	\
  {								\
    switch(_elm)						\
      {								\
      case __emkZELM_triangle:					\
	{							\
	  __mmkS_lagrangebubble_localspl_tria( _degree,		\
					       _coo,		\
					       _off);		\
	  break;						\
	}							\
      case __emkZELM_quadrangle:				\
	{								\
	  mk_msg("__mmkS_lagrangebubble_localspl:__emkZELM_quadrangle not available"); \
	  break;							\
	}								\
      case __emkZELM_tetrahedra:					\
	{								\
	  mk_msg("__mmkS_lagrangebubble_localspl:__emkZELM_quadrangle not available"); \
	  break;							\
	}								\
      case __emkZELM_interval:						\
	{								\
	  mk_msg("__mmkS_lagrangebubble_localspl:__emkZELM_quadrangle not available"); \
	  break;						\
	}							\
      case __emkZELM_n:						\
      case __emkZELM_error:					\
	{							\
	  break;						\
	}							\
      }								\
  }

#endif




void __mmkS_lagrange_localspl_tria( cst_pI 	_degree,
				 pR 	_coo,
				 cst_pI	_off)
  {

#define __mmkZELM_triangle_refcoo_x0 0.0
#define __mmkZELM_triangle_refcoo_x1 1.0
#define __mmkZELM_triangle_refcoo_x2 0.0

#define __mmkZELM_triangle_refcoo_y0 0.0
#define __mmkZELM_triangle_refcoo_y1 0.0
#define __mmkZELM_triangle_refcoo_y2 1.0

    
    const R __mmkS_lagrange_localspl_triax = (R)1.0/((R)(_degree)[0]);	\
    switch((_degree)[0])						\
      {									\
      case 0:								\
	{								\
	  (_coo)[0] = (__mmkZELM_triangle_refcoo_x0+__mmkZELM_triangle_refcoo_x1+__mmkZELM_triangle_refcoo_x2)/((R)3.0); \
	  (_coo)[(_off)[0]] = ((__mmkZELM_triangle_refcoo_x0+__mmkZELM_triangle_refcoo_x1+__mmkZELM_triangle_refcoo_x2))/((R)3.0); \
	  break;							\
	}								\
      case 1:								\
	{								\
	  (_coo)[0] =__mmkZELM_triangle_refcoo_x0;			\
	  (_coo)[1] =__mmkZELM_triangle_refcoo_x1;			\
	  (_coo)[2] =__mmkZELM_triangle_refcoo_x2;			\
	  (_coo)[(_off)[0]+0] =__mmkZELM_triangle_refcoo_y0;		\
	  (_coo)[(_off)[0]+1] =__mmkZELM_triangle_refcoo_y1;		\
	  (_coo)[(_off)[0]+2] =__mmkZELM_triangle_refcoo_y2;		\
	  break;							\
	}								\
      case 2:								\
	{								\
	  (_coo)[0] =__mmkZELM_triangle_refcoo_x0;			\
	  (_coo)[1] =__mmkZELM_triangle_refcoo_x1;			\
	  (_coo)[2] =__mmkZELM_triangle_refcoo_x2;			\
	  (_coo)[3] =(__mmkZELM_triangle_refcoo_x0+__mmkZELM_triangle_refcoo_x1)*((R)0.5); \
	  (_coo)[4] =(__mmkZELM_triangle_refcoo_x1+__mmkZELM_triangle_refcoo_x2)*((R)0.5);	\
	  (_coo)[5] =(__mmkZELM_triangle_refcoo_x2+__mmkZELM_triangle_refcoo_x0)*((R)0.5);	\
	  (_coo)[(_off)[0]+0] =__mmkZELM_triangle_refcoo_y0;		\
	  (_coo)[(_off)[0]+1] =__mmkZELM_triangle_refcoo_y1;		\
	  (_coo)[(_off)[0]+2] =__mmkZELM_triangle_refcoo_y2;		\
	  (_coo)[(_off)[0]+3] =(__mmkZELM_triangle_refcoo_y0+__mmkZELM_triangle_refcoo_y1)*((R)0.5); \
	  (_coo)[(_off)[0]+4] =(__mmkZELM_triangle_refcoo_y1+__mmkZELM_triangle_refcoo_y2)*((R)0.5); \
	  (_coo)[(_off)[0]+5] =(__mmkZELM_triangle_refcoo_y2+__mmkZELM_triangle_refcoo_y0)*((R)0.5); \
	  break;							\
	}								\
      case 3:								\
	{								\
	  (_coo)[0] =__mmkZELM_triangle_refcoo_x0;			\
	  (_coo)[1] =__mmkZELM_triangle_refcoo_x1;			\
	  (_coo)[2] =__mmkZELM_triangle_refcoo_x2;			\
	  (_coo)[3] = __mmkZELM_triangle_refcoo_x0*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x1*(((R)1.0)/((R)3.0)); \
	  (_coo)[4] = __mmkZELM_triangle_refcoo_x0*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x1*(((R)2.0)/((R)3.0)); \
	  (_coo)[5] = __mmkZELM_triangle_refcoo_x1*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x2*(((R)1.0)/((R)3.0)); \
	  (_coo)[6] = __mmkZELM_triangle_refcoo_x1*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x2*(((R)2.0)/((R)3.0)); \
	  (_coo)[7] = __mmkZELM_triangle_refcoo_x2*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x0*(((R)1.0)/((R)3.0)); \
	  (_coo)[8] = __mmkZELM_triangle_refcoo_x2*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_x0*(((R)2.0)/((R)3.0)); \
	  (_coo)[9] = (__mmkZELM_triangle_refcoo_x0+__mmkZELM_triangle_refcoo_x1+__mmkZELM_triangle_refcoo_x2)/((R)3.0); \
	  (_coo)[(_off)[0]+0] =__mmkZELM_triangle_refcoo_y0;		\
	  (_coo)[(_off)[0]+1] =__mmkZELM_triangle_refcoo_y1;		\
	  (_coo)[(_off)[0]+2] =__mmkZELM_triangle_refcoo_y2;		\
	  (_coo)[(_off)[0]+3] = __mmkZELM_triangle_refcoo_y0*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y1*(((R)1.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+4] = __mmkZELM_triangle_refcoo_y0*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y1*(((R)2.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+5] = __mmkZELM_triangle_refcoo_y1*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y2*(((R)1.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+6] = __mmkZELM_triangle_refcoo_y1*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y2*(((R)2.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+7] = __mmkZELM_triangle_refcoo_y2*(((R)2.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y0*(((R)1.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+8] = __mmkZELM_triangle_refcoo_y2*(((R)1.0)/((R)3.0))+__mmkZELM_triangle_refcoo_y0*(((R)2.0)/((R)3.0)); \
	  (_coo)[(_off)[0]+9] = (__mmkZELM_triangle_refcoo_y0+__mmkZELM_triangle_refcoo_y1+__mmkZELM_triangle_refcoo_y2)/((R)3.0); \
	  break;							\
	}								\
      default:								\
	{								\
	  I __mmkS_lagrange_localspl_triak=(I)0;		\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= __mmkZELM_triangle_refcoo_x0;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= __mmkZELM_triangle_refcoo_y0; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= __mmkZELM_triangle_refcoo_x1;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= __mmkZELM_triangle_refcoo_y1; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  (_coo)[__mmkS_lagrange_localspl_triak] 	= __mmkZELM_triangle_refcoo_x2;	\
	  (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak]	= __mmkZELM_triangle_refcoo_y2; \
	  ++__mmkS_lagrange_localspl_triak;				\
	  { I __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak] = ((R)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((R)0.0); \
	      }}							\
	  								\
	  { I __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak] = ((R)((_degree)[0]-__mmkS_lagrange_localspl_triaj))*__mmkS_lagrange_localspl_triax; \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((R)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
	      }}							\
	  								\
	  { I __mmkS_lagrange_localspl_triaj;			\
	    for (__mmkS_lagrange_localspl_triaj=1;			\
		 __mmkS_lagrange_localspl_triaj<(_degree)[0];		\
		 ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
	      {								\
		(_coo)[__mmkS_lagrange_localspl_triak]      	= ((R)0.0); \
		(_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] 	= ((R)((_degree)[0]-__mmkS_lagrange_localspl_triaj))*__mmkS_lagrange_localspl_triax; \
	      } }							\
									\
	  { I __mmkS_lagrange_localspl_triai;			\
	    for (__mmkS_lagrange_localspl_triai=1;			\
		 __mmkS_lagrange_localspl_triai<(_degree)[0]-1;		\
		 ++__mmkS_lagrange_localspl_triai)			\
	      { I __mmkS_lagrange_localspl_triaj;			\
		for (__mmkS_lagrange_localspl_triaj=1;			\
		     __mmkS_lagrange_localspl_triaj<(_degree)[0]-__mmkS_lagrange_localspl_triai; \
		     ++__mmkS_lagrange_localspl_triaj,++__mmkS_lagrange_localspl_triak) \
		  {							\
		    (_coo)[__mmkS_lagrange_localspl_triak]      = ((R)__mmkS_lagrange_localspl_triai)*__mmkS_lagrange_localspl_triax; \
		    (_coo)[(_off)[0]+__mmkS_lagrange_localspl_triak] = ((R)__mmkS_lagrange_localspl_triaj)*__mmkS_lagrange_localspl_triax; \
		  } }							\
	  } }								\
	break;								\
      }									\
  }								       


#define __mmkZELM_interval_refcoo_x0 -1.0
#define __mmkZELM_interval_refcoo_x1 1.0

#define __mmkS_lagrange_localspl_interval(_degree,_coo,_off)		\
  {									\
    switch((_degree)[0])						\
      {									\
      case 0:								\
	{								\
	  (_coo)[0] = (__mmkZELM_interval_refcoo_x0+__mmkZELM_interval_refcoo_x1)*( (R)0.5 ); \
	  break;							\
	}								\
      case 1:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  break;							\
	}								\
      case 2:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] =  __mmkZELM_interval_refcoo_x0*((R)0.5)+__mmkZELM_interval_refcoo_x1*((R)0.5); \
	  break;							\
	}								\
      case 3:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((R)2.0)/((R)3.0))+__mmkZELM_interval_refcoo_x1*(((R)1.0)/((R)3.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((R)1.0)/((R)3.0))+__mmkZELM_interval_refcoo_x1*(((R)2.0)/((R)3.0)); \
	  break;							\
	}								\
      case 4:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((R)3.0)/((R)4.0))+__mmkZELM_interval_refcoo_x1*(((R)1.0)/((R)4.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*((R)0.5)+__mmkZELM_interval_refcoo_x1*((R)0.5); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*(((R)1.0)/((R)4.0))+__mmkZELM_interval_refcoo_x1*(((R)3.0)/((R)4.0)); \
	  break;							\
	}								\
      case 5:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((R)4.0)/((R)5.0))+__mmkZELM_interval_refcoo_x1*(((R)1.0)/((R)5.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((R)3.0)/((R)5.0))+__mmkZELM_interval_refcoo_x1*(((R)2.0)/((R)5.0)); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*(((R)2.0)/((R)5.0))+__mmkZELM_interval_refcoo_x1*(((R)3.0)/((R)5.0)); \
	  (_coo)[5] = __mmkZELM_interval_refcoo_x0*(((R)1.0)/((R)5.0))+__mmkZELM_interval_refcoo_x1*(((R)4.0)/((R)5.0)); \
	  break;							\
	}								\
      case 6:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  (_coo)[2] = __mmkZELM_interval_refcoo_x0*(((R)5.0)/((R)6.0))+__mmkZELM_interval_refcoo_x1*(((R)1.0)/((R)6.0)); \
	  (_coo)[3] = __mmkZELM_interval_refcoo_x0*(((R)4.0)/((R)6.0))+__mmkZELM_interval_refcoo_x1*(((R)2.0)/((R)6.0)); \
	  (_coo)[4] = __mmkZELM_interval_refcoo_x0*((R)0.5)+__mmkZELM_interval_refcoo_x1*((R)0.5); \
	  (_coo)[5] = __mmkZELM_interval_refcoo_x0*(((R)2.0)/((R)6.0))+__mmkZELM_interval_refcoo_x1*(((R)4.0)/((R)6.0)); \
	  (_coo)[6] = __mmkZELM_interval_refcoo_x0*(((R)1.0)/((R)6.0))+__mmkZELM_interval_refcoo_x1*(((R)5.0)/((R)6.0)); \
	  break;							\
	}								\
      default:								\
	{								\
	  (_coo)[0] = __mmkZELM_interval_refcoo_x0;			\
	  (_coo)[1] = __mmkZELM_interval_refcoo_x1;			\
	  I _i;							\
	  const R _dt = (__mmkZELM_interval_refcoo_x1-__mmkZELM_interval_refcoo_x0)/((R)((_degree)[0])); \
	  for (_i=1;_i<(_degree)[0];++_i)				\
	    (_coo)[_i+1] = __mmkZELM_interval_refcoo_x0 + _dt*((R)_i); \
	  break;							\
	}								\
      }									\
  }



#if 0
void mkS_lagrange_localspl( eTopology elm_,
			    cst_pI 	degree_,
			    pR 	coo_,
			    cst_pI	off_)
{
  __mmkS_lagrange_localspl(elm_,degree_,coo_,off_);
}
#endif

void mkS_lagrange_localspl_interval( cst_pI 	degree_,
				     pR		coo_,
				     cst_pI		off_)
{
  __mmkS_lagrange_localspl_interval(degree_,coo_,off_);
}

void mkS_lagrange_localspl_tria( cst_pI 	degree_,
				 pR 	coo_,
				 cst_pI	off_)
{
  __mmkS_lagrange_localspl_tria(degree_,coo_,off_);
}
