#ifndef __ns_config__
#define __ns_config__

#define _shape_u 		__ensBASIS_LAGRANGE_2
#define _shape_p 		__ensBASIS_LAGRANGE_1
#define _shape_dg 		__ensBASIS_L2ORTHONORMAL_1

#define	_ndg 			3
#define	_dim 			__eDim_2
#define	_nu 			6
#define	_np 		        3
#define	_total_nddlelm  	(_nu*_dim+_np)
#define	_ngface 		7
#define	_ng 			36
#define	_ngfacex3 		(3*_ngface)

#define	_ju			0
#define	_jv			(_nu)
#define	_jw			(_nu*2)
#define	_jp			(_nu*_dim)

#ifdef __cplusplus
extern "C"
{
#endif
void ns_basename(const char * 	t_,
		 char		o_[512]);

#ifdef __cplusplus
}
#endif

#endif
