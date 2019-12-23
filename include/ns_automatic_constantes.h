#ifndef __header_ns_automatic_constantes_h__
#define __header_ns_automatic_constantes_h__
/* AUTOMATIC MACRO */
#define 	_dimxdim 	(_dim*_dim)
#define 	_dimp1   	(_dim+1)
#define 	_dimp1   	(_dim+1)
#define 	_ngfacexnface 		((_dim+1)*_ngface)
#define 	_nuxdim			(_nu*_dim)
#define 	_nuxnu 			(_nu*_nu)
#define 	_nuxnuxnu 		(_nu*_nuxnu)
#define 	_npxnu 			(_np*_nu)
#define 	_ndgxnu 		(_ndg*_nu)
#define 	_ndgxndg 		(_ndg*_ndg)

/* AUTOMATIC CONSTANTES */

static const I 	dim		= ((I)_dim);
static const I 	dimp1		= ((I)_dimp1);
static const I 	dimxdim		= ((I)_dimxdim);

static const I   	np 		= (I)_np;
static const I   	ndg 		= (I)_ndg;
static const I   	nu 		= (I)_nu;

static const I   	nuxnu 		= (I)_nuxnu;
static const I   	nuxnp 		= (I)_npxnu;
static const I   	nuxndg 		= (I)_ndgxnu;
static const I   	ndgxndg 		= (I)_ndgxndg;
static const I   	nuxnuxnu	= (I)_nuxnuxnu;

static const I   	offu 		= (I)_ju;
static const I   	offv 		= (I)_ju+_nu;
static const I   	offw 		= (I)_ju+_nu*2;
static const I   	offp 		= (I)_jp;
static const I   	off  		= (I)_total_nddlelm;

static const I   	ng		= (I)_ng;
static const I 	ngfacex3	= (I)_ngfacexnface;
static const I 	ngface 		= (I)_ngface;
static const I 	total_nddlelmxtotal_nddlelm 	= (I)(_total_nddlelm*_total_nddlelm);

#endif
