#ifndef __header_mkS_h__
#define __header_mkS_h__

#include "Type.h"
#include "eTopology.h"
#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif

enum __emkCONTINUITY	{ __emk_discontinuous = 0,
			  __emk_continuous };


#if 0
void 		mkS_basis	(cst_mkS S_,
				 cst_mkI 	n_,
				 mkR 	r_,
				 cst_mkI 	roff_,
				 cst_mkR	p_,
				 cst_mkI 	poff_,
				 mkR 	rwork_,
				 cst_mkI 	rwork_n_,
				 mkI  	ferr_);
#endif

typedef enum __emkCONTINUITY emkCONTINUITY;


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
#if 0
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
#endif

enum  __emkS_FAMILY { __emkS_FAMILY_error=0,
		       __emkS_FAMILY_lagrange,
		       __emkS_FAMILY_canonic,
		       __emkS_FAMILY_l2orthonormal,
		       __emkS_FAMILY_lagrangebubble,
		       __emkS_FAMILY_n };

typedef enum __emkS_FAMILY emkS_FAMILY;

enum  __emk_s_diff { 
  __emk_s_diff_no=0,
  __emk_s_diff_dx,
  __emk_s_diff_dy,
  __emk_s_diff_dz,
  __emk_s_diff_n };

typedef enum __emk_s_diff emk_s_diff;

typedef struct
{
  eTopology 		element;
  emkS_FAMILY 		family;
  emk_s_diff       	diff;
  I 			nshape;
  I 			degree;
  emk_s            	shape;
  emkCONTINUITY       	continu;

#if 0
  //  emkZELM 		element;
  


#endif
} mkS_st,*mkS;


void mkS_def	(mkS 			shape_,
		 cst_emk_s 		kind_,
		 Err * 		err_);

void 			mkS_definit			(mkS			shape_,
							 const eTopology    	element_,
							 emkS_FAMILY    	family_,
							 const I 		degree_,
							 emkCONTINUITY   	continuity_,
							 Err * 		err_);


#if 0
#define			mkS_lagrange_localspl_tria     __mmkS_lagrange_localspl_tria
#define			mkS_lagrange_localspl_quad     __mmkS_lagrange_localspl_quad
#define			mkS_lagrange_localspl_tetra    __mmkS_lagrange_localspl_tetra
#define			mkS_lagrange_localspl_interval __mmkS_lagrange_localspl_interval
//#define			mkS_lagrange_localspl          __mmkS_lagrange_localspl
#endif


void __mmkS_lagrange_localspl_tria( cst_pI 	_degree,
				 pR 	_coo,
				    cst_pI	_off);

typedef const mkS_st 	cst_mkS_st;
typedef cst_mkS_st* 	cst_mkS;

emk_s_diff 		mkS_diff	(cst_mkS S_);
I 			mkS_k		(cst_mkS S_);
I 			mkS_n		(cst_mkS S_);

void 			mkS_copy	(mkS 	 dest_,
					 cst_mkS src_);
emkS_FAMILY 		mkS_family	(cst_mkS S_);
cst_mkS 		mkS_b		(mkS shape_);
cst_mkS 		mkS_dx		(mkS shape_);
cst_mkS 		mkS_dy		(mkS shape_);
cst_mkS 		mkS_dz		(mkS shape_);

  cst_mkS 		mkS_derivative	(mkS shape_,const I idim);
  
  void mkS_lagrange_localspl_tria( cst_pI 	_degree,
				   pR 	_coo,
				   cst_pI	_off);
  

  //#define			mkS_lagrange_localspl_tria     __mmkS_lagrange_localspl_tria

#if 0



  
#if 0
#include "emkS.h"
#include "mkS_macros.h"
#include "emkCONTINUITY.h"

#include "mok_m.h"

#if __mk_debug__
#define __mkS_debug(_P,_where) ( (_P) ? (void)0 : (fflush(stdout),fflush(stderr),fprintf(stderr,"*** mkS-debug-info\n\troutine '%s'\n\tfile '%s'\n\tline '%d'\n\tparameter '%s' is null\n",#_where,__FILE__,__LINE__,#_P)) )
#else
#define __mkS_debug(_P,_where) ((void)0)
#endif
void mkS_errmsg(const char * msg,...);

#endif




  
mkS 			mkS_new_transformation		(cst_emkZELM		elm_,
							 emk_err * 		err_);

void 			mkS_def_transformation		(mkS 			shape_,
							 cst_emkZELM		elm_,
							 emk_err * 		err_);

mkS 			mkS_kill			(mkS 			shape_);
void 			mkS_def				(mkS 			shape_,
							 cst_emk_s    		kind_,
							 emk_err * 		err_);
mkS 			mkS_new				(cst_emk_s 		kind_);

mkS 			mkS_newinit			(cst_emkZELM    	element_,
							 cst_emkS_FAMILY    	family_,
							 const I 		degree_,
							 cst_emkCONTINUITY   	continuity_,
							 emk_err * 		err_);

mkS 			mkS_dup				(cst_mkS		A_);
emk_err 		mkS_fread			(FILE * fich_,mkS shape_);
emk_err 		mkS_fwrite			(FILE * fich_,cst_mkS S_);

I			mkS_nshapevertex		(cst_mkS S_);
I 			mkS_nshapedge			(cst_mkS S_);
I 			mkS_nshapelm			(cst_mkS S_);
emk_s              	mkS_enum			(cst_mkS S_);
emk_logic 		mkS_cmp				(cst_mkS a_,
							 cst_mkS b_);
mok_m 			mkS_sample_new_mok_m		(cst_mkS S_);

  
#if __mk_debug__

void 			mkS_lagrange_localspl_tria	( cst_mkI 	degree,
							  mkR 	coo_,
							  cst_mkI	off_);
void 			mkS_lagrange_localspl_tetra	( cst_mkI 	degree,
							  mkR 	coo_,
							  cst_mkI	off_);
void 			mkS_lagrange_localspl_interval	( cst_mkI 	degree,
							  mkR 	coo_,
							  cst_mkI	off_);
void 			mkS_lagrange_localspl_quad	( cst_mkI 	degree,
							  mkR 	coo_,
							  cst_mkI	off_);

void 			mkS_lagrange_localspl		( cst_emkZELM 	elm_,
							  cst_mkI 	degree,
							  mkR 	coo_,
							  cst_mkI	off_);






emkCONTINUITY 		mkS_continuity	(cst_mkS S_);

emkZELM        		mkS_emkZELM	(cst_mkS shape_);

#else


#define 		mkS_emkZELM(_S) (_S)->element

#define			mkS_lagrange_localspl_tria     __mmkS_lagrange_localspl_tria
#define			mkS_lagrange_localspl_quad     __mmkS_lagrange_localspl_quad
#define			mkS_lagrange_localspl_tetra    __mmkS_lagrange_localspl_tetra
#define			mkS_lagrange_localspl_interval __mmkS_lagrange_localspl_interval
#define			mkS_lagrange_localspl          __mmkS_lagrange_localspl
#define 		mkS_family(_A) 		(_A)->family
#define 		mkS_k(_A) 		(_A)->degree
#define 		mkS_n(_A) 		(_A)->nshape
#define 		mkS_continuity(_A) 	(_A)->continu
#define 		mkS_diff(_A) 		(_A)->diff

#define 		mkS_copy(_A,_B)			\
  (_A)->element = (_B)->element;	\
  (_A)->family 	= (_B)->family;				\
  (_A)->diff 	= (_B)->diff;				\
  (_A)->shape 	= (_B)->shape;				\
      (_A)->degree 	= (_B)->degree;			\
	  (_A)->nshape 	= (_B)->nshape;			\
	      (_A)->continu  = (_B)->continu


#define 		mkS_b(_A) 	( (_A)->diff=__emk_s_diff_no , (_A) )
#define 		mkS_dx(_A) 	( (_A)->diff=__emk_s_diff_dx , (_A) )
#define 		mkS_dy(_A) 	( (_A)->diff=__emk_s_diff_dy , (_A) )
#define 		mkS_dz(_A) 	( (_A)->diff=__emk_s_diff_dz , (_A) )

#endif
cst_mkS mkS_continuous(mkS shape_);
cst_mkS mkS_discontinuous(mkS shape_);
mok_m mkS_shape_new_mok_m(cst_mkS S_,
			  cst_mok_m lc_,
			  const char * Tr);
void   mkS_shape_def_mok_m(cst_mkS 	S_,
			   cst_mok_m 	lc_,
			   mok_m 		m_,
			   const char * 	Tr,
			   nsREAL*		rwork,
			   const I * 	rwork_n,
			   emk_err * 	err_);



//#define mkS_basis 	mkS_basis_

void 		mkS_basis_tria_	(cst_mkS S_,
				 cst_mkI 	n_,
				 mkR 	r_,
				 cst_mkI 	roff_,
				 cst_mkR	p_,
				 cst_mkI 	poff_,
				 mkR 	rwork_,
				 cst_mkI 	rwork_n_,
				 mkI  	ferr_);

#define mkS_basis_tria 	mkS_basis_tria_

void 		mkS_basis_quad_	(cst_mkS S_,
				 cst_mkI 	n_,
				 mkR 	r_,
				 cst_mkI 	roff_,
				 cst_mkR	p_,
				 cst_mkI 	poff_,
				 mkR 	rwork_,
				 cst_mkI 	rwork_n_,
				 mkI  	ferr_);

#define mkS_basis_quad 	mkS_basis_quad_

void 	mkS_basis_tetra_		(cst_mkS S_,
					 cst_mkI 	n_,
					 mkR 	r_,
					 cst_mkI 	roff_,
					 cst_mkR	p_,
					 cst_mkI 	poff_,
					 mkR 	rwork_,
					 cst_mkI 	rwork_n_,
					 mkI  	ferr_);

#define mkS_basis_tetra 	mkS_basis_tetra_

void 	mkS_lagr_spl_ 		(cst_mkI 		degree,
				 mkR		coo,
				 cst_mkI 		off_);



void 	mkS_bwji_nei		(cst_mkI  		qface_n_,
				 cst_mkR 		qface_p_,
				 cst_mkR 		qface_w_,
				 cst_mkS	trial_,
				 cst_mkS	test_,
				 mkR 		x_,
				 const I 		xoff_,
				 mkI 		ferr_);

void 	mkS_spl_lc2gl_		(cst_mkI 		degree_,
				 mkR 		L_,
				 cst_mkI 		Loff_);

#define mkS_lagr_spl 		mkS_lagr_spl_





/* HERMITE ########################################################################################################### */

#define __mmkS_hermite_interval(_n,_r,_roff,_p,_poff) { I _i;for (_i=0;_i<(_n);++_i) { const nsREAL s  = ((_p)[_i*(_poff)]+1.0)*0.5; const nsREAL s2 = s*s; (_r)[(_roff)*_i+0] = (s*2.0-3.0)*s2+1.0; (_r)[(_roff)*_i+1] = (s*(-2.0)+3.0)*s2;  (_r)[(_roff)*_i+2] = s*(s-1.0)*(s-1.0);  (_r)[(_roff)*_i+3] = s2*(s-1.0);  }}
#define __mmkS_hermite_interval_dt(_n,_r,_roff,_p,_poff) { I _i;for (_i=0;_i<(_n);++_i) { const nsREAL s  = ((_p)[_i*(_poff)]+1.0)*0.5;  (_r)[(_roff)*_i+0] = (s-1.0)*s*6.0;(_r)[(_roff)*_i+1] = (1.0-s)*s*6.0;(_r)[(_roff)*_i+2] = s*(s*3.0-4.0)+1.0;(_r)[(_roff)*_i+3] = s*(s*3.0-2.0); }}
#define __mmkS_hermite_interval_dtt(_n,_r,_roff,_p,_poff) { I _i;for (_i=0;_i<(_n);++_i) { const nsREAL s  = ((_p)[_i*(_poff)]+1.0)*0.5; (_r)[(_roff)*_i+0] = (s*12.0-6.0)*0.25; (_r)[(_roff)*_i+1] = (6.0-s*12.0)*0.25; (_r)[(_roff)*_i+2] = (s*6.0-4.0)*0.25; (_r)[(_roff)*_i+3] = (s*6.0-2.0)*0.25; }}

/*
  x  = x0 * (2.0*t*t*t-3.0*t*t+1.0) + tx0 * (t*t*t-2*t*t+t) + x1 * (-2.0*t*t*t+3.0*t*t)  + tx1 * (t*t*t-t*t);      
  y  = y0 * (2.0*t*t*t-3.0*t*t+1.0) + ty0 * (t*t*t-2*t*t+t) + y1 * (-2.0*t*t*t+3.0*t*t)  + ty1 * (t*t*t-t*t);
*/

#if __mk_debug__
void mkS_hermite_interval_	(cst_mkI 	n,
				 mkR		r,
				 cst_mkI 	roff_,
				 cst_mkR 	p,
				 cst_mkI 	poff_);

void mkS_hermite_interval_dt_	(cst_mkI 	n,
				 mkR		r,
				 cst_mkI 	roff_,
				 cst_mkR 	p,
				 cst_mkI 	poff_);

void mkS_hermite_interval_dtt_	(cst_mkI 	n,
				 mkR		r,
				 cst_mkI 	roff_,
				 cst_mkR 	p,
				 cst_mkI 	poff_);
#else

#define mkS_hermite_interval_(_n,_r,_roff,_p,_poff) 	__mmkS_hermite_interval((_n)[0],_r,(_roff)[0],_p,(_poff)[0])
#define mkS_hermite_interval_dt_(_n,_r,_roff,_p,_poff) 	__mmkS_hermite_interval_dt((_n)[0],_r,(_roff)[0],_p,(_poff)[0])
#define mkS_hermite_interval_dtt_(_n,_r,_roff,_p,_poff) __mmkS_hermite_interval_dtt((_n)[0],_r,(_roff)[0],_p,(_poff)[0])

#endif

#define mkS_hermite_interval(_n,_r,_roff,_p,_poff) 	mkS_hermite_interval_(_n,_r,_roff,_p,_poff)
#define mkS_hermite_interval_dt(_n,_r,_roff,_p,_poff) 	mkS_hermite_interval_dt_(_n,_r,_roff,_p,_poff)
#define mkS_hermite_interval_dtt(_n,_r,_roff,_p,_poff) 	mkS_hermite_interval_dtt_(_n,_r,_roff,_p,_poff)

#endif



void 	mkS_kji			(cst_mkS   	teta_,
				 cst_mkS 	trial_,
				 cst_mkS 	test_,
				 pR 		x_,
				 cst_pI		xoff_,
				 cst_pI 		qelm_n_,
				 cst_pR 		qelm_w_,
				 cst_pR		qelm_pos_,
				 cst_pI 		qelm_posoff_,
				 cst_pI 		rwork_n_,
				 pR 		rwork_,
				 pI 		err_);

  
void 	mkS_bwji			(cst_pI  		qface_n_,
					 cst_pR 		qface_p_,
					 cst_pR 		qface_w_,
					 cst_mkS	trial_,
					 cst_mkS	test_,
					 pR 		x_,
					 const I 		xoff_,
					 pI 		ferr_);


  void 		mkS_bmapping		(const I 	n,
					 pR 		r,
					 cst_pI 	roff,
					 cst_pR 	p);

  void 		mkS_basis	(cst_mkS S_,
				 cst_pI 	n_,
				 pR 	r_,
				 cst_pI 	roff_,
				 cst_pR	p_,
				 cst_pI 	poff_,
				 pR 	rwork_,
				 cst_pI 	rwork_n_,
				 pI  	ferr_);
void 	mkS_bwji_nei		(cst_pI  		qface_n_,
				 cst_pR 		qface_p_,
				 cst_pR 		qface_w_,
				 cst_mkS	trial_,
				 cst_mkS	test_,
				 pR 		x_,
				 const I 		xoff_,
				 pI 		ferr_);


  //  #define 	mkQ_collapse 		mkQ_collapse_
  //#define 	mkQ_legendre_interval 	mkQ_legendre_interval_

void  		mkQ_collapse_		(cst_pI  		qinterval_n_,
					 cst_pR 		qinterval_p_,
					 cst_pI  		qinterval_poff_,
					 cst_pR 		qinterval_w_,
					 cst_pI  		qinterval_woff_,
					 pI 			qelm_n_,
					 pR 			qelm_p_,
					 cst_pI 		qelm_poff_,
					 pR 			qelm_w_,
					 cst_pI 		qelm_woff_,
					 pI 			ferr_);

void 		mkQ_legendre_interval_	(cst_pI 		nspl_,
					 pR 			p_,
					 cst_pI 		p_off_, 
					 pR 			w_,
					 cst_pI 		w_off_, 
					 pR 			work_,
					 cst_pI 		work_n_,
					 pI 			ferr_);

#ifdef __cplusplus
}
#endif

#endif
