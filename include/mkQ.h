#ifndef __HEADER_mkQ_H__
#define __HEADER_mkQ_H__

#include "eTopology.h"
#include "Err.h"

enum __emkQ { __emkQ_error = 0, 			
	      __emkQ_HAMMER_TRIANGLE_1,
	      __emkQ_HAMMER_TRIANGLE_3,
	      __emkQ_HAMMER_TRIANGLE_4,
	      __emkQ_HAMMER_TRIANGLE_6,
	      __emkQ_HAMMER_TRIANGLE_12,
	      __emkQ_HAMMER_TRIANGLE_42,
	      __emkQ_HAMMER_INTERVAL_1,
	      __emkQ_HAMMER_INTERVAL_2,
	      __emkQ_HAMMER_INTERVAL_3,
	      __emkQ_HAMMER_INTERVAL_4,
	      __emkQ_HAMMER_INTERVAL_5,
	      __emkQ_HAMMER_TETRAHEDRA_1,
	      __emkQ_HAMMER_TETRAHEDRA_4,
	      __emkQ_HAMMER_TETRAHEDRA_16,
	      __emkQ_GAUSS_INTERVAL_1,
	      __emkQ_GAUSS_INTERVAL_2,
	      __emkQ_GAUSS_INTERVAL_3,
	      __emkQ_GAUSS_INTERVAL_4,
	      __emkQ_GAUSS_INTERVAL_5,
	      __emkQ_GAUSS_INTERVAL_6,
	      __emkQ_GAUSS_INTERVAL_7,
	      __emkQ_GAUSS_INTERVAL_8,
	      __emkQ_GAUSS_INTERVAL_9,
	      __emkQ_GAUSS_INTERVAL_10,
	      __emkQ_GAUSS_TRIANGLE_4,
	      __emkQ_GAUSS_TRIANGLE_9,
	      __emkQ_GAUSS_TRIANGLE_16,
	      __emkQ_GAUSS_TRIANGLE_25,
	      __emkQ_GAUSS_TRIANGLE_36,
	      __emkQ_GAUSS_TRIANGLE_49,
	      __emkQ_GAUSS_TRIANGLE_64,
	      __emkQ_GAUSS_TRIANGLE_81,
	      __emkQ_GAUSS_TRIANGLE_100,
	      __emkQ_GAUSS_TRIANGLE_121,			 
	      __emkQ_GAUSS_QUADRANGLE_4,
	      __emkQ_GAUSS_QUADRANGLE_9,
	      __emkQ_GAUSS_QUADRANGLE_16,
	      __emkQ_GAUSS_QUADRANGLE_25,
	      __emkQ_GAUSS_QUADRANGLE_36,
	      __emkQ_GAUSS_QUADRANGLE_49,
	      __emkQ_GAUSS_QUADRANGLE_64,
	      __emkQ_GAUSS_QUADRANGLE_81,
	      __emkQ_GAUSS_QUADRANGLE_100,
	      __emkQ_GAUSS_QUADRANGLE_121,			 
	      __emkQ_GAUSS_TETRAHEDRA_8,
	      __emkQ_GAUSS_TETRAHEDRA_27,
	      __emkQ_GAUSS_TETRAHEDRA_64,
	      __emkQ_GAUSS_TETRAHEDRA_125,
	      __emkQ_GAUSS_TETRAHEDRA_216,
	      __emkQ_GAUSS_TETRAHEDRA_343,
	      __emkQ_HAMMER_QUADRANGLE_1,
	      __emkQ_HAMMER_QUADRANGLE_2,
	      __emkQ_HAMMER_QUADRANGLE_3,
	      __emkQ_HAMMER_QUADRANGLE_4,
	      __emkQ_HAMMER_QUADRANGLE_5,
	      __emkQ_DUNAVANT_TRIANGLE_1,
	      __emkQ_DUNAVANT_TRIANGLE_3,
	      __emkQ_DUNAVANT_TRIANGLE_4,
	      __emkQ_DUNAVANT_TRIANGLE_6,
	      __emkQ_DUNAVANT_TRIANGLE_7,
	      __emkQ_DUNAVANT_TRIANGLE_12,
	      __emkQ_DUNAVANT_TRIANGLE_13,
	      __emkQ_DUNAVANT_TRIANGLE_16,
	      __emkQ_DUNAVANT_TRIANGLE_19,
	      __emkQ_adaptative,
	      __emkQ_n }; 

typedef enum __emkQ emkQ;

enum __emkQ_family 	{ __emkQ_family_error = 0, 			
			  __emkQ_family_gauss, 
			  __emkQ_family_hammer, 
			  __emkQ_family_dunavant, 
			  __emkQ_family_n}; 

typedef enum __emkQ_family emkQ_family;

typedef const emkQ_family cst_emkQ_family;
typedef const emkQ cst_emkQ;

#if 0
#include "emk_template.h"

__emk_template_headers(emkQ_family);
__emk_template_headers(emkQ);
#endif

typedef struct
{
  emkQ 		kind;
  eTopology     elm;

  cst_pI       qn;
  cst_pR       qp;
  cst_pR       qw;

  I 	own_qn;
  pR           own_qp;
  pR           own_qw;

  cst_pR       lc2gl;
  cst_pR       lc2gl_dx;
  cst_pR       lc2gl_dy;
  cst_pR       lc2gl_dz;

  pR 		shape_data;
} mkQ_ST,*mkQ;

typedef const mkQ_ST cst_mkQ_ST;
typedef cst_mkQ_ST * cst_mkQ;

#ifndef NDEBUG

#define 	mkQ_w(_Q)   		(_Q)->qw
#define 	mkQ_p(_Q)   		(_Q)->qp
#define 	mkQ_n(_Q) 		(_Q)->qn[0]
#define 	cst_mkQ_w(_Q)   	(_Q)->qw
#define 	cst_mkQ_p(_Q)   	(_Q)->qp
#define 	mkQ_enum(_Q)   		(_Q)->kind
#define 	mkQ_exact_for(_Q,_k) 	((_k)+2)/2+((_k)+2)%2

#else

I 		mkQ_n			(cst_mkQ 		Q_);
cst_pR 	cst_mkQ_w		(cst_mkQ 		Q_);
cst_pR 	cst_mkQ_p		(cst_mkQ 		Q_);
pR 		mkQ_w			(mkQ 			Q_);
pR 		mkQ_p			(mkQ 			Q_);
emkQ 		mkQ_enum   		(cst_mkQ 		Q_);
I 		mkQ_exact_for		(cst_emkQ_family 	family_,
					 const I 		k);

#endif

mkQ  		mkQ_kill		(mkQ 			Q_);

void    	mkQ_def			(mkQ 			Q_,
					 cst_emkQ_family 	family_,
					 const eTopology		elm_,
					 cst_pI 		degree_,
					 Err*		err_);

mkQ		mkQ_new			(cst_emkQ_family 	family_,
					 const eTopology 		elm_,
					 cst_pI 		degree_);

emkQ 		mkQ_compute_enum	(cst_emkQ_family 	family_,
					 const eTopology		elm_,
					 cst_pI 		degree_);


void 		mkQ_def_gauss		(mkQ 			q_,
					 const eTopology 		elm_,
					 const I 		nspl_,
					 Err*		err_);


void 		mkQ_def_cadyf		(mkQ 			q_,
			      		Err*		err_);


void    	mkQ_def_adaptative	(mkQ 			q_,
					 const I 		N_,
					 const eTopology		elm_,
					 Err*		err_);

mkQ		mkQ_new_adaptative	(const I 		N,
					 const eTopology 		elm_);

#define 	mkQ_collapse 		mkQ_collapse_
#define 	mkQ_legendre_interval 	mkQ_legendre_interval_

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

void 		mkQ_write_eps		(const char * 		name_,
					 cst_mkQ 		q_,
					 const L 		q_show_,
					 const L 		e_show_);

#endif
