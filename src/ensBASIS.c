#include "Type.h"
#include "ns_config_lapack.h"
#include "ensBASIS.h"

typedef void 	 (*__eval_ensbasis_n)	(cst_ensBASIS );
typedef void 	 (*__eval_ensbasis_eval)(cst_pS ,cst_pS ,cst_pI ,cst_pR ,cst_pI ,pR ,cst_pI );
extern void 	 ensbasis_edge_lagrbubble0_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagrbubble0_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_edge_lagrbubble1_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagrbubble1_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_edge_lagrbubble2_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagrbubble2_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_lagr0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagr0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_lagr1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagr1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_lagr2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagr2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_lagr3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_lagr3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_ortho5_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_ortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);


extern void 	 ensbasis_edge_drlagrbubble0_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagrbubble0_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_edge_drlagrbubble1_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagrbubble1_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_edge_drlagrbubble2_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagrbubble2_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drlagr0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagr0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drlagr1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagr1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drlagr2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagr2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drlagr3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drlagr3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_edge_drortho5_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_drortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);






extern void 	 ensbasis_tria_dslagrbubble0_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dslagrbubble1_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dslagrbubble2_	(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_tria_dslagr0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dslagr1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dslagr2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dslagr3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);

extern void 	 ensbasis_tria_dsortho0_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dsortho1_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dsortho2_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dsortho3_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dsortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);
extern void 	 ensbasis_tria_dsortho4_		(cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_);


__eval_ensbasis_eval ensbasis_edge_eval_table [__ensBASIS_ALL] = {NULL,
								  ensbasis_edge_lagr0_,
								  ensbasis_edge_lagr1_,
								  ensbasis_edge_lagr2_,
								  ensbasis_edge_lagr3_,
								  ensbasis_edge_ortho0_,
								  ensbasis_edge_ortho1_,
								  ensbasis_edge_ortho2_,
								  ensbasis_edge_ortho3_,
								  NULL,
								  NULL,
								  NULL};
  

__eval_ensbasis_eval ensbasis_tria_eval_table [__ensBASIS_ALL] = {NULL,
								ensbasis_tria_lagr0_,
								ensbasis_tria_lagr1_,
								ensbasis_tria_lagr2_,
								ensbasis_tria_lagr3_,
								ensbasis_tria_ortho0_,
								ensbasis_tria_ortho1_,
								ensbasis_tria_ortho2_,
								ensbasis_tria_ortho3_,
								ensbasis_tria_lagrbubble0_,
								ensbasis_tria_lagrbubble1_,
								ensbasis_tria_lagrbubble2_};




__eval_ensbasis_eval ensbasis_edge_dreval_table [__ensBASIS_ALL] = {NULL,
								ensbasis_edge_drlagr0_,
								ensbasis_edge_drlagr1_,
								ensbasis_edge_drlagr2_,
								ensbasis_edge_drlagr3_,
								ensbasis_edge_drortho0_,
								ensbasis_edge_drortho1_,
								ensbasis_edge_drortho2_,
								ensbasis_edge_drortho3_,
								  NULL,
								NULL,
								NULL};
  

__eval_ensbasis_eval ensbasis_tria_dreval_table [__ensBASIS_ALL] = {NULL,
								ensbasis_tria_drlagr0_,
								ensbasis_tria_drlagr1_,
								ensbasis_tria_drlagr2_,
								ensbasis_tria_drlagr3_,
								ensbasis_tria_drortho0_,
								ensbasis_tria_drortho1_,
								ensbasis_tria_drortho2_,
								ensbasis_tria_drortho3_,
								ensbasis_tria_drlagrbubble0_,
								ensbasis_tria_drlagrbubble1_,
								ensbasis_tria_drlagrbubble2_};


  

__eval_ensbasis_eval ensbasis_tria_dseval_table [__ensBASIS_ALL] = {NULL,
								  ensbasis_tria_dslagr0_,
								  ensbasis_tria_dslagr1_,
								  ensbasis_tria_dslagr2_,
								  ensbasis_tria_dslagr3_,
								  ensbasis_tria_dsortho0_,
								  ensbasis_tria_dsortho1_,
								  ensbasis_tria_dsortho2_,
								  ensbasis_tria_dsortho3_,
								  ensbasis_tria_dslagrbubble0_,
								  ensbasis_tria_dslagrbubble1_,
								  ensbasis_tria_dslagrbubble2_};

								   
static __eval_ensbasis_eval * ensbasis_eval_table[__eFace_ALL] 	= {NULL,
								   ensbasis_tria_eval_table,
								   NULL,
								   NULL};
static __eval_ensbasis_eval * ensbasis_dreval_table[__eFace_ALL] = 	{NULL,
									 ensbasis_tria_dreval_table,
									 NULL,
									 NULL};

static __eval_ensbasis_eval * ensbasis_dseval_table[__eFace_ALL] = 	 {NULL,
									  ensbasis_tria_dseval_table,
									  NULL,
									  NULL};

static __eval_ensbasis_eval * ensbasis_dteval_table[__eFace_ALL] = 	{NULL,
									 NULL,
									 NULL,
									 NULL};
  
static const I 		ensbasis_tria_n_table[__ensBASIS_ALL] 	= {0,1,3,6,10,1,3,6,10,2,4,7};
static const I 		ensbasis_edge_n_table[__ensBASIS_ALL] 	= {0,1,2,3,4,1,2,3,4,2,3,4};
static const I *	ensbasis_n_table[__eFace_ALL] 		= {NULL,
								   ensbasis_tria_n_table,
								   NULL,
								   NULL};

static const I 		ensbasis_degree_table[__ensBASIS_ALL] = {-1,0,1,2,3,0,1,2,3,0,1,2};


I	 ensBASIS_n		(cst_ensBASIS b_,cst_eFace elm_)	
{ 
  return ensbasis_n_table[elm_][b_] ; 
};

I	 ensBASIS_degree	(cst_ensBASIS b_)			{ return ensbasis_degree_table[b_] ; };



void 	 ensBASIS_b		(cst_ensBASIS basis_,cst_eFace elm_,cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_)
{
  ensbasis_eval_table[elm_][basis_] (ptr_,ptb_,n_,p_,poff_,b_,boff_);
}


void 	 ensBASIS_dr		(cst_ensBASIS basis_,cst_eFace elm_,cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_)
{
  ensbasis_dreval_table[elm_][basis_] (ptr_,ptb_,n_,p_,poff_,b_,boff_);
}


void 	 ensBASIS_ds		(cst_ensBASIS basis_,cst_eFace elm_,cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_)
{
  ensbasis_dseval_table[elm_][basis_] (ptr_,ptb_,n_,p_,poff_,b_,boff_);
}

void 	 ensBASIS_dt		(cst_ensBASIS basis_,cst_eFace elm_,cst_pS ptr_,cst_pS ptb_,cst_pI n_,cst_pR p_,cst_pI poff_,pR b_,cst_pI boff_)
{
  ensbasis_dteval_table[elm_][basis_] (ptr_,ptb_,n_,p_,poff_,b_,boff_);
}

static const I 		ensbasis_tria_n_onface_table[__ensBASIS_ALL] 	= {0,0,1,2,3,0,1,2,3,0,0,1};
static const I 		ensbasis_edge_n_onface_table[__ensBASIS_ALL] 	= {0,1,1,1,1,1,1,1,1,1,1,1};
static const I * 	ensbasis_n_onface_table[__eFace_ALL] 		= {NULL,
									   ensbasis_tria_n_onface_table,
									   NULL,
									   NULL};

I ensBASIS_n_onface(cst_ensBASIS basis_,cst_eFace elm_)
{
  return (ensbasis_n_onface_table[elm_])[basis_];
}


static const R __cootreilli_TRIA_LAGRANGE_0[2] 	= {(R)3.333333333333333e-01,(R)3.333333333333333e-01};
static const R __cootreilli_TRIA_LAGRANGE_1[3*2] 	= {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,
							   (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00};

static const R __cootreilli_TRIA_LAGRANGE_2[6*2] 	= {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)0.000000000000000e+00,
							   (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01};

static const R __cootreilli_TRIA_LAGRANGE_3[10*2] = {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,(R)6.666666666666666e-01,(R)6.666666666666666e-01,(R)3.333333333333333e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
							  (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,(R)6.666666666666666e-01,(R)6.666666666666666e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01};


static const R __cootreilli_TRIA_LAGRANGE_4[15*2] = {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)7.500000000000000e-01,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)2.500000000000000e-01,(R)5.000000000000000e-01,
						      (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)7.500000000000000e-01,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01};

static const R __cootreilli_TRIA_LAGRANGE_5[21*2] = {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)8.000000000000000e-01,(R)8.000000000000000e-01,(R)6.000000000000001e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,
						      (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)8.000000000000000e-01,(R)8.000000000000000e-01,(R)6.000000000000001e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01};

static const R __cootreilli_TRIA_LAGRANGE_6[28*2] = {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)8.333333333333333e-01,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,
							  
							  (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)8.333333333333333e-01,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01};



static const R __cootreilli_TRIA_LAGRANGEBUBBLE_0[4] 	= {(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01};
static const R __cootreilli_TRIA_LAGRANGEBUBBLE_1[4*2] 	= {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
								   (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)3.333333333333333e-01};

static const R __cootreilli_TRIA_LAGRANGEBUBBLE_2[7*2] = {(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
						     (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01};



static cst_pR __cootreilli_TRIA[__ensBASIS_ALL] = {
  NULL,
  __cootreilli_TRIA_LAGRANGE_0,
  __cootreilli_TRIA_LAGRANGE_1,
  __cootreilli_TRIA_LAGRANGE_2,
  __cootreilli_TRIA_LAGRANGE_3,
  __cootreilli_TRIA_LAGRANGE_0,
  __cootreilli_TRIA_LAGRANGE_1,
  __cootreilli_TRIA_LAGRANGE_2,
  __cootreilli_TRIA_LAGRANGE_3,
  __cootreilli_TRIA_LAGRANGEBUBBLE_0,
  __cootreilli_TRIA_LAGRANGEBUBBLE_1,
  __cootreilli_TRIA_LAGRANGEBUBBLE_2};


static cst_pR *  __cootreilli[__eFace_ALL] = {NULL,
					      __cootreilli_TRIA,
					      NULL,
					      NULL};


cst_pR 	ensBASIS_cootreilli		(cst_ensBASIS shape_,cst_eFace elm_)
{
  cst_pR * __cootreilli_elm = __cootreilli[elm_];
  return __cootreilli_elm[shape_];
}


static const I __cnctreilli_TRIA_LAGRANGEBUBBLE_0[9] 	= {(I)0,(I)1,(I)3,
							   (I)1,(I)2,(I)3,
							   (I)2,(I)0,(I)3};

static const I __cnctreilli_TRIA_LAGRANGEBUBBLE_1[9] 	= {(I)0,(I)1,(I)3,
							   (I)1,(I)2,(I)3,
							   (I)2,(I)0,(I)3};

static const I __cnctreilli_TRIA_LAGRANGEBUBBLE_2[18] 	= {(I)0,(I)3,(I)5,
							   (I)3,(I)1,(I)4,
							   (I)4,(I)2,(I)5,
							   (I)3,(I)4,(I)7,
							   (I)4,(I)5,(I)7,
							   (I)5,(I)3,(I)7};


static const I __cnctreilli_TRIA_LAGRANGE_1[3] 		= {(I)0,(I)1,(I)2};
static const I __cnctreilli_TRIA_LAGRANGE_2[12] 	= {(I)0,(I)3,(I)5,
							   (I)3,(I)1,(I)4,
							   (I)4,(I)2,(I)5,
							   (I)3,(I)4,(I)5};

static const I __cnctreilli_TRIA_LAGRANGE_3[18] 	= {(I)0,(I)3,(I)6,
							   (I)3,(I)1,(I)6,
							   (I)1,(I)4,(I)6,
							   (I)4,(I)2,(I)6,
							   (I)2,(I)5,(I)6,
							   (I)5,(I)0,(I)6};


static const I __cnctreilli_TRIA_LAGRANGE_4[27]  	= {(I)0,(I)3,(I)8,
							   (I)3,(I)9,(I)8,
							   (I)3,(I)4,(I)9,
							   (I)4,(I)5,(I)9,
							   (I)4,(I)1,(I)5,
							   (I)9,(I)5,(I)6,
							   (I)9,(I)6,(I)7,
							   (I)9,(I)7,(I)8,
							   (I)6,(I)2,(I)7};

static const I __cnctreilli_TRIA_LAGRANGE_5[48] 	= {(I)0,(I)3,(I)11,
							   (I)3,(I)12,(I)11,
							   (I)11,(I)12,(I)10,
							   (I)12,(I)13,(I)10,
							   (I)10,(I)13,(I)9,
							   (I)13,(I)8,(I)9,
							   (I)9,(I)8,(I)2,
							   (I)3,(I)4,(I)12,
							   (I)4,(I)14,(I)12,
							   (I)12,(I)14,(I)13,
							   (I)14,(I)7,(I)13,
							   (I)7,(I)8,(I)13,
							   (I)4,(I)5,(I)14,
							   (I)5,(I)6,(I)14,
							   (I)6,(I)7,(I)14,
							   (I)5,(I)1,(I)6};

/*
  2
  11 10
  12 17  9 
  13 16 19  8
  14 15 18 20 7
  0   3  4  5 6 1
*/
static const I __cnctreilli_TRIA_LAGRANGE_6[75] 	= {(I)0,(I)3,(I)14,
							   (I)3,(I)15,(I)14,
							   (I)3,(I)4,(I)15,
							   (I)4,(I)18,(I)15,
							   (I)4,(I)5,(I)18,
							   (I)5,(I)20,(I)18,
							   (I)5,(I)6,(I)20,
							   (I)6,(I)7,(I)20,
							   (I)6,(I)1,(I)7,
							   (I)14,(I)15,(I)13,
							   (I)15,(I)16,(I)13,
							   (I)15,(I)18,(I)16,
							   (I)18,(I)19,(I)16,
							   (I)18,(I)20,(I)19,
							   (I)20,(I)8,(I)19,
							   (I)20,(I)7,(I)8,
							   (I)13,(I)16,(I)12,
							   (I)16,(I)17,(I)12,
							   (I)16,(I)19,(I)17,
							   (I)19,(I)9,(I)17,
							   (I)19,(I)8,(I)9,							   
							   (I)12,(I)17,(I)11,
							   (I)17,(I)10,(I)11,
							   (I)17,(I)9,(I)10,
							   (I)11,(I)10,(I)2};
/*
  2
  13 12 
  14 21 11
  15 20 24 10
  16 19 23 26  9
  17 18 22 25 27 8
  0   3  4  5  6 7 1
*/
static const I __cnctreilli_TRIA_LAGRANGE_7[108] 		= {
  (I)0,(I)3,(I)17,
  (I)3,(I)18,(I)17,
  (I)3,(I)4,(I)18,
  (I)4,(I)22,(I)18,
  (I)4,(I)5,(I)22,
  (I)5,(I)25,(I)22,
  (I)5,(I)6,(I)25,
  (I)6,(I)27,(I)25,
  (I)6,(I)7,(I)27,
  (I)7,(I)8,(I)27,
  (I)7,(I)1,(I)8,

  (I)17,(I)18,(I)16,
  (I)18,(I)19,(I)16,
  (I)18,(I)22,(I)19,
  (I)22,(I)23,(I)19,
  (I)22,(I)25,(I)23,
  (I)25,(I)26,(I)23,
  (I)25,(I)27,(I)26,
  (I)27,(I)9,(I)26,
  (I)27,(I)8,(I)9,
    
  (I)16,(I)19,(I)15,
  (I)19,(I)20,(I)15,
  (I)19,(I)23,(I)20,
  (I)23,(I)24,(I)20,
  (I)23,(I)26,(I)24,
  (I)26,(I)10,(I)24,
  (I)26,(I)9,(I)10,

  (I)15,(I)20,(I)14,
  (I)20,(I)21,(I)14,
  (I)20,(I)24,(I)21,
  (I)24,(I)11,(I)21,
  (I)24,(I)10,(I)11,

  (I)14,(I)21,(I)13,
  (I)21,(I)12,(I)13,
  (I)21,(I)11,(I)12,

  (I)13,(I)12,(I)2};


static cst_pI __cnctreilli_TRIA[__ensBASIS_ALL] = {
  NULL,
  __cnctreilli_TRIA_LAGRANGE_1,
  __cnctreilli_TRIA_LAGRANGE_1,
  __cnctreilli_TRIA_LAGRANGE_2,
  __cnctreilli_TRIA_LAGRANGE_3,
  __cnctreilli_TRIA_LAGRANGE_1,
  __cnctreilli_TRIA_LAGRANGE_1,
  __cnctreilli_TRIA_LAGRANGE_2,
  __cnctreilli_TRIA_LAGRANGE_3,
  __cnctreilli_TRIA_LAGRANGEBUBBLE_0,
  __cnctreilli_TRIA_LAGRANGEBUBBLE_1,
  __cnctreilli_TRIA_LAGRANGEBUBBLE_2};


static cst_pI *  __cnctreilli[__eFace_ALL] = { NULL,
					       __cnctreilli_TRIA,
					       NULL,
					       NULL};

cst_pI 	ensBASIS_cnctreilli		(cst_ensBASIS shape_,cst_eFace elm_)
{
  cst_pI * __cnctreilli_elm = __cnctreilli[elm_];
  return __cnctreilli_elm[shape_];
}

static const I __nelmtreilli_TRIA[__ensBASIS_ALL] = {
  0,
  1,
  1,
  4,
  6,
  1,
  1,
  4,
  6,
  3,
  3,
  6};

static cst_pI __nelmtreilli[__eFace_ALL] = { NULL,
					     __nelmtreilli_TRIA,
					     NULL,
					     NULL};


I  		ensBASIS_nelmtreilli		(cst_ensBASIS shape_,cst_eFace elm_)
{
  return __nelmtreilli[elm_][shape_];
}



static const R __lc2gltreilli_TRIA_LAGRANGE_0[3] = {(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01};

static const R __lc2gltreilli_TRIA_LAGRANGE_1[3*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,
							     (R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,
							     (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00};

static const R __lc2gltreilli_TRIA_LAGRANGE_2[6*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)0.000000000000000e+00,(R)5.000000000000000e-01,
							       (R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)0.000000000000000e+00,
							       (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01};

static const R __lc2gltreilli_TRIA_LAGRANGE_3[10*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)6.666666666666666e-01,(R)3.333333333333333e-01,(R)0.0,(R)0.0,(R)3.333333333333333e-01,(R)6.666666666666666e-01,(R)3.333333333333333e-01,
								(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,(R)6.666666666666666e-01,(R)6.666666666666666e-01,(R)3.333333333333333e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
								(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,(R)6.666666666666666e-01,(R)6.666666666666666e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01};

static const R __lc2gltreilli_TRIA_LAGRANGE_4[15*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)2.500000000000000e-01,
								(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)7.500000000000000e-01,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)2.500000000000000e-01,(R)5.000000000000000e-01,
								(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)7.500000000000000e-01,(R)7.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01,(R)2.500000000000000e-01,(R)5.000000000000000e-01,(R)2.500000000000000e-01};

static const R __lc2gltreilli_TRIA_LAGRANGE_5[21*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)8.000000000000000e-01,(R)6.000000000000000e-01,(R)3.999999999999999e-01,(R)2.000000000000000e-01,(R)0.0,(R)0.0,(R)0.0,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)3.999999999999999e-01,(R)6.000000000000000e-01,(R)8.000000000000000e-01,(R)6.000000000000001e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)1.999999999999999e-01,
								(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)8.000000000000000e-01,(R)8.000000000000000e-01,(R)6.000000000000001e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,
								(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)8.000000000000000e-01,(R)8.000000000000000e-01,(R)6.000000000000001e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)6.000000000000001e-01,(R)2.000000000000000e-01,(R)4.000000000000000e-01,(R)2.000000000000000e-01};

static const R __lc2gltreilli_TRIA_LAGRANGE_6[28*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)0.0,(R)0.0,(R)0.000000000000000e+00,(R)0.0,(R)0.0,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)5.000000000000001e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,
								(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)8.333333333333333e-01,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,
								(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)8.333333333333333e-01,(R)8.333333333333333e-01,(R)6.666666666666666e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)6.666666666666666e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)5.000000000000000e-01,(R)1.666666666666666e-01,(R)3.333333333333333e-01,(R)1.666666666666666e-01};

static const R __lc2gltreilli_TRIA_LAGRANGEBUBBLE_0[2*3] = {(R)3.333333333333333e-01,(R)3.333333333333333e-01,
								   (R)3.333333333333333e-01,(R)3.333333333333333e-01,
								   (R)3.333333333333333e-01,(R)3.333333333333333e-01};


static const R __lc2gltreilli_TRIA_LAGRANGEBUBBLE_1[4*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
								   (R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
								   (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)3.333333333333333e-01};

static const R __lc2gltreilli_TRIA_LAGRANGEBUBBLE_2[7*3] = {(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)3.333333333333333e-01,
							       (R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)0.000000000000000e+00,(R)3.333333333333333e-01,
							       (R)0.000000000000000e+00,(R)0.000000000000000e+00,(R)1.000000000000000e+00,(R)0.000000000000000e+00,(R)5.000000000000000e-01,(R)5.000000000000000e-01,(R)3.333333333333333e-01};




static cst_pR __lc2gltreilli_TRIA[__ensBASIS_ALL] = {
  NULL,
  __lc2gltreilli_TRIA_LAGRANGE_0,
  __lc2gltreilli_TRIA_LAGRANGE_1,
  __lc2gltreilli_TRIA_LAGRANGE_2,
  __lc2gltreilli_TRIA_LAGRANGE_3,
  __lc2gltreilli_TRIA_LAGRANGE_0,
  __lc2gltreilli_TRIA_LAGRANGE_1,
  __lc2gltreilli_TRIA_LAGRANGE_2,
  __lc2gltreilli_TRIA_LAGRANGE_3,
  __lc2gltreilli_TRIA_LAGRANGEBUBBLE_0,
  __lc2gltreilli_TRIA_LAGRANGEBUBBLE_1,
  __lc2gltreilli_TRIA_LAGRANGEBUBBLE_2};


static cst_pR *  __lc2gltreilli[__eFace_ALL] = {NULL,
						__lc2gltreilli_TRIA,
						NULL,
						NULL};


cst_pR ensBASIS_lc2gltreilli_ptr	(cst_ensBASIS 	shape_,
					 cst_eFace 	elm_)
{
  return __lc2gltreilli[elm_][shape_];
}

void ensBASIS_lc2gltreilli	(cst_ensBASIS 	shape_,
				 cst_eFace 	elm_,
				 cst_pR 	cooelm_,
				 cst_pI 	cooelmoff_,
				 pR 		gl_,
				 cst_pI 	gloff_)
{
  const I nshape = ensBASIS_n(shape_,elm_);
  static const I negal2=2;
  static const R regal1=1.0;
  static const R regal0=0.0;
  static const I negal3=3;
  nsblas_dgemm("N",
	 "N",
	 &nshape,
	 &negal2,
	 &negal3,
	 &regal1,
	 __lc2gltreilli[elm_][shape_],
	 &nshape,  
	 cooelm_,
	 cooelmoff_,
	 &regal0,
	 gl_,
	 gloff_);
}



/* 
   sub-linear connectivity 

*/

/*
  2


  5     4

     6   

  0   3    1
*/
