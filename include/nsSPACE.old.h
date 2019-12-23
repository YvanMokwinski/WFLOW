#ifndef __header_nsSPACE_h__
#define __header_nsSPACE_h__

#include "ns_mesh.h"
#include "Sparse.h"
/**\brief structure to define finite element space */
typedef struct
{
  /**\brief only symbolic sparse matrix  */
  pSparse		endomorphism; /* uniquement le symbolic, chaque matrice construite pointe (en const) sur les tableaux d'entiers de l'endomorphism */
  /**\brief pointer to the mesh structure */
  cst_nsMESHZONE 	mesh;
  /**\brief number of dof */
  I 		nddl;
  /**\brief number of dof perelm */
  I 		nddlelm;
  /**\brief number of dof peredge */
  I 		nddledge;
  /**\brief kind shape */
  ensBASIS 		shape;
  /**\brief shift for element connectivity */
  I 		ddlcncoff;
  /**\brief dof coordinates */
  cst_pR		ddlcoo;
  /**\brief dof codes */
  cst_pI		ddlcod;
  /**\brief element connectivity */
  cst_pI		ddlcnc;
  /**\brief (modifiable) dof coordinates */
  pR			own_ddlcoo; 
  /**\brief (modifiable) dof codes */
  pI			own_ddlcod; 
  /**\brief (modifiable) element connectivity */
  pI			own_ddlcnc; 
  /**\brief pointer to exact local mass matrices  */
} nsSPACE_ST,*nsSPACE;

typedef const nsSPACE_ST cst_nsSPACE_ST;
typedef cst_nsSPACE_ST* cst_nsSPACE;

cst_nsSPACE	ns_mesh_get_space		(cst_nsMESHZONE 	mesh_,
						 const ensBASIS 	shape_);
cst_nsSPACE	ns_mesh_addspace		(nsMESHZONE 		mesh_,
						 const ensBASIS 	shape_);

void  		nsSPACE_get_ddlcoo		(cst_nsSPACE		s_,
						 cst_pI 		nddl_,
						 cst_pI 		ddl_,
						 pR 			pos_,
						 cst_pI 		posoff_);

I 		nsSPACE_iddl2ddlP6		(cst_nsSPACE 		SF_,
						 const I 		iddl_);

#ifndef NDEBUG

cst_pR 	nsSPACE_mass		(cst_nsSPACE const S_,const ensBASIS shape_);
cst_pR 	nsSPACE_convx		(cst_nsSPACE const S_,const ensBASIS shape_);
cst_pR 	nsSPACE_convy		(cst_nsSPACE const S_,const ensBASIS shape_);
const ns_mesh * nsSPACE_mesh		(cst_nsSPACE const S_);
I 		nsSPACE_nddl		(cst_nsSPACE const S_);
I 		nsSPACE_nddlelm		(cst_nsSPACE const S_);
I 		nsSPACE_nddledge	(cst_nsSPACE const S_);
ensBASIS 	nsSPACE_shape		(cst_nsSPACE const S_);
cst_pR 	nsSPACE_ddlcoo		(cst_nsSPACE const S_);
cst_pI 	nsSPACE_ddlcnc		(cst_nsSPACE const S_);
cst_pI 	nsSPACE_ddlcod		(cst_nsSPACE const S_);
cst_eElement 	nsSPACE_elm		(cst_nsSPACE const S_);
#else

#define 	nsSPACE_mass(_s,_shape)  	(_s)->MASS	[_shape]
#define 	nsSPACE_convx(_s,_shape)  	(_s)->CONVX	[_shape]
#define 	nsSPACE_convy(_s,_shape)  	(_s)->CONVY	[_shape]
#define 	nsSPACE_mesh(_s) 		(_s)->mesh
#define 	nsSPACE_elm(_s) 		ns_mesh_elm((_s)->mesh)
#define 	nsSPACE_nddl(_s) 		(_s)->nddl
#define 	nsSPACE_nddlelm(_s) 		(_s)->nddlelm
#define 	nsSPACE_nddledge(_s) 		(_s)->nddledge
#define 	nsSPACE_shape(_s) 		(_s)->shape
#define 	nsSPACE_ddlcoo(_s) 		(_s)->ddlcoo
#define 	nsSPACE_ddlcod(_s) 		(_s)->ddlcod
#define 	nsSPACE_ddlcnc(_s) 		(_s)->ddlcnc

#endif


void 		nsSPACE_free			(nsSPACE_ST*const		space_);
nsSPACE_ST*	nsSPACE_kill			(nsSPACE_ST*const		space_);

void 		nsSPACE_def 			(nsSPACE_ST*		space_,  
						 const ns_mesh * 	mesh_,
						 cst_ensBASIS 	shape_);

nsSPACE_ST* 	nsSPACE_new 			(const ns_mesh * 	mesh_,
						 cst_ensBASIS 	shape_);

pSparse 	nsSPACE_get_endomorphism		(cst_nsSPACE const S_);


void 		nsSPACE_leastsquare_galerkin_def(cst_nsSPACE		J_,
						 cst_nsSPACE		I_,
						 cst_pI 			level_,
						 cst_pI 			level_boundary_,
						 cst_pI			is_boundary_,
						 pSparse F_);

pSparse 	nsSPACE_galerkin_new		(cst_nsSPACE	S_);

void 		nsSPACE_get_ddlelm		(cst_nsSPACE  S_,cst_pI ielm_,cst_pR x_,pR ddlelm_);


void 		nsSPACE_codelm			(cst_nsSPACE	space_,
						 cst_pI 		ielm_,
						 pI   			codelm_);
void 		nsSPACE_cncelm_dec		(cst_nsSPACE	space_,
						 cst_pI 		ielm_,
						 pI   			ddlelm_,
						 const I 		dec_);
void 		nsSPACE_cncelm			(cst_nsSPACE	space_,
						 cst_pI 		ielm_,
						 pI   			ddlelm_);

void 		nsSPACE_write_medit		(cst_nsSPACE	s_,
						 const char * 		name_,
						 ...);

void 		nsSPACE_dg_write_medit		(cst_nsSPACE	s_,
						 const char * 		name_,
						 ...);

void 		nsSPACE_infopatch		(cst_nsSPACE	space_,
						 cst_pI 		ivertex_,
						 cst_pI 		level_,
						 pI 			patchvertex_n_,
						 pI 			patchvertex_,
						 pI 			patchelm_n_,
						 pI 			patchelm_,
						 pI 			patchddl_n_,
						 pI 			patchddl_,
						 pI 			blank_);


void nsSPACE_K_CONVECTION_matrix	(cst_nsSPACE 		space_,
					 cst_ensBASIS 		shape_velocity_,
					 cst_pR			belm_,
					 cst_pR 			xcoeff_,
					 cst_pR			coeff_,
					 cst_pI			coeffoff_,
					 pR 				matrix_,
					 cst_pI 			matrixoff_);
void nsSPACE_K_MASS_matrix		(cst_nsSPACE 		space_,
					 cst_ensBASIS         	shape_density_,
					 cst_pR			jacelm_,
					 cst_pR 			xcoeff_,
					 cst_pR			coeff_,
					 cst_pI			coeffoff_,
					 pR 				matrix_,
					 cst_pI 			matrixoff_);


#endif
