#ifndef __header_mkS_eval_h__
#define __header_mkS_eval_h__


void mkS_lagrange_interval(cst_pI 	degree,
			   cst_pI 	n,
			   pR 	r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   pR 	rwork,
			   cst_pI 	rwork_n,
			   pI		err_);

void mkS_dx_lagrange_interval(cst_pI 	degree,
			      cst_pI 	n,
			      pR		r,
			      cst_pI 	roff_,
			      cst_pR 	p,
			      cst_pI 	poff_,
			      pR		rwork,
			      cst_pI 	rwork_n,
			      pI 		err_);

void mkS_canonic_interval(cst_pI 	degree,
			  cst_pI 	n,
			  pR 	r,
			  cst_pI 	roff_,
			  cst_pR 	p,
			  cst_pI 	poff_,
			  pR 	rwork,
			  cst_pI 	rwork_n,
			  pI 		err_);

void mkS_dx_canonic_interval(cst_pI 	degree,
			     cst_pI 	n,
			     pR	r,
			     cst_pI 	roff_,
			     cst_pR 	p,
			     cst_pI 	poff_,
			     pR	rwork,
			     cst_pI 	rwork_n,
			     pI 	err_);


void mkS_canonic_tria(cst_pI 	degree,
		      cst_pI 	n,
		      pR 		r,
		      cst_pI 	roff_,
		      cst_pR 	p,
		      cst_pI 	poff_,
		      pR 		rwork,
		      cst_pI 	rwork_n,
		      pI 		err_);

void mkS_dx_canonic_tria(cst_pI 	degree,
			 cst_pI 	n,
			 pR	r,
			 cst_pI 	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 pR	rwork,
			 cst_pI 	rwork_n,
			 pI 	err_);


void mkS_dy_canonic_tria(cst_pI 	degree,
			 cst_pI 	n,
			 pR 		r,
			 cst_pI 	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 R*		rwork,
			 cst_pI 	rwork_n,
			 pI 		err_);



void mkS_canonic_quad(cst_pI 	degree,
		      cst_pI 	n,
		      pR 		r,
		      cst_pI 	roff_,
		      cst_pR 	p,
		      cst_pI 	poff_,
		      pR 		rwork,
		      cst_pI 	rwork_n,
		      pI 		err_);

void mkS_dx_canonic_quad(cst_pI 	degree,
			 cst_pI 	n,
			 pR	r,
			 cst_pI 	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 pR	rwork,
			 cst_pI 	rwork_n,
			 pI 	err_);


void mkS_dy_canonic_quad(cst_pI 	degree,
			 cst_pI 	n,
			 pR 		r,
			 cst_pI 	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 R*		rwork,
			 cst_pI 	rwork_n,
			 pI 		err_);



void mkS_canonic_tetra(cst_pI 	degree,
		       cst_pI 	n,
		       pR 		r,
		       cst_pI 	roff_,
		       cst_pR 	p,
		       cst_pI 	poff_,
		       pR 		rwork,
		       cst_pI 	rwork_n,
		       pI 		err_);

void mkS_dx_canonic_tetra(cst_pI 	degree,
			  cst_pI 	n,
			  pR	r,
			  cst_pI 	roff_,
			  cst_pR 	p,
			  cst_pI 	poff_,
			  pR	rwork,
			  cst_pI 	rwork_n,
			  pI 	err_);


void mkS_dy_canonic_tetra(cst_pI 	degree,
			  cst_pI 	n,
			  pR 		r,
			  cst_pI 	roff_,
			  cst_pR 	p,
			  cst_pI 	poff_,
			  R*		rwork,
			  cst_pI 	rwork_n,
			  pI 		err_);

void mkS_dz_canonic_tetra(cst_pI 	degree,
			  cst_pI 	n,
			  pR 		r,
			  cst_pI 	roff_,
			  cst_pR 	p,
			  cst_pI 	poff_,
			  R*		rwork,
			  cst_pI 	rwork_n,
			  pI 		err_);



void mkS_lagrange_tetra(cst_pI 	degree,
			cst_pI 	n,
			pR 		r,
			cst_pI 	roff_,
			cst_pR 	p,
			cst_pI 	poff_,
			pR 		rwork,
			cst_pI 	rwork_n,
			pI 		err_);

void mkS_dx_lagrange_tetra(cst_pI 	degree,
			   cst_pI 	n,
			   pR	r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   pR	rwork,
			   cst_pI 	rwork_n,
			   pI 	err_);


void mkS_dy_lagrange_tetra(cst_pI 	degree,
			   cst_pI 	n,
			   pR 		r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   R*		rwork,
			   cst_pI 	rwork_n,
			   pI 		err_);

void mkS_dz_lagrange_tetra(cst_pI 	degree,
			   cst_pI 	n,
			   pR 		r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   R*		rwork,
			   cst_pI 	rwork_n,
			   pI 		err_);


void mkS_lagrange_quad(cst_pI  	degree,
		       cst_pI 	n,
		       pR 	r,
		       cst_pI 	roff_,
		       cst_pR 	p,
		       cst_pI 	poff_,
		       pR 	rwork,
		       cst_pI 	rwork_n,
		       pI 	err_);


void mkS_dx_lagrange_quad(	cst_pI degree,
				cst_pI n,
				pR r,
				cst_pI roff_,
				cst_pR p,
				cst_pI  poff_,
				pR rwork,
				cst_pI rwork_n,
				pI err_);

void mkS_dy_lagrange_quad(cst_pI degree,
			  cst_pI n,
			  pR r,
			  cst_pI roff_,
			  cst_pR p,
			  cst_pI poff_,
			  pR rwork,
			  cst_pI rwork_n,
			  pI err_);


void mkS_lagrange_tria(cst_pI  	degree,
		       cst_pI 	n,
		       pR 	r,
		       cst_pI 	roff_,
		       cst_pR 	p,
		       cst_pI 	poff_,
		       pR 	rwork,
		       cst_pI 	rwork_n,
		       pI 	err_);


void mkS_dx_lagrange_tria(	cst_pI degree,
				cst_pI n,
				pR r,
				cst_pI roff_,
				cst_pR p,
				cst_pI  poff_,
				pR rwork,
				cst_pI rwork_n,
				pI err_);

void mkS_dy_lagrange_tria(cst_pI degree,
			  cst_pI n,
			  pR r,
			  cst_pI roff_,
			  cst_pR p,
			  cst_pI poff_,
			  pR rwork,
			  cst_pI rwork_n,
			  pI err_);


void mkS_l2ortho_tria(cst_pI  	degree,
		      cst_pI 	n,
		      pR 	r,
		      cst_pI 	roff_,
		      cst_pR 	p,
		      cst_pI 	poff_,
		      pR 	rwork,
		      cst_pI 	rwork_n,
		      pI 	err_);


void mkS_dx_l2ortho_tria(	cst_pI degree,
				cst_pI n,
				pR r,
				cst_pI roff_,
				cst_pR p,
				cst_pI  poff_,
				pR rwork,
				cst_pI rwork_n,
				pI err_);

void mkS_dy_l2ortho_tria(cst_pI degree,
			 cst_pI n,
			 pR r,
			 cst_pI roff_,
			 cst_pR p,
			 cst_pI poff_,
			 pR rwork,
			 cst_pI rwork_n,
			 pI err_);
#endif
