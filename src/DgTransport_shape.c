#include "mkS_basis_tria.c"
#include "ns_constantes.h"


emkS_FAMILY 	mkS_family	(cst_mkS shape_) {return shape_->family;}
I		mkS_k		(cst_mkS shape_) {return shape_->degree;}
I		mkS_n		(cst_mkS shape_) {return shape_->nshape;}

emk_s_diff 	mkS_diff	(cst_mkS shape_) {return shape_->diff;}
void 		mkS_copy	(mkS 	 dest_,
				 cst_mkS src_)
{
  dest_->element	= src_->element;	
  dest_->family 	= src_->family;		
  dest_->diff 		= src_->diff;		
  dest_->degree 	= src_->degree;	
  dest_->nshape 	= src_->nshape;	
  //  dest_->continu  	= src_->continu;
}

cst_mkS mkS_b	(mkS shape_) { shape_->diff = __emk_s_diff_no; return shape_;}
cst_mkS mkS_dx	(mkS shape_) { shape_->diff = __emk_s_diff_dx; return shape_;}
cst_mkS mkS_dy	(mkS shape_) { shape_->diff = __emk_s_diff_dy; return shape_;}
cst_mkS mkS_dz	(mkS shape_) { shape_->diff = __emk_s_diff_dz; return shape_;}



#include "mkS_def.c"



#include "mkS_lagrange_localspl.c"

#include "mkS_eval_canonic_interval.c"
#include "mkS_eval_lagrange_interval.c"


#include "mkS_eval_canonic_tria.c"
#include "mkS_eval_lagrange_tria.c"

#include "mkS_l2ortho_28.c"
#include "mkS_l2ortho_21.c"
#include "mkS_l2ortho_15.c"
#include "mkS_l2ortho_10.c"
#include "mkS_l2ortho_6.c"
#include "mkS_l2ortho_3.c"
#include "mkS_l2ortho_1.c"
#include "mkS_eval_l2ortho_tria.c"



void mkS_basis(cst_mkS 		shape_,
	       cst_pI 		n_,
	       pR 		r_,
	       cst_pI 		roff_,
	       cst_pR		p_,
	       cst_pI 		poff_,
	       pR 		rwork_,
	       cst_pI 		rwork_n_,
	       pI  		ferr_)
{
  emkS_FAMILY family 	= mkS_family(shape_);
  emk_s_diff   diff   	= mkS_diff(shape_);
  const I degree 		= mkS_k(shape_);
  eTopology kindelm 		= shape_->element;

  if (__eTopology_TRIANGLE==kindelm)
    {
      mkS_basis_tria(shape_,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);
      return;
    }

  switch(family)
    {
    case __emkS_FAMILY_lagrangebubble:
      {
	ferr_[0]=(I)1;
	break;
      }
    case __emkS_FAMILY_lagrange:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
		    mkS_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    printf("bb\n");
		    break;
		  }
		case __eTopology_SEGMENT:
		  {		    
#if 0		    
		    mkS_lagrange_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_lagrange_interval(&degree,n_,r_,roff_,p_,&negal1,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_dx_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		  {
#if 0
		    mkS_dx_lagrange_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_dx_lagrange_interval(&degree,n_,r_,roff_,p_,&negal1,rwork_,rwork_n_,ferr_);  


		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dx_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dx_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_dy_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dy_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dy_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dz_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
		case __eTopology_SEGMENT:
		case __eTopology_QUADRILATERAL:
		case __eTopology_TRIANGLE:

		default:

		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }

	  default:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }
	  }
	break;
      }


    case __emkS_FAMILY_canonic:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		  {
#if 0
		    mkS_canonic_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_canonic_interval(&degree,n_,r_,roff_,p_,&negal1,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_dx_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		  {
#if 0
		    mkS_dx_canonic_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_dx_canonic_interval(&degree,n_,r_,roff_,p_,&negal1,rwork_,rwork_n_,ferr_);
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dx_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dx_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    //		    mkS_dy_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dy_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dy_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __eTopology_TETRAHEDRON:
		  {
		    //		    mkS_dz_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
		case __eTopology_QUADRILATERAL:
		case __eTopology_TRIANGLE:
		case __eTopology_SEGMENT:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }
#if NOT __mk_debug__
	  default:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }
#endif
	  }
	break;
      }


    case __emkS_FAMILY_l2orthonormal:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_SEGMENT:
		  {
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    break;
		  }
		case __eTopology_SEGMENT:
		  {
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dx_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __eTopology_QUADRILATERAL:
		  {
		    break;
		  }
		case __eTopology_TRIANGLE:
		  {
		    mkS_dy_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __eTopology_TETRAHEDRON:
		  {
		    break;
		  }
		case __eTopology_SEGMENT:
		case __eTopology_ALL:
		case __eTopology_ERROR:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __eTopology_TETRAHEDRON:
		  {
		    break;
		  }
		case __eTopology_ALL:
		case __eTopology_ERROR:
		case __eTopology_SEGMENT:
		case __eTopology_QUADRILATERAL:

		case __eTopology_TRIANGLE:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    fprintf(stderr,"mkS_basis:switch on kindelm failed\n");
		    exit(1);

		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }
#if NOT __mk_debug__
	  default:
	    {
	      ferr_[0]=(I)1;
	      break;
	    }
#endif
	  }
	break;
      }

    case __emkS_FAMILY_error:
    case __emkS_FAMILY_n:
      {
	ferr_[0]=(I)1;
	break;
      }
      

    default:
      {
	ferr_[0]=(I)1;
	break;
      }

    }
}



void mkS_kji(cst_mkS 	teta_,
	     cst_mkS  	trial_,
	     cst_mkS  	test_,
	     pR 		x_,
	     cst_pI		xoff_,
	     cst_pI 	qelm_n_,
	     cst_pR 	qelm_w_,
	     cst_pR	qelm_pos_,
	     cst_pI 	qelm_posoff_,
	     cst_pI 	rwork_n_,
	     pR 		rwork_,
	     pI 		err_)

{
  err_[0] 			= (I)0;
  const I trial_n 		= mkS_n(trial_);
  const I teta_n  		= mkS_n(teta_);
  const I test_n  		= mkS_n(test_);  
  const I mx_ttt_n          = MAX(MAX(teta_n,test_n),trial_n);
  const I rwork_size 	= qelm_n_[0]*(trial_n+teta_n+test_n) + mx_ttt_n*(mx_ttt_n+3);
  const I rwork_n  	 	= (rwork_n_[0] > rwork_size ) ? (rwork_n_[0] - qelm_n_[0]*(trial_n+teta_n+test_n)):(I)0;
  if (rwork_n<((I)1))
    {
      err_[0] = (I)1;
      
      fprintf(stderr,"not enough memory\n");
      return;
    }
  if (xoff_[0]<test_n*trial_n)
    {
      err_[0] = (I)1;
      fprintf(stderr,"mkS_kji wrong offset\n");  
      return;
    }
  R * trial_rwork 	= &rwork_[0];
  R * test_rwork  	= &rwork_[qelm_n_[0]*trial_n];
  R * teta_rwork  	= &rwork_[qelm_n_[0]*trial_n+qelm_n_[0]*test_n];
  R * rwork 	= &rwork_[qelm_n_[0]*(trial_n+teta_n+test_n)];
  I i,j;
  R w;

  mkS_basis(trial_,qelm_n_,trial_rwork,&trial_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);
  mkS_basis(test_,qelm_n_,test_rwork,&test_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);
  mkS_basis(teta_,qelm_n_,teta_rwork,&teta_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);          
  for (i=0;i<teta_n;++i)
    {
      for (j=0;j<test_n*trial_n;++j)   
	x_[i*xoff_[0]+j]=(R)0.0;   
      for (j=0;j<qelm_n_[0];++j)
	{
	  w = qelm_w_[j] * teta_rwork[teta_n*j+i];
	  dger((I*)&test_n,
		      &trial_n,
		      &w,
		      &test_rwork[j*test_n],
		      &negal1,
		      &trial_rwork[j*trial_n],
		      &negal1,
		      &x_[i*xoff_[0]],
		      &test_n);
	}
    }
  return;
}


void mkS_bmapping(const I 	n,
		  pR 		r,
		  cst_pI 	roff,
		  cst_pR 	p)
{
  I i;
  for (i=0;i<n;++i) r[i] 		= ( ((R)1.0) + p[i] ) *((R)0.5);
  for (i=0;i<n;++i) r[n+i] 		= ( ((R)1.0) - p[i] ) *((R)0.5);
  for (i=0;i<n;++i) r[2*n+i] 		= (R)0.0;
  for (i=0;i<n;++i) r[roff[0]+i] 	= (R)0.0;
  for (i=0;i<n;++i) r[roff[0]+n+i] 	= ( ((R)1.0) + p[i] ) *((R)0.5);
  for (i=0;i<n;++i) r[roff[0]+2*n+i] 	= ( ((R)1.0) - p[i] ) *((R)0.5);  
}







void mkS_bwji(cst_pI  	qface_n_,
	      cst_pR 	qface_p_,
	      cst_pR 	qface_w_,
	      cst_mkS	trial_,
	      cst_mkS 	test_,
	      pR 		x_,
	      const I 	xoff_,
	      pI 		ferr_)
{    
  I rwork_n_new,
    i,j,n,N;
  R rwork[2048];
  I rwork_n=2048;
#if __mk_debug__
  __mkS_debug(qface_n_,mkS_bwji);
  __mkS_debug(qface_p_,mkS_bwji);
  __mkS_debug(qface_w_,mkS_bwji);
  __mkS_debug(trial_,mkS_bwji);
  __mkS_debug(test_,mkS_bwji);
  __mkS_debug(x_,mkS_bwji);
#endif
  const I ntrial = mkS_n(trial_);
  const I ntest  = mkS_n(test_);
  N   		= 3*qface_n_[0];  
  n   		= 2*qface_n_[0]*3;
  
  mkS_bmapping(qface_n_[0],&rwork[0],&N,qface_p_); 
  
  rwork_n_new = rwork_n - (n + ntrial * N + ntest * N);
  mkS_basis(trial_,
	    &N,
	    &rwork[n],
	    &ntrial,
	    &rwork[0],
	    &N,
	    &rwork[n+ntrial*N + ntest * N],
	    &rwork_n_new,
	    ferr_); 
  mkS_basis(test_,
	    &N,
	    &rwork[n + ntrial * N ],
	    &ntest,
	    &rwork[0],
	    &N,
	    &rwork[n + ntrial * N + ntest * N],
	    &rwork_n_new,
	    ferr_);
  for (j=0;j<3*qface_n_[0];++j)
    for (i=0;i<ntest*ntrial;++i)	   
      x_[j*xoff_+i]=(R)0.0;

  for (j=0;j<3;++j)
    for (i=0;i<qface_n_[0];++i)
      dger((I*)&ntest,
		  &ntrial,
		  &qface_w_[i],
		  &rwork[n + ntrial * N + j*ntest*qface_n_[0]+ i*ntest],
		  &negal1,
		  &rwork[n+j*ntrial*qface_n_[0]+ i*ntrial],
		  &negal1,
		  &x_[( j * qface_n_[0] + i)*xoff_],
		  &ntest);      
}


void mkS_bwji_nei(cst_pI  	qface_n_,
		  cst_pR 	qface_p_,
		  cst_pR 	qface_w_,
		  cst_mkS trial_,
		  cst_mkS test_,
		  pR 	x_,
		  const I 	xoff_,
		  pI 	ferr_)
{

  I rwork_n_new,i,j,k,n,p,N;
  R rwork[2048];
  I rwork_n=2048;
#if __mk_debug__
  __mkS_debug(trial_,mkS_bwji_nei);
  __mkS_debug(test_,mkS_bwji_nei);
  __mkS_debug(qface_n_,mkS_bwji_nei);
  __mkS_debug(qface_p_,mkS_bwji_nei);
  __mkS_debug(qface_w_,mkS_bwji_nei);
  __mkS_debug(x_,mkS_bwji_nei);
#endif

  const I ntrial		= mkS_n(trial_);
  const I ntest 		= mkS_n(test_);
  N   		= 3*qface_n_[0];  
  n   		= 2*qface_n_[0]*3;  

  mkS_bmapping(qface_n_[0],&rwork[0],&N,qface_p_); 

  rwork_n_new = rwork_n - (n+ntrial*N+ntest*N);
  mkS_basis(trial_,
	    &N,
	    &rwork[n],
	    &ntrial,
	    &rwork[0],
	    &N,
	    &rwork[n+ntrial*N+ntest*N],
	    &rwork_n_new,
	    ferr_); 
  mkS_basis(test_,
	    &N,
	    &rwork[n+ntrial*N],
	    &ntest,
	    &rwork[0],
	    &N,
	    &rwork[n+ntrial*N+ntest*N],
	    &rwork_n_new,
	    ferr_); 
  for (k=0;k<3;++k)
    for (j=0;j<3;++j)
      for (i=0;i<qface_n_[0];++i)
	{
	  for (p=0;p<ntrial*ntest;++p)
	    x_[(3*qface_n_[0]*k+j*qface_n_[0]+i)*xoff_+p]=(R)0.0;
	  dger((I*)&ntest,
		      &ntrial,
		      &qface_w_[i],
		      &rwork[n+ntrial*N+k*ntest*qface_n_[0] + i*ntest],
		      &negal1,
		      &rwork[n+j*ntrial*qface_n_[0] + (qface_n_[0]-1-i)*ntrial],
		      &negal1,
		      &x_[(3*qface_n_[0]*k+j*qface_n_[0]+i)*xoff_],
		      &ntest);	      
	}
  
}
