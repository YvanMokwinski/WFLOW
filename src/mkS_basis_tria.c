
#include "mkS.h"
#include "mkS_eval.h"

void mkS_basis_tria(cst_mkS 		shape_,
		    cst_pI 		n_,
		    pR 		r_,
		    cst_pI 		roff_,
		    cst_pR		p_,
		    cst_pI 		poff_,
		    pR			rwork_,
		    cst_pI 		rwork_n_,
		    pI			ferr_)
{
#if 0
#if __mk_debug__
  __mkS_debug(n_,mkS);
  __mkS_debug(r_,mkS);
  __mkS_debug(roff_,mkS);
  __mkS_debug(p_,mkS);
  __mkS_debug(poff_,mkS);
  __mkS_debug(rwork_,mkS);
  __mkS_debug(rwork_n_,mkS);
  __mkS_debug(poff_[0]>=n_[0],mkS);
#endif  
#endif

  emkS_FAMILY family 	= mkS_family(shape_);
  emk_s_diff   diff   	= mkS_diff(shape_);
  const I  degree	= mkS_k(shape_);
  switch(family)
    {
    case __emkS_FAMILY_lagrange:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      mkS_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      mkS_dx_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      mkS_dy_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
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


    case __emkS_FAMILY_canonic:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      mkS_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      mkS_dx_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      mkS_dy_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
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
	      mkS_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      mkS_dx_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      mkS_dy_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
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
    case __emkS_FAMILY_lagrangebubble:
      {
	switch(diff)
	  {
	  case __emk_s_diff_no:
	    {
	      { I i;
		for (i=0;i<n_[0];++i)
		  {
		    const R x   	= p_[i];
		    const R y   	= p_[poff_[0]+i];
		    const R lam 	= ((R)1.0)-(x+y);
		    r_[roff_[0]*i+0] 	= lam*(((R)2.0)*lam-((R)1.0));
		    r_[roff_[0]*i+1] 	= x*(((R)2.0)*x-((R)1.0));	
		    r_[roff_[0]*i+2] 	= y*(((R)2.0)*y-((R)1.0));	
		    r_[roff_[0]*i+3] 	= ((R)4.0)*x*lam;
		    r_[roff_[0]*i+4] 	= ((R)4.0)*x*y;	  
		    r_[roff_[0]*i+5] 	= ((R)4.0)*y*lam;
		    r_[roff_[0]*i+6] 	= ((R)0.5)*((R)54.0)*lam*x*y; 
		  } }
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      { I i;
		for (i=0;i<n_[0];++i)
		  {
		    const R x   	= p_[i];
		    const R y   	= p_[poff_[0]+i];
		    const R lam 	= ((R)1.0)-(x+y);
		    r_[roff_[0]*i+0] = ((R)1.0)-((R)4.0)*lam;
		    r_[roff_[0]*i+1] = ((R)-1.0)+((R)4.0)*x;
		    r_[roff_[0]*i+2] = (R)0.0;
		    r_[roff_[0]*i+3] = ((R)4.0)*(lam-x);
		    r_[roff_[0]*i+4] = ((R)4.0)*y;  
		    r_[roff_[0]*i+5] = ((R)-4.0)*y;
		    r_[roff_[0]*i+6]= ((R)0.5)*((R)54.0)*y*(lam-x);
		  } }
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      
	      { I i;
		for (i=0;i<n_[0];++i)
		  {
		    const R x   	= p_[i];
		    const R y   	= p_[poff_[0]+i];
		    const R lam 	= ((R)1.0)-(x+y);
		    r_[roff_[0]*i+0] = ((R)1.0)-((R)4.0)*lam;
		    r_[roff_[0]*i+1] = (R)0.0;
		    r_[roff_[0]*i+2] = ((R)-1.0)+((R)4.0)*y;
		    r_[roff_[0]*i+3] = ((R)-4.0)*x;
		    r_[roff_[0]*i+4] = ((R)4.0)*x;
		    r_[roff_[0]*i+5] = ((R)4.0)*(lam-y);
		    r_[roff_[0]*i+6] = ((R)0.5)*((R)54.0)*x*(lam-y);
		  } }
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
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
      
#if NOT __mk_debug__
    default:
      {
	ferr_[0]=(I)1;
	break;
      }
#endif
    }
}
