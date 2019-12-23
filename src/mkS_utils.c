#include "mkSYS.h"

#include "mkS.h"
#include "mkS_eval.h"
#include "mok_m.h"


#if 0
mok_m mkS_allbasis_new_mok_m(mkS 	shape_,
			     cst_mok_m 	q_pos_,
			     mkR 		rwork_,
			     const nsINT 	rwork_n)
{
  const nsINT shape_n 	= mkS_n(shape_);
  cst_emkDIM dim 	= mkS_dim(shape_);
  const nsINT q_n 	= mok_m_n(q_pos_);
  mok_m A 		= mok_m_new(shape_n,q_n*(1+dim));
  const nsINT q_pos_off = mok_m_off(q_pos_);
  nsINT ferr;
  mkS_basis(mkS_b(shape_),
	    &q_n, 
	    mok_m_x(A),
	    &A->off,
	    cst_mok_m_x(q_pos_),
	    &q_pos_off,
	    rwork_,
	    &rwork_n,
	    &ferr);
  switch(dim)
    {
    case __emkDIM_1:
      {
	mkS_basis(mkS_dx(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);
	break;
      }
      
    case __emkDIM_2:
      {
	mkS_basis(mkS_dx(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);
	mkS_basis(mkS_dy(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n*2),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);
	break;
      }
    case __emkDIM_3:
      {
	mkS_basis(mkS_dx(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);
	mkS_basis(mkS_dy(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n*2),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);
	mkS_basis(mkS_dz(shape_),
		  &q_n, 
		  mok_m_xcol(A,q_n*3),
		  &A->off,
		  cst_mok_m_x(q_pos_),
		  &q_pos_off,
		  rwork_,
		  &rwork_n,
		  &ferr);      
	break;
      }
    case __emkDIM_n:
    case __emkDIM_error:
    default:
      {
	ferr=1;
	mk_msg("mkS_allbasis_new_mok_m:switch failed");
	break;
      }
    }
  if (ferr)
    {
      return mok_m_kill(A);
    }
  return A;
}
#endif

emk_logic mkS_cmp(cst_mkS a_,
		  cst_mkS b_)
{
  emk_logic cmp_elm = (a_->element==b_->element)?__emk_yes:__emk_no;
  if (cmp_elm)
    {
      if ( (a_->degree==b_->degree)AND(a_->nshape==b_->nshape))
	{
	  return __emk_yes;
	}
      else return __emk_no;
    }
  else return __emk_no;

}

nsINT	mkS_nshapevertex(cst_mkS shape_)
{
  nsINT n=0;
  const nsINT degree 	= mkS_k(shape_);
  cst_emkZELM elm 	= mkS_emkZELM(shape_);
  switch(elm)
    {
    case __emkZELM_interval:
      {
	n =(degree==0)?0:2;
	break;
      }
    case __emkZELM_triangle:
      {
	n =(degree==0)?0:3;
	break;
      }
    case __emkZELM_quadrangle:
      {
	n =(degree==0)?0:4;
	break;
      }
    case __emkZELM_tetrahedra:
      {
	n =(degree==0)?0:4;
	break;
      }
    case __emkZELM_error:
    case __emkZELM_n:
      {
	mk_err("mkS_nshapevertex:switch failed on __emkZELM");
	break;
      }
    }
  return n;
}

nsINT 			mkS_nshapedge	(cst_mkS shape_)
{
  const nsINT degree 	= mkS_k(shape_);
  cst_emkZELM elm 	= mkS_emkZELM(shape_);
  nsINT n=0;
  switch(elm)
    {
    case __emkZELM_interval:
      {
	n=0;
	break;
      }
    case __emkZELM_triangle:
    case __emkZELM_quadrangle:
    case __emkZELM_tetrahedra:
      {
	n =(degree>1)?degree-1:0;
	break;
      }
    case __emkZELM_error:
    case __emkZELM_n:
      {

	mk_err("mkS_nshapedge:switch failed on __emkZELM");

	break;
      }
    }
  return n;
}


nsINT mkS_nshapelm	(cst_mkS shape_)
{
  const nsINT degree 	= mkS_k(shape_);
  cst_emkZELM elm 	= mkS_emkZELM(shape_);
  cst_emkS_FAMILY family = mkS_family(shape_);
  nsINT n=0;
  switch(elm)
    {
    case __emkZELM_interval:
      {
	n=degree-1;
	break;
      }
    case __emkZELM_quadrangle:
      {
	switch(family)
	  {
	  case __emkS_FAMILY_lagrangebubble:
	    {
	      n=(degree==1)?1:((degree==2)?1:0);
	      break;
	    }
	  case __emkS_FAMILY_canonic:
	  case __emkS_FAMILY_l2orthonormal:
	  case __emkS_FAMILY_lagrange:
	    {
	      n=(degree>1)?((degree-1)*(degree-1)):0;
	      break;
	    }
	  case __emkS_FAMILY_n:
	  case __emkS_FAMILY_error:
	    {
	      mk_err("mkS_nshapelm:switch failed on __emkZELM");
	      break;
	    }
	  }	
	break;
      }
    case __emkZELM_triangle:
    case __emkZELM_tetrahedra:
      {
	switch(family)
	  {
	  case __emkS_FAMILY_lagrangebubble:
	    {
	      n=(degree==1)?1:((degree==2)?1:0);
	      break;
	    }
	  case __emkS_FAMILY_canonic:
	  case __emkS_FAMILY_l2orthonormal:
	  case __emkS_FAMILY_lagrange:
	    {
	      n=(degree>2)?((degree-1)*(degree-2))/2:0;
	      break;
	    }
	  case __emkS_FAMILY_n:
	  case __emkS_FAMILY_error:
	    {
	      mk_err("mkS_nshapelm:switch failed on __emkZELM");
	      break;
	    }
	  }	
	break;
      }
    case __emkZELM_error:
    case __emkZELM_n:
      {
	mk_err("mkS_nshapelm:switch failed on __emkZELM");
	break;
      }
    }
  return n;
}

#if 0

emkDIM           mkS_dim		(cst_mkS shape_){return mkZELM_dim(&shape_->element);}
#endif

#if __mk_debug__

emkZELM        mkS_emkZELM	(cst_mkS shape_){return shape_->element;}

emkS_FAMILY 	mkS_family	(cst_mkS shape_) {return shape_->family;}
nsINT 			mkS_k		(cst_mkS shape_) {return shape_->degree;}
nsINT 			mkS_n		(cst_mkS shape_) {return shape_->nshape;}
emkCONTINUITY 	mkS_continuity	(cst_mkS shape_) {return shape_->continu;}
emk_s_diff 	mkS_diff		(cst_mkS shape_) {return shape_->diff;}
void 			mkS_copy		(mkS 	 dest_,
						 cst_mkS src_)
{
  dest_->element	= src_->element;	
  dest_->family 	= src_->family;		
  dest_->diff 		= src_->diff;		
  dest_->degree 	= src_->degree;	
  dest_->nshape 	= src_->nshape;	
  dest_->continu  	= src_->continu;
}

cst_mkS mkS_b	(mkS shape_) { shape_->diff = __emk_s_diff_no; return shape_;}
cst_mkS mkS_dx	(mkS shape_) { shape_->diff = __emk_s_diff_dx; return shape_;}
cst_mkS mkS_dy	(mkS shape_) { shape_->diff = __emk_s_diff_dy; return shape_;}
cst_mkS mkS_dz	(mkS shape_) { shape_->diff = __emk_s_diff_dz; return shape_;}
#endif
cst_mkS mkS_continuous(mkS shape_)    { shape_->continu = __emk_continuous;    return shape_;}
cst_mkS mkS_discontinuous(mkS shape_) { shape_->continu = __emk_discontinuous; return shape_;}
mkS mkS_kill(mkS A_)
{
  if (A_)
    {
      free(A_);
    }
  return NULL;
}
emk_err mkS_fwrite(FILE * fich_,cst_mkS shape_)
{
  size_t rc;
  emk_err err = __emkZELM_fwrite(fich_,shape_->element);
  if (err)
    {
      mk_msg("mkIerp_fwrite:__emkZELM_fwrite failed");
      return err;
    }  
#if 0
  rc = fwrite(&shape_->shape, sizeof(emk_s), (size_t)1, fich_);
#endif
  err = __emk_s_fwrite(fich_,shape_->shape);
  if (err)
    {
      mk_msg("mkIerp_fwrite:__emkIerp_fwrite failed");
      return err;
    }  
  rc = fwrite(&shape_->family, sizeof(emkS_FAMILY), (size_t)1, fich_);
  if (err)
    {
      mk_msg("mkIerp_fwrite:__emkIerp_fwrite failed");
      return err;
    }  

  rc = fwrite(&shape_->continu, sizeof(emkCONTINUITY), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mkIerp_fwrite:can not write continuity enumeration "nsPRINT_SIZET"",rc);
      return err;
    }
  rc = mok_fwrite_integer(&shape_->degree, sizeof(nsINT), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mpsolzone_fwrite:can not read version "nsPRINT_SIZET"",rc);
      return err;
    }
  rc = mok_fwrite_integer(&shape_->nshape, sizeof(nsINT), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mpsolzone_fwrite:can not read version "nsPRINT_SIZET"",rc);
      return err;
    }
  return err;
}


emk_err mkS_fread(FILE * fich_,mkS shape_)
{
  size_t rc;
  emk_err err = __emkZELM_fread(fich_,&shape_->element);
  if (err)
    {
      mk_msg("mkIerp_fwrite:mpselm_fwrite failed");
      return err;
    }  
  rc = fread(&shape_->shape, sizeof(emk_s), (size_t)1, fich_);
  if (err)
    {
      mk_msg("mkIerp_fwrite:__emkIerp_fwrite failed");
      return err;
    }  
  rc = fread(&shape_->family, sizeof(emkS_FAMILY), (size_t)1, fich_);
  if (err)
    {
      mk_msg("mkIerp_fwrite:__emkIerp_fwrite failed");
      return err;
    }  

  rc = fread(&shape_->continu, sizeof(emkCONTINUITY), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mkIerp_fwrite:can not write continuity enumeration "nsPRINT_SIZET"",rc);
      return err;
    }
  rc = fread(&shape_->degree, sizeof(nsINT), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mpsolzone_fwrite:can not read version "nsPRINT_SIZET"",rc);
      return err;
    }
  rc = fread(&shape_->nshape, sizeof(nsINT), (size_t)1, fich_);
  if (rc != ((size_t)1) ) 
    {
      err = __emk_err_file;
      mk_msg("mpsolzone_fwrite:can not read version "nsPRINT_SIZET"",rc);
      return err;
    }
  return err;
}


mkS mkS_dup(cst_mkS A_)
{
#if __mk_debug__
  __mkS_debug(A_,mkS_dup);
#endif
  mkS s = calloc(1,sizeof(mkS_st));
  memcpy(s,A_,sizeof(mkS_st));
  return s;
}

emk_s mkS_enum			(cst_mkS S_)
{
  return S_->shape;
}


#if 0
mkS mkS_newinit(cst_emkS_FAMILY    	family_,
		cst_emkZELM    		element_,
		const nsINT		degree_)
{
  emk_err err;
  mkS s = calloc(1,sizeof(mkS_st));
#if 0
  if (NOT s)
    {
      mk_msg("mkS_newinit:calloc failed");
    }
#endif
  mkS_def(s,element_,family_,degree_,__emk_discontinuous,&err);
#if 0
  if (err)
    {
      mk_msg("mkS_newinit:mkS_def failed");
    }
#endif
  return s;
}
#endif




void mkS_basis(cst_mkS 		shape_,
	       cst_mkI 		n_,
	       mkR 		r_,
	       cst_mkI 		roff_,
	       cst_mkR		p_,
	       cst_mkI 		poff_,
	       mkR 		rwork_,
	       cst_mkI 		rwork_n_,
	       mkI  		ferr_)
{

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
  cst_emkS_FAMILY family 	= mkS_family(shape_);
  cst_emk_s_diff   diff   	= mkS_diff(shape_);
  const nsINT degree 		= mkS_k(shape_);
  cst_emkZELM kindelm 		= mkS_emkZELM(shape_);

  if (__emkZELM_triangle==kindelm)
    {
      mkS_basis_tria(shape_,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);
      return;
    }

  switch(family)
    {
    case __emkS_FAMILY_lagrangebubble:
      {
	ferr_[0]=(nsINT)1;
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
		case __emkZELM_quadrangle:
		  {
		    mkS_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
		    mkS_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    printf("bb\n");
		    break;
		  }
		case __emkZELM_interval:
		  {		    
#if 0		    
		    mkS_lagrange_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_lagrange_interval(&degree,n_,r_,roff_,p_,&__vmps_blas_negal1,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    mkS_dx_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		  {
#if 0
		    mkS_dx_lagrange_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_dx_lagrange_interval(&degree,n_,r_,roff_,p_,&__vmps_blas_negal1,rwork_,rwork_n_,ferr_);  


		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dx_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_dx_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    mkS_dy_lagrange_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dy_lagrange_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_dy_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __emkZELM_tetrahedra:
		  {
		    mkS_dz_lagrange_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
		case __emkZELM_interval:
		case __emkZELM_quadrangle:
		case __emkZELM_triangle:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(nsINT)1;
	      break;
	    }
#if NOT __mk_debug__
	  default:
	    {
	      ferr_[0]=(nsINT)1;
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
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    mkS_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		  {
#if 0
		    mkS_canonic_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_canonic_interval(&degree,n_,r_,roff_,p_,&__vmps_blas_negal1,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    mkS_dx_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		  {
#if 0
		    mkS_dx_canonic_interval(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
#endif
		    mkS_dx_canonic_interval(&degree,n_,r_,roff_,p_,&__vmps_blas_negal1,rwork_,rwork_n_,ferr_);
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dx_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_dx_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    mkS_dy_canonic_quad(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dy_canonic_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    mkS_dy_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __emkZELM_tetrahedra:
		  {
		    mkS_dz_canonic_tetra(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
		case __emkZELM_quadrangle:
		case __emkZELM_triangle:
		case __emkZELM_interval:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(nsINT)1;
	      break;
	    }
#if NOT __mk_debug__
	  default:
	    {
	      ferr_[0]=(nsINT)1;
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
		case __emkZELM_quadrangle:
		  {
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_interval:
		  {
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dx:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    break;
		  }
		case __emkZELM_interval:
		  {
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dx_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dy:
	    {
	      switch(kindelm)
		{
		case __emkZELM_quadrangle:
		  {
		    break;
		  }
		case __emkZELM_triangle:
		  {
		    mkS_dy_l2ortho_tria(&degree,n_,r_,roff_,p_,poff_,rwork_,rwork_n_,ferr_);  
		    break;
		  }
		case __emkZELM_tetrahedra:
		  {
		    break;
		  }
		case __emkZELM_interval:
		case __emkZELM_n:
		case __emkZELM_error:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_dz:
	    {
	      switch(kindelm)
		{
		case __emkZELM_tetrahedra:
		  {
		    break;
		  }
		case __emkZELM_n:
		case __emkZELM_error:
		case __emkZELM_interval:
		case __emkZELM_quadrangle:
		case __emkZELM_triangle:
#if NOT __mk_debug__
		default:
#endif
		  {
		    ferr_[0] = 1;
		    mk_msg("mkS_basis:switch on kindelm failed");
		    break;		  
		  }
		}
	      break;
	    }
	  case __emk_s_diff_n:
	    {
	      ferr_[0]=(nsINT)1;
	      break;
	    }
#if NOT __mk_debug__
	  default:
	    {
	      ferr_[0]=(nsINT)1;
	      break;
	    }
#endif
	  }
	break;
      }

    case __emkS_FAMILY_error:
    case __emkS_FAMILY_n:
      {
	ferr_[0]=(nsINT)1;
	break;
      }
      
#if NOT __mk_debug__
    default:
      {
	ferr_[0]=(nsINT)1;
	break;
      }
#endif
    }
}








mok_m mkS_shape_new_mok_m(cst_mkS shape_,
			  cst_mok_m lc_,
			  const char * Tr)
{
  nsINT 
    rwork_n,
    n   	  	= mok_m_n(lc_),
    noff   	  	= mok_m_off(lc_),
    nshape 	  	= mkS_n(shape_);  
  const nsREAL *t 	= cst_mok_m_x(lc_);
  nsREAL * rwork 	= NULL;
  mok_m m 		= mok_m_new(nshape,n);
  nsINT dim  = 2;
  rwork_n 		= MAX ( (nshape*(nshape+dim)*10),nshape*n);
  rwork 		= calloc(rwork_n,sizeof(nsREAL));
#if 0
  const nsINT moff 	= mok_m_off(m);
#endif
  if (NOT rwork)
    {
      mk_msg("mkIerp_shape_new_mpsmat:calloc failed");
      return mok_m_kill(m);
    }  
  nsINT ferr=1;
  mkS_basis(shape_,
	    &n,
	    mok_m_x(m),
	    &nshape,
	    t,
	    &noff,
	    rwork,
	    &rwork_n,
	    &ferr);    
  if ( (Tr[0]=='T') OR (Tr[0]=='t') )
    {
      nsINT i,j;
      for (j=0;j<n;++j)
	for (i=0;i<nshape;++i)
	  rwork[n*i+j] = m->x[m->off*j+i];
      for (i=0;i<nshape;++i)
	for (j=0;j<n;++j)
	  m->x[n*i+j] = rwork[n*i+j];
      m->n = n;
      m->m = nshape;
      m->off=n;
    }
  free(rwork);
  return m;

}

void   mkS_shape_def_mok_m(cst_mkS 	S_,
			   cst_mok_m 	lc_,
			   mok_m 		m_,
			   const char * 	Tr,
			   nsREAL*		rwork,
			   const nsINT * 	rwork_n,
			   emk_err * 	err_)
{
  const nsINT offm  	= mok_m_off(m_); 
  const nsINT offs 	= mok_m_off(lc_);
  const nsINT n   	= mok_m_n(lc_); 
  err_[0] 	= __emk_err_no;
  nsINT ferr=0;
  mkS_basis(S_,
	    &n,
	    mok_m_x(m_),
	    &offm,
	    cst_mok_m_x(lc_),
	    &offs,
	    rwork,
	    rwork_n,
	    &ferr);    
#if 0
  const nsINT nshape 	= mkS_n(S_);
  if (rwork_n[0]<nshape*(nshape+mkIerp_dim(S_)))
    {      
      mps_errmsg("mkIerp_shape_def_mpsmat failed rwork_n="ifmt" should be greater than\n",rwork_n[0]);
      err_[0] 	= __emk_err_internal;
      return;
    }
  mkIerp_shape(interpolant,
	       trans,
	       &n,
	       mpsmat_x(m_),
	       &offm,
	       cst_mpsmat_x(lc_),
	       &offs,	   			   
	       rwork,
	       rwork_n,
	       err_);   
  if (err_[0])
    {
      mps_errmsg("mkIerp_shape_def_mpsmat:mkIerp_shape failed");
      err_[0] 	= __emk_err_internal;
      return;
    }
#endif
}





void mkS_bmapping(const nsINT 	n,
		  mkR 		r,
		  cst_mkI 	roff,
		  cst_mkR 	p)
{
  nsINT i;
  for (i=0;i<n;++i) r[i] 		= ( ((nsREAL)1.0) + p[i] ) *((nsREAL)0.5);
  for (i=0;i<n;++i) r[n+i] 		= ( ((nsREAL)1.0) - p[i] ) *((nsREAL)0.5);
  for (i=0;i<n;++i) r[2*n+i] 		= (nsREAL)0.0;
  for (i=0;i<n;++i) r[roff[0]+i] 	= (nsREAL)0.0;
  for (i=0;i<n;++i) r[roff[0]+n+i] 	= ( ((nsREAL)1.0) + p[i] ) *((nsREAL)0.5);
  for (i=0;i<n;++i) r[roff[0]+2*n+i] 	= ( ((nsREAL)1.0) - p[i] ) *((nsREAL)0.5);  
}



void mkS_kji(cst_mkS 	teta_,
	     cst_mkS  	trial_,
	     cst_mkS  	test_,
	     mkR 		x_,
	     cst_mkI		xoff_,
	     cst_mkI 	qelm_n_,
	     cst_mkR 	qelm_w_,
	     cst_mkR	qelm_pos_,
	     cst_mkI 	qelm_posoff_,
	     cst_mkI 	rwork_n_,
	     mkR 		rwork_,
	     mkI 		err_)

{
  err_[0] 			= (nsINT)0;
  const nsINT trial_n 		= mkS_n(trial_);
  const nsINT teta_n  		= mkS_n(teta_);
  const nsINT test_n  		= mkS_n(test_);  
  const nsINT mx_ttt_n          = MAX(MAX(teta_n,test_n),trial_n);
  const nsINT rwork_size 	= qelm_n_[0]*(trial_n+teta_n+test_n) + mx_ttt_n*(mx_ttt_n+3);
  const nsINT rwork_n  	 	= (rwork_n_[0] > rwork_size ) ? (rwork_n_[0] - qelm_n_[0]*(trial_n+teta_n+test_n)):(nsINT)0;
  if (rwork_n<((nsINT)1))
    {
      err_[0] = (nsINT)1;
      
      fprintf(stderr,"not enough memory\n");
      return;
    }
  if (xoff_[0]<test_n*trial_n)
    {
      err_[0] = (nsINT)1;
      fprintf(stderr,"mkS_kji wrong offset\n");  
      return;
    }
  nsREAL * trial_rwork 	= &rwork_[0];
  nsREAL * test_rwork  	= &rwork_[qelm_n_[0]*trial_n];
  nsREAL * teta_rwork  	= &rwork_[qelm_n_[0]*trial_n+qelm_n_[0]*test_n];
  nsREAL * rwork 	= &rwork_[qelm_n_[0]*(trial_n+teta_n+test_n)];
  nsINT i,j;
  nsREAL w;

  mkS_basis(trial_,qelm_n_,trial_rwork,&trial_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);
  mkS_basis(test_,qelm_n_,test_rwork,&test_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);
  mkS_basis(teta_,qelm_n_,teta_rwork,&teta_n,qelm_pos_,qelm_posoff_,rwork,&rwork_n,err_);          
  for (i=0;i<teta_n;++i)
    {
      for (j=0;j<test_n*trial_n;++j)   
	x_[i*xoff_[0]+j]=(nsREAL)0.0;   
      for (j=0;j<qelm_n_[0];++j)
	{
	  w = qelm_w_[j] * teta_rwork[teta_n*j+i];
	  nsblas_dger((nsINT*)&test_n,
		      &trial_n,
		      &w,
		      &test_rwork[j*test_n],
		      &__vmps_blas_negal1,
		      &trial_rwork[j*trial_n],
		      &__vmps_blas_negal1,
		      &x_[i*xoff_[0]],
		      &test_n);
	}
    }
  return;
}
#if 0
void mkS_ji( cst_mkS  	trial_,
	     cst_mkS  	test_,
	     mkR 		x_,
	     cst_mkI		xoff_,
	     cst_mkI 	qelm_n_,
	     cst_mkR 	qelm_w_,
	     cst_mkR	qelm_pos_,
	     cst_mkI 	qelm_posoff_,
	     cst_mkI 	rwork_n_,
	     mkR 		rwork_,
	     emk_err* 		err_)
{
  const nsINT ntrial		= mkS_n(trial_);
  const nsINT ntest		= mkS_n(test_);
  const nsINT ntestxntrial	= ntest*ntrial;
  nsINT rwork_n			= (ntrial+ntest)*qelm_n_[0]+MAX(ntrial,ntest)*(MAX(ntrial,ntest)+2);
  nsINT rwork_n_new 		= (rwork_n_[0]>rwork_n)?rwork_n_[0]-rwork_n:(nsINT)0;
  err_[0]=__emk_err_no;
  if (NOT rwork_n_new)
    {
      err_[0]=__emk_err_memory;
      mk_msg("mkS_ji:rwork_n_ = "ifmt" <(ntrial+ntest)*qelm_n_[0]+MAX(ntrial,ntest)*(MAX(ntrial,ntest)+2) = "ifmt"",rwork_n_[0],rwork_n);
      return;
    }
  
  rwork_n_new = rwork_n_[0] - (ntrial * q_n + ntest * q_n);
  mkS_basis(trial_,
	    qelm_n_,
	    &rwork[0],
	    &ntrial,
	    qelm_p_,
	    qelm_poff_,
	    &rwork[ntrial*qelm_n_[0] + ntest * qelm_n_[0]],
	    &rwork_n_new,
	    err_);
  if (err_[0])
    {
      mk_msg("mkS_ji:mkS_basis failed");
      return;
    }
  
  mkS_basis(test_,
	    qelm_n_,
	    &rwork[ntrial * qelm_n_[0]],
	    &ntest,
	    qelm_p_,
	    qelm_poff_,
	    &rwork[ntrial * qelm_n_[0] + ntest * qelm_n_[0]],
	    &rwork_n_new,
	    err_);
  if (err_[0])
    {
      mk_msg("mkS_ji:mkS_basis failed");
      return;
    }
  { nsINT j;
    for (j=0;j<ntestxntrial;++j)	   
      {
	x_[j]=(nsREAL)0.0;
      } }  
  { nsINT i;
    for (i=0;i<qelm_n_[0];++i)
      {
	nsblas_dger(&ntest,
		    &ntrial,
		    &qelm_w_[i],
		    &rwork[ntrial * qelm_n_[0] + ntest*i],
		    &negal1,
		    &rwork[ntrial*i],
		    &negal1,
		    &x_[0],
		    &ntest);   
      } }     
}
#endif

#if 0
void mkIerp_wji(cst_mpsquadrature q_,
		cst_mkIerp trial_,
		cst_mkIerp test_,
		nsREAL*x_,
		const nsINT xoff_,
		enum __emps_error*err_)
{
  const nsREAL* w   	= cst_mpsquadrature_w(q_);  
  nsINT q_n 	= mpsquadrature_n(q_),ntrial=mkIerp_n(trial_),ntest=mkIerp_n(test_);
  nsINT i,j,negal1=(nsINT)1,ntestxntrial=ntest*ntrial,rwork_n_new;
  nsINT rwork_n=(ntrial+ntest)*q_n+MAX(ntrial,ntest)*(MAX(ntrial,ntest)+2);
  nsREAL * rwork = malloc(sizeof(nsREAL)*rwork_n);
  err_[0]=__emps_error_no;
  if (NOT rwork)
    {
      err_[0]=__emps_error_memory;
      mps_errmsg("mkIerp_wji:malloc failed");
      return;
    }
  rwork_n_new = rwork_n - (ntrial * q_n + ntest * q_n);
  (void)mkIerp_shape(trial_,
		     "No transpose",
		     &q_n,
		     &rwork[0],
		     &ntrial,
		     cst_mpsquadrature_p(q_),
		     &q_n,
		     &rwork[ntrial*q_n + ntest * q_n],
		     &rwork_n_new,
		     err_); __mps_checkerror("mkIerp_wji:mkIerp_shape",err_);
      
  (void)mkIerp_shape(test_,
		     "No transpose",
		     &q_n,
		     &rwork[ntrial * q_n ],
		     &ntest,
		     cst_mpsquadrature_p(q_),
		     &q_n,
		     &rwork[ntrial * q_n + ntest * q_n],
		     &rwork_n_new,
		     err_); __mps_checkerror("mkIerp_wji:mkIerp_shape",err_);
      
  for (i=0;i<q_n;++i)
    for (j=0;j<ntestxntrial;++j)	   
      x_[i*xoff_+j]=(nsREAL)0.0;
      
  for (i=0;i<q_n;++i)
    (void)nsblas_dger(&ntest,
		      &ntrial,
		      (nsREAL*)&w[i],
		      &rwork[ ntrial * q_n + ntest*i],
		      &negal1,
		      &rwork[ ntrial*i],
		      &negal1,
		      &x_[i*xoff_],
		      &ntest);      

  free(rwork);
}
#endif
void mkS_bwji(cst_mkI  	qface_n_,
	      cst_mkR 	qface_p_,
	      cst_mkR 	qface_w_,
	      cst_mkS	trial_,
	      cst_mkS 	test_,
	      mkR 		x_,
	      const nsINT 	xoff_,
	      mkI 		ferr_)
{    
  nsINT rwork_n_new,
    i,j,n,N;
  nsREAL rwork[2048];
  nsINT rwork_n=2048;
#if __mk_debug__
  __mkS_debug(qface_n_,mkS_bwji);
  __mkS_debug(qface_p_,mkS_bwji);
  __mkS_debug(qface_w_,mkS_bwji);
  __mkS_debug(trial_,mkS_bwji);
  __mkS_debug(test_,mkS_bwji);
  __mkS_debug(x_,mkS_bwji);
#endif
  const nsINT ntrial = mkS_n(trial_);
  const nsINT ntest  = mkS_n(test_);
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
      x_[j*xoff_+i]=(nsREAL)0.0;

  for (j=0;j<3;++j)
    for (i=0;i<qface_n_[0];++i)
      nsblas_dger((nsINT*)&ntest,
		  &ntrial,
		  &qface_w_[i],
		  &rwork[n + ntrial * N + j*ntest*qface_n_[0]+ i*ntest],
		  &__vmps_blas_negal1,
		  &rwork[n+j*ntrial*qface_n_[0]+ i*ntrial],
		  &__vmps_blas_negal1,
		  &x_[( j * qface_n_[0] + i)*xoff_],
		  &ntest);      
}


void mkS_bwji_nei(cst_mkI  	qface_n_,
		  cst_mkR 	qface_p_,
		  cst_mkR 	qface_w_,
		  cst_mkS trial_,
		  cst_mkS test_,
		  mkR 	x_,
		  const nsINT 	xoff_,
		  mkI 	ferr_)
{

  nsINT rwork_n_new,i,j,k,n,p,N;
  nsREAL rwork[2048];
  nsINT rwork_n=2048;
#if __mk_debug__
  __mkS_debug(trial_,mkS_bwji_nei);
  __mkS_debug(test_,mkS_bwji_nei);
  __mkS_debug(qface_n_,mkS_bwji_nei);
  __mkS_debug(qface_p_,mkS_bwji_nei);
  __mkS_debug(qface_w_,mkS_bwji_nei);
  __mkS_debug(x_,mkS_bwji_nei);
#endif

  const nsINT ntrial		= mkS_n(trial_);
  const nsINT ntest 		= mkS_n(test_);
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
	    x_[(3*qface_n_[0]*k+j*qface_n_[0]+i)*xoff_+p]=(nsREAL)0.0;
	  nsblas_dger((nsINT*)&ntest,
		      &ntrial,
		      &qface_w_[i],
		      &rwork[n+ntrial*N+k*ntest*qface_n_[0] + i*ntest],
		      &__vmps_blas_negal1,
		      &rwork[n+j*ntrial*qface_n_[0] + (qface_n_[0]-1-i)*ntrial],
		      &__vmps_blas_negal1,
		      &x_[(3*qface_n_[0]*k+j*qface_n_[0]+i)*xoff_],
		      &ntest);	      
	}
  
}

void mkS_spl_lc2gl_(cst_mkI 	degree_,
		    mkR 	L_,
		    cst_mkI 	Loff_)
{
  nsINT i,j,k;
  if (degree_[0]==0)
    {
      L_[0] = (nsREAL)1.0/(nsREAL)3.0;
      L_[1] = (nsREAL)1.0/(nsREAL)3.0;
      L_[2] = (nsREAL)1.0/(nsREAL)3.0;
    }
  else
    {      
      L_[0]	= (nsREAL)1.0;      
      L_[1]  	= (nsREAL)0.0;
      L_[2]	= (nsREAL)0.0;      
      
      L_[Loff_[0]+0]	= (nsREAL)0.0;      
      L_[Loff_[0]+1]  	= (nsREAL)1.0;
      L_[Loff_[0]+2]	= (nsREAL)0.0;      

      L_[2*Loff_[0]+0]	= (nsREAL)0.0;      
      L_[2*Loff_[0]+1] = (nsREAL)0.0;
      L_[2*Loff_[0]+2]	= (nsREAL)1.0;      

      k=3;	
      for (j=1;j<degree_[0];++j)
	{
	  L_[k*Loff_[0]+0]   	= ((nsREAL)1.0)-((nsREAL)j)/((nsREAL)degree_[0]);
	  L_[k*Loff_[0]+1]   	= ((nsREAL)j)/((nsREAL)degree_[0]);
	  L_[k*Loff_[0]+2] 	= ((nsREAL)0.0);
	  k = k+1;
	}	
      for (j=1;j<degree_[0];++j)
	{
	  L_[k*Loff_[0]+0] 	= ((nsREAL)0.0);
	  L_[k*Loff_[0]+1] 	= ((nsREAL)(degree_[0]-j))/((nsREAL)degree_[0]);
	  L_[k*Loff_[0]+2] 	= ((nsREAL)j)/((nsREAL)degree_[0]);
	  k = k+1;
	}	
      for (j=1;j<degree_[0];++j)
	{
	  L_[k*Loff_[0]+0] 	= ((nsREAL)1.0)-((nsREAL)(degree_[0]-j))/((nsREAL)degree_[0]);
	  L_[k*Loff_[0]+1] 	= ((nsREAL)0.0);
	  L_[k*Loff_[0]+2] 	= ((nsREAL)(degree_[0]-j))/((nsREAL)degree_[0]);	  
	  k = k+1;
	}	
      for (i=1;i<degree_[0]-1;++i)
	{
	  for (j=1;j<degree_[0]-i;++j)
	    {
	      L_[k*Loff_[0]+0] 	= ((nsREAL)1.0)-((nsREAL)i)/((nsREAL)degree_[0])-((nsREAL)j)/((nsREAL)degree_[0]);
	      L_[k*Loff_[0]+1] 	= ((nsREAL)i)/((nsREAL)degree_[0]);	 
	      L_[k*Loff_[0]+2] 	= ((nsREAL)j)/((nsREAL)degree_[0]);	  

	      k = k+1;
	    }
	}
    }
}



      
#if 0

static const nsREAL __vmkS_lagrange_quadrangle_spl0[2] = {((nsREAL)0.0),
							  ((nsREAL)0.0)};
static const nsREAL __vmkS_lagrange_quadrangle_spl1[8] = {((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)1.0),
							  ((nsREAL)1.0),((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0)};
static const nsREAL __vmkS_lagrange_quadrangle_spl2[18] = {((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)1.0),((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)1.0),((nsREAL)1.0),
							   ((nsREAL)1.0),((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)-1.0),((nsREAL)1.0),((nsREAL)-1.0),((nsREAL)-1.0),};

static const nsREAL __vmkS_lagrange_interval_spl0[1] = {((nsREAL)0.0)};
static const nsREAL __vmkS_lagrange_interval_spl1[2] = {((nsREAL)-1.0),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl2[3] = {((nsREAL)-1.0),((nsREAL)0.0),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl3[4] = {((nsREAL)-1.0),((nsREAL)-1.0)/((nsREAL)3.0),((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl4[5] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};

static const nsREAL __vmkS_lagrange_interval_spl5[6] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl6[7] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl7[8] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl8[9] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl9[10] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_interval_spl10[11] = {((nsREAL)-1.0),((nsREAL)-0.5),((nsREAL)0.0),((nsREAL)0,5),((nsREAL)1.0)};


static const nsREAL __vmkS_lagrange_triangle_spl0[2] = {((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)1.0)/((nsREAL)3.0)};
static const nsREAL __vmkS_lagrange_triangle_spl1[6] = {((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),
							((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0)};
static const nsREAL __vmkS_lagrange_triangle_spl2[12] = {((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)0.5),((nsREAL)0.5),((nsREAL)0.0),
							 ((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)0.5),((nsREAL)0.5)};
static const nsREAL __vmkS_lagrange_triangle_spl3[20] = {((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)2.0)/((nsREAL)3.0),((nsREAL)2.0)/((nsREAL)3.0),((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0)/((nsREAL)3.0),
							 ((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)2.0)/((nsREAL)3.0),((nsREAL)2.0)/((nsREAL)3.0),((nsREAL)1.0)/((nsREAL)3.0),((nsREAL)1.0)/((nsREAL)3.0)};

static const nsREAL __vmkS_lagrangebubble_triangle_spl1[8] = {((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)1.0)/((nsREAL)3.0),
							      ((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0),((nsREAL)1.0)/((nsREAL)3.0)};
static const nsREAL __vmkS_lagrangebubble_triangle_spl2[14] = {((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)0.5),((nsREAL)0.5),((nsREAL)0.0),((nsREAL)1.0)/((nsREAL)3.0),
							       ((nsREAL)0.0),((nsREAL)0.0),((nsREAL)1.0),((nsREAL)0.0),((nsREAL)0.5),((nsREAL)0.5),((nsREAL)1.0)/((nsREAL)3.0)};

#endif



mok_m mkS_sample_new_mok_m(cst_mkS r_interp)
{
  const nsINT degree 	= mkS_k(r_interp);
  const nsINT n 	= mkS_n(r_interp);
  mok_m coo 		= mok_m_new(n,2);
  mkR xcoo = mok_m_x(coo);
  const nsINT off 	= mok_m_off(coo);
  cst_emkZELM elm 	= mkS_emkZELM(r_interp);
  __mmkS_lagrange_localspl(elm, &degree,xcoo,&off);
  return coo;
}



