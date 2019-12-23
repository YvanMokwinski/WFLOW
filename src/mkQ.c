#include "Type.h"
#include <math.h>

#include "mkQ.h"
#include "Blas.h"
#include "mkQ_data_tria_dunavant.c"
#include "mkQ_data_interval_hammer.c"
#include "mkQ_data_tria_hammer.c"
#include "mkQ_data_tetra_hammer.c"
#include "ns_constantes.h"

#if 0
#ifndef NDEBUG


emkQ 	mkQ_enum     	(cst_mkQ q_)	{  return Q_->kind; }
I 	mkQ_n		(cst_mkQ Q_)	{  return Q_->qn[0]; }
cst_pR	cst_mkQ_p	(cst_mkQ Q_) 	{  return Q_->qp; }
cst_pR	cst_mkQ_w	(cst_mkQ Q_) 	{  return Q_->qw; }

I 	mkQ_exact_for	(cst_emkQ_family __emkQ_family,	 
			 const I 	 k) {  return (k+2)/2+(k+2)%2;}
#endif
#endif

void  mkQ_collapse_(cst_pI  		qinterval_n_,
		    cst_pR 		qinterval_p_,
		    cst_pI  		qinterval_poff_,
		    cst_pR 		qinterval_w_,
		    cst_pI  		qinterval_woff_,
		    pI 		qelm_n_,
		    pR 		qelm_p_,
		    cst_pI 		qelm_poff_,
		    pR 		qelm_w_,
		    cst_pI 		qelm_woff_,
		    pI 		ferr_)
{
  I i,j;
  const I n = qinterval_n_[0];
  ferr_[0]	= (I)0;
  qelm_n_[0] 	= n*n;
  if (qelm_poff_[0]<qelm_n_[0])
    {
      fprintf(stderr,"__rmpsquadrature_collapse:error:not enough work memory qelm_poff="ifmt" < qelm_n = "ifmt "\n",
	      qelm_poff_[0],
	      qelm_n_[0]);
      ferr_[0]	= (I)1;
      return;
    }
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      {
	qelm_w_[(j*n+i)*qelm_woff_[0]] = ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * qinterval_w_[i*qinterval_woff_[0]] * qinterval_w_[j*qinterval_woff_[0]] * ((R)0.125);
	qelm_p_[j*n+i]                 = ( ((R)1.0) + qinterval_p_[i*qinterval_poff_[0]] ) * ((R)0.5);
	qelm_p_[qelm_poff_[0]+j*n+i]   = ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * ( ((R)1.0) + qinterval_p_[j*qinterval_poff_[0]] ) *  ((R)0.25);
      }
  return;
}

void  mkQ_collapse_triangle_(cst_pI  	qinterval_n_,
			     cst_pR 	qinterval_p_,
			     cst_pI  	qinterval_poff_,
			     cst_pR 	qinterval_w_,
			     cst_pI  	qinterval_woff_,
			     pI 		qelm_n_,
			     pR 		qelm_p_,
			     cst_pI 	qelm_poff_,
			     pR 		qelm_w_,
			     cst_pI 	qelm_woff_,
			     pI 		ferr_)
{
  I i,j;
  const I n = qinterval_n_[0];
  ferr_[0]	= (I)0;
  qelm_n_[0] 	= n*n;
  if (qelm_poff_[0]<qelm_n_[0])
    {
      fprintf(stderr,"__rmpsquadrature_collapse:error:not enough work memory qelm_poff="ifmt" < qelm_n = " ifmt "\n",
	      qelm_poff_[0],
	      qelm_n_[0]);
      ferr_[0]	= (I)1;
      return;
    }
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      {
	qelm_w_[(j*n+i)*qelm_woff_[0]] = ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * qinterval_w_[i*qinterval_woff_[0]] * qinterval_w_[j*qinterval_woff_[0]] * ((R)0.125);
	qelm_p_[j*n+i]                 = ( ((R)1.0) + qinterval_p_[i*qinterval_poff_[0]] ) * ((R)0.5);
	qelm_p_[qelm_poff_[0]+j*n+i]   = ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * ( ((R)1.0) + qinterval_p_[j*qinterval_poff_[0]] ) *  ((R)0.25);
      }
  return;
}

void  mkQ_collapse_tetrahedra_(cst_pI  	qinterval_n_,
			       cst_pR 		qinterval_p_,
			       cst_pI  	qinterval_poff_,
			       cst_pR 		qinterval_w_,
			       cst_pI  	qinterval_woff_,
			       pI 		qelm_n_,
			       pR 		qelm_p_,
			       cst_pI 		qelm_poff_,
			       pR 		qelm_w_,
			       cst_pI 		qelm_woff_,
			       pI 		ferr_)
{
  const I n = qinterval_n_[0];
  ferr_[0]	= (I)0;
  qelm_n_[0] 	= n*n*n;
  if (qelm_poff_[0]<qelm_n_[0])
    {
      fprintf(stderr,"__rmpsquadrature_collapse:error:not enough work memory qelm_poff="ifmt" < qelm_n = " ifmt "\n",
	      qelm_poff_[0],
	      qelm_n_[0]);
      ferr_[0]	= (I)1;
      return;
    }

  { I i;
    for (i=0;i<n;++i)
      {
	{ I j;
	  for (j=0;j<n;++j)
	    {
	      { I k;
		for (k=0;k<n;++k)
		  {
		    qelm_w_[(j*n+i)*qelm_woff_[0]] 	= ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * qinterval_w_[i*qinterval_woff_[0]] * qinterval_w_[j*qinterval_woff_[0]] * ((R)0.125);
		    qelm_p_[j*n+i]                 	= ( ((R)1.0) + qinterval_p_[i*qinterval_poff_[0]] ) * ((R)0.5);
		    qelm_p_[qelm_poff_[0]+j*n+i]   	= ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * ( ((R)1.0) + qinterval_p_[j*qinterval_poff_[0]] ) *  ((R)0.25);
		    qelm_p_[2*qelm_poff_[0]+j*n+i]   	= ( ((R)1.0) - qinterval_p_[i*qinterval_poff_[0]] ) * ( ((R)1.0) + qinterval_p_[j*qinterval_poff_[0]] ) *  ((R)0.25);
		  } }
	    } }
      } }
  return;
}



static void mkQ_legendre_c(cst_pI n,pR b)
{
  I k;
  for (k=0;k<n[0];++k)
    b[k]=(R)0.0;
  b[0]=(R)2.0;
  b[1]=(R)0.0;
  for (k=2;k<n[0];++k)
    {
      b[k]=((R)1.0)/(((R)4.0) - ((R)1.0)/( (  (R)(k-1) ) * ((R)(k-1)) ) );
    }
#if 0
  if (k==1)printf("aaaaaaaaaa %e\n",((R)1.0)/( (  (R)(k-1) ) * ((R)(k-1)) ));
#endif
}


static void sort_pour_divector_q(pR C,
				 pR II,
				 const I g,
				 const I d)
{
  if ( g>=d ) return;
  else
    {
      I i,der; R tmpv,tmpi;      
      tmpv 		=  C[g];      C[g] 		=  C[( (g+d)/2 )];      C[( (g+d)/2 )]	= tmpv;
      tmpi 		=  II[g];      II[g] 		=  II[( (g+d)/2 )];      II[( (g+d)/2 )]	= tmpi;      
      der=g;
      for(i=g+1; i<=d; i++)
	if ( C[i] < C[g] )
	  {
	    der=der+1;
	    tmpv 	=  C[der];	      C[der] 	=  C[i];	      C[i]	=  tmpv;
	    tmpi 	=  II[der];	      II[der] 	=  II[i];	      II[i]	=  tmpi;	      
	  }      
      tmpv 	=  C[der];      C[der] 	=  C[g];      C[g]	=  tmpv;
      tmpi 	=  II[der];      II[der] 	=  II[g];      II[g]	=  tmpi;
      sort_pour_divector_q(C,II,g,der-1);      sort_pour_divector_q(C,II,der+1,d);
    }
}/*-- --*/

void mkQ_legendre_interval_(cst_pI 	nspl_,
			    pR 		p_,
			    cst_pI 	p_off_, 
			    pR 		w_,
			    cst_pI 	w_off_, 
			    pR 		work_,
			    cst_pI 	work_n_,
			    pI 		ferr_)
{
  ferr_[0]	= (I)NO;
  I workn	= work_n_[0];
  I d 	= nspl_[0]+1;
  I n    	= d;
  I nn 	= n+1;
  char trN[1] 	= {'N'};
  char trV[1] 	= {'V'};
  I i,j,ldvl=1,info;
#if 0
  if (lwork<4*(p->degree+1)+2*n+n*n+2*n+n*n)
      * a    = &work_[0],
#endif
    R 
      * vl = NULL,

      * b    = &work_[nn],
      * jacM = &work_[2*nn],
      * wr   = &work_[2*nn+nn*nn],
      * wi   = &work_[2*nn+nn*nn+nn],
      * vr   = &work_[2*nn+nn*nn+2*nn];  
  if (work_n_[0] >= (d+1)*(5 + d) + 4*d )
    {

      mkQ_legendre_c(&n,b);
      
      for (i=0;i<n;++i)
	for (j=0;j<n;++j)
	  jacM[i+j*n]=(R)0.0;
      for (i=0;i<n;++i)
	jacM[i+i*n]=(R)0.0;
      for (i=0;i<n-1;++i)
	jacM[i+i*n+1]=nsSQRT(b[1+i]);
      for (i=1;i<n;++i)
	jacM[i+i*n-1]=nsSQRT(b[i]);
      
      dgeev(trN,
	    trV,
	    &n,
	    jacM,
	    &n, 
	    wr, 
	    wi, 
	    vl, 
	    &ldvl, 
	    vr, 
	    &n, 
	    &work_[2*nn+2*nn*nn+2*nn], 
	    &workn,
	    &info);
      for (i=0;i<n-1;++i) 
	{
	  p_[i * p_off_[0]] = wr[i];
	  w_[i * w_off_[0]] = vr[1+i*n]*vr[1+i*n]*((R)2.0);
	}
      sort_pour_divector_q(p_,w_,0,n-2);
#if 0
      for (i=0;i<n-1;++i) 
	{
	  printf("%e\n",w_[i]);
	}
#endif
    }
  else
    {

      fprintf(stderr,"not enough memory cquadrature_legendre_interval_ work_n_[0] = "ifmt" < "ifmt"= (d+1)*(5 + d) + 4*d with d = (nspl+1) ",work_n_[0], (d+1)*(5 + d) + 4*d );
      ferr_[0]=(I)1;
    }
}

void 	mkQ_write_eps	(const char * 	name_,
			 cst_mkQ 	q_,
			 const L 	q_show_,
			 const L 	e_show_)
{
#if 0
  const I elmn 	=  __eTo_nvertex(q_->elm);
  cst_pR elmcoo 	= __emkZELM_refcoo(q_->elm);
  const I q_show 	= (I)q_show_;
  const I e_show 	= (I)e_show_;
  mkEPS_FILE_mkQ 	(name_,
			 q_->elm,
			 &elmn,
			 elmcoo,
			 &elmn,  
			 &e_show, 
			 &q_->qn,
			 q_->qp, 
			 &q_->qn,
			 &q_show);
#endif
  fprintf(stderr,"come and correct\n");
  exit(1);
}


void mkQ_free(mkQ	Q_)
{
  if (Q_)
    {
      if (Q_->own_qp)
	free(Q_->own_qp);
      if (Q_->own_qw)
	free(Q_->own_qw);
      memset(Q_,0,sizeof(mkQ_ST));
    }
}

mkQ	mkQ_kill(mkQ	Q_)
{
  if (Q_)
    {
      mkQ_free(Q_);
      free(Q_);
    }
  return NULL;
}

void     mkQ_def_adaptative	(mkQ 		Q_,
				 const I 	N_,
				 const eTopology	elm_,
				 Err*	err_)
{
#if __mk_debug__
  __mk_debug(Q_,mkQ_def_adaptative);
  __mk_debug(err_,mkQ_def_adaptative);
  __mk_debug(((N_+1)>0),mkQ_def_adaptative);
#endif
  memset(Q_,0,sizeof(mkQ_ST));
  err_[0] 	= __eErr_no;
  Q_->kind	= __emkQ_adaptative;
  switch(elm_)
    {
    case __eTopology_QUADRILATERAL:
      {
	Q_->own_qp = malloc(N_*2*sizeof(R));
	Q_->own_qw = malloc(N_*sizeof(R));
	break;
      }
    case __eTopology_SEGMENT:
      {
	Q_->own_qp = malloc(N_*sizeof(R));
	Q_->own_qw = malloc(N_*sizeof(R));
	break;
      }
    case __eTopology_TRIANGLE:
      {
	Q_->own_qp = malloc(N_*2*sizeof(R));
	Q_->own_qw = malloc(N_*sizeof(R));
	break;
      }
    case __eTopology_TETRAHEDRON:
      {
	Q_->own_qp = malloc(N_*3*sizeof(R));
	Q_->own_qw = malloc(N_*sizeof(R));	
	break;
      }
    case __eTopology_ERROR:
    case __eTopology_ALL:
#if NOT  __mok_pure_enum__
    default:
#endif
      {
	err_[0] = __eErr_switch;
	fprintf(stderr,"mkQ_def_adaptative:switch failed on __emkZELM");
	break;
      }
    }
  
  Q_->qp = Q_->own_qp;
  Q_->qw = Q_->own_qw;

  if (err_[0])
    return;
  return;
}


mkQ	mkQ_new_adaptative	(const I 	N_,
				 const eTopology 	elm_)
{
#if __mk_debug__
  __mk_debug(((N_+1)>0),mkQ_new_adaptative);
#endif
  Err err 	= __eErr_no;
  mkQ Q 	= calloc(1,sizeof(mkQ_ST));
  if (NOT Q)
    {
      fprintf(stderr,"mkQ_new_adaptative:calloc failed\n");
      return NULL;
    }
  mkQ_def_adaptative(Q,
		     N_,
		     elm_,
		     &err);
  if (err)
    {
      fprintf(stderr,"mkQ_new_adaptative:mkQ_def_adaptative failed\n");
      return mkQ_kill(Q);
    }
  return Q;  
}


emkQ mkQ_compute_enum(cst_emkQ_family 	family_,
		      const eTopology	elm_,
		      cst_pI		degree_)
{
  emkQ e 		= __emkQ_error;
  const I degree 	= degree_[0];
  switch(family_)
    {
    case __emkQ_family_dunavant:
      {
	switch(elm_)
	  {
	  case __eTopology_TRIANGLE:
	    {
	      switch(degree)
		{
		case 1:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_1;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_3;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_4;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_6;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_7;
		    break;
		  }
		case 12:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_12;
		    break;
		  }
		case 13:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_13;
		    break;
		  }
		case 16:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_16;
		    break;
		  }
		case 19:
		  {
		    e=__emkQ_DUNAVANT_TRIANGLE_19;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle dunavant nshape = (1,3,4,6,7,12,13,16,19) but not nshape = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_QUADRILATERAL:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_QUADRILATERAL dunavant family not available");
	      break;
	    }
	  case __eTopology_POLYGON:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_POLYGON dunavant family not available");
	      break;
	    }
	  case __eTopology_TETRAHEDRON:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_TETRAHEDRON dunavant family not available");
	      break;
	    }

	  case __eTopology_HEXAHEDRON:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_HEXAHEDRON dunavant family not available");
	      break;
	    }

	  case __eTopology_WEDGE:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_WEDGE dunavant family not available");
	      break;
	    }
	  case __eTopology_PYRAMID:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_PYRAMID dunavant family not available");
	      break;
	    }
	  case __eTopology_POLYHEDRON:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_POLYHEDRON dunavant family not available");
	      break;
	    }

	  case __eTopology_SEGMENT:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_SEGMENT dunavant family not available");
	      break;
	    }
	  case __eTopology_NODE:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_NODE dunavant family not available");
	      break;
	    }
	  case __eTopology_POLYLINE:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:__eTopology_POLYLINE dunavant family not available");
	      break;
	    }
	  case __eTopology_ALL:
	  case __eTopology_ERROR:
	    {
	      e=__emkQ_error;
	      fprintf(stderr,"mkQ_compute_enum:switch failed on __emkZELM");
	      break;
	    }
	  }
	break;
      }
    case __emkQ_family_gauss:
      {
	switch(elm_)
	  {
	  case __eTopology_TRIANGLE:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_4;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_9;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_16;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_25;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_36;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_49;
		    break;
		  }
		case 8:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_64;
		    break;
		  }
		case 9:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_81;
		    break;
		  }
		case 10:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_100;
		    break;
		  }
		case 11:
		  {
		    e=__emkQ_GAUSS_TRIANGLE_121;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_QUADRILATERAL:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_4;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_9;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_16;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_25;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_36;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_49;
		    break;
		  }
		case 8:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_64;
		    break;
		  }
		case 9:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_81;
		    break;
		  }
		case 10:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_100;
		    break;
		  }
		case 11:
		  {
		    e=__emkQ_GAUSS_QUADRANGLE_121;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}
	      break;
	    }
	  case __eTopology_TETRAHEDRON:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_8;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_27;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_64;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_125;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_216;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_GAUSS_TETRAHEDRA_343;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}
	      break;
	    }
	  case __eTopology_SEGMENT:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_GAUSS_INTERVAL_2;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_GAUSS_INTERVAL_3;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_GAUSS_INTERVAL_4;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_GAUSS_INTERVAL_5;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_GAUSS_INTERVAL_6;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_GAUSS_INTERVAL_7;
		    break;
		  }
		case 8:
		  {
		    e=__emkQ_GAUSS_INTERVAL_8;
		    break;
		  }
		case 9:
		  {
		    e=__emkQ_GAUSS_INTERVAL_9;
		    break;
		  }
		case 10:
		  {
		    e=__emkQ_GAUSS_INTERVAL_10;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}
	      break;
	    }
	  case __eTopology_ALL:
	  case __eTopology_ERROR:
	    {
	      break;
	    }
	  default:
	    {
	      break;
	    }

	  }
	break;
      }
    case __emkQ_family_hammer:
      {
	switch(elm_)
	  {
	  case __eTopology_TRIANGLE:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_1;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_3;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_4;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_6;
		    break;
		  }
		case 6:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_12;
		    break;
		  }
		case 7:
		  {
		    e=__emkQ_HAMMER_TRIANGLE_42;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:hammer triangle  degree>7 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_QUADRILATERAL:
	    {
	      switch(degree)
		{
		case 1:
		  {
		    e=__emkQ_HAMMER_QUADRANGLE_1;
		    break;
		  }
		case 2:
		  {
		    e=__emkQ_HAMMER_QUADRANGLE_2;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_HAMMER_QUADRANGLE_3;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_HAMMER_QUADRANGLE_4;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_HAMMER_QUADRANGLE_5;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_TETRAHEDRON:
	    {
	      switch(degree)
		{
		case 2:
		  {
		    e=__emkQ_HAMMER_TETRAHEDRA_1;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_HAMMER_TETRAHEDRA_4;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_HAMMER_TETRAHEDRA_16;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_SEGMENT:
	    {
	      switch(degree)
		{
		case 1:
		  {
		    e=__emkQ_HAMMER_INTERVAL_1;
		    break;
		  }
		case 2:
		  {
		    e=__emkQ_HAMMER_INTERVAL_2;
		    break;
		  }
		case 3:
		  {
		    e=__emkQ_HAMMER_INTERVAL_3;
		    break;
		  }
		case 4:
		  {
		    e=__emkQ_HAMMER_INTERVAL_4;
		    break;
		  }
		case 5:
		  {
		    e=__emkQ_HAMMER_INTERVAL_5;
		    break;
		  }
		default:
		  {
		    fprintf(stderr,"mkQ_compute_enum:triangle gauss degree>11 not available degree = "ifmt"",degree);
		    e=__emkQ_error;
		    break;
		  }
		}	      
	      break;
	    }
	  case __eTopology_ALL:
	  case __eTopology_ERROR:
	    {
	      break;
	    }
	  default:
	    {
	      break;
	    }
	  }
	break;
      }
    case __emkQ_family_error:
    case __emkQ_family_n:
#if NOT __mok_pure_enum__
    default:
#endif
      {
	e = __emkQ_error;
	fprintf(stderr,"mkQ_compute_enum:switch failed on __emkQ_family");
	break;
      }
    }
  return e;
}

mkQ	mkQ_new	(cst_emkQ_family 	family_,
		 const eTopology 	elm_,
		 cst_pI			degree_)
{
  Err err 	= __eErr_no;
  mkQ Q 	= calloc(1,sizeof(mkQ_ST));
  if (NOT Q)
    {
      fprintf(stderr,"mkQ_new:calloc failed\n");
      return NULL;
    }
  mkQ_def(Q,
	  family_,
	  elm_,
	  degree_,
	  &err);
  if (err)
    {
      fprintf(stderr,"mkQ_new:mkQ_def failed\n");
      return mkQ_kill(Q);
    }
  return Q;
}






void mkQ_def_gauss(mkQ 	q_,
		   const eTopology elm_,
		   const I nspl_,
		   Err*		err_)
{
  err_[0] 	= __eErr_no;
  I ferr;
  switch(elm_)
    {
    case __eTopology_TRIANGLE:
      {
	const I n1d		= nspl_+1;	      
	q_->own_qp 		= malloc(n1d*n1d*2*sizeof(R));
	q_->own_qw 		= malloc(n1d*n1d*sizeof(R));
	q_->qp = q_->own_qp;
	q_->qw = q_->own_qw;
	q_->own_qn = n1d*n1d;
	q_->qn = &q_->own_qn;
	const I rwork_n 	= n1d*n1d+12*n1d+16; 
	pR rwork 	  	= malloc((n1d*2+rwork_n*2)*sizeof(R));
	const I Q_off   	= n1d*n1d;
	I q_n 		= (I)0;
	mkQ_legendre_interval	(&n1d,
				 &rwork[n1d],
				 &negal1, 
				 &rwork[0],
				 &negal1,
				 &rwork[2*n1d],
				 &rwork_n,
				 &ferr);
	mkQ_collapse		(&n1d,
				 &rwork[n1d],
				 &negal1,
				 &rwork[0],
				 &negal1,
				 &q_n,
				 q_->own_qp,
				 &Q_off,
				 q_->own_qw,
				 &negal1,
				 &ferr);            	
	free(rwork);
	break;
      }
    case __eTopology_SEGMENT:
      {
	const I n1d		= nspl_;	      
	q_->own_qp 		= malloc(n1d*sizeof(R));
	q_->own_qw 		= malloc(n1d*sizeof(R));
	q_->qp = q_->own_qp;
	q_->qw = q_->own_qw;
	q_->own_qn = n1d;
	q_->qn = &q_->own_qn;

	const I rwork_n 	= n1d*n1d+12*n1d+16;
	pR rwork 		= malloc(rwork_n*2*sizeof(R));
	mkQ_legendre_interval(&nspl_,
			      q_->own_qp,
			      &negal1, 
			      q_->own_qw,
			      &negal1,
			      rwork,
			      &rwork_n,
			      &ferr);
	free(rwork);
	break;
      }
    case __eTopology_QUADRILATERAL:
      {
	const I n1d		= nspl_;	      
	q_->own_qp 		= malloc(n1d*n1d*2*sizeof(R));
	q_->own_qw 		= malloc(n1d*n1d*sizeof(R));
	q_->qp = q_->own_qp;
	q_->qw = q_->own_qw;
	q_->own_qn = n1d*n1d;
	q_->qn = &q_->own_qn;
	const I rwork_n 	= n1d*n1d+12*n1d+16;
	pR rwork 		= malloc((n1d*2+2*rwork_n)*sizeof(R));
	mkQ_legendre_interval	(&nspl_,
				 &rwork[0],
				 &negal1, 
				 &rwork[n1d],
				 &negal1,
				 &rwork[2*n1d],
				 &rwork_n,
				 &ferr);
	Blas_dger(&n1d,&n1d,&regal1,&rwork[n1d],&negal1,&rwork[n1d],&negal1,q_->own_qw,&n1d);
	{  I i;for (i=0;i<n1d;++i) { I j;for (j=0;j<n1d;++j){ q_->own_qp[i*n1d+j] = rwork[i]; } } } 
	{  I i;for (i=0;i<n1d;++i) { Blas_dcopy(&n1d,rwork,&negal1,&q_->own_qp[n1d+i*n1d],&negal1); } } 
	free(rwork);
	break;
      }
    case __eTopology_TETRAHEDRON:
      {
	const I n1d		= nspl_;	      
	q_->own_qp 		= malloc(n1d*n1d*n1d*3*sizeof(R));
	q_->own_qw 		= malloc(n1d*n1d*n1d*sizeof(R));
	q_->qp = q_->own_qp;
	q_->qw = q_->own_qw;
	q_->own_qn = n1d*n1d*n1d;
	q_->qn = &q_->own_qn;
	const I rwork_n 	= n1d*n1d+12*n1d+16;
	pR rwork 		= malloc((n1d*2+2*rwork_n)*sizeof(R));
	mkQ_legendre_interval	(&nspl_,
				 &rwork[0],
				 &negal1, 
				 &rwork[n1d],
				 &negal1,
				 &rwork[2*n1d],
				 &rwork_n,
				 &ferr);
	free(rwork);
	break;
      }
    case __eTopology_ALL:
    case __eTopology_ERROR:
#if NOT __mok_pure_enum__
    default:
#endif
      {
	break;
      }
    }
  if (err_[0])
    return;
}



void mkQ_def_cadyf(mkQ 	q_,
		   Err*err_)
{
  q_->own_qp 		= malloc(42*2*sizeof(R));
  q_->own_qw 		= malloc(42*sizeof(R));
  q_->qp = q_->own_qp;
  q_->qw = q_->own_qw;
  q_->own_qn = 42;
  q_->qn = &q_->own_qn;

  R * 	rgauss = q_->own_qp;
  R * 	sgauss = &q_->own_qp[42];
  R*	wgauss = q_->own_qw;

#define w1 0.021883581369429 / 2.0
#define a1    0.022072179275643                            
#define b1    0.488963910362179 

rgauss[0-1+1]=a1;
sgauss[0-1+1]=b1;
wgauss[0-1+1]=w1;

rgauss[0-1+2]=b1;
sgauss[0-1+2]=a1;
wgauss[0-1+2]=w1;

rgauss[0-1+3]=b1;
sgauss[0-1+3]=b1;
wgauss[0-1+3]=w1;

#define w2 0.032788353544125/2.0
#define a2 0.164710561319092
#define b2 0.417644719340454

rgauss[0-1+4]=a2;
sgauss[0-1+4]=b2;
wgauss[0-1+4]=w2;

rgauss[0-1+5]=b2;
sgauss[0-1+5]=a2;
wgauss[0-1+5]=w2;

rgauss[0-1+6]=b2;
sgauss[0-1+6]=b2;
wgauss[0-1+6]=w2;



#define w3 0.051774104507292/2.0
#define a3 0.453044943382323
#define b3 0.273477528308839

rgauss[0-1+7]=a3;
sgauss[0-1+7]=b3;
wgauss[0-1+7]=w3;

rgauss[0-1+8]=b3;
sgauss[0-1+8]=a3;
wgauss[0-1+8]=w3;

rgauss[0-1+9]=b3;
sgauss[0-1+9]=b3;
wgauss[0-1+9]=w3;



#define w4 0.042162588736993/2.0
#define a4 0.645588935174913
#define b4 0.177205532412543

rgauss[0-1+10]=a4;
sgauss[0-1+10]=b4;
wgauss[0-1+10]=w4;

rgauss[0-1+11]=b4;
sgauss[0-1+11]=a4;
wgauss[0-1+11]=w4;

rgauss[0-1+12]=b4;
sgauss[0-1+12]=b4;
wgauss[0-1+12]=w4;



#define w5 0.014433699669777/2.0
#define a5 0.876400233818255
#define b5 0.061799883090873

rgauss[0-1+13]=a5;
sgauss[0-1+13]=b5;
wgauss[0-1+13]=w5;

rgauss[0-1+14]=b5;
sgauss[0-1+14]=a5;
wgauss[0-1+14]=w5;

rgauss[0-1+15]=b5;
sgauss[0-1+15]=b5;
wgauss[0-1+15]=w5;



#define w6 0.004923403602400/2.0
#define a6 0.961218077502598
#define b6 0.019390961248701

rgauss[0-1+16]=a6;
sgauss[0-1+16]=b6;
wgauss[0-1+16]=w6;

rgauss[0-1+17]=b6;
sgauss[0-1+17]=a6;
wgauss[0-1+17]=w6;

rgauss[0-1+18]=b6;
sgauss[0-1+18]=b6;
wgauss[0-1+18]=w6;



#define w7 0.024665753212564/2.0
#define a7 0.057124757403648
#define b7 0.172266687821356
#define c7 0.770608554774996

rgauss[0-1+19]=a7;
sgauss[0-1+19]=b7;
wgauss[0-1+19]=w7;

rgauss[0-1+20]=a7;
sgauss[0-1+20]=c7;
wgauss[0-1+20]=w7;

rgauss[0-1+21]=b7;
sgauss[0-1+21]=a7;
wgauss[0-1+21]=w7;

rgauss[0-1+22]=b7;
sgauss[0-1+22]=c7;
wgauss[0-1+22]=w7;

rgauss[0-1+23]=c7;
sgauss[0-1+23]=a7;
wgauss[0-1+23]=w7;

rgauss[0-1+24]=c7;
sgauss[0-1+24]=b7;
wgauss[0-1+24]=w7;



#define w8 0.038571510787061/2.0
#define a8 0.092916249356972
#define b8 0.336861459796345
#define c8 0.570222290846683

rgauss[0-1+25]=a8;
sgauss[0-1+25]=b8;
wgauss[0-1+25]=w8;

rgauss[0-1+26]=a8;
sgauss[0-1+26]=c8;
wgauss[0-1+26]=w8;

rgauss[0-1+27]=b8;
sgauss[0-1+27]=a8;
wgauss[0-1+27]=w8;

rgauss[0-1+28]=b8;
sgauss[0-1+28]=c8;
wgauss[0-1+28]=w8;

rgauss[0-1+29]=c8;
sgauss[0-1+29]=a8;
wgauss[0-1+29]=w8;

rgauss[0-1+30]=c8;
sgauss[0-1+30]=b8;
wgauss[0-1+30]=w8;



#define w9 0.014436308113534/2.0
#define a9 0.014646950055654
#define b9 0.298372882136258
#define c9 0.686980167808088

rgauss[0-1+31]=a9;
sgauss[0-1+31]=b9;
wgauss[0-1+31]=w9;

rgauss[0-1+32]=a9;
sgauss[0-1+32]=c9;
wgauss[0-1+32]=w9;

rgauss[0-1+33]=b9;
sgauss[0-1+33]=a9;
wgauss[0-1+33]=w9;

rgauss[0-1+34]=b9;
sgauss[0-1+34]=c9;
wgauss[0-1+34]=w9;

rgauss[0-1+35]=c9;
sgauss[0-1+35]=a9;
wgauss[0-1+35]=w9;

rgauss[0-1+36]=c9;
sgauss[0-1+36]=b9;
wgauss[0-1+36]=w9;



#define w10 0.005010228838501/2.0
#define a10 0.001268330932872
#define b10 0.118974497696957
#define c10 0.879757171370171

rgauss[0-1+37]=a10;
sgauss[0-1+37]=b10;
wgauss[0-1+37]=w10;

rgauss[0-1+38]=a10;
sgauss[0-1+38]=c10;
wgauss[0-1+38]=w10;

rgauss[0-1+39]=b10;
sgauss[0-1+39]=a10;
wgauss[0-1+39]=w10;

rgauss[0-1+40]=b10;
sgauss[0-1+40]=c10;
wgauss[0-1+40]=w10;

rgauss[0-1+41]=c10;
sgauss[0-1+41]=a10;
wgauss[0-1+41]=w10;

rgauss[0-1+42]=c10;
sgauss[0-1+42]=b10;
wgauss[0-1+42]=w10;

}


void mkQ_def_dunavant(mkQ 		q_,
		      const eTopology 	elm_,
		      const I 	nspl_,
		      Err*		err_)
{
  err_[0] 	= __eErr_no;
  switch(elm_)
    {
    case __eTopology_TRIANGLE:
      {
	{ I i;
	  for (i=0;i<9;++i)
	    if (nspl_<=n_tria_dunavant[i])
	      break;
	  if (i>9)
	    {
	      i= 9;
	      err_[0] = __eErr_user;
	      fprintf(stderr,"mkQ_def_dunavant:2D-rule with "ifmt" points is not available",nspl_);
exit(1);
		      break;
	    }

	  q_->own_qp = NULL;
	  q_->own_qw = NULL;
	  q_->qp     = pos_tria_dunavant[i];
	  q_->qw     = wei_tria_dunavant[i];
	  q_->qn     = &n_tria_dunavant[i]; }
	break;
      }
    case __eTopology_SEGMENT:
      {
	err_[0] = __eErr_user;
	fprintf(stderr,"mkQ_def_dunavant:__eTopology_SEGMENT rule not available");
exit(1);
		break;
      }
    case __eTopology_QUADRILATERAL:
      {
	err_[0] = __eErr_user;
	fprintf(stderr,"mkQ_def_dunavant:__eTopology_QUADRILATERAL rule not available");
	exit(1);
	break;
      }
    case __eTopology_TETRAHEDRON:
      {
	err_[0] = __eErr_user;
	fprintf(stderr,"mkQ_def_dunavant:__eTopology_TETRAHEDRON rule not available");
	exit(1);
	break;
      }
    case __eTopology_ALL:
    case __eTopology_ERROR:
#if NOT __mok_pure_enum__
    default:
#endif
      {
	break;
      }
    }
  if (err_[0])
    return;
}



/* HARD DATA */

void mkQ_def_hammer	(mkQ 		q_,
			 const eTopology 	elm_,
			 const I 	nspl_,
			 Err*	err_)
{
  err_[0] 	= __eErr_no;
  switch(elm_)
    {
    case __eTopology_SEGMENT:
      {
	{ I i;
	  for (i=0;i<5;++i)
	    if (nspl_<=n_interval_hammer[i])
	      break;
	  if (i>4)
	    {
	      i=4;
	      err_[0] = __eErr_user;
	      fprintf(stderr,"mkQ_def_hammer:1D-rule with "ifmt" points is not available",nspl_);
	exit(1);
	      break;
	    }
	  q_->own_qp = NULL;
	  q_->own_qw = NULL;
	  q_->qp     = pos_interval_hammer[i];
	  q_->qw     = wei_interval_hammer[i];
	  q_->qn     = &n_interval_hammer[i]; }
	break;
      }
    case __eTopology_TRIANGLE:
      {
	{ I i;
	  for (i=0;i<6;++i)
	    if (nspl_<=n_tria_hammer[i])
	      break;
	  if (i>5)
	    {
	      i= 5;
	      err_[0] = __eErr_user;
	      fprintf(stderr,"mkQ_def_hammer:2D-rule with "ifmt" points is not available",nspl_);
	exit(1);
	      break;
	    }

	  q_->own_qp = NULL;
	  q_->own_qw = NULL;
	  q_->qp     = pos_tria_hammer[i];
	  q_->qw     = wei_tria_hammer[i];
	  q_->qn     = &n_tria_hammer[i]; }
	break;
      }

    case __eTopology_TETRAHEDRON:
      {
	{ I i;
	  for (i=0;i<3;++i)
	    if (nspl_<=n_tetra_hammer[i])
	      break;
	  if (i>2)
	    {
	      i=2;
	      err_[0] = __eErr_user;
	      fprintf(stderr,"mkQ_def_hammer:3D-rule with "ifmt" points is not available",nspl_);
exit(1);
		      break;
	    }
	  q_->own_qp = NULL;
	  q_->own_qw = NULL;
	  q_->qp     = pos_tetra_hammer[i];
	  q_->qw     = wei_tetra_hammer[i];
	  q_->qn     = &n_tetra_hammer[i]; }
	break;
      }
      

    case __eTopology_QUADRILATERAL:
      {
#if 0
	for (i=0;i<5;++i)
	  if (nspl_<=n_interval_hammer[i])
	    break;
	if (i>4)
	  {
	    i=4;
	    mk_warn("mkQ_def_hammer:2D-rule with "ifmt" points is not available, set N to "ifmt"",nspl_,level_interval[i]);
	  }
	N 	= n_interval_hammer[i];
	q_->own_qp = malloc(N*N*2*sizeof(R));
	q_->own_qw = malloc(N*sizeof(R));
	q_->qp = q_->own_qp;
	q_->qw = q_->own_qw;
	q_->own_qn = N*N;       

	{
	  const I poff  = mok_m_off(q->lc);
	  pR rwork = malloc(N*sizeof(R));
	  Blas_dger(&N,&N,&regal1,&wei_interval_hammer[i],&negal1,&wei_interval_hammer[i],&negal1,mok_m_x(q->w),&N);
	  { pR px = q->own_qp; { I k;for (k=0;k<N;++k) { I j;for (j=0;j<N;++j){ px[k*poff+j] = pos_interval[i][k]; } } } }
	  { pR py = &q->own_qp[N]; { I k;for (k=0;k<N;++k) { Blas_dcopy(&N,pos_interval[i],&regal1,&py[i*poff],&regal1); } } }
	  free(rwork);
	}
#endif
	break;
      }
    case __eTopology_ALL:
    case __eTopology_ERROR:
#if NOT __mok_pure_enum__
    default:
#endif
      {
	err_[0] = __eErr_switch;
	fprintf(stderr,"mkQ_def_hammer:switch failed on __emkZELM");	
exit(1);
		break;
      }
    }
  if (err_[0])
    return;
}




void  mkQ_def(mkQ 		q_,
	      cst_emkQ_family 	family_,
	      const eTopology	elm_,
	      cst_pI 		degree,
	      Err*		err_)
{
  memset(q_,0,sizeof(mkQ_ST));
  err_[0] 	= __eErr_no;
  q_->kind	= mkQ_compute_enum(family_,elm_,degree);
  if (q_->kind==__emkQ_error)
    {
      err_[0] = __eErr_user;
      fprintf(stderr,"mkQ_compute_enum failed\n");
      return;
    }
  switch(family_)
    {
    case __emkQ_family_hammer:
      {
	mkQ_def_hammer(q_,
		       elm_,
		       degree[0],
		       err_);
	if (err_[0])
	  {
	    fprintf(stderr,"mkQ_def:mkQ_def_hammer failed\n");
	    break;
	  }
	break;
      }
    case __emkQ_family_gauss:
      {
	mkQ_def_gauss(q_,
		      elm_,
		      degree[0],
		      err_);
	if (err_[0])
	  {
	    fprintf(stderr,"mkQ_def:mkQ_def_gauss failed\n");
	    break;
	  }	
	break;
      }
    case __emkQ_family_dunavant:
      {
	mkQ_def_dunavant(q_,
			 elm_,
			 degree[0],
			 err_);	
	if (err_[0])
	  {
	    fprintf(stderr,"mkQ_def:mkQ_def_dunavant failed\n");
	    exit(1);
	    break;
	  }	
	break;
      }
    case __emkQ_family_error:
    case __emkQ_family_n:
      {
	err_[0] = __eErr_switch;
	fprintf(stderr,"mkQ_def:switch failed on __emkQ_family failed\n");
	exit(1);
	break;
      }
    }
  if (err_[0])
    return;
  return;
}
