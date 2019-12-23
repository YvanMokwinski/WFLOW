#include "mkS.h"
#include "Blas.h"
#include "ns_constantes.h"
#include <math.h>

#if 0
#include "mkSYS.h"

#include "mkS_eval.h"
#include "mok_m.h"
#endif


static void __rmkS_pascal(cst_pI  	pascal_n_,
			  pR 		pascal_x_,
			  cst_pI  	pascal_x_n_,
			  pI  		err)
{
  static I i,j;
  err[0] = (I)0;
  if ( pascal_x_n_[0] >= pascal_n_[0] * pascal_n_[0] )
    {
      pascal_x_[0]=(R)1.0;
      pascal_x_[1]=(R)1.0;
      pascal_x_[2]=(R)1.0;
      for (i=0;i<pascal_n_[0];++i)
	for (j=0;j<pascal_n_[0];++j)
	  pascal_x_[i+j*pascal_n_[0]] = (R)0.0;
      for (i=0;i<pascal_n_[0];++i)
	pascal_x_[i+i*pascal_n_[0]] = (R)1.0;
      for (i=0;i<pascal_n_[0];++i)
	pascal_x_[i] = (R)1.0;
      for (i=2;i<pascal_n_[0];++i)
	for (j=1;j<pascal_n_[0];++j)
	    pascal_x_[i+j*pascal_n_[0]] = pascal_x_[(i-1)+(j-1)*pascal_n_[0]] + pascal_x_[(i-1)+j*pascal_n_[0]];
    }
  else
    {
      fprintf(stderr,"__rmkS_pascal:not enough memory\n");
      err[0] = (I)1;
    }
}

static void __rmkS_canonic_l2sp_monom(const I    	p,
				      const I    	q,
				      pR  		m,
				      pR  		triangle_pascal_,
				      cst_pI   		N_)
{
  R x1,x2;
  I k;
  x1 = (R)0.0;
  x2 = (R)0.0;
  for (k=0;k<=q+1;k+=2)
    x1+=triangle_pascal_[q+1+k*N_[0]]/( ((R)p) + ((R)1.0) + ((R)k)  );
  for (k=1;k<=q+1;k+=2)
    x2+=triangle_pascal_[q+1+k*N_[0]]/( ((R)p) + ((R)1.0) + ((R)k)  );
  m[0] = (x1-x2)/( ((R)q) + ((R)1.0));
}


void __rmkS_L2Gram(pR 		m_,
		   cst_pI 	degree_,
		   cst_pI 	nshapes_,
		   cst_pI 	rwork_n_,
		   pR 		rwork_,
		   pI 		err)
{
  I i,p,q,k = nshapes_[0];
  const I pascal_n = (degree_[0]+1)*2;
  const I pascal_len = pascal_n*pascal_n;
  if (rwork_n_[0]<pascal_len)
    {
      fprintf(stderr,"*** DGERR mkS_L2Gram not enough rwork memory (rwork_n="ifmt" < ask_n="ifmt")\n",rwork_n_[0],pascal_len);
      err[0] = (I)1;
      return;
    }
  __rmkS_pascal(&pascal_n,rwork_,&pascal_len,err);
  for (p=0;p<=degree_[0];++p)
    for (q=0;q<=p;++q)
      for (k=0;k<=degree_[0];++k)
	for (i=0;i<=k;++i)
	  __rmkS_canonic_l2sp_monom(p-q+k-i,
				    q+i,
				    &m_[((p*(p+1))/2+q)*(((degree_[0]+1)*(degree_[0]+2))/2 )+ (k*(k+1))/2+i],
				    rwork_,
				    &pascal_n);
}




void mok_l2ortho(cst_pI 	n_,
		 pR 		m_,
		 pR 		RR, 
		 cst_pI 	rwork_n_,
		 pR 		rwork_)
{ 
  I i,j;
  for (i=0;i<n_[0];++i)
    for (j=0;j<n_[0];++j)
      RR[i*n_[0]+j]=(R)0.0;
  RR[0] = nsSQRT((R)2.0);
  for (i=1;i<n_[0];++i)
    RR[i*n_[0]+i]=(R)1.0;
  for (i=1;i<n_[0];++i)
    {
      /* produit scalaire avec les precedents */
      Blas_dgemv("T",n_,&i,&regal1,&RR[0],n_,&m_[i*n_[0]],&negal1,&regal0,&rwork_[0],&negal1);      
      Blas_dgemv("N",n_,&i,&mregal1,&RR[0],n_,&rwork_[0],&negal1,&regal1,&RR[i*n_[0]],&negal1);      
      Blas_dgemv("N",n_,n_,&regal1,&m_[0],n_,&RR[i*n_[0]],&negal1,&regal0,&rwork_[0],&negal1);      
      const R xy = Blas_ddot(n_,&rwork_[0],&negal1,&RR[i*n_[0]],&negal1);
      const R xy1 = nsSQRT((R)1.0/xy);
      Blas_dscal(n_,&xy1,&RR[i*n_[0]],&negal1);
    }
}

void mkS_l2ortho_tria(cst_pI 	degree,
		      cst_pI 	n,
		      pR 		r,
		      cst_pI 	roff_,
		      cst_pR 	p,
		      cst_pI 	poff_,
		      pR 		rwork_,
		      cst_pI 	rwork_n_,
		      pI 	err_)
{
  err_[0] = (I)0;
  const I nshapes 		= ( (degree[0]+1)*(degree[0]+2) )/2;
  const I i_evalL2Ortho 	= 0;
  const I i_L2Gram      	= nshapes*nshapes;
  I i_newrwork    		= MAX(nshapes*nshapes + nshapes*n[0],2*nshapes*nshapes);
  if (rwork_n_[0]<i_newrwork)
    {
      fprintf(stderr,"*** DGERR mkS_l2ortho not enough rwork memory (rwork_n="ifmt" < ask_n="ifmt"), rule  is  MAX(n*n + n*nsamples,2*n*n)\n",rwork_n_[0],i_newrwork);
      err_[0] = (I)1;
      return;
    }
  I new_rwork_n  = rwork_n_[0] - i_newrwork;  
  __rmkS_L2Gram(&rwork_[i_L2Gram],
		    degree,
		    &nshapes,
		    &new_rwork_n,
		    &rwork_[i_newrwork],
		    err_);
  
  mok_l2ortho(&nshapes,&rwork_[i_L2Gram],&rwork_[i_evalL2Ortho],&new_rwork_n,&rwork_[i_newrwork]);
  /* maintenant on evalue les fonctions de bases canoniques */
  
  new_rwork_n  = new_rwork_n + nshapes*nshapes;
  i_newrwork   = i_newrwork  - nshapes*nshapes;
  
  if (i_newrwork<0)
    {
      fprintf(stderr,"mkS_eval_l2ortho_tria:not enough space\n");
      exit(1);
    }
  
  mkS_canonic_tria( degree,
		    n,
		    &rwork_[i_newrwork],
		    &nshapes,
		    p,
		    poff_,
		    NULL,
		    NULL,
		    err_);


  Blas_dgemm("T",
	       "N",
	       &nshapes,
	       n,
	       &nshapes,
	       &regal1,
	       &rwork_[i_evalL2Ortho],
	       &nshapes,
	       &rwork_[i_newrwork],
	       &nshapes,
	       &regal0,
	       r,
	       roff_);


}

void mkS_dx_l2ortho_tria(cst_pI 	degree,
			 cst_pI 	n,
			 pR 	r,
			 cst_pI 	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 pR 	rwork_,
			 cst_pI 	rwork_n_,
			 pI 	err_)
{
  err_[0] = (I)0;
  const I nshapes = ( (degree[0]+1)*(degree[0]+2) )/2;
  const I i_evalL2Ortho = 0;
  const I i_L2Gram      = nshapes*nshapes;
  I i_newrwork    = MAX(nshapes*nshapes + nshapes*n[0],2*nshapes*nshapes);
  if (rwork_n_[0]<i_newrwork)
    {
      fprintf(stderr,"*** DGERR mkS_l2ortho not enough rwork memory (rwork_n="ifmt" < ask_n="ifmt"), rule  is  MAX(n*n + n*nsamples,2*n*n)\n",rwork_n_[0],i_newrwork);
      err_[0] = (I)1;
      return;
    }
  I new_rwork_n  = rwork_n_[0] - i_newrwork;
  
  __rmkS_L2Gram(&rwork_[i_L2Gram],
		    degree,
		    &nshapes,
		    &new_rwork_n,
		    &rwork_[i_newrwork],
		    err_);
  mok_l2ortho(&nshapes,&rwork_[i_L2Gram],&rwork_[i_evalL2Ortho],&new_rwork_n,&rwork_[i_newrwork]);

  /* maintenant on evalue les fonctions de bases canoniques */
  
  new_rwork_n  = new_rwork_n + nshapes*nshapes;
  i_newrwork   = i_newrwork  - nshapes*nshapes;
  mkS_dx_canonic_tria( degree,
	       n,
	       &rwork_[i_newrwork],
	       &nshapes,
	       p,
	       poff_,
	       NULL,
	       NULL,
	       err_);


  Blas_dgemm("T",
	       "N",
	       &nshapes,
	       n,
	       &nshapes,
	       &regal1,
	       &rwork_[i_evalL2Ortho],
	       &nshapes,
	       &rwork_[i_newrwork],
	       &nshapes,
	       &regal0,
	       r,
	       roff_);
}

void mkS_dy_l2ortho_tria(cst_pI  	degree,
			 cst_pI  	n,
			 pR 	 	r,
			 cst_pI  	roff_,
			 cst_pR 	p,
			 cst_pI 	poff_,
			 pR 		rwork_,
			 cst_pI 	rwork_n_,
			 pI 		err_)
{
  err_[0] = (I)0;
  const I nshapes 		= ( (degree[0]+1)*(degree[0]+2) )/2;
  const I i_evalL2Ortho 	= 0;
  const I i_L2Gram      	= nshapes*nshapes;
  I i_newrwork    		= MAX(nshapes*nshapes + nshapes*n[0],2*nshapes*nshapes);
  if (rwork_n_[0]<i_newrwork)
    {
      err_[0] = (I)1;
      fprintf(stderr,"*** DGERR mkS_l2ortho not enough rwork memory (rwork_n="ifmt" < ask_n="ifmt"), rule  is  MAX(n*n + n*nsamples,2*n*n)\n",rwork_n_[0],i_newrwork);
      return;
    }
  I new_rwork_n  = rwork_n_[0] - i_newrwork;  
  __rmkS_L2Gram(&rwork_[i_L2Gram],degree,&nshapes,&new_rwork_n,&rwork_[i_newrwork],err_);
  mok_l2ortho(&nshapes,&rwork_[i_L2Gram],&rwork_[i_evalL2Ortho],&new_rwork_n,&rwork_[i_newrwork]);
  /* maintenant on evalue les fonctions de bases canoniques */
  new_rwork_n  = new_rwork_n + nshapes*nshapes;
  i_newrwork   = i_newrwork  - nshapes*nshapes;
  mkS_dy_canonic_tria( degree,n,&rwork_[i_newrwork],&nshapes,p,poff_,NULL,NULL,err_);   
  Blas_dgemm("T",
	       "N",
	       &nshapes,
	       n,
	       &nshapes,
	       &regal1,
	       &rwork_[i_evalL2Ortho],
	       &nshapes,
	       &rwork_[i_newrwork],
	       &nshapes,
	       &regal0,
	       r,
	       roff_);
}

