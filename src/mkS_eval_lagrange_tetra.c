#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"
#include "mkS_eval.h"

nsINT mkS_tetra_n(cst_mkI degree)
{
  if (degree[0]==0)
    {
      return 1;
    }
  else if (degree[0]==1)
    {
      return 4;
    }
  else if (degree[0]==2)
    {
      return 10;
    }
  else if (degree[0]==3)
    {
      return 20;
    }
  else
    {
      mk_msg("mkS_tetra_n failed");
      return 0;
    }
}

#if __mk_debug__
void mkS_lagrange_localspl_tetra( cst_mkI 	degree,
				  mkR		coo_,
				  cst_mkI	off_)
{
  __mmkS_lagrange_localspl_tetra(degree,coo_,off_);
}
#endif



void mkS_dy_lagrange_tetra(cst_mkI degree,
			   cst_mkI n,
			   mkR r,
			   cst_mkI roff_,
			   cst_mkR p,
			   cst_mkI poff_,
			   mkR rwork,
			   cst_mkI rwork_n,
			   mkI err_)
{
  err_[0] = (nsINT)0;
  nsINT info,nshapes = mkS_tetra_n(degree);
  nsINT perm[nshapes];
  mkS_lagrange_localspl_tetra(degree,&rwork[0],&nshapes);
  mkS_canonic_tetra(degree,
		    &nshapes,
		    &rwork[3*nshapes],
		    &nshapes,
		    &rwork[0],
		    &nshapes,
		    NULL,
		    NULL,
		    err_);  
  mkS_dy_canonic_tetra(degree,
		       n,
		       r,
		       roff_,
		       p,
		       poff_,
		       NULL,
		       NULL,
		       err_);  
#if 1
  nslapack_dgesv(&nshapes,n,&rwork[3*nshapes],&nshapes,perm,r,roff_,&info);
#else
  char A1[2] = {'A','\0'};
  char A2[2] = {'A','\0'};
  nsREAL AS[36],AU[36*36],AV[36*36],rwork2[28*5];nsINT rwork2_n=28*5;
  nslapack_dgesvd(A1,A2,&nshapes,&nshapes,&rwork[2*nshapes],&nshapes,AS,AU,&nshapes,AV,&nshapes,rwork2,&rwork2_n,&info);
  { nsINT i,j;for (i=0;i<nshapes;++i){ nsREAL y =AS[i];if (nsFABS(AS[i])>((nsREAL)1.0e-15)) { y=((nsREAL)1.0)/AS[i]; for (j=0;j<nshapes;++j) {AV[i*nshapes+j] *=y; }   }else {for (j=0;j<nshapes;++j) {AV[i*nshapes+j] =(nsREAL)0.0;} } } }
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AU,&nshapes,r,roff_,&__vmps_blas_regal0,rwork2,&nshapes);
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AV,&nshapes,rwork2,&nshapes,&__vmps_blas_regal0,r,roff_);
#endif
}




void mkS_dx_lagrange_tetra(	cst_mkI degree,
				cst_mkI n,
				mkR r,
				cst_mkI roff_,
				cst_mkR p,
				cst_mkI  poff_,
				mkR rwork,
				cst_mkI rwork_n,
				mkI err_)
{
  nsINT info,nshapes = mkS_tetra_n(degree);
  nsINT perm[nshapes];
  mkS_lagrange_localspl_tetra		(degree,&rwork[0],&nshapes);
  mkS_canonic_tetra	(degree,&nshapes,&rwork[3*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
  mkS_dx_canonic_tetra(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
#if 1
  nslapack_dgesv(&nshapes,n,&rwork[3*nshapes],&nshapes,perm,r,roff_,&info);
#else
  char A1[2] = {'A','\0'};
  char A2[2] = {'A','\0'};
  nsREAL AS[36],AU[36*36],AV[36*36],rwork2[28*5];nsINT rwork2_n=28*5;
  nslapack_dgesvd(A1,A2,&nshapes,&nshapes,&rwork[2*nshapes],&nshapes,AS,AU,&nshapes,AV,&nshapes,rwork2,&rwork2_n,&info);
  { nsINT i,j;for (i=0;i<nshapes;++i){ nsREAL y =AS[i];if (nsFABS(AS[i])>((nsREAL)1.0e-15)) { y=((nsREAL)1.0)/AS[i]; for (j=0;j<nshapes;++j) {AV[i*nshapes+j] *=y; }   }else {for (j=0;j<nshapes;++j) {AV[i*nshapes+j] =(nsREAL)0.0;} } } }
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AU,&nshapes,r,roff_,&__vmps_blas_regal0,rwork2,&nshapes);
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AV,&nshapes,rwork2,&nshapes,&__vmps_blas_regal0,r,roff_);
#endif
}


void mkS_dz_lagrange_tetra(cst_mkI degree,
			   cst_mkI n,
			   mkR r,
			   cst_mkI roff_,
			   cst_mkR p,
			   cst_mkI  poff_,
			   mkR rwork,
			   cst_mkI rwork_n,
			   mkI err_)
{
  nsINT info,nshapes = mkS_tetra_n(degree);
  nsINT perm[nshapes];
  mkS_lagrange_localspl_tetra		(degree,&rwork[0],&nshapes);
  mkS_canonic_tetra(degree,&nshapes,&rwork[3*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
  mkS_dz_canonic_tetra(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
#if 1
  nslapack_dgesv(&nshapes,n,&rwork[3*nshapes],&nshapes,perm,r,roff_,&info);
#else
  char A1[2] = {'A','\0'};
  char A2[2] = {'A','\0'};
  nsREAL AS[36],AU[36*36],AV[36*36],rwork2[28*5];nsINT rwork2_n=28*5;
  nslapack_dgesvd(A1,A2,&nshapes,&nshapes,&rwork[2*nshapes],&nshapes,AS,AU,&nshapes,AV,&nshapes,rwork2,&rwork2_n,&info);
  { nsINT i,j;for (i=0;i<nshapes;++i){ nsREAL y =AS[i];if (nsFABS(AS[i])>((nsREAL)1.0e-15)) { y=((nsREAL)1.0)/AS[i]; for (j=0;j<nshapes;++j) {AV[i*nshapes+j] *=y; }   }else {for (j=0;j<nshapes;++j) {AV[i*nshapes+j] =(nsREAL)0.0;} } } }
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AU,&nshapes,r,roff_,&__vmps_blas_regal0,rwork2,&nshapes);
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AV,&nshapes,rwork2,&nshapes,&__vmps_blas_regal0,r,roff_);
#endif
}

void mkS_lagrange_tetra(cst_mkI  	degree,
			cst_mkI 	n,
			mkR 		r,
			cst_mkI 	roff_,
			cst_mkR 	p,
			cst_mkI 	poff_,
			mkR 		rwork,
			cst_mkI 	rwork_n,
			mkI 		err_)
{
  err_[0] = (nsINT)0;
  nsINT info,nshapes = mkS_tetra_n(degree);
  nsINT perm[nshapes];
  mkS_lagrange_localspl_tetra(degree,&rwork[0],&nshapes);
  mkS_canonic_tetra(degree,&nshapes,&rwork[3*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
  mkS_canonic_tetra(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
#if 1
  nslapack_dgesv(&nshapes,n,&rwork[3*nshapes],&nshapes,perm,r,roff_,&info);  
#else
  char A1[2] = {'A','\0'};
  char A2[2] = {'A','\0'};
  nsREAL AS[64],AU[64*64],AV[64*64],rwork2[2048];nsINT rwork2_n=2048;
  nslapack_dgesvd(A1,A2,&nshapes,&nshapes,&rwork[2*nshapes],&nshapes,AS,AU,&nshapes,AV,&nshapes,rwork2,&rwork2_n,&info);
  { nsINT i,j;for (i=0;i<nshapes;++i){ nsREAL y =AS[i];if (nsFABS(AS[i])>((nsREAL)1.0e-15)) { y=((nsREAL)1.0)/AS[i]; for (j=0;j<nshapes;++j) {AV[j*nshapes+i] *=y; }   }else {for (j=0;j<nshapes;++j) {AV[j*nshapes+i] =(nsREAL)0.0;} } } }
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AU,&nshapes,r,roff_,&__vmps_blas_regal0,rwork2,&nshapes);
  nsblas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AV,&nshapes,rwork2,&nshapes,&__vmps_blas_regal0,r,roff_);
#endif

}
