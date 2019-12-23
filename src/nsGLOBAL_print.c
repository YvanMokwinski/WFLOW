#include "nsGLOBAL.h"
#include "ns_constantes.h"
#include "ns_config_lapack.h"
#include "Monitor.h"
#if 0
Err nsGLOBAL_print_dgsol	(const nsGLOBAL_ST * const G_,
				 const char * name_,
				 cst_ensBASIS shape_,
				 cst_pR x_,
				 cst_pI xoff_)
{
  ns_mesh * mesh = (ns_mesh*)G_[0].mesh_usrptr;
  cst_eElement elm  =ns_mesh_elm(mesh);
  FILE * fich_ 	= fopen(name_,"w");
  const I subnelm=ensBASIS_nelmtreilli(shape_,elm);
  const I nelm = ns_mesh_nelm(mesh);
  fprintf(fich_,"2 1 "ifmt" 2\n",nelm*3*subnelm);
  R tmpdg_basis[4096];
  R lbasis[128];
  const I n = ensBASIS_n(shape_,elm);
  { 
#if 0
    I degree=ensBASIS_degree(shape_);
    I 	locrwork_n = 2048;
    R 	locrwork[2048];   
    Err err = __eErr_no;    
    cst_pR lcoo = nsTREILLI_TRIA_get_lc(shape_);
    nsSHAPE_l2ortho_tria(&degree,
			 &n,
			 tmpdg_basis,
			 &n,
			 lcoo,
			 &n,
			 locrwork,
			 &locrwork_n,
			 &G_[0].ferr);
#else
    Monitor_errmsg(G_->iproc,"nsGLOBAL_print_dgsol:SHAPE");
#endif
    
  }
  { I i;
    cst_pI split = ensBASIS_cnctreilli(shape_,elm);
    for (i=0;i<nelm;++i)
      { 
	nsblas_dgemv(transT,
		     &n,
		     &n,
		     &regal1,
		     tmpdg_basis,
		     &n,
		     &x_[i*xoff_[0]],
		     &negal1,
		     &regal0,
		     lbasis,
		     &negal1);
	{ I k;
	  for (k=0;k<subnelm;++k)
	    {
	      fprintf(fich_,"" rfmt "\n" rfmt "\n" rfmt "\n",
		      lbasis[split[3*k+0]],
		      lbasis[split[3*k+1]],
		      lbasis[split[3*k+2]]);	    
	    } }



      } }
  fclose(fich_);
  return __eErr_no;
}



void nsGLOBAL_dg_print(const nsGLOBAL_ST * const G_,
		       const char * 	name_,		 
		       cst_ensBASIS shape_,
		       cst_pR 	x_,
		       cst_pI 	xoff_)
{

  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(G_);
  char ctmp[512];
  /*  Err err 	= __eErr_no;*/
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".mesh",name_,gTimeInfo->itime);
  /*  err 		= */nsGLOBAL_print_dgspace(G_,ctmp,shape_);
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".bb",name_,gTimeInfo->itime);	
  /*err 		= */nsGLOBAL_print_dgsol(G_,ctmp,shape_,x_,xoff_);
}



void nsGLOBAL_dg_print_sol(const nsGLOBAL_ST * const G_,
			   const char * 	name_,		 	 
			   cst_ensBASIS 	shape_,
			   cst_pR 		x_,
			   cst_pI 		xoff_)
{
  char ctmp[512];
  /*  Err err = __eErr_no;*/
  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(G_);
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".bb",name_,gTimeInfo->itime);	
  /*err = */nsGLOBAL_print_dgsol(G_,ctmp,shape_,x_,xoff_);
}

void nsGLOBAL_dg_print_sol1(const nsGLOBAL_ST * const G_,
		      const char * 	name_,		 
		      const ensBASIS shape_,
		     cst_pR 		x_,
		     cst_pI 		xoff_)
{
  char ctmp[512];
  /*  Err err = __eErr_no;*/
  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(G_);
  sprintf(ctmp,"integrator.%s."ifmt".bb",name_,gTimeInfo->itime);	
  /*err = */nsGLOBAL_print_dgsol(G_,ctmp,shape_,x_,xoff_);
}









#endif

Err nsGLOBAL_print_dgspace(const nsGLOBAL_ST * const G_,
			   const char * 	name_,
			   const ensBASIS shape_)
{

  FILE * fich_ = fopen(name_,"w");
  R glcoo[32];
  ns_mesh * mesh = (ns_mesh*)G_[0].mesh_usrptr;
  cst_eFace elm  =ns_mesh_elm(mesh);
  Err err	= __eErr_no;
  /*  const I degree=ensBASIS_degree(shape_);*/
  const I subnelm=ensBASIS_nelmtreilli(shape_,elm);
  const I n = ensBASIS_n(shape_,elm);
  cst_pR lcoo = ensBASIS_cootreilli(shape_,elm);
  fprintf(fich_,"MeshVersionFormatted 1\nDimension\n2\n");
  const I nelm = ns_mesh_nelm( ((ns_mesh*)G_[0].mesh_usrptr) );
  fprintf(fich_,"\nVertices\n" ifmt "\n",subnelm*nelm*3);
  cst_pI split = ensBASIS_cnctreilli(shape_,elm);
  R elmcoo[6];
  { I j; 
    for (j=0;j<nelm;++j)
      {	
	elmcoo[0] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+0]-1)+0];
	elmcoo[3] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+0]-1)+1];
	elmcoo[1] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+1]-1)+0];
	elmcoo[4] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+1]-1)+1];
	elmcoo[2] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+2]-1)+0];
	elmcoo[5] 	= G_[0].ddlcoo[2*(G_[0].ddlcnc[6*j+2]-1)+1];
	{ I k;
	  for (k=0;k<n;++k)
	    {
	      glcoo[k] 		= (((R)1.0)-lcoo[k]-lcoo[n+k])*elmcoo[0]+lcoo[k]*elmcoo[1]+lcoo[k+n]*elmcoo[2];
	      glcoo[n+k] 	= (((R)1.0)-lcoo[k]-lcoo[n+k])*elmcoo[3]+lcoo[k]*elmcoo[4]+lcoo[k+n]*elmcoo[5];
	    } }       
	{ I k;
	  for (k=0;k<subnelm;++k)
	    {
	      fprintf(fich_,"%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n%8.15"nsFORMAT_REAL" %8.15"nsFORMAT_REAL" 0\n",
		      glcoo[split[3*k+0]],
		      glcoo[n+split[3*k+0]],
		      glcoo[split[3*k+1]],
		      glcoo[n+split[3*k+1]],
		      glcoo[split[3*k+2]],
		      glcoo[n+split[3*k+2]]);		      
	    } }	
      } }
  fprintf(fich_,"\nTriangles\n" ifmt "\n",subnelm*nelm);      
  { I j;
    for (j=0;j<nelm*subnelm;++j)
      {	
	fprintf(fich_,"" ifmt " " ifmt " " ifmt " 0\n",j*3+1,j*3+2,j*3+3);	    
       } }
  fprintf(fich_,"\nEnd\n");
  fclose(fich_);
  return  err;
}



void nsGLOBAL_dg_print_mesh_interval(const nsGLOBAL_ST * const 	G_,
				     const char * 		name_,		 
				     const ensBASIS 		shape_)
{
  char ctmp[512];
  /*  Err err = __eErr_no;*/
  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(G_);
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".mesh",name_,gTimeInfo->itime_interval);
  /*err = */nsGLOBAL_print_dgspace(G_,ctmp,shape_);
}
