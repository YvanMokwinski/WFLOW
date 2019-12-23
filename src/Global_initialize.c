#include "nsGLOBAL.h"
#include "ns_constantes.h"
#include "ns_config.h"
#include "Monitor.h"

#if 0
#include "nsDG.h"
#include "nsNONLINEAR_STOKES_GLOBAL.h"
#endif

#if 0

typedef struct
{
  I 		iinfo	[DG_I_n];
  R		rres	[DG_rres_n];
  I			ires	[DG_ires_n];
  I 		rinfo_n;
  pR 			rinfo;
  I 		rworkdg_n;
  I           	iworkdg_n;
  pR			rworkdg;
  pI			iworkdg;
  R 		boundary_basis[21*10*3];
  R 		qface_w[256];
} ns_dg_op;

ns_dg_op dg_op[__NS_MAXNUM_THREADS__];

#define ns_dg_dcpl_iinfo 	dg_op[iproc_].iinfo
#define ns_dg_dcpl_rres 	dg_op[iproc_].rres
#define ns_dg_dcpl_ires 	dg_op[iproc_].ires
#define ns_dg_dcpl_rinfo_n 	dg_op[iproc_].rinfo_n
#define ns_dg_dcpl_rinfo 	dg_op[iproc_].rinfo

#define ns_dg_dcpl_rworkdg_n 	dg_op[iproc_].rworkdg_n
#define ns_dg_dcpl_rworkdg 	dg_op[iproc_].rworkdg
#define ns_dg_dcpl_iworkdg_n 	dg_op[iproc_].iworkdg_n
#define ns_dg_dcpl_iworkdg 	dg_op[iproc_].iworkdg

I ns_dg_fluxelm_boundary_get_coo(nsGLOBAL_ST*const  	G_,const I jelm,const I jadj,pR pos_)
{
  cst_pR lc 	= &dg_op[0].rinfo[dg_op[0].iinfo[DG_IA_lc_face]];
  const I N = dg_op[0].iinfo[DG_I_QFACE_N];
  memcpy(pos_,&lc[N*jadj],sizeof(R)*N);
  memcpy(pos_+N,&lc[dg_op[0].iinfo[DG_I_nTot]+N*jadj],sizeof(R)*N);
  R cc[32];
  ns_mesh_cooelm(G_[0].mesh_usrptr,&jelm,cc);
  { I i;
    for (i=0;i<N;++i)
      {
	const R r = pos_[i];
	const R s = pos_[N+i];
	pos_[i] = cc[0] * (((R)1.0)-(r+s)) + cc[1]*r+cc[2]*s;
	pos_[N+i] = cc[3] * (((R)1.0)-(r+s)) + cc[4]*r+cc[5]*s;
      } }
#if 0
  { I i;
    for (i=0;i<N;++i)
      {
	printf("%e %e\n",pos_[i],pos_[N+i]);
      } }
  exit(1);
#endif
  return N;
}

void ns_dg_fluxelm_boundary(nsGLOBAL_ST*const  	G_,
			    cst_pR 		xa_,
			    const I 	jelm,
			    const I 	jadj,
			    pR 		rhs_,
			    cst_pR 		loc_,
			    cst_pR 		uvw_)
{

  const I jedge 	= G_[0].pointeur_maillage->cnc[6*jelm+3+jadj]-G_[0].pointeur_maillage->nvertex-1;  
  cst_pR lc 		= &dg_op[0].rinfo[dg_op[0].iinfo[DG_IA_lc]];
  const I N 	= dg_op[0].iinfo[DG_I_QFACE_N];
  const R nx 	= G_[0].pointeur_maillage->normaledge[2*jedge+0];
  const R ny 	= G_[0].pointeur_maillage->normaledge[2*jedge+1];
  R rwork[128];
  { I i;
    for (i=0;i<N;++i)
      {
	rwork[i] = nx*uvw_[i]+ny*uvw_[N+i];
	rwork[i] = -MIN(0.0,rwork[i]);
	rwork[i] *= loc_[i];	
#if 0
	printf("ddddddddddddd %e %e %e\n",uvw_[i],uvw_[N+i],rwork[i]);
#endif
      } }

  const R xa  =G_[0].pointeur_maillage->jacedge[jedge];
  nsblas_dgemv(transN,&G_[0].nddlelmdg,&N,&xa,&dg_op[0].boundary_basis[jadj*G_[0].nddlelmdg*N],&G_[0].nddlelmdg,rwork,&negal1,&regal1,&rhs_[jelm*G_[0].nddlelmdg],&negal1);

#if 0
  printf("%e %e\n",rhs_[jelm*G_[0].nddlelmdg+0],xa);
  printf("%e\n",rhs_[jelm*G_[0].nddlelmdg+1]);
  printf("%e\n",rhs_[jelm*G_[0].nddlelmdg+2]);
  exit(1);
#endif
#if 0
  nsSHAPE_bwji(qface_n_,qface_p_,qface_w_,nsSHAPE_b(&s_trial),nsSHAPE_b(&s_test),&bmat_x[(nA+2*nU)*bmat_n],bmat_n,ferr_);     
  nsSHAPE_bwji(qface_n_,qface_p_,qface_w_,nsSHAPE_b(&s_teta),nsSHAPE_b(&s_test),&bmat_x[(nA+2*nU)*bmat_n+nxn],bmat_n,ferr_);
#endif
}


void ns_dg_decoupled_init(nsGLOBAL_ST*const  	G_,
			  cst_pI dg_n_)
{
  const I iproc_ = G_[0].iproc;
  { I degreef=1;

    if (dg_n_[0]==1) 	
      degreef=0;
    else if (dg_n_[0]==3) 	degreef=1;
    else if (dg_n_[0]==6) 	degreef=2;
    else if (dg_n_[0]==10) 	degreef=3;
    else if (dg_n_[0]==15) 	degreef=4;
    else if (dg_n_[0]==21) 	degreef=5;
    else ns_err("ns_dg_decoupled_init:wrong dg_n_");
    nsSHAPE_ST teta_a,shape_F,shape_U;
    Err err	= __eErr_no;
    nsSHAPE_definit	(&teta_a,__eFace_TRIANGLE,__ensSHAPE_FAMILY_lagrange,0,__ens_discontinuous,&err);		        
    nsSHAPE_definit	(&shape_U,__eFace_TRIANGLE,__ensSHAPE_FAMILY_lagrange,2,__ens_discontinuous,&err);
    nsSHAPE_definit	(&shape_F,__eFace_TRIANGLE,__ensSHAPE_FAMILY_l2orthonormal,degreef,__ens_discontinuous,&err);
    ns_dg_dcpl_iinfo[DG_I_QELM_N]	= 7;
    ns_dg_dcpl_iinfo[DG_I_QFACE_N]  	= 7;  
    ns_dg_dcpl_rinfo_n 			= (I)64000;
    ns_dg_dcpl_rinfo 			= rmalloc(ns_dg_dcpl_rinfo_n);
    ns_dg_dcpl_iworkdg_n 		= (I)128000;
    ns_dg_dcpl_rworkdg_n 		= (I)128000;
    ns_dg_dcpl_rworkdg 			= rmalloc(ns_dg_dcpl_rworkdg_n);
    ns_dg_dcpl_iworkdg 			= imalloc(ns_dg_dcpl_iworkdg_n);
    ns_dg_init_	(&teta_a,
		 &shape_U,
		 &shape_F,
		 &shape_F,
		 &shape_F,
		 ns_dg_dcpl_iinfo,
		 &ns_dg_dcpl_rinfo_n,
		 ns_dg_dcpl_rinfo,
		 &ns_dg_dcpl_rworkdg_n,
		 ns_dg_dcpl_rworkdg,
		 dg_op[G_[0].iproc].qface_w); 
#if 0
    

    I fferr;
    const I N 	= dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N]*3;
    cst_pR lc 		= &dg_op[G_[0].iproc].rinfo[dg_op[G_[0].iproc].iinfo[DG_IA_lc_face]];
#if 0
    I i;
    for (i=0;i<N;++i)
      {
	printf("%e %e\n",lc[dg_op[0].iinfo[DG_IA_lc_face]+i],lc[dg_op[0].iinfo[DG_IA_lc_face]+i+dg_op[0].iinfo[DG_I_nTot]]);
      }
    exit(1);
#endif
    I i;
    nsSHAPE_basis	(nsSHAPE_b(&shape_F),
			 &N,
			 dg_op[G_[0].iproc].boundary_basis,
			 &G_[0].nddlelmdg,
			 lc,
			 &dg_op[G_[0].iproc].iinfo[DG_I_nTot],			   			   
			 ns_dg_dcpl_rworkdg,
			 &ns_dg_dcpl_rworkdg_n,
			 &fferr); 
    for (i=0;i<dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N];++i)
      dscal(&G_[0].nddlelmdg,&dg_op[G_[0].iproc].qface_w[i],&dg_op[G_[0].iproc].boundary_basis[ (  i   ) *G_[0].nddlelmdg ],&negal1);
    for (i=0;i<dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N];++i)
      dscal(&G_[0].nddlelmdg,&dg_op[G_[0].iproc].qface_w[i],&dg_op[G_[0].iproc].boundary_basis[ (  i +dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N]  ) *G_[0].nddlelmdg ],&negal1);
    for (i=0;i<dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N];++i)
      dscal(&G_[0].nddlelmdg,&dg_op[G_[0].iproc].qface_w[i],&dg_op[G_[0].iproc].boundary_basis[ (  i+2*dg_op[G_[0].iproc].iinfo[DG_I_QFACE_N]   ) *G_[0].nddlelmdg ],&negal1);
#endif
  }
}

void ns_dg_decoupled	(nsGLOBAL_ST*const  	G_,
			 cst_ensMETHOD_TRANSIENT scheme_,
			 cst_pR t_,
			 cst_pR ti_,
			 cst_pR tii_,
			 cst_pR u_,
			 cst_pR ui_,
			 cst_pR uii_,
			 cst_pR v_,
			 cst_pR vi_,
			 cst_pR vii_,
			 cst_pI dg_n_,
			 pR     dg_f_,
			 cst_pI dg_f_off_,
			 pR     dg_fi_,
			 cst_pI dg_fi_off_,
			 pR     dg_fii_,
			 cst_pI dg_fii_off_,
			 cst_pR dg_rhs_,
			 cst_pI dg_rhsoff_,
			 pR     dg_c_,
			 cst_pI dg_c_off_,
			 pI     transition_detection_)
{ 
  nsGLOBAL_ST*constZ=G_;
  const I iproc_ = Z->iproc;

  extern const R 	regal0;
  extern const R 	regal1;
  extern const I 	negal1;


  const R aii 	= ti_[0]-tii_[0];
  const R ai 	= t_[0]-ti_[0];
  const I N 	= dg_n_[0] * Z->pointeur_maillage->nelm;
  R coeff_a 	= ((R)0.0);
  R coeff_u 	= ((R)1.0);
  R coeff_a_rhs_i 	= ((R)0.0);
  R coeff_a_rhs_ii = ((R)0.0);
#if 0
  { I i;
    for (i=0;i<Z->pointeur_maillage->nelm;++i)
      {
	{ I j;
	  for (j=0;j<dg_n_[0];++j)
	    {	
	      dg_rhs_[i*dg_rhsoff_[0]+j] = ((R)0.0);
	    } }
      } }
#endif

  switch(scheme_)
    {
    case __ensMETHOD_TRANSIENT_euler:
      {
	coeff_a 		= ((R)1.0);
	coeff_u 		= ai;
	if( (dg_fi_off_[0]==dg_rhsoff_[0]) AND (dg_fi_off_[0]==dg_n_[0]) )
	  {
	    nsblas_daxpy(&N,&regal1,dg_fi_,&negal1,dg_rhs_,&negal1);
	    { I i;
	      for (i=0;i<Z->pointeur_maillage->nelm;++i)
		{
		  nsblas_dscal(dg_n_,&Z->jacelm[i],&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
		} }      	  
	  }
	else
	  {
	    { I i;
	      for (i=0;i<Z->pointeur_maillage->nelm;++i)
		{
		  nsblas_daxpy(dg_n_,&Z->jacelm[i],&dg_fi_[i*dg_fi_off_[0]],&negal1,&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
		} }      
	  }     
	break;
      }

    case __ensMETHOD_TRANSIENT_imr:
      {
	coeff_a    	= ((R)1.0);
      coeff_u 		= ai*((R)0.5);
      if( (dg_fi_off_[0]==dg_rhsoff_[0])AND (dg_fi_off_[0]==dg_n_[0]) )
	{
	  nsblas_daxpy(&N,&regal1,dg_fi_,&negal1,dg_rhs_,&negal1);
	  { const I nelm = Z->nelm;
	    { I i;
	      for (i=0;i<nelm;++i)
		{
		  nsblas_dscal(dg_n_,&Z->jacelm[i],&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
		} }  }    	  
	}
      else
	{
	  {const I nelm = Z->nelm;
	    { I i;
	      for (i=0;i<nelm;++i)
		{
		  nsblas_daxpy(dg_n_,&Z->jacelm[i],&dg_fi_[i*dg_fi_off_[0]],&negal1,&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
		} }  }    
	}     
      break;
    }

    case __ensMETHOD_TRANSIENT_gear_euler:
    case __ensMETHOD_TRANSIENT_gear_imr:
      {
	
	const R coeff_a_rhs_i 	= (((R)2.0)+aii/ai+ai/aii)/(((R)2.0)+aii/ai);
	const R coeff_a_rhs_ii 	= -(ai/aii)/(((R)2.0)+aii/ai);
      coeff_a 				= ((R)1.0);
      coeff_u 				= (ai+aii)/(((R)2.0)+aii/ai);
      nsblas_daxpy(&N,&coeff_a_rhs_i,dg_fi_,&negal1,dg_rhs_,&negal1);
      nsblas_daxpy(&N,&coeff_a_rhs_ii,dg_fii_,&negal1,dg_rhs_,&negal1);
      { I i;
	for (i=0;i<Z->pointeur_maillage->nelm;++i)
	  {
	    nsblas_dscal(dg_n_,&Z->jacelm[i],&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	  } }      	  
#if 0
      coeff_a 				= ((R)1.0);
      coeff_u 				= (ai+aii)/(((R)2.0)+aii/ai);
      const R coeff_a_rhs_i 	= (((R)2.0)+aii/ai+ai/aii)/(((R)2.0)+aii/ai);
      const R coeff_a_rhs_ii 	= (ai/aii)/(((R)2.0)+aii/ai);
      if( (dg_fi_off_[0]==dg_rhsoff_[0])AND(dg_fii_off_[0]==dg_rhsoff_[0])AND (dg_fi_off_[0]==dg_n_[0]) )
	{
	  nsblas_daxpy(&N,&coeff_a_rhs_i,dg_fi_,&negal1,dg_rhs_,&negal1);
	  nsblas_daxpy(&N,&coeff_a_rhs_ii,dg_fii_,&negal1,dg_rhs_,&negal1);
	  { I i;
	    for (i=0;i<Z->pointeur_maillage->nelm;++i)
	      {
		nsblas_dscal(dg_n_,&Z->jacelm[i],&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	      } }      	  
	}
      else
	{
	  { I i;
	    for (i=0;i<Z->pointeur_maillage->nelm;++i)
	      {
		const R a = coeff_a_rhs_i * Z->jacelm[i];
		nsblas_daxpy(dg_n_,&a,&dg_fi_[i*dg_fi_off_[0]],&negal1,&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	      } }
	  { I i;
	    for (i=0;i<Z->pointeur_maillage->nelm;++i)
	      {
		const R a = coeff_a_rhs_ii * Z->jacelm[i];
		nsblas_daxpy(dg_n_,&a,&dg_fii_[i*dg_fii_off_[0]],&negal1,&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	      } }          
	} 
#endif
      break;
      }
      case __ensMETHOD_TRANSIENT_trapeze:
    {

      coeff_a 		= ((R)2.0);
      coeff_u 		= ai*((R)0.5);
      if( (dg_fi_off_[0]==dg_rhsoff_[0])AND (dg_fi_off_[0]==dg_n_[0]) )
	{
	  const I a  = ((R)2.0);
	  nsblas_daxpy(&N,&a,dg_fi_,&negal1,dg_rhs_,&negal1);
	  { I i;
	    for (i=0;i<Z->pointeur_maillage->nelm;++i)
	      {
		nsblas_dscal(dg_n_,&Z->jacelm[i],&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	      } }      	  
	}
      else
	{
	  { I i;
	    for (i=0;i<Z->pointeur_maillage->nelm;++i)
	      {
		const R a = ((R)2.0) * Z->jacelm[i];
		nsblas_daxpy(dg_n_,&a,&dg_fi_[i*dg_fi_off_[0]],&negal1,&dg_rhs_[i*dg_rhsoff_[0]],&negal1);
	      } }
	}
      /* ON RAJOUTE DANS LE RESIDU LE TERME DE CONVECTION AU TEMPS PRECEDENT *1/2 */
      /* ET LE TERME DE SAUT AU TEMPS PRECEDENT *1/2 */
      break;
    }

    case __ensMETHOD_TRANSIENT_no:
      {
	coeff_a		= ((R)0.0);
	coeff_u 		= ((R)1.0);
	coeff_a_rhs_i 	= ((R)0.0);
	coeff_a_rhs_ii 	= ((R)0.0);
	break;
      }
    case __ensMETHOD_TRANSIENT_n:
    case __ensMETHOD_TRANSIENT_error:
      {
	Monitor_errmsg(Z->iproc,"ns_dg_decoupled:wrong case");
	break;
      }  
    }
  /* compute correction */
#if 0
  printf("coeff a %e coeff u %e\n",coeff_a,coeff_u);
#endif
  const I dimp1=3;
  const I dim=2;
  ns_dg_routine_	(&coeff_a,
			 &coeff_u,
			 dg_rhs_,
			 dg_rhsoff_,
			 Z->ddlcnc,
			 &Z->nddlelmu,
			 u_,	
			 v_,	
			 dg_f_,	
			 dg_f_off_,
			 dg_c_,
			 dg_c_off_,
			 t_,
			 &Z->pointeur_maillage->nelm,
			 Z->ddlcoo,
			 &dim,
			 Z->ddlcnc,
			 &Z->nddlelmu,
			 Z->adj,
			 &dimp1,
			 Z->ddlcod,
			 &Z->params.iinfo[__ens_iinfo_noboundary_vcod],
			 &ns_dg_dcpl_rworkdg_n,
			 ns_dg_dcpl_rworkdg,
			 &ns_dg_dcpl_iworkdg_n,
			 ns_dg_dcpl_iworkdg,
			 ns_dg_dcpl_rinfo,
			 ns_dg_dcpl_iinfo,
			 ns_dg_dcpl_rres,
			 ns_dg_dcpl_ires);  

  
  if (NOT ns_dg_dcpl_ires[DG_ires_convergence])
    {
      Monitor_warn(Z->iproc,"DG:no convergence");
    }
  else
    {
      Monitor_msg(Z->iproc,"DG:niter "ifmt", max value %e",ns_dg_dcpl_ires[DG_ires_iter_gauss_seidel],ns_dg_dcpl_rres[DG_rres_max]);
    }


 /* add correction */  
  if (dg_f_off_[0]==dg_n_[0])
    {
      nsblas_daxpy(&N,
		   &regal1,
		   dg_c_,
		   &negal1,
		   dg_f_,
		   &negal1);    
    }
  else
    {
      { I i;
	for (i=0;i<Z->pointeur_maillage->nelm;++i)
	  {
	    nsblas_daxpy(dg_n_,
			 &regal1,
			 &dg_c_[dg_c_off_[0]*i],
			 &negal1,
			 &dg_f_[dg_f_off_[0]*i],
			 &negal1);    
	  } }
    }




}
#endif




void nsGLOBAL_linsys_clear(nsGLOBAL_ST*const self_)
{
  pParameters const parameters 	= Global_get_Parameters(self_);
  const L pressure_uncoupled 	= Parameters_getl(parameters,__ens_linfo_pressure_uncoupled);
  const L pspg 			= Parameters_getl(parameters,__ens_linfo_pspg);
  SparseBlock_clear(&self_->sparseStokes->block_UU);  
  if (NOT pressure_uncoupled)
    {
      if (pspg)
	{
#ifndef NDEBUG
	  Monitor_msg(self_->iproc,"debug:clear sparse_A_block_BB");
	  Monitor_msg(self_->iproc,"debug:clear sparse_A_block_BU");
#endif
	  SparseBlock_clear(&self_->sparseStokes->block_BB);
	  SparseBlock_clear(&self_->sparseStokes->block_BU);
	}        
    }
  else
    {
      /* 
	 uncoupled pressure
	 we do nothing
      */
    }
}

