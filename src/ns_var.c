#include "nsSPACE.h"
#include "ns_sys.h"
#include "ns_var.h"


#include "ns_constantes.h"
#include "ns_config_lapack.h"
#ifndef NDEBUG
#define wrong_param_ns_var(_cond,_i,_s) if ( (_cond) ) { fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); exit(1);} ((void)0)
#else
#define wrong_param_ns_var(_cond,_i,_s) ((void)0)
#endif
#ifndef NDEBUG

cst_pR 	cst_ns_var_x	(const ns_var*const F_) { wrong_param_ns_var(NOT F_,1,cst_ns_var_x); return F_->x;}
cst_pR 	cst_ns_var_y	(const ns_var*const F_) { wrong_param_ns_var(NOT F_,1,cst_ns_var_y); return F_->y;}

#endif

void ns_var_def(ns_var * 		v_,
		const ns_mesh * 	mesh_,
		const char * 		name_,
		const ensBASIS 	shape_,
		const L            is_discont_,
		cst_pI                 ncomp_,
		cst_pR 		t_,
		cst_pR 		ti_,
		cst_pR 		tii_,
		pR 			x_,
		pR 			xi_,
		pR 			xii_)
{
#ifndef NDEBUG
  wrong_param_ns_var(NOT v_,1,ns_var_def);
  wrong_param_ns_var(NOT mesh_,2,ns_var_def);
  wrong_param_ns_var(NOT name_,3,ns_var_def);
  wrong_param_ns_var(NOT ncomp_,6,ns_var_def);
#endif
  memset(v_,0,sizeof(ns_var));
  v_->lagr_id		= shape_;
  v_->mesh      	= mesh_;
  v_->t 		= t_;
  v_->ti 		= ti_;
  v_->tii 		= tii_;
  v_->ncomp     	= ncomp_[0];
  v_->nddlelm   	= ensBASIS_n(shape_,mesh_->elm);  
  v_->cnddl  		= ns_mesh_nddlspace(mesh_,shape_);
  v_->nddl 		= (is_discont_)?(v_->nddlelm*mesh_->nelm):v_->cnddl;
#if 0
  I ncomp_from_tii = 1;
  if (ti_)
    ncomp_from_tii=2;
  if (tii_)
    ncomp_from_tii=3;
#endif
  v_->x 			= x_;
  if (ncomp_[0]>1)
    {
      v_->y                 	= &x_[v_->nddl];
      if (ncomp_[0]>2)
	{
	  v_->z             	= &x_[v_->nddl*2];
	}


#if 0
      if (ncomp_[0]>3)
	{
	  ns_err("ns_var_def_malloc max comp is reached (>3)");
	}
#endif

    }
  if (ti_)
    {
      v_->xi 			= xi_;
      if (ncomp_[0]>1)
	{
	  v_->yi                 	= &xi_[v_->nddl];
	  if (ncomp_[0]>2)
	    {
	      v_->zi             	= &xi_[v_->nddl*2];
	    }
#if 0
	  if (ncomp_[0]>3)
	    {
	      ns_err("ns_var_def_malloc max comp is reached (>3)");
	    }
#endif
	}
    }
  if (tii_)
    {
      v_->xii 				= xii_;
      if (ncomp_[0]>1)
	{
	  v_->yii                 	= &xii_[v_->nddl];
	  if (ncomp_[0]>2)
	    {
	      v_->zii             	= &xii_[v_->nddl*2];
	    }
#if 0
	  if (ncomp_[0]>3)
	    {
	      ns_err("ns_var_def_malloc max comp is reached (>3)");
	    }
#endif
	}  
    }
  v_->is_discont 	= is_discont_;
  v_->name[0] 		= name_[0];
  v_->space 		= ns_mesh_get_space(v_->mesh,v_->lagr_id);
#if 0
#ifndef NDEBUG
  if (NOT v_->space)
    {
      ns_err("space not available %d",v_->lagr_id);
    }
#endif
#endif
}


void ns_var_welm(const ns_var*	self_,
		 cst_pI 	ielm_,
		 pR 		w_,
		 cst_pI 	woff_)
{
#ifndef NDEBUG
  wrong_param_ns_var(woff_[0]<self_->nddlelm,4,ns_var_welm);
  DebugVerif(ielm_[0]<self_->mesh->nelm);
#endif
  if (self_->is_discont)
    {
      { I i;
	for (i=0;i<self_->ncomp;++i)
	  {
	    { I j;
	      for (j=0;j<self_->nddlelm;++j)
		{
		  w_[woff_[0]*i+j] = self_->x[i*self_->nddl+ielm_[0]*self_->nddlelm+j];
		} }
	  } }
    }
  else
    {
      I cncelm[128];
      nsSPACE_cncelm(self_->space,ielm_,cncelm);
      { I i;
	for (i=0;i<self_->ncomp;++i)
	  {
	    { I j;
	      for (j=0;j<self_->nddlelm;++j)
		{
		  w_[woff_[0]*i+j] = self_->x[i*self_->nddl+cncelm[j]];
		} }
	  } }
    }
}


void ns_var_zeroi(ns_var * v_)
{
#ifndef NDEBUG
  wrong_param_ns_var(NOT v_,1,ns_var_zeroi);
  wrong_param_ns_var(NOT v_->xi,1,ns_var_zeroi);
#endif
  const I N = v_->nddl;
  { I i = 0;
    for (i=0;i<N;++i)
      {
	v_->xi[i]=((R)0.0);
      } }
  if (v_->ncomp>1)
    {
#ifndef NDEBUG
      wrong_param_ns_var(NOT v_->yi,1,ns_var_zeroi);
#endif
  { I i = 0;
    for (i=0;i<N;++i)
      {
	v_->yi[i]=((R)0.0);
      } }
    }
  if (v_->ncomp>2)
    {
#ifndef NDEBUG
      wrong_param_ns_var(NOT v_->zi,1,ns_var_zeroi);
#endif
      { I i = 0;
	for (i=0;i<N;++i)
	  {
	    v_->zi[i]=((R)0.0);
	  } }
    }
}

void ns_var_print(const ns_var * 	v_,
		  const char * 		name_,
		  cst_pI 		itime_interval_,
		  cst_pI 		itime_)
{
  char ctmp[512];
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".bb",name_,itime_[0]);
  FILE*fil = fopen(ctmp,"w");				
  fprintf(fil,"2 "ifmt" "ifmt" 2\n",v_->ncomp,v_->nddl);		
  if (v_->ncomp==1)
    {
      { I i;
	for (i=0;i<v_->nddl;++i)
	  { 	
	    fprintf(fil,""rfmt"\n",v_->x[i]); 
	  } }   
    }
  else if (v_->ncomp>1)
    {
      if (v_->ncomp!=2)
	{
	  fprintf(stderr,"*** ERREUR ns_var_print ncomp!=dimension\n");
	}
      { I i;
	for (i=0;i<v_->nddl;++i)
	  { 	
	    fprintf(fil,""rfmt" "rfmt"\n",v_->x[i],v_->y[i]); 
#if 0
	    else if (dimension==3)
	      fprintf(fil,""rfmt" "rfmt" "rfmt"\n",v_->x[i],v_->y[i],v_->z[i]); 
#endif
	  } }   
    }
  else
    {
      fprintf(stderr,"ns_var_print(%s):ncomp=  "ifmt"\n",name_,v_->ncomp);
      exit(1);
    }
  fclose(fil);							
  sprintf(ctmp,"integrator.%s.%.5"nsFORMAT_INTEGER".timeinfo.txt",name_,itime_[0]);	
  fil = fopen(ctmp,"w");
  fprintf(fil,""ifmt" "ifmt" "rfmt"\n",itime_interval_[0],itime_[0],v_->t[0]);
  fclose(fil);
}

void ns_var_cpy_xi2x(ns_var * self_)
{  
  nsblas_dcopy(&self_->nddl,self_->xi,&negal1,self_->x,&negal1);
  if (self_->ncomp>1)
    {
      nsblas_dcopy(&self_->nddl,self_->yi,&negal1,self_->y,&negal1);
      if (self_->ncomp>2)	
	nsblas_dcopy(&self_->nddl,self_->zi,&negal1,self_->z,&negal1);
    }
}

void nsISO_INFO_def(nsISO_INFO 			iso_info_,
		    cst_ens_method_capturing 	method_,		    
		    const R 			iso_,
		    const R 			tol_,
		    const I 			lev_) 
{
  
  iso_info_ 		= (nsISO_INFO)memset(iso_info_,0,sizeof(nsISO_INFO_ST));
  iso_info_->method	= method_;
  iso_info_->iso 	= iso_;
  iso_info_->tol	= tol_;
  iso_info_->lev	= (lev_>4)?4:lev_;
#if 0
  if (lev_>4)
    ns_warn("nsISO_INFO_def:max lev_ is 4");
#endif
}

void nsOP_INFO_def(nsOP_INFO 			operator_,
		    cst_eOperator 		op_,
		    cst_eSmoothingMethod 	qh_method_,
		    const I 		qh_degree_)
{
  operator_=(nsOP_INFO)memset(operator_,0,sizeof(nsOP_INFO_ST));  
  operator_->op 		= op_;
  operator_->smoothing_method 	= qh_method_;
  operator_->qh_degree 		= qh_degree_;
}

void CsfInfo_def(pCsfInfo 			csf_,
		 const L			color_,
		 cst_eOperator		f_op_,
		 cst_eOperator		nabla_op_,
		 cst_eOperator		normal_op_,
		 const R 			normal_tol_,
		 cst_eOperator		kappa_op_,
		 const R 			kappa_tol_,
		 cst_eOperator		tension_op_,
		 cst_eOperator		distance_op_,
		 const R 			distance_iso_,
		 cst_ens_method_capturing	kind_curve_,
		 const I 		distance_level_,
		 cst_eSmoothingMethod 	qh_method_,
		 const I 		qh_degree_,
		 const eHeaviside	heaviside_,
		 const eDirac		dirac_,
		 cst_pR			eps_,
		 pErr 			err_)
{
  csf_=(pCsfInfo)memset(csf_,0,sizeof(CsfInfo));
  csf_->rinfo[__ensCSF_rinfo_distance_iso] 	= distance_iso_;
  csf_->rinfo[__ensCSF_rinfo_normal_tol] 	= normal_tol_;
  csf_->rinfo[__ensCSF_rinfo_kappa_tol] 	= kappa_tol_;
  csf_->iinfo[__ensCSF_iinfo_kappa_filtering] 	= 1;
  csf_->color 					= color_;
  nsOP_INFO_def		(&csf_->op_f,
			 f_op_,
			 qh_method_,
			 qh_degree_);
  nsOP_INFO_def		(&csf_->op_kappa,
			 kappa_op_,
			 qh_method_,
			 qh_degree_);
  nsOP_INFO_def		(&csf_->op_normal,
			 normal_op_,
			 qh_method_,
			 qh_degree_);
  nsOP_INFO_def		(&csf_->op_nabla,
			 nabla_op_,
			 qh_method_,
			 qh_degree_);
  nsOP_INFO_def		(&csf_->op_tension,
			 tension_op_,
			 qh_method_,
			 qh_degree_);
  nsOP_INFO_def		(&csf_->op_distance,
			 distance_op_,
			 qh_method_,
			 qh_degree_);
  nsISO_INFO_def	(&csf_->iso_info,
			 kind_curve_,
			 distance_iso_,
			 ((R)0.0001),
			 distance_level_);

  
  SmoothedDirac_def	(&csf_->m_smoothedDirac,
			 dirac_,
			 eps_,
			 err_);

  
  SmoothedHeaviside_def	(&csf_->m_smoothedHeaviside,
			 heaviside_,
			 eps_,
			 err_);

}


void ns_mesh_free(ns_mesh * mesh_)
{
  if (mesh_)
    {
      if (mesh_->time_x)
	free(mesh_->time_x);
      if (mesh_->own_cnc)
	free(mesh_->own_cnc);
      if (mesh_->own_cod)
	free(mesh_->own_cod);
      if (mesh_->own_adj)
	free(mesh_->own_adj);
      if (mesh_->own_bpatch)
	free(mesh_->own_bpatch);
      if (mesh_->own_patch)
	free(mesh_->own_patch);
      if (mesh_->own_coo)
	free(mesh_->own_coo);
      if (mesh_->own_jacelm)
	free(mesh_->own_jacelm);
      if (mesh_->own_jacedge)
	free(mesh_->own_jacedge);
      if (mesh_->own_normaledge)
	free(mesh_->own_normaledge);
      if (mesh_->own_trelm)
	free(mesh_->own_trelm);
      {	ensBASIS i = __ensBASIS_ERROR;
	for (++i;i<__ensBASIS_ALL;++i)
	  {
	    mesh_->spaces[i]=nsSPACE_kill((pSpace)mesh_->spaces[i]);
	  } }
      {	I j;
	for (j=0;j<64;++j)
	  {
	    mesh_->var_scalar[j] = ns_var_kill((ns_var*)mesh_->var_scalar[j]);
	  } }
      {	I j;
	for (j=0;j<64;++j)
	  {
	    mesh_->var_tensor[j] = ns_var_kill((ns_var*)mesh_->var_tensor[j]);
	  } }
      {	I j;
	for (j=0;j<64;++j)
	  {
	    mesh_->var_vector[j] = ns_var_kill((ns_var*)mesh_->var_vector[j]);
	  } }
#if 0
      {	I j;
	for (j=0;j<64;++j)
	  {
	    mesh_->metrics[j] = nsMETRIC_kill(mesh_->metrics[j]);
	  } }      
#endif
#if 0
      {	I j;
	for (j=0;j<16;++j)
	  {
	    mesh_->instances_cubature_triangle[j] = nsCUBA_TRIA_kill(mesh_->instances_cubature_triangle[j]);
	  } }
      {	I j;
	for (j=0;j<16;++j)
	  {
	    mesh_->instances_cubature_triangle_interval[j] = nsCUBA_TRIA_INTERVAL_kill(mesh_->instances_cubature_triangle_interval[j]);
	  } }
#endif
      memset(mesh_,0,sizeof(ns_mesh));
    }
}



void ns_var_free(ns_var*v_)
{
  if (v_)
    {
      if (v_->x)
	{
	  free(v_->x);
	}
      memset(v_,0,sizeof(ns_var));
    }
}


ns_var*ns_var_kill(ns_var*v_)
{
  if (v_)
    {
      ns_var_free(v_);
      free(v_);
    }
  return NULL;
}
