#if 0
#include "DgTransport.h"

void mok_dg_compute_(cst_pR 	xu_,
		     cst_pR	rhs_,
		     cst_pI  	rhsoff_,
		     cst_pR 	data_u_,
		     cst_pI 	data_uoff_,
		     cst_pR	data_v_,
		     cst_pI 	data_voff_,
		     cst_pR	sol_,
		     cst_pI 	soloff_,
		     pR 		corr_,
		     cst_pI 	corroff_ ,

		     cst_pR 	t_,
#if 0		     
		     cst_mkZ   zone_,
#endif
#if 1
		     cst_pI 	nvertex_,
		     cst_pI 	nelm_,
		     cst_pR 	coo_,
		     cst_pI 	cooff_,
		     cst_pI 	cnc_,
		     cst_pI 	cncoff_,
		     cst_pI 	adj_,
		     cst_pI 	adjoff_,
		     cst_pI 	vcod_,
		     cst_pI 	noboundary_cod_,
#endif

		     cst_pI 	rwork_n_,
		     pR  		rwork_,
		     cst_pI 	iwork_n_,
		     pI		iwork_,
		     cst_pR 	rinfo_,
		     const I  	iinfo_[DG_I_n],
		     R 		rres_[DG_rres_n],
		     I 		ires_[DG_ires_n])
{  

#if 0
  const I 	nelm_[1] 		= {mkZ_nelm(zone_)};
  cst_pR 	coo_ 			= cst_mok_m_x(cst_mkZ_coo(zone_));
  cst_pI 	cooff_ 			= mok_m_ptroff(cst_mkZ_coo(zone_));
  cst_pI 	cnc_ 			= cst_pIM_x(cst_mkZ_cnc(zone_));
  cst_pI 	cncoff_ 		= pIM_ptroff(cst_mkZ_cnc(zone_));
  cst_pI 	adj_ 			= cst_pIM_x(cst_mkZ_adj(zone_));
  cst_pI 	adjoff_ 		= pIM_ptroff(cst_mkZ_adj(zone_));
  cst_pI 	vcod_ 			= cst_pIM_x(cst_mkZ_vertex_cod(zone_));
  const I 	noboundary_cod_ [1] 	= {mkZ_noboundary_vcod(zone_)};
#endif

  cst_pI	trial_n = &iinfo_[DG_I_TRIAL_NBASIS];
  cst_pI	test_n  = &iinfo_[DG_I_TEST_NBASIS];
  cst_pI	teta_n  = &iinfo_[DG_I_TETA_NBASIS];
  cst_pI	qface_n = &iinfo_[DG_I_QFACE_N];
  cst_pI	nu      = &iinfo_[DG_I_TETA_U_NBASIS];
  /*_____________________________________________________*/
  const R mr1=(R)-1.0,r1=(R)1.0,r0=(R)0.0;
  const I n1=(I)1, n2 = (I)2;
  /*_____________________________________________________*/
#if DG_CONSERVATIVE_WEAK_FORM
  const R mxu=xu_[0]*((R)-1.0);
#else
  const R mxu=xu_[0]*((R)1.0);
#endif

  /*_____________________________________________________*/
  pI lcperm  	= NULL;
  pI perm  	        = NULL;
  pI graph 		= NULL;
  pI marker 		= NULL;
  pI marker_begin 	= NULL;
  pI marker_flag 	= NULL;
  I marker_size	= 0;
  pR brhs         = NULL;
  pR brhs_uelm    = NULL;  
  pR brhs_a       = NULL;
  pR uface0 	= NULL;
  pR uface1 	= NULL;
  pR uface2 	= NULL;
  R mx;  
  I  id;
  I  ielm;
  I info_lapack;
  I iter_gauss_seidel;
  I b1,b2,pp,n;
  R jacelm[1],x;  
  R longueurs[3];
  I vcod[3*1],cnc[3*1];
  I neids[8*1],neids_face[8*1],codface[8*1];
  R jacface[3*1],nrmelm[6*1],cooelm[6*1];
  R belm[4*1],lcsol[21];
  const I q_nx3 		= qface_n[0]*3;
  const I bmat_n 		= iinfo_[DG_I_bmat_n];
  const I bmat_m 		= iinfo_[DG_I_bmat_m];  
  cst_pR eval_u_ 	= &rinfo_[iinfo_[DG_RA_EVAL_TETA_U]];
  const I eval_uoff_        = iinfo_[DG_I_TETA_U_NBASIS];
  cst_pR bmat 	 	= &rinfo_[iinfo_[DG_RA_bmat]];
  cst_pR bmatflux 	= &rinfo_[iinfo_[DG_RA_bmatx]];
  cst_pR fpart[3][3];
  { I i,j;
    for (i=0;i<3;++i)
      for (j=0;j<3;++j)
	fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
  }

#if __mk_debug__
  if (bmat_m!=1 + 2*nu[0]+3*qface_n[0]) 
    {
      fprintf(stderr,"*** DGERR "ifmt" "ifmt"",1 + 2*nu[0]+3*qface_n[0],bmat_m);exit(1);
    }
#endif
  I N_pts_u = bmat_m-1;/*nTot-nA*/
  
  R 	matflux[21*(21+21)];/*trial_n*(trial_n+teta_n)*/
  pR 		matflux_corr 	= &matflux[0];
  pR 		matflux_residu 	= &matflux[trial_n[0]*trial_n[0]];
  
  R 	lcmat[21*21*2];
  R 	lcrhs[21*2];
  pR 		lcmatc = &lcmat[0];
  pR 		lcmatr = &lcmat[trial_n[0]*trial_n[0]];

  const I renum_iwork_n = 3*nelm_[0] + 4*(nelm_[0]+1) + trial_n[0];
  ires_[DG_ires_err] = (I)0;
  
  /*_____________________________________________________*/
  if (rwork_n_[0]<bmat_m + 2*bmat_m)
    {      
      ires_[DG_ires_err] = (I)1;
      ires_[DG_ires_required_rw_n] = 3*bmat_m;
      fprintf(stderr,"too small rwork_n_ "ifmt" "ifmt"\n",3*bmat_m,rwork_n_[0]);
      return;
    }  

  brhs      	= &rwork_[0];
  brhs_a    	= &rwork_[0];  
  brhs_uelm 	= &rwork_[1];  
  uface0 	= &rwork_[iinfo_[DG_IA_lc_face]];
  uface1 	= &rwork_[iinfo_[DG_IA_lc_face]+qface_n[0]];
  uface2 	= &rwork_[iinfo_[DG_IA_lc_face]+2*qface_n[0]];
  
  const I tmpbrhs_ufaceoff 	= bmat_m;
  pR tmpbrhs        	= &rwork_[bmat_m]; 
  pR tmpbrhs_ufaces 	= &tmpbrhs[1+2*nu[0]]; 
  pR tmpbrhs_uface0 	= &tmpbrhs[1+2*nu[0]]; 
  pR tmpbrhs_uface1 	= &tmpbrhs[1+2*nu[0]+qface_n[0]];
  pR tmpbrhs_uface2 	= &tmpbrhs[1+2*nu[0]+2*qface_n[0]];
  pR tmpbrhs_uelm   	= &tmpbrhs[1]; 
  pR tmpbrhs_u 		= tmpbrhs_uelm;

  /*_____________________________________________________*/
  if (renum_iwork_n>iwork_n_[0])
    {
      ires_[DG_ires_err] = (I)1;
      ires_[DG_ires_required_iw_n] = renum_iwork_n;
      fprintf(stderr,"too small iwork_n_ "ifmt" "ifmt"\n",renum_iwork_n,iwork_n_[0]);
      return;
    }

  lcperm  	= &iwork_[0];
  perm  	= &iwork_[trial_n[0]];
  graph 	= &iwork_[(nelm_[0]+1) + trial_n[0] ];
  marker 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]];
  marker_begin 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ (nelm_[0]+1)];
  marker_flag 	= &iwork_[(nelm_[0]+1) + trial_n[0] + 3*nelm_[0]+ 2*(nelm_[0]+1)];
  /*_____________________________________________________*/

  goto __state_renum;

  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  

  
 __state_renum:
  b1=1;
  b2=0;    
  { I j;for (j=trial_n[0];j<renum_iwork_n;++j) iwork_[j]=(I)0;}  
  /* */
  /* on donne un sens a l adjacence */
  { I jelm;
    for (jelm=0;jelm<nelm_[0];++jelm)
      {        
	/* compute jaface and nrmelm */
	neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];
	cnc[0]     = cnc_[jelm*cncoff_[0]+0];cnc[1]     = cnc_[jelm*cncoff_[0]+1];cnc[2]     = cnc_[jelm*cncoff_[0]+2];
	cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
	cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
	jacface[0] = nsSQRT( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
	jacface[1] = nsSQRT( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
	jacface[2] = nsSQRT( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
	longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];  
	nrmelm[0] = cooelm[4] - cooelm[3];nrmelm[1] = cooelm[0] - cooelm[1];				
	x = ((R)1.0)/jacface[0];nrmelm[0] *= x;nrmelm[1] *= x;						  
	nrmelm[2] = cooelm[5] - cooelm[4];nrmelm[3] = cooelm[1] - cooelm[2];
	x = ((R)1.0)/jacface[1];nrmelm[2] *= x;nrmelm[3] *= x;  
	nrmelm[4] = cooelm[3] - cooelm[5];nrmelm[5] = cooelm[2] - cooelm[0] ;
	x = ((R)1.0)/jacface[2];nrmelm[4] *= x;nrmelm[5] *= x;
	jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	/*____ INFO ELEMENT DONE */ 	
	/*___________________________________________________________________________________________________________________*/  
	/* 
	   eval u faces
	*/  
	/*___________________________________________________________________________________________________________________*/
	Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*iinfo_[DG_IA_lc_face]],&eval_uoff_,&data_u_[jelm*data_uoff_[0]],&n1,&r0,tmpbrhs_ufaces,&n1);
	Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*iinfo_[DG_IA_lc_face]],&eval_uoff_,&data_v_[jelm*data_voff_[0]],&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);
	Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&nrmelm[2*0],&n1,&r0,uface0,&n1);
	Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&nrmelm[2*1],&n1,&r0,uface1,&n1);
	Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&nrmelm[2*2],&n1,&r0,uface2,&n1);	
	
	{ I k;
	  for (k=0;k<qface_n[0];++k) if (uface0[k]<((R)0.0)) break;
	  graph[jelm*3+0] = (neids[0]) ? ( (k<qface_n[0])?-neids[0]:neids[0] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	{ I k;
	  for (k=0;k<qface_n[0];++k)if (uface1[k]<((R)0.0))break;
	  graph[jelm*3+1] = (neids[1]) ? ( (k<qface_n[0])?-neids[1]:neids[1] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }
	{ I k;
	  for (k=0;k<qface_n[0];++k)if (uface2[k]<((R)0.0))break;
	  graph[jelm*3+2] = (neids[2]) ? ( (k<qface_n[0])?-neids[2]:neids[2] ) : ( (k<qface_n[0]) ? -nelm_[0]-1 : (I)0 ); }	
      } }
  /* on met tous ceux qui n'ont pas d antecedent en premier*/
  { I jelm;
    for (jelm=0;jelm<nelm_[0];++jelm)
      { 
	if (graph[jelm*3]<-nelm_[0]) {perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	if (graph[jelm*3+1]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
	if (graph[jelm*3+2]<-nelm_[0]){perm[++b2]=jelm+1; marker_flag[jelm+1]=1;}
      } }
  /*________________________________________________________________________________________________________________*/  
  pp	= b2;
  n	= (I)1;
 __state_renum_do1:
  marker_size=0;
  /* INJECTION */
  n=0;
  { I i;	
    for (i=b1;i<=b2;++i)
      {
	{ I j;
	  for (j=0;j<3;++j)
	    {
	      const I s = graph[(perm[i]-1)*3+j];
	      if ( (s>0) AND ( (NOT marker[s]) AND (NOT marker_flag[s]) ))
		{
		  marker_size=marker_size+1; marker_begin[marker_size] = s;
		  marker[s] = marker_size;n=n+1;
		}
	    } }
      } }
    
  /* FIN INJECTION */
  { I i;	    
    for (i=1;i<=marker_size;++i)
      {
	marker_flag[marker_begin[i]]=1;perm[b2+i]=marker_begin[i];
	marker[marker_begin[i]]=0;marker_begin[i]=0;
      } }

  
  pp+=n;b1 = b2+1;b2 = b2+n;
  if (b1>b2)
    {
      if (pp==nelm_[0])
	goto __state_renum_do1_done;
      else
	{
	  /* on cherche un element*/
	  { I jelm;
	    for (jelm=1;jelm<=nelm_[0];++jelm)
	      {
		if (NOT marker_flag[jelm])
		  {
		    marker_flag[jelm]=(I)1;
		    perm[b1] = jelm;
		    ++n;b2 = b1;
		    break;
		  }
	      } }
	}
    }
  if (n>0)
    {
      /* on continue a injecter */
      goto __state_renum_do1;
    }
    
    
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
 __state_renum_do1_done:
  if (b2!=nelm_[0])
    {      
      fprintf(stderr,"*** DGERR RENUM FAILED "ifmt"!="ifmt"",b2,nelm_[0]);
      goto __state_error;
    }
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  iter_gauss_seidel = 0;
 __state_next_iter_gauss_seidel:
#if __mk_debug__
  fprintf(stdout,"next iter gauss seidel \n");
#endif
  mx = (R)0.0;
  /*___________________________________________________________________________________________________________________*/  
  ielm=0;      
 __state_next_elm:
  id 		= perm[ielm+1]-1;
  neids[0]	= adj_[id*adjoff_[0]];
  neids[1]	= adj_[id*adjoff_[0]+1];
  neids[2]	= adj_[id*adjoff_[0]+2];
  vcod[0] 	= vcod_[cnc_[id*cncoff_[0]+0]-1];
  vcod[1] 	= vcod_[cnc_[id*cncoff_[0]+1]-1];
  vcod[2] 	= vcod_[cnc_[id*cncoff_[0]+2]-1];
  cnc[0]        = cnc_[id*cncoff_[0]+0];
  cnc[1]        = cnc_[id*cncoff_[0]+1];
  cnc[2]        = cnc_[id*cncoff_[0]+2];
  codface[0] 	= (neids[0])? noboundary_cod_[0] :  MAX(vcod_[0],vcod_[1]);
  codface[1] 	= (neids[1])? noboundary_cod_[0] :  MAX(vcod_[1],vcod_[2]);
  codface[2] 	= (neids[2])? noboundary_cod_[0] :  MAX(vcod_[2],vcod_[0]);
  neids_face[0]=0;neids_face[1]=0;neids_face[2]=0;
  /*---*/
  if (neids[0])
    {	if (neids[0]<0) neids[0]=-neids[0];
      if ( (adj_[(neids[0]-1)*adjoff_[0]+0]==id+1) ) neids_face[0]=0;
      else if ( (adj_[(neids[0]-1)*adjoff_[0]+1]==id+1) ) neids_face[0]=1;
      else if ( (adj_[(neids[0]-1)*adjoff_[0]+2]==id+1) ) neids_face[0]=2; }
  if (neids[1])
    {	if (neids[1]<0) neids[1]=-neids[1];
      if ( (adj_[(neids[1]-1)*adjoff_[0]+0]==id+1) ) neids_face[1]=0;
      else if ( (adj_[(neids[1]-1)*adjoff_[0]+1]==id+1) ) neids_face[1]=1;
      else if ( (adj_[(neids[1]-1)*adjoff_[0]+2]==id+1) ) neids_face[1]=2; }
  if (neids[2])
    {	if (neids[2]<0) neids[2]=-neids[2];
      if ( (adj_[(neids[2]-1)*adjoff_[0]+0]==id+1) ) neids_face[2]=0;
      else if ( (adj_[(neids[2]-1)*adjoff_[0]+1]==id+1) ) neids_face[2]=1;
      else if ( (adj_[(neids[2]-1)*adjoff_[0]+2]==id+1) ) neids_face[2]=2; }
  /*---*/
  cooelm[0]  = coo_[cooff_[0]*(cnc[0]-1)+0];
  cooelm[1]  = coo_[cooff_[0]*(cnc[1]-1)+0];
  cooelm[2]  = coo_[cooff_[0]*(cnc[2]-1)+0];
  cooelm[3]  = coo_[cooff_[0]*(cnc[0]-1)+1];
  cooelm[4]  = coo_[cooff_[0]*(cnc[1]-1)+1];
  cooelm[5]  = coo_[cooff_[0]*(cnc[2]-1)+1];
  jacface[0] = nsSQRT( (cooelm[0]-cooelm[1])*(cooelm[0]-cooelm[1]) + (cooelm[3]-cooelm[4])*(cooelm[3]-cooelm[4]) ); 
  jacface[1] = nsSQRT( (cooelm[1]-cooelm[2])*(cooelm[1]-cooelm[2]) + (cooelm[4]-cooelm[5])*(cooelm[4]-cooelm[5]) ); 
  jacface[2] = nsSQRT( (cooelm[2]-cooelm[0])*(cooelm[2]-cooelm[0]) + (cooelm[5]-cooelm[3])*(cooelm[5]-cooelm[3]) ); 	  
  longueurs[0] = jacface[0];longueurs[1] = jacface[1];longueurs[2] = jacface[2];	  
  nrmelm[0] = cooelm[4] - cooelm[3];nrmelm[1] = cooelm[0] - cooelm[1];					
  x = ((R)1.0)/jacface[0];nrmelm[0] *= x;nrmelm[1] *= x;							  
  nrmelm[2] = cooelm[5] - cooelm[4];nrmelm[3] = cooelm[1] - cooelm[2];
  x = ((R)1.0)/jacface[1];nrmelm[2] *= x;nrmelm[3] *= x;	  
  nrmelm[4] = cooelm[3] - cooelm[5];nrmelm[5] = cooelm[2] - cooelm[0] ;
  x = ((R)1.0)/jacface[2];nrmelm[4] *= x;nrmelm[5] *= x;
  if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
  if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
  if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
  jacelm[0]  = nsSQRT((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  
  belm[0] = cooelm[5]-cooelm[3];belm[1] = cooelm[0]-cooelm[2];belm[2] = cooelm[3]-cooelm[4];belm[3] = cooelm[1]-cooelm[0];  	  
  x = belm[2];belm[2]=belm[1];belm[1]=x;	  
  jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
  /*____ INFO ELEMENT DONE */ 
  /*____ COPY LOCAL RHS */ 
  {I j;for (j=0;j<test_n[0];++j) lcrhs[j] = rhs_[rhsoff_[0]*id+j];}	
  brhs_a[0] = jacelm[0];	
  /*___________________________________________________________________________________________________________________*/  
  /* 
     eval u : elm and faces
  */  
  /*___________________________________________________________________________________________________________________*/  
  Blas_dgemv("T",nu,&N_pts_u,&r1,eval_u_,&eval_uoff_,&data_u_[id*data_uoff_[0]],&n1,&r0,tmpbrhs_u,&n1);
  Blas_dgemv("T",nu,&N_pts_u,&r1,eval_u_,&eval_uoff_,&data_v_[id*data_voff_[0]],&n1,&r0,&tmpbrhs_u[tmpbrhs_ufaceoff],&n1);	  	  
  /*____ TRANSFORM ON UELM */   
  Blas_dgemm("N","T",nu,&n2,&n2,&mxu,tmpbrhs_uelm,&tmpbrhs_ufaceoff/*=bmat_m*/,belm,&n2,&r0,brhs_uelm,nu);
  /*____ COMPUTE MAX(-U.N,0) */ 
  Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&nrmelm[2*0],&n1,&r0,uface0,&n1);
  { I j;for (j=0;j<qface_n[0];++j) uface0[j] = (uface0[j]<0.0) ? uface0[j]*(-xu_[0]) : (R)0.0;}       
  Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&nrmelm[2*1],&n1,&r0,uface1,&n1);
  { I j;for (j=0;j<qface_n[0];++j) uface1[j] = (uface1[j]<0.0) ? uface1[j]*(-xu_[0]) : (R)0.0;}	
  Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&nrmelm[2*2],&n1,&r0,uface2,&n1);
  { I j;for (j=0;j<qface_n[0];++j) uface2[j] = (uface2[j]<0.0) ? uface2[j]*(-xu_[0]) : (R)0.0;}
  /*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
  if ( neids[0] )
    {	Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[0][neids_face[0]],&bmat_n,uface0,&n1,&r0,matflux,&n1);		  		  
      Blas_dgemv("N",trial_n,trial_n,&r1,matflux_corr,trial_n,&corr_[(neids[0]-1)*corroff_[0]],&n1,&r1,lcrhs,&n1);		  
      Blas_dgemv("N",trial_n,teta_n,&r1,matflux_residu,trial_n,&sol_[(neids[0]-1)*soloff_[0]],&n1,&r1,lcrhs,&n1); }
  if ( neids[1] )
    {	Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[1][neids_face[1]],&bmat_n,uface1,&n1,&r0,matflux,&n1);		  		  
      Blas_dgemv("N",trial_n,trial_n,&r1,matflux_corr,trial_n,&corr_[(neids[1]-1)*corroff_[0]],&n1,&r1,lcrhs,&n1);		  
      Blas_dgemv("N",trial_n,teta_n,&r1,matflux_residu,trial_n,&sol_[(neids[1]-1)*soloff_[0]],&n1,&r1,lcrhs,&n1); }	
  if ( neids[2] )
    { Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[2][neids_face[2]],&bmat_n,uface2,&n1,&r0,matflux,&n1);		  		  
      Blas_dgemv("N",trial_n,trial_n,&r1,matflux_corr,trial_n,&corr_[(neids[2]-1)*corroff_[0]],&n1,&r1,lcrhs,&n1);		  
      Blas_dgemv("N",trial_n,teta_n,&r1,matflux_residu,trial_n,&sol_[(neids[2]-1)*soloff_[0]],&n1,&r1,lcrhs,&n1); }
#if DG_CONSERVATIVE_WEAK_FORM
  /*____ COMPUTE MAX(U.N,0), APPLY FLUX MATRIX ON DOWNWIND STREAM */ 
  Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&nrmelm[2*0],&n1,&r0,uface0,&n1);
  {I j;for (j=0;j<qface_n[0];++j)uface0[j] = (uface0[j]>0.0) ? uface0[j]*(xu_[0]) : (R)0.0;}
  Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&nrmelm[2*1],&n1,&r0,uface1,&n1);
  {I j;for (j=0;j<qface_n[0];++j)uface1[j] = (uface1[j]>0.0) ? uface1[j]*(xu_[0]) : (R)0.0;}
  Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&nrmelm[2*2],&n1,&r0,uface2,&n1);
  {I j;for (j=0;j<qface_n[0];++j)uface2[j] = (uface2[j]>0.0) ? uface2[j]*(xu_[0]) : (R)0.0;}
#endif  
  Blas_dgemv("N",&bmat_n,&bmat_m,&r1,bmat,&bmat_n,brhs,&n1,&r0,lcmat,&n1);	  
  /*--- COMPUTE LOCAL RESIDUAL  */
  Blas_dgemv("N",trial_n,teta_n,&mr1,lcmatr,trial_n,&sol_[id*soloff_[0]],&n1,&r1,lcrhs,&n1);
  /*--- COMPUTE LOCAL CORRECTION   */
  { I j; for (j=0;j<test_n[0];++j) {lcsol[j] = corr_[corroff_[0]*id+j];}}

  dgesv(test_n,
	&n1,
	lcmatc,
	test_n,
	lcperm,
	lcrhs,
	test_n,
	&info_lapack);

  

   
  /*--- COMPUTE NRMS OF LOCAL CORRECTION */
  { I j;for (j=0;j<test_n[0];++j) corr_[corroff_[0]*id+j]=lcrhs[j];}
  { I j;for (j=0;j<test_n[0];++j) lcrhs[j]-=lcsol[j];}
  { R xx=(R)0.0; {I j;for (j=0;j<test_n[0];++j) xx+=lcrhs[j]*lcrhs[j];}
    xx=nsSQRT(xx);if (mx<xx) mx = xx; }		    


  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  ++ielm;
  if (ielm<nelm_[0])
    {
      goto  __state_next_elm;
    }  
  else
    {
      goto  __state_next_elm_done;
    }


  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
 __state_next_elm_done:
  {
    fprintf(stdout,"*** DGINFO time %.5f iter "ifmt"/"ifmt" tol = %8.15e\n",t_[0],iter_gauss_seidel,nelm_[0],mx);  
    if ( (mx>((R)1.0e-12)) AND (iter_gauss_seidel<nelm_[0]) )
      {
	++iter_gauss_seidel;
	goto __state_next_iter_gauss_seidel;
      }
    else
      {
	
	/* CALCUL DE NORMES */  
	
	
	goto __state_computation_done;
      }
  }
  
  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
  /*___________________________________________________________________________________________________________________*/  
 __state_computation_done:
  {    
    ires_[DG_ires_iter_gauss_seidel] 	= iter_gauss_seidel;
    rres_[DG_rres_max] 			= mx;
    if (iter_gauss_seidel<nelm_[0])
      {
	ires_[DG_ires_convergence] = (I)1;
      }
    else
      {
	ires_[DG_ires_convergence] = (I)0;
      }
    return;
  }
  
  /*_________________________________________________________________________________________________________________*/    
 __state_error:
  {
    printf("out\n");
    ires_[DG_ires_err] = (I)1;
    return;
  }
}













































#endif
