


#if 0
  void dgadvection_local_system(cst_pR 	xa_,
				cst_pR 	xu_,
				struct matrix_handle *    matelm,
				struct vector_handle *    cooelm,
				struct vector_handle *    rhselm,
				struct vector_handle *    aelm,
				struct vector_handle *    uelm,
				struct vector_handle *    velm,
				struct vector_handle *    solelm,
				struct vector_handle *    neisolelm[],
				struct vector_handle *    correlm,
				struct vector_handle *    neicorrelm[],
				struct matrix_handle *    u_eval,
				struct matrix_handle *    belm,
				struct matrix_handle *    nrmelm,
				I 			neids[],
				I 			neids_face[],
				cst_pI 			rwork_n_,
				pR 			rwork_,			      
				cst_pR 			rinfo_,
				const I   		iinfo_[DG::I_n])
  {    
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    cst_pI	nu      	= &iinfo_[DG::I_TETA_U_NBASIS];
    I N_pts_u2 			= bmat_m-1;/*nTot-nA*/

  
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    //    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    cst_pR bmat 	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
  
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }


    struct matrix_handle mat_tmpbrhs_uface0;
    struct matrix_handle mat_tmpbrhs_uface1;
    struct matrix_handle mat_tmpbrhs_uface2;
    struct matrix_handle mat_matflux_corr;
    struct matrix_handle mat_matflux_residu;
    struct matrix_handle hbmat;
 
    struct vector_handle vec_uface0;
    struct vector_handle vec_uface1;
    struct vector_handle vec_uface2;
    struct vector_handle vec_nrmelm0;
    struct vector_handle vec_nrmelm1;
    struct vector_handle vec_nrmelm2;
    struct vector_handle hbrhs;
    struct vector_handle hsol;
  
    struct matrix_handle mat_tmpbrhs_uelm;
    struct matrix_handle mat_brhs_uelm;
    R
      jacface[3*1];

    R 	matflux[21*(21+21)];/*trial_n*(trial_n+teta_n)*/
    pR	matflux_corr 	= &matflux[0];
    pR	matflux_residu 	= &matflux[trial_n[0]*trial_n[0]];
  
    I n1 = 1;
    R r1 = 1.0;

    R r0 = 0.0;
#if DG_CONSERVATIVE_WEAK_FORM
    const R mxu=xu_[0]*((R)-1.0);
#else
    const R mxu=xu_[0]*((R)1.0);
#endif
    pR ufaces[3] = {nullptr,nullptr,nullptr};

    /*---*/

    const I tmpbrhs_ufaceoff 	= bmat_m;
  
    pR brhs      	= &rwork_[0];
    pR brhs_uelm 	= &rwork_[1];
    ufaces[0] 	= &rwork_[iinfo_[DG::IA_lc_face]];
    ufaces[1] 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    ufaces[2] 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

  
    pR tmpbrhs        	= &rwork_[bmat_m]; 
    pR tmpbrhs_uface0 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface1 	= &tmpbrhs[1+2*nu[0]+qface_n[0]];
    pR tmpbrhs_uface2 	= &tmpbrhs[1+2*nu[0]+2*qface_n[0]];
    pR tmpbrhs_uelm   	= &tmpbrhs[1]; 
    pR tmpbrhs_u 		= tmpbrhs_uelm;
  
    matrix_handle_def(&mat_tmpbrhs_uelm,nu[0],2,tmpbrhs_uelm,tmpbrhs_ufaceoff);
    matrix_handle_def(&mat_brhs_uelm,nu[0],2,brhs_uelm,nu[0]);  
    matrix_handle_def(&hbmat,bmat_n,bmat_m,(pR)bmat,bmat_n);
    vector_handle_def(&hbrhs,bmat_m,brhs,n1);

    struct vector_handle built_matrices;
    struct matrix_handle built_matrix_residu;

    struct matrix_handle mat_lcmatc;  
    R 	        lcmat[21*21*2];

    //    pR 		lcmatc = &lcmat[0];
    pR 		lcmatr = &lcmat[trial_n[0]*trial_n[0]];

    vector_handle_def(&built_matrices,bmat_n,lcmat,n1);
    matrix_handle_def(&built_matrix_residu,trial_n[0],teta_n[0], lcmatr, trial_n[0]);
  
    //  matrix_handle_def(&mat_lcmatc,test_n[0],test_n[0],lcmatc,test_n[0]);

  
    //  vector_handle_def(&hlcrhs,trial_n[0],lcrhs,1);
    //  vector_handle_def(&hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
    matrix_handle_def(&mat_matflux_corr,trial_n[0],trial_n[0],matflux_corr,trial_n[0]);
    matrix_handle_def(&mat_matflux_residu,trial_n[0],teta_n[0],matflux_residu,trial_n[0]);

  
    matrix_handle_def(&mat_tmpbrhs_uface0,qface_n[0],2,tmpbrhs_uface0,tmpbrhs_ufaceoff);
    matrix_handle_def(&mat_tmpbrhs_uface1,qface_n[0],2,tmpbrhs_uface1,tmpbrhs_ufaceoff);
    matrix_handle_def(&mat_tmpbrhs_uface2,qface_n[0],2,tmpbrhs_uface2,tmpbrhs_ufaceoff);
  
    vector_handle_def(&vec_uface0,qface_n[0],ufaces[0],1);
    vector_handle_def(&vec_uface1,qface_n[0],ufaces[1],1);
    vector_handle_def(&vec_uface2,qface_n[0],ufaces[2],1);
  
    vector_handle_def(&vec_nrmelm0,2,nrmelm->x,1);
    vector_handle_def(&vec_nrmelm1,2,nrmelm->x+nrmelm->ld,1);
    vector_handle_def(&vec_nrmelm2,2,nrmelm->x+2*nrmelm->ld,1);
  
    struct vector_handle brhs_u;
    struct vector_handle brhs_v;

    vector_handle_def(&brhs_u,N_pts_u2,tmpbrhs_u,1);
    vector_handle_def(&brhs_v, N_pts_u2,tmpbrhs_u + tmpbrhs_ufaceoff,1);

    //
    // Evaluate u (only needed if non-nodal).
    //
    matrix_handle_gemv(u_eval,
		       "T",
		       &r1,
		       uelm,
		       &r0,
		       &brhs_u);  

    //
    // Evaluate u (only needed if non-nodal).
    //
    matrix_handle_gemv(u_eval,
		       "T",
		       &r1,
		       velm,
		       &r0,
		       &brhs_v);

    //
    // TRANSFORM ON UELM
    //
    matrix_handle_gemm(&mat_tmpbrhs_uelm,
		       "N",
		       "T",
		       &mxu,
		       belm,
		       &r0,
		       &mat_brhs_uelm);

    vec_uface0 = jacface[0] * mat_tmpbrhs_uface0 * vec_nrmelm0;
    vec_uface1 = jacface[1] * mat_tmpbrhs_uface1 * vec_nrmelm1;
    vec_uface2 = jacface[2] * mat_tmpbrhs_uface2 * vec_nrmelm2;

  
#if 0
    matrix_handle_gemv(&mat_tmpbrhs_uface0,"N",&jacface[0],&vec_nrmelm0,&r0,&vec_uface0);
    matrix_handle_gemv(&mat_tmpbrhs_uface1,"N",&jacface[1],&vec_nrmelm1,&r0,&vec_uface1);
    matrix_handle_gemv(&mat_tmpbrhs_uface2,"N",&jacface[2],&vec_nrmelm2,&r0,&vec_uface2);
#endif  
    //
    // Apply min mod.
    //
    for (I k=0;k<3;++k)
      {
	pR uface = ufaces[k];
	for (I j=0;j<qface_n[0];++j)
	  {
	    uface[j] = (uface[j]<0.0) ? uface[j]*(-xu_[0]) : (R)0.0;
	  }
      }

    /*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
    for (int localFaceIndex=0;localFaceIndex<3;++localFaceIndex)
      {
	if (neids[localFaceIndex])
	  {	 
	    Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,matflux,&n1);		  		  
	    //
	    // Compute the flux matrix for the correction.
	    //
	    matrix_handle_gemv(&mat_matflux_corr,"N",&r1,neicorrelm[localFaceIndex],&r1,rhselm);
	    //
	    // Compute the flux matrix for the residu.
	    //
	    matrix_handle_gemv(&mat_matflux_residu,"N",&r1,neisolelm[localFaceIndex],&r1,rhselm);
	  }
      }

#if 0
#if DG_CONSERVATIVE_WEAK_FORM
    /*____ COMPUTE MAX(U.N,0), APPLY FLUX MATRIX ON DOWNWIND STREAM */ 
    Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&nrmelm[2*0],&n1,&r0,uface0,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface0[j] = (uface0[j]>0.0) ? uface0[j]*(xu_[0]) : (R)0.0;}
    Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&nrmelm[2*1],&n1,&r0,uface1,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface1[j] = (uface1[j]>0.0) ? uface1[j]*(xu_[0]) : (R)0.0;}
    Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&nrmelm[2*2],&n1,&r0,uface2,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface2[j] = (uface2[j]>0.0) ? uface2[j]*(xu_[0]) : (R)0.0;}
#endif
#endif
 
    built_matrices = hbmat * hbrhs;
    *rhselm += -1.0 *  built_matrix_residu * hsol;
 
    // matrix_handle_gemv(&built_matrix_residu,"N",&mr1,&hsol,&r1,rhselm);

    //
    // copy
    //
  
  };
#endif  
