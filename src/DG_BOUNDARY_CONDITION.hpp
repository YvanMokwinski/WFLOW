#pragma once
struct DG_BOUNDARY_CONDITION
{
  double m_mem[16000];
  I m_qn;
  I m_nu;
  I m_nx;
  I m_ntest;
  //
  // Matrix to build the evaluation of u.
  //
  matrix_handle m_beval_uvw;
  matrix_handle m_beval_uvw_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  matrix_handle m_beval_xyz;
  matrix_handle m_beval_xyz_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  matrix_handle m_beval_test;
  matrix_handle m_beval_test_part[nfaceinelm];

  
  matrix_handle m_qrst;
  matrix_handle m_qrst_part[nfaceinelm];
  
  
  DG_BOUNDARY_CONDITION(){};
  DG_BOUNDARY_CONDITION(mkS s_test_,
			mkS s_teta_u_,
			mkS s_teta_x_,
			I 	qn,
			cst_pR 	qw,
			cst_pR 	qp)
  {
    const I qnXnfaceinelm = qn * nfaceinelm;
    const I nu = mkS_n(s_teta_u_);
    const I nx = mkS_n(s_teta_x_);
    const I ntest = mkS_n(s_test_);
    m_qn = qn;
    m_nu = nu;
    m_nx = nx;
    m_ntest = ntest;
    
    I at = 0;

    //
    // Reference quadrature points 
    //
    {
      I qrst_size = qnXnfaceinelm*dim;
      matrix_handle_def(&this->m_qrst,
			qnXnfaceinelm,
			dim,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    this->m_qrst.x + localFaceIndex * qn,
			    this->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   this->m_qrst.x,
		   &this->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      matrix_handle_def(&this->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&this->m_mem[at],
			nx);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    this->m_beval_xyz.x + (localFaceIndex * qn) * this->m_beval_xyz.ld,
			    this->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      matrix_handle_def(&this->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&this->m_mem[at],
			ntest);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    this->m_beval_test.x + (localFaceIndex * qn) * this->m_beval_test.ld,
			    this->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      matrix_handle_def(&this->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&this->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    this->m_beval_uvw.x + localFaceIndex * qn * this->m_beval_uvw.ld,
			    this->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    I err;
    R rwork[4096];
    I rwork_n = 4096;

    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      this->m_beval_xyz.x,
	      &this->m_beval_xyz.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      this->m_beval_uvw.x,
	      &this->m_beval_uvw.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      this->m_beval_test.x,
	      &this->m_beval_test.ld,
	      this->m_qrst.x,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    for (I l=0;l<nfaceinelm;++l)
      {
	for (I j=0;j<qn;++j)
	  {
	    for (I i=0;i<ntest;++i)
	      {
		m_beval_test_part[l].x[m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }



  static void define(DG_BOUNDARY_CONDITION*self,
		     mkS s_test_,
		     mkS s_teta_u_,
		     mkS s_teta_x_,
		     I 	qn,
		     cst_pR 	qw,
		     cst_pR 	qp)
  {
    const I qnXnfaceinelm = qn * nfaceinelm;
    const I nu = mkS_n(s_teta_u_);
    const I nx = mkS_n(s_teta_x_);
    const I ntest = mkS_n(s_test_);
    self->m_qn = qn;
    self->m_nu = nu;
    self->m_nx = nx;
    self->m_ntest = ntest;
    
    I at = 0;

    //
    // Reference quadrature points 
    //
    {
      I qrst_size = qnXnfaceinelm*dim;
      matrix_handle_def(&self->m_qrst,
			qnXnfaceinelm,
			dim,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    self->m_qrst.x + localFaceIndex * qn,
			    self->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   self->m_qrst.x,
		   &self->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      matrix_handle_def(&self->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&self->m_mem[at],
			nx);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    self->m_beval_xyz.x + (localFaceIndex * qn) * self->m_beval_xyz.ld,
			    self->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      matrix_handle_def(&self->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&self->m_mem[at],
			ntest);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    self->m_beval_test.x + (localFaceIndex * qn) * self->m_beval_test.ld,
			    self->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      matrix_handle_def(&self->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  matrix_handle_def(&self->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    self->m_beval_uvw.x + localFaceIndex * qn * self->m_beval_uvw.ld,
			    self->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    I err;
    R rwork[4096];
    I rwork_n = 4096;

    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      self->m_beval_xyz.x,
	      &self->m_beval_xyz.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      self->m_beval_uvw.x,
	      &self->m_beval_uvw.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      self->m_beval_test.x,
	      &self->m_beval_test.ld,
	      self->m_qrst.x,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    for (I l=0;l<nfaceinelm;++l)
      {
	for (I j=0;j<qn;++j)
	  {
	    for (I i=0;i<ntest;++i)
	      {
		self->m_beval_test_part[l].x[self->m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }

struct DATA
{
  double m_mem[512];
  
  matrix_handle eval_xyz;
    matrix_handle eval_uvw;    

    vector_handle eval_udotn;
    vector_handle eval_f;
    vector_handle lrhs;

    DATA(I qn,I ntest)
    {
      // qn * (2 * dim + 3)
      I at = 0;
      matrix_handle_def(&eval_xyz,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      matrix_handle_def(&eval_uvw,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      vector_handle_def(&eval_udotn,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&eval_f,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&lrhs,ntest,&m_mem[at],1);
      at+=ntest;
    };
};    

DATA * CreateData()
{
  return new DATA(this->m_qn,
		  this->m_ntest);
};

  R m_eps;
  
  SmoothedHeaviside m_hea;
  template <typename userfct_t>
  void boundary_condition(const I 		localFaceIndex,
			  cst_pR xu_,
			  const vector_handle&	normal,
			  const matrix_handle&	xyz,
			  const matrix_handle&	uvw,
			  DATA*			data,
			  userfct_t             userfct)
  {
    m_eps = 0.25;
    {
      Err e;
      SmoothedHeaviside_def	( &this->m_hea,
				  __eHeaviside_m4p3,
				  &m_eps,
				  &e);
    }
    
    
    //
    // Evaluate xyz.
    //
    data->eval_xyz = this->m_beval_xyz_part[localFaceIndex].transpose() * xyz;
#if 0    
    fprintf(stdout,"xyz " ifmt " " ifmt " " ifmt " " ifmt "\n",data->eval_xyz.n,data->eval_xyz.m,data->eval_xyz.ld,data->eval_f.ld);
    matrix_handle_print(&data->eval_xyz,stdout);
#endif

    //
    // User function.
    //


    //
    // eval_f pouvait ne pas etre initialise dans userfct
    //
    userfct(data->eval_xyz.n,
	    data->eval_xyz.x,
	    data->eval_xyz.ld,
	    data->eval_f.x,
	    data->eval_f.ld);
    
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif
    //
    // Evaluate uvw.
    //
#if 0
    fprintf(stdout,"uvw\n");
    matrix_handle_print(&uvw,stdout);
#endif
    
    data->eval_uvw = this->m_beval_uvw_part[localFaceIndex].transpose() * uvw;
    
#if 0
    fprintf(stdout,"eval_uvw\n");
    matrix_handle_print(&data->eval_uvw,stdout);
#endif    
    //
    // Compute the dot product with the normal
    //
    data->eval_udotn = data->eval_uvw * normal;
#if 0
    fprintf(stdout,"eval_udotn\n");
    vector_handle_print(&data->eval_udotn,stdout);
#endif

    //
    //
    //
    
    //
    // Form the flux functions.
    //
    
    for (I i=0;i<data->eval_f.n;++i)
      {	//data->eval_f.x[i*data->eval_f.ld]=1.0;
	data->eval_f.x[i*data->eval_f.ld]
	  *= ( data->eval_udotn.x[i*data->eval_udotn.ld] < ((R)0.0) ) ? -data->eval_udotn.x[i*data->eval_udotn.ld] * xu_[0] : ((R)0.0) ;
      }
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif	      
    //
    // Form the residu.
    //
    data->lrhs = this->m_beval_test_part[localFaceIndex] * data->eval_f;
#if 0
    fprintf(stdout,"lrhs\n");
    vector_handle_print(&data->lrhs,stdout);
#endif
  }

  

  void boundary_condition(ns_mesh * 	mesh,
			  cst_pI  	cnc_u_,
			  cst_pI  	cncoff_u_,
			  cst_pR 	data_u_,		       
			  cst_pR	data_v_,
			  pR 		rhs_,
			  I 		rhsoff_)
  {    

    vector_handle normal;
    matrix_handle xyz;
    matrix_handle uvw;
    double normal_values[6];
    double cooelm[32];
    double uvw_values[32];
    matrix_handle_def(&xyz,m_nx,dim,cooelm,m_nx);
    matrix_handle_def(&uvw,m_nu,dim,uvw_values,m_nu);
    vector_handle_def(&normal,2,normal_values,1);
    
    DATA* data = this->CreateData();
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {

		//
		//
		//
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const I jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		normal_values[0]	= mesh->normaledge[2*jedge+0];
		normal_values[1] 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];
		//

		//
		//
		//
		
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];;
		  }
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[m_nu+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		  }
		//
		
		boundary_condition(jadj,
				   &xjac,
				   normal,
				   xyz,
				   uvw,
				   data,
				   [this](I n_,cst_pR pos,I posoff,pR f,I foff)
				   {
				     for (I k=0;k<n_;++k)
				       {
#if 1
					 
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];
					 if (x<1.0e-13)
					   {
					     f[foff*k] = sin(12.0*y);
					   }
					 
					 
#else
					 
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];

					 
					 if ( (x < 1.0e-13) && (y>0.5 && y<=1.0) )
					   {
					     f[foff*k] = 0.0;
				       }
					 else if ( (x < 1.0e-13) && (y<=0.5 && y >0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x>0.5 && x<=1.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x<=0.5 && x >=0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (x > 1.0-1.0e-13) )
					   {
					 double z  = sin(acos(-1.0)*y);
					 f[foff*k] = z*z*z;
				       }
					 else
					   {
					 f[foff*k] = 0.0;
				       }
#endif
				       }
				     
				     
#if 0
				     
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = y - 0.5; // exp(-sin(32.0*y)*y);
				       }
				     Err e;
				     SmoothedHeaviside_eval	(&this->m_hea,
								 n_,
								 f,
								 &e);
#endif
#if 0
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = exp(-sin(32.0*y)*y);
				       }
#endif
				   });
		
		for (I k=0;k<m_ntest;++k)
		  {
		    rhs_[jelm * rhsoff_ + k] += data->lrhs.x[k];
		  }		
	      }
	  }
      }
    
#if 0
    printf("ddddddddddddddddddd\n");
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (I k=0;k<m_ntest;++k)
	  {
	    std::cout << rhs_[jelm * rhsoff_ + k] << std::endl;
	  }
	printf("\n");
      }
    printf("ddddddddddddddddddd\n");
#endif    

    
#if 0
    double mem[2048];
    matrix_handle beval_uvw[nfaceinelm];
    matrix_handle eval_uvw;
    matrix_handle beval_xyz[nfaceinelm];
    matrix_handle nrmelm;
    matrix_handle cooelm;
    matrix_handle uvwelm;
    vector_handle udotn;
    matrix_handle qxyz;
    vector_handle f;
  
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const I jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		const R nx	= mesh->normaledge[2*jedge+0];
		const R ny 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];	    
		
	      //
	      // Interpolation of u * nx + v*ny + w * nz 
	      //

	      
	      //
	      // Evaluation of u on the edge.
	      //
	      eval_uvw = beval_uvw[jadj] * dofs_uvw;

	      //
	      // Evaluation of u on the edge.
	      //
	      udotn = eval_uvw * nrmedge;

	      //
	      //
	      // Cartesian coordinates of the quadrature points.
	      //
	      //
	      qxyz = build_qxyz * cooelm;

	      //
	      // Evaluation of the user function.
	      //
	      

	      //
	      // Form the flux functions.
	      //
	      for (I i=0;i<nq;++i)
		{
		  f.x[i] *= ( udotn.x[i] < ((R)0.0) ) ? -udotn.x[i] * xu_[0] : ((R)0.0) ;
		}
	      
	      //
	      // Form the residu.
	      //
	      lrhs = qbasis * f;

	    }
	}
    }
#endif
  };

  
  void boundary_condition(CG_VAR& 	velocity_,
			  DG_VAR&       residual_)
  {    
    ns_mesh * mesh = velocity_.m_mesh;
    vector_handle normal;
    matrix_handle xyz;
    matrix_handle uvw;
    double normal_values[6];
    double cooelm[32];
    double uvw_values[32];
    matrix_handle_def(&xyz,m_nx,dim,cooelm,m_nx);
    matrix_handle_def(&uvw,m_nu,dim,uvw_values,m_nu);
    vector_handle_def(&normal,2,normal_values,1);
    
    DATA* data = this->CreateData();
    for (I jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {

		//
		//
		//
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const I jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		normal_values[0]	= mesh->normaledge[2*jedge+0];
		normal_values[1] 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];
		//

		//
		//
		//
		velocity_.dofselm(jelm, 0, uvw_values, 1);
		velocity_.dofselm(jelm, 1, &uvw_values[m_nu], 1);
#if 0
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];;
		  }
		for (I k=0;k<m_nu;++k)
		  {
		    uvw_values[m_nu+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		  }
#endif
		//
		
		boundary_condition(jadj,
				   &xjac,
				   normal,
				   xyz,
				   uvw,
				   data,
				   [this](I n_,cst_pR pos,I posoff,pR f,I foff)
				   {
				     // std::cout << "----" << std::endl;
				     for (I k=0;k<n_;++k)
				       {
#if 1
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];
#if 0
					 if (x<1.0e-10)
					   {
					     f[foff*k] = sin(12.0*y);
					   }
#endif
#if 1
					 //					 std::cout << " " << x << " " << y << std::endl;
					 if ( (x < 1.0e-13) && (y>0.5 && y<=1.0) )
					   {
					 f[foff*k] = 0.0;
				       }
					 else if ( (x < 1.0e-13) && (y<=0.5 && y >0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x>0.5 && x<=1.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x<=0.5 && x >=0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (x > 1.0-1.0e-13) )
					   {
					 double z  = sin(acos(-1.0)*y);
					 f[foff*k] = z*z*z;
				       }
					 else
					   {
					     f[foff*k] = 0.0;
				       }
#endif
#endif

				       }
#if 0
				     
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = y - 0.5; // exp(-sin(32.0*y)*y);
				       }
				     Err e;
				     SmoothedHeaviside_eval	(&this->m_hea,
								 n_,
								 f,
								 &e);
#endif
				     
#if 0
				     for (I k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = exp(-sin(32.0*y)*y);
				       }
#endif
				   });

		//		residual_.setdofselm(jelm,0,data->lrhs.x,1);
		WLA::matrix_h r = residual_.matrix(); 
		for (I k=0;k<m_ntest;++k)
		  {
		    //		    residual_.m_values3.x[jelm * residual_.m_values3.ld + k] += data->lrhs.x[k];
		    r.x[jelm * r.ld + k] += data->lrhs.x[k];
		  }		
	      }
	  }
      }
    
  }

};
