#pragma once

struct DG_Jacobian
{
  I  m_nelm;
  I  m_nfaceinelm;
  I  m_ni;
  I  m_nj;
  I  m_niXnj;
  pR m_x;

  
  I  m_n;
  I  m_m;
  I  m_nc;
  pI m_begin;
  pI m_index;
  pR m_values;

  
  void spy(const char * filename)
  {
    FILE * f = fopen(filename,"w");
    for (I i = 0;i<m_n;++i)
      {
	for (I at = m_begin[i];at<m_begin[i+1];++at)
	  {
	      fprintf(f," " ifmt " " ifmt "\n",m_n - i, m_index[at] + 1);
	      //  fprintf(f," " ifmt " " ifmt "\n",i, m_index[at]);
	  }
	fprintf(f,"\n");
      }
    fclose(f);
  };
  
  DG_Jacobian(I 	nelm_,
	      I 	nfaceinelm_,
	      I 	ni_,
	      I 	nj_,
	      cst_pI 	adj_,
	      I 	adjoff_)
  {
    this->m_ni 		= ni_;
    this->m_nj 		= nj_;
    this->m_niXnj 	= ni_*nj_;
    
    this->m_nelm			= nelm_;
    this->m_nfaceinelm			= nfaceinelm_;
    this->m_x 				= (pR)malloc(sizeof(R) * this->m_niXnj * (this->m_nfaceinelm + 1) * this->m_nelm);
    
    this->m_n 				= this->m_ni * this->m_nelm;
    this->m_m 				= this->m_nj * this->m_nelm;
    this->m_nc 				= this->m_ni * this->m_nj * (this->m_nfaceinelm + 1) * this->m_nelm;
    this->m_begin 			= (pI)calloc(this->m_n+1,sizeof(I));
    this->m_index 			= (pI)malloc(this->m_nc*sizeof(I));
    this->m_values 			= (pR)malloc(this->m_nc*sizeof(R));
    
    for (I ielm=0;ielm<nelm_;++ielm)
      {
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_ni;++k)
	  {
	    m_begin[ielm * m_ni + k + 1] += m_nj;
	  }
#if 1
	//
	// Extra diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_*ielm+localFaceIndex];
	    if (nei)
	      {
		for (I k = 0;k <  this->m_ni;++k)
		  {
		    m_begin[ielm * m_ni + k + 1] += m_nj;
		  }
	      }
	  }
#endif
      }
    
    for (I i=2;i<=this->m_n;++i)
      {
	m_begin[i]+=m_begin[i-1];
      }
#if 0
    fprintf(stdout," begin[" ifmt "] = " ifmt "\n",i,m_begin[i]);
    fprintf(stdout," "ifmt" \n",m_begin[this->m_n]);
    exit(1);
#endif
    for (I ielm=0;ielm<nelm_;++ielm)
      {
#if 1
	//
	// Before diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei-1 < ielm)
		  {
		    for (I k = 0;k <  this->m_ni;++k)
		      {
			for (I j = 0;j <  this->m_nj;++j)
			  {					    
			    m_index[m_begin[ielm * m_ni + k]] = (nei-1) * m_nj + j;
			    m_begin[ielm * m_ni + k]+=1;
			  }
		      }
		  }
	      }
	  }	
#endif	
	//
	// Diagonal.
	//
	for (I k = 0;k <  this->m_ni;++k)
	  {
	    for (I j = 0;j <  this->m_nj;++j)
	      {		
		m_index[m_begin[ielm * m_ni + k]] = ielm * m_nj + j;
		m_begin[ielm * m_ni + k]+=1;
	      }
	  }
#if 1
	//
	// After diagonal.
	//
	for (I localFaceIndex=0;localFaceIndex<this->m_nfaceinelm;++localFaceIndex)
	  {
	    I nei = adj_[adjoff_*ielm+localFaceIndex];
	    if (nei)
	      {
		if (nei - 1 > ielm)
		  {
		    for (I k = 0;k <  this->m_ni;++k)
		      {
			for (I j = 0;j <  this->m_nj;++j)
			  {					    
			    m_index[m_begin[ielm * m_ni + k]] = (nei-1) * m_nj + j;
			    m_begin[ielm * m_ni + k]+=1;
			  }
		      }
		  }
	      }
	  }
#endif	
      }

    for (I i=this->m_n;i>0;--i)
      {
	m_begin[i] = m_begin[i-1];
      }
    m_begin[0] = 0;

    for (I i=0;i<this->m_n;++i)
      {
	qsort(&m_index[m_begin[i]],
	      m_begin[i+1]-m_begin[i],
	      sizeof(I),
	      [](const void * a ,const void * b)
	      {
		cst_pI a_ = (cst_pI)a;
		cst_pI b_ = (cst_pI)b;
		if (a_[0] < b_[0] )
		  {
		    return -1;
		  }
		else if (a_[0] > b_[0] )
		  {
		    return 1;
		  }
		else
		  {
		    return 0;
		  }
	      });
      }
    this->spy("dg.txt");
  };
  
  void addelm(I 	ielm_,
	      I 	jelm_,
	      R         s_,
	      pR 	block_x_,
	      I  	block_off_)
  {
    for (I i = 0;i<m_ni;++i)
      {
	I bound = m_begin[ielm_*m_ni + i + 1];
	for (I at = m_begin[ielm_*m_ni + i];at<bound;++at)
	  {
	    if (m_index[at] == jelm_ * m_nj)
	      {
		for (I k=0;k<m_nj;++k)
		  {
		    m_values[at+k] = block_x_[block_off_*i+k] + s_ * m_values[at+k];
		  }
		break;
	      }
	  }
      }
  };

};
