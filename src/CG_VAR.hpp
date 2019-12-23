#pragma once
#include "VAR.hpp"

struct CG_VAR : public VAR<CG_VAR>
{
private: WLA::matrix_h 	m_values3;
  pR 			m_memory3;
  
public:
  ns_mesh * 	m_mesh;
  I		m_ncomponents;
  I		m_ndofsPerComponent;
  mkS_st	m_shape;
  bool 		m_owner;
  
  template <typename fct_t>
  void ih(fct_t f)
  {
    R uvw[3];
    I nddlu;
    if (mkS_k(&m_shape) ==1)
      {
	nddlu = m_mesh->m_numEntities[0];
      }
    else if (mkS_k(&m_shape) ==2)
      {
	nddlu = m_mesh->nddlP6;
      }
    else
      {
	fprintf(stderr," ih: non supported\n");
	exit(1);
      }
    I nddlv	= nddlu;
    for (int i=0;i<nddlu;++i)
      {
	double x= m_mesh->coo[2*i+0];
	double y= m_mesh->coo[2*i+1];
	f(x,y,uvw);
	m_memory3[i] = uvw[0];
	m_memory3[nddlu+i] = uvw[1];
      }
  };
  
  inline void 		clear()
  {
    this->m_values3.clear();
  };
  
  inline cst_mkS shape() const
  {
    return &this->m_shape;
  };

  inline void dofselm(I id, I icomp, pR dofs, I inc) const
  {
    I n = ndofselm();
    for (I i=0;i<n;++i)
      {
	dofs[i*inc] = m_values3.x[m_values3.ld * icomp + m_mesh->cnc[6*id+i]-1];
      }
  };
  
  inline void setdofselm(I id, I icomp, cst_pR dofs, I inc) 
  {
    I n = ndofselm();
    for (I i=0;i<n;++i)
      {
	m_values3.x[m_values3.ld * icomp + m_mesh->cnc[6*id+i]-1] = dofs[i*inc];
      }
  };
  
  inline I ndofselm() const 	{ return mkS_n(&m_shape); }
  inline I ndofs() const	{ return this->m_ndofsPerComponent * this->m_ncomponents; }
  
  inline CG_VAR(ns_mesh * mesh_, I ncomponents, mkS shape) : m_mesh(mesh_)
  {
    if (mkS_k(shape)==2)
      {
	this->m_ndofsPerComponent = mesh_->nddlP6;
      }
    else if (mkS_k(shape)==1)
      {
	this->m_ndofsPerComponent = mesh_->m_numEntities[0];
      }
    else
      {
	fprintf(stderr,"invalid shape for CG_VAR\n");
	exit(1);
      }
    this->m_ncomponents = ncomponents;
    memcpy(&this->m_shape,shape,sizeof(mkS_st));
    this->m_memory3 = (pR)calloc(this->m_ndofsPerComponent * this->m_ncomponents,sizeof(R));
    this->m_owner = true;
    //    this->m_values3 = WLA::matrix_h(this->m_ndofsPerComponent,this->m_ncomponents,this->m_memory3,this->m_ndofsPerComponent);
    this->m_values3.x = this->m_memory3;
    this->m_values3.m = this->m_ndofsPerComponent;
    this->m_values3.n = this->m_ncomponents;
    this->m_values3.ld = this->m_ndofsPerComponent;
    //    matrix_handle_def(&this->m_values3,this->m_ndofsPerComponent,this->m_ncomponents,this->m_memory3,this->m_ndofsPerComponent);
  };
  
  inline CG_VAR(ns_mesh * mesh_,I ncomponents, mkS shape,pR mem_,pR memory_,I ld) : m_mesh(mesh_)
  {
    if (mkS_k(shape)==2)
      {
	m_ndofsPerComponent = mesh_->nddlP6;
      }
    else if (mkS_k(shape)==1)
      {
	m_ndofsPerComponent = mesh_->m_numEntities[0];
      }
    else
      {
	fprintf(stderr,"invalid shape for CG_VAR\n");
	exit(1);
      }
    m_ncomponents = ncomponents;
    if (ld < this->m_ndofsPerComponent)
      {
	std::cerr << "Invalid ld(= " << ld << ") < ndofsPerComponent (= " << this->m_ndofsPerComponent << " )" <<  std::endl;
      }
    memcpy(&this->m_shape,shape,sizeof(mkS_st));
    this->m_owner = false;
    this->m_memory3 = memory_;
    this->m_values3.x = this->m_memory3;
    this->m_values3.m = this->m_ndofsPerComponent;
    this->m_values3.n = this->m_ncomponents;
    this->m_values3.ld = ld;
    //    matrix_handle_def(&this->m_values3,this->m_ndofsPerComponent,this->m_ncomponents,this->m_memory3,ld);
  };

  inline ~CG_VAR()
  {
    if (this->m_owner && this->m_memory3)
      {
	free(this->m_memory3);
	this->m_memory3 = NULL;
      }
  };
  
};
