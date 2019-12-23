#pragma once
#include "VAR.hpp"

struct DG_VAR : public VAR<DG_VAR>
{

private:
  WLA::matrix_h m_values3;
  pR 		m_memory;
  mkS_st	m_shape;
  bool 		m_owner;
  ns_mesh * 	m_mesh;

public:

  inline const ns_mesh* 	mesh()		const;
  inline const WLA::matrix_h& 	matrix() 	const;
  inline cst_mkS 		shape() 	const;
  inline I 			ndofselm() 	const;
  inline I 			ndofs() 	const;
  
  inline ns_mesh* 		mesh();
  inline WLA::matrix_h& 	matrix();
  inline void 			clear();

  inline DG_VAR& operator *= 	(double alpha_);
  inline DG_VAR& operator += 	(const DG_VAR&that);
  inline DG_VAR& operator -= 	(const DG_VAR&that);
  inline DG_VAR& operator = 	(const temp_jacvar&that);;

  inline DG_VAR(ns_mesh * 	mesh_,
		mkS 		shape);

  inline DG_VAR(ns_mesh * 	mesh_,
		mkS 		shape,
		pR 		mem_,
		pR 		memory_,
		I 		ld);
  
  inline ~DG_VAR();
  
  inline void 	dofselm		(I id, I icomp, pR dofs, I inc) const;
  inline I 	compute_nrms	(double nrms[]) const;
  inline void 	setdofselm	(I id, I icomp, cst_pR dofs, I inc);
  
};

inline const ns_mesh* 		DG_VAR::mesh()	 const 	{ return this->m_mesh; };
inline ns_mesh* 		DG_VAR::mesh() 		{ return this->m_mesh; };
inline cst_mkS 			DG_VAR::shape()  const 	{ return &this->m_shape; };
inline WLA::matrix_h& 		DG_VAR::matrix() 	{ return this->m_values3; };
inline const WLA::matrix_h& 	DG_VAR::matrix() const 	{ return this->m_values3; };

inline DG_VAR& DG_VAR::operator *= (double alpha_)
{
  this->m_values3 *= alpha_;
  return *this;
};

inline void DG_VAR::clear()
{
  this->m_values3.clear();
};

inline DG_VAR& DG_VAR::operator += (const DG_VAR&that)
{
  this->m_values3 += that.m_values3;
  return *this;
};

inline DG_VAR& DG_VAR::operator -= (const DG_VAR&that)
{
  this->m_values3 -= that.m_values3;
  return *this;
};
  
inline I DG_VAR::ndofselm() 	const 	{ return mkS_n(&m_shape); }
inline I DG_VAR::ndofs() 	const	{ return this->ndofselm() * m_mesh->nelm; }

inline DG_VAR::DG_VAR(ns_mesh * mesh_, mkS shape) : m_mesh(mesh_)
{
  memcpy(&this->m_shape,shape,sizeof(mkS_st));
  this->m_memory = (pR)calloc(m_mesh->nelm*mkS_n(shape),sizeof(R));
  this->m_owner = true;
  WLA::matrix_h::define(this->m_values3,mkS_n(&this->m_shape),m_mesh->nelm,this->m_memory,mkS_n(&this->m_shape));
};
  
inline DG_VAR::DG_VAR(ns_mesh * mesh_,mkS shape,pR mem_,pR memory_,I ld) : m_mesh(mesh_)
{
  if (ld < mkS_n(shape))
    {
      std::cerr << "Invalid ld(= " << ld << ") < nshape (= " << mkS_n(shape) << " )" <<  std::endl;
    }
    
  memcpy(&this->m_shape,shape,sizeof(mkS_st));
  this->m_owner = false;
  this->m_memory = memory_;
  WLA::matrix_h::define(this->m_values3,mkS_n(&this->m_shape),this->m_mesh->nelm,this->m_memory,ld);
};

inline DG_VAR::~DG_VAR()
{
  if (this->m_owner && this->m_memory)
    {
      free(this->m_memory);
      this->m_memory = NULL;
    }
};

inline void DG_VAR::dofselm(I id, I icomp, pR dofs, I inc) const
{
  I n = ndofselm();
  for (I i=0;i<n;++i)
    {
      dofs[i*inc] = m_values3.x[m_values3.ld * id + i];
    }
};
  
inline void DG_VAR::setdofselm(I id, I icomp, cst_pR dofs, I inc) 
{
  I n = ndofselm();
  for (I i=0;i<n;++i)
    {
      m_values3.x[m_values3.ld * id + i] = dofs[i*inc];
    }
};
  
inline I DG_VAR::compute_nrms(double nrms[]) const
{

  I degree = mkS_k(&this->m_shape);
  I n1 = 1;
  for (I ideg=0;ideg<=degree;++ideg)
    {
      double nrm=0.0;
      for (I ielm=0;ielm<m_mesh->nelm;++ielm)
	{
	  double h1 = 0.0;
	  I start = ((ideg) * (ideg+1) )/2;
	  I bound = ((ideg+1) * (ideg+2) )/2;
	  for (I j = start;j<bound;++j)
	    {
	      double x = m_values3.x[ielm*m_values3.ld+j];
	      h1 += x*x;
	    }
	  h1 *= m_mesh->jacelm[ielm];
	  nrm+=h1;
	}
      nrms[1+ideg] = nrm;
    }

  nrms[0] = 0.0;
  for (I ideg=0;ideg<=degree;++ideg)
    {
      nrms[0] += nrms[1+ideg];
    }

  nrms[0] = sqrt(nrms[0]);
  for (I ideg=0;ideg<=degree;++ideg)
    {
      nrms[1+ideg] = sqrt(nrms[1+ideg]);
    }
    
  return 2+degree;
};
  
  

