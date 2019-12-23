#pragma once
template <typename impl_t>
struct VAR
{  

  inline I 		ndofselm() const
  {
    return static_cast<const impl_t&>(*this).ndofselm();
  };
  
  inline I 		ndofs() const
  {
    return static_cast<const impl_t&>(*this).ndofs();
  };
  
  inline void 		dofselm(I 	id,
				I 	icomp,
				pR 	dofs,
				I 	inc) 	const
  {
    return static_cast<const impl_t&>(*this).dofselm(id,icomp,dofs,inc);
  };
  
  inline cst_mkS 	shape() const
  {
    return static_cast<const impl_t&>(*this).shape();
  };

  inline void 		clear()
  {
    return static_cast<impl_t&>(*this).clear();
  };
  
  inline void 		setdofselm(I 		id,
				   I 		icomp,
				   cst_pR 	dofs,
				   I 		inc)
  {
    return static_cast<impl_t&>(*this).setdofselm(id,icomp,dofs,inc);
  };
  
};
