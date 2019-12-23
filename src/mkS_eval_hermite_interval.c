#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"


#if __mk_debug__

void mkS_hermite_interval_(cst_mkI 	n,
			   mkR		r,
			   cst_mkI 	roff_,
			   cst_mkR 	p,
			   cst_mkI 	poff_)
{
  __mmkS_hermite_interval(n[0],r,roff_[0],p,poff_[0]);
}

void mkS_hermite_interval_dt_(cst_mkI 	n,
			      mkR	r,
			      cst_mkI 	roff_,
			      cst_mkR 	p,
			      cst_mkI 	poff_)
{
  __mmkS_hermite_interval_dt(n[0],r,roff_[0],p,poff_[0]);
}

void mkS_hermite_interval_dtt_(cst_mkI 	n,
			       mkR	r,
			       cst_mkI 	roff_,
			       cst_mkR 	p,
			       cst_mkI 	poff_)
{
  __mmkS_hermite_interval_dtt(n[0],r,roff_[0],p,poff_[0]);
}

#endif

