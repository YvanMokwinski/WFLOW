#include <signal.h>
#include <pthread.h>
#include <stdio.h>
#include <math.h>

#include "Cmdline.h"
#include "Monitor.h"
#include "Blas.h"
#include "Type.h"
#include "ns_sys.h"
#include "ns_mesh.h"
#include "mkS.h"
#include "ns_config_lapack.h"
#include "ns_constantes.h"
ns_mesh * mesh;
#define DG_CONSERVATIVE_WEAK_FORM 	NO

#define DGERR_MEMORY  2
#define DGERR_USER    3

#define DG_r_lc        			0
#define DG_r_n         			1

#define DG_ires_err 			0 
#define DG_ires_convergence 		1
#define DG_ires_iter_gauss_seidel 	2
#define DG_ires_required_iw_n 		3
#define DG_ires_required_rw_n 		4
#define DG_ires_n 			5

#define DG_rres_max 			0 
#define DG_rres_nrmL2 			1 
#define DG_rres_nrmLInf			2 
#define DG_rres_areaL1 			3 
#define DG_rres_jumpL2 			4 
#define DG_rres_johnson			5 
#define DG_rres_n 			6
cst_mkS 		mkS_derivative	(mkS shape_,const I idim)
{
  if (idim==0)
    {
      return mkS_dx(shape_);
    }
  else if (idim==1)
    {
      return mkS_dy(shape_);
    }
  else if (idim==2)
    {
      return mkS_dz(shape_);
    }
  else
    {
      fprintf(stderr,"mkS_derivative error idim = " ifmt "\n",idim);
      exit(1);
    }
  return NULL;
}



// #include "DgTransport.h"


void dg_print_sol(I N,cst_pR sol,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.bb",ctmp2); }   
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"2 1 " ifmt " 2\n",N);
    { I i;
      for (i=0;i<N;++i)
	{
	  fprintf(fil,"" rfmt "\n",sol[i]);
	} } 
    fclose(fil); }
}

void dg_print_mesh(const ns_mesh*s_,const char * name_,...)
{
  const I numNodes = ns_mesh_get_numNodes(s_);
  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",s_->nelm*3);
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  ns_mesh_get_cellToNodes(s_,&i,cncelm);
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	} } 

    fprintf(fil,"Triangles\n" ifmt "\n",s_->nelm); 
    { I i;
      for (i=0;i<s_->nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,((I)0));
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}

struct vector_handle;
struct matrix_handle;

struct temp_gemm
{
  const R a;
  const char transA;
  const char transB;
  const matrix_handle&A;
  const matrix_handle&B;
};

struct temp_gemv
{
  const R a;
  const char trans;
  const matrix_handle&A;
  const vector_handle&b;
};
struct temp_scal;
struct temp_transpose
{
  const char trans;
  const matrix_handle&A;
  inline temp_scal operator * (const R&v_) const;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
};

struct temp_scal
{
  const R a;
  const char trans;
  const matrix_handle&A;
  inline temp_gemv operator * (const vector_handle&v_) const;
  inline temp_gemm operator * (const matrix_handle&v_) const;
  inline temp_gemm operator * (const temp_transpose&v_) const;
};

inline temp_gemm temp_scal::operator * (const temp_transpose&v_) const
{
  return {a,trans,v_.trans,A,v_.A};
}

inline temp_gemv temp_scal::operator * (const vector_handle&v_) const
{
  return {a,trans,A, v_};
};
inline temp_gemm temp_scal::operator * (const matrix_handle&v_) const
{
  return {a, trans,'N',A, v_};
};

inline temp_scal temp_transpose::operator * (const R&v_) const
{
  return {v_,trans,A};
};
inline temp_gemv temp_transpose::operator * (const vector_handle&v_) const
{
  return {1.0,trans,A, v_};
};
  
inline temp_gemm temp_transpose::operator * (const matrix_handle&v_) const
{
  return {1.0, trans,'N', A, v_};
};

struct matrix_handle
{
  pR  x{};
  I   n{};
  I   m{};
  I   ld{};
  matrix_handle();
  matrix_handle(I n_, I m_,pR  x_,I ld_);
  inline temp_gemm operator * (const matrix_handle&v_)const;  
  inline temp_gemv operator * (const vector_handle&v_)const;  
  inline temp_scal operator * (const R&v_)const;  
  inline temp_transpose transpose() const
  {
    return {'T',*this};
  };


  inline matrix_handle& operator=(const temp_gemm&temp)
  {
    const R r0 = 0.0;
    I k = (temp.transA == 'N') ? temp.A.m : temp.A.n;
    Blas_dgemm(&temp.transA,
	       &temp.transB,
	       &n,
	       &m,
	       &k,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.B.x,
	       &temp.B.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };

};


temp_scal operator * (const R&v_,const matrix_handle &m) 
{
  return {v_,'N',m};
};

matrix_handle::matrix_handle(){};
matrix_handle::matrix_handle(I n_, I m_,pR  x_,I ld_) : x(x_), n(n_), m(m_), ld(ld_)
  {
  };

  inline temp_gemv matrix_handle::operator * (const vector_handle&v_) const
  {
    return {1.0,'N',*this, v_};
  };

inline temp_gemm matrix_handle::operator * (const matrix_handle&v_) const
  {
    return {1.0,'N','N',*this, v_};
  };

inline temp_scal matrix_handle::operator * (const R&v_) const
  {
    return {v_,'N',*this};
  };
  

struct vector_handle
{
  pR  x{};
  I   n{};
  I   ld{};
  vector_handle(){};
  vector_handle(I n_, pR  x_,I ld_)
    : x(x_),n(n_),ld(ld_)
  {};

  template <typename F>
  inline void apply(F f)
  {
    for (I j=0;j<n;++j)
      {
	pR e = x + j*ld;
	*e = f(*e);
      }    
  };
  
  inline vector_handle& operator=(const temp_gemv&temp)
  {
    const R r0 = 0.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r0,
	       x,
	       &ld);
    return *this;
  };
  
  inline vector_handle& operator+=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &temp.a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };

  inline vector_handle& operator-=(const temp_gemv&temp)
  {
    const R r1 = 1.0;
    R a = -temp.a;
    Blas_dgemv(&temp.trans,
	       &temp.A.n,
	       &temp.A.m,
	       &a,
	       temp.A.x,
	       &temp.A.ld,
	       temp.b.x,
	       &temp.b.ld,
	       &r1,
	       x,
	       &ld);
    return *this;

  };


};

void vector_handle_def(struct vector_handle * h,I n_, pR  x_,I ld_)
{
  h->n  = n_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_def(struct matrix_handle * h,I n_, I m_,pR  x_,I ld_)
{
  h->n  = n_;
  h->m  = m_;
  h->x  = x_;
  h->ld = ld_;
};

void matrix_handle_gemv_low(const struct matrix_handle * h,const char * trans,cst_pR a,cst_pR x,cst_pI xoff,cst_pR ry,pR y,cst_pI yoff)
{

  Blas_dgemv(trans,
	     &h->n,
	     &h->m,
	     a,
	     h->x,
	     &h->ld,
	     x,
	     xoff,
	     ry,
	     y,
	     yoff);

}

void matrix_handle_gemv(const struct matrix_handle * h,const char * trans,cst_pR a,const struct vector_handle * x,cst_pR ry,struct vector_handle * y)
{
  matrix_handle_gemv_low(h,trans,a,x->x,&x->ld,ry,y->x,&y->ld);
}


void matrix_handle_gemm(const struct matrix_handle * h,const char * transA,const char * transB,cst_pR a,const struct matrix_handle * x,cst_pR ry,struct matrix_handle * y)
{
  const I s = transA[0]=='N' ? h->m : h->n;
  Blas_dgemm(transA,transB,&y->n,&y->m,&s,a,h->x,&h->ld,x->x,&x->ld,ry,y->x,&y->ld);

}

void matrix_handle_gesv(const struct matrix_handle * h,struct vector_handle * rhs,pI lcperm)
{
  I n1=1;
  I info_lapack;
  dgesv(&h->n,
	&n1,
	h->x,
	&h->ld,
	lcperm,
	rhs->x,
	&h->n,
	&info_lapack);
}




static const double tria_L1_L1[9] = {
((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)

,((double)4.1666666666666664E-2)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)

};

static const double tria_L2_L2[36] = {
  ((double)1.6666666666666666E-2)
  ,((double)-2.7777777777777779E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)0.)
  ,((double)-1.1111111111111112E-2)
  ,((double)0.)
,((double)-2.7777777777777779E-3)
,((double)1.6666666666666666E-2)
,((double)-2.7777777777777779E-3)
,((double)0.)
,((double)0.)
,((double)-1.1111111111111112E-2)

,((double)-2.7777777777777779E-3)
,((double)-2.7777777777777779E-3)
,((double)1.6666666666666666E-2)
,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)0.)

,((double)0.)
,((double)0.)
,((double)-1.1111111111111112E-2)
,((double)8.8888888888888892E-2)
,((double)4.4444444444444446E-2)
,((double)4.4444444444444446E-2)

,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)0.)
,((double)4.4444444444444446E-2)
,((double)8.8888888888888892E-2)
,((double)4.4444444444444446E-2)

,((double)0.)
,((double)-1.1111111111111112E-2)
,((double)0.)
,((double)4.4444444444444446E-2)
,((double)4.4444444444444446E-2)
,((double)8.8888888888888892E-2)

};

static const double tria_L2p_L2p[49] = {
((double)0.2388888888888889)
,((double)0.10833333333333334)
,((double)0.10833333333333334)
,((double)-0.12222222222222222)
,((double)0.12222222222222222)
,((double)-0.12222222222222222)
,((double)0.16071428571428573)
,((double)0.10833333333333334)
,((double)8.3333333333333329E-2)
,((double)4.1666666666666664E-2)
,((double)-6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)7.4999999999999997E-2)

,((double)0.10833333333333334)
,((double)4.1666666666666664E-2)
,((double)8.3333333333333329E-2)
,((double)-6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)7.4999999999999997E-2)

,((double)-0.12222222222222222)
,((double)-6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.5714285714285715E-2)

,((double)0.12222222222222222)
,((double)6.6666666666666666E-2)
,((double)6.6666666666666666E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.5714285714285715E-2)

,((double)-0.12222222222222222)
,((double)-6.6666666666666666E-2)
,((double)-6.6666666666666666E-2)
,((double)8.8888888888888892E-2)
,((double)-8.8888888888888892E-2)
,((double)8.8888888888888892E-2)
,((double)-8.5714285714285715E-2)

,((double)0.16071428571428573)
,((double)7.4999999999999997E-2)
,((double)7.4999999999999997E-2)
,((double)-8.5714285714285715E-2)
,((double)8.5714285714285715E-2)
,((double)-8.5714285714285715E-2)
,((double)0.14464285714285716)

};

static const double tria_L3_L3[100] = {
((double)5.6547619047619046E-3)
,((double)8.1845238095238097E-4)
,((double)8.1845238095238097E-4)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)8.1845238095238097E-4)
,((double)5.6547619047619046E-3)
,((double)8.1845238095238097E-4)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)2.6785714285714286E-3)

,((double)8.1845238095238097E-4)
,((double)8.1845238095238097E-4)
,((double)5.6547619047619046E-3)
,((double)2.0089285714285712E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.6785714285714286E-3)

,((double)1.3392857142857143E-3)
,((double)0.)
,((double)2.0089285714285712E-3)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)1.2053571428571429E-2)

,((double)0.)
,((double)1.3392857142857143E-3)
,((double)2.0089285714285712E-3)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)1.2053571428571429E-2)

,((double)2.0089285714285712E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)1.2053571428571429E-2)

,((double)2.0089285714285712E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)1.2053571428571429E-2)

,((double)0.)
,((double)2.0089285714285712E-3)
,((double)1.3392857142857143E-3)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)2.0089285714285716E-2)
,((double)4.0178571428571432E-2)
,((double)-1.40625E-2)
,((double)1.2053571428571429E-2)

,((double)1.3392857142857143E-3)
,((double)2.0089285714285712E-3)
,((double)0.)
,((double)2.0089285714285716E-2)
,((double)-1.0044642857142858E-2)
,((double)-4.0178571428571425E-3)
,((double)-1.0044642857142858E-2)
,((double)-1.40625E-2)
,((double)4.0178571428571432E-2)
,((double)1.2053571428571429E-2)

,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.14464285714285716)

};

static const double tria_L1_L1_L1[27] = {
((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)

,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)

};

static const double tria_L2_L2_L2[216] = {
((double)7.1428571428571426E-3)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)4.7619047619047623E-3)
,((double)1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)7.1428571428571426E-3)
,((double)-7.9365079365079365E-4)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)

,((double)-7.9365079365079365E-4)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-7.9365079365079365E-4)
,((double)-7.9365079365079365E-4)
,((double)7.1428571428571426E-3)
,((double)1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)

,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-6.3492063492063492E-3)
,((double)5.7142857142857141E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)

,((double)1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)-6.3492063492063492E-3)
,((double)0.)
,((double)0.)
,((double)1.9047619047619049E-2)
,((double)5.7142857142857141E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)

,((double)4.7619047619047623E-3)
,((double)0.)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)1.5873015873015873E-3)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)-6.3492063492063492E-3)
,((double)-1.5873015873015873E-3)
,((double)0.)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-3.1746031746031746E-3)
,((double)0.)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)0.)
,((double)-6.3492063492063492E-3)
,((double)0.)
,((double)1.9047619047619049E-2)
,((double)1.9047619047619049E-2)
,((double)5.7142857142857141E-2)

};

static const double tria_L2p_L2p_L2p[343] = {
((double)0.17857142857142858)
,((double)7.6984126984126988E-2)
,((double)7.6984126984126988E-2)
,((double)-9.3650793650793651E-2)
,((double)9.3650793650793651E-2)
,((double)-9.3650793650793651E-2)
,((double)0.11785714285714285)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)9.3650793650793651E-2)
,((double)4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.4285714285714279E-2)
,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)0.11785714285714285)
,((double)5.3571428571428568E-2)
,((double)5.3571428571428568E-2)
,((double)-6.4285714285714279E-2)
,((double)6.4285714285714279E-2)
,((double)-6.4285714285714279E-2)
,((double)0.10607142857142857)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)5.0000000000000003E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)

,((double)9.3650793650793651E-2)
,((double)4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.4285714285714279E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.4285714285714279E-2)
,((double)3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)5.7857142857142857E-2)

,((double)-9.3650793650793651E-2)
,((double)-4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)6.9841269841269843E-2)
,((double)-6.9841269841269843E-2)
,((double)6.9841269841269843E-2)
,((double)-6.4285714285714279E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.9841269841269843E-2)
,((double)-3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-4.2857142857142858E-2)
,((double)6.9841269841269843E-2)
,((double)3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-5.7142857142857141E-2)
,((double)5.7142857142857141E-2)
,((double)-5.7142857142857141E-2)
,((double)4.2857142857142858E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)

,((double)0.11785714285714285)
,((double)5.3571428571428568E-2)
,((double)5.3571428571428568E-2)
,((double)-6.4285714285714279E-2)
,((double)6.4285714285714279E-2)
,((double)-6.4285714285714279E-2)
,((double)0.10607142857142857)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)
,((double)6.4285714285714279E-2)
,((double)3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)5.7857142857142857E-2)
,((double)-6.4285714285714279E-2)
,((double)-3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.2857142857142858E-2)
,((double)-4.2857142857142858E-2)
,((double)4.2857142857142858E-2)
,((double)-5.7857142857142857E-2)
,((double)0.10607142857142857)
,((double)4.8214285714285716E-2)
,((double)4.8214285714285716E-2)
,((double)-5.7857142857142857E-2)
,((double)5.7857142857142857E-2)
,((double)-5.7857142857142857E-2)
,((double)0.10650974025974026)

};

static const double tria_L3_L3_L3[1000] = {
((double)2.5162337662337662E-3)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)7.3051948051948055E-4)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)-1.2175324675324675E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)9.8620129870129864E-4)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)9.8620129870129864E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)6.5746753246753243E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)-1.0957792207792207E-3)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)8.7662337662337668E-4)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)2.5162337662337662E-3)
,((double)1.8939393939393939E-4)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)7.3051948051948055E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.2175324675324675E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)6.5746753246753243E-4)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.0957792207792207E-3)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)6.5746753246753243E-4)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)9.8620129870129864E-4)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)9.8620129870129864E-4)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)8.7662337662337668E-4)

,((double)1.8939393939393939E-4)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)-1.2175324675324675E-5)
,((double)5.4112554112554113E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-1.2175324675324675E-5)
,((double)1.8939393939393939E-4)
,((double)1.8939393939393939E-4)
,((double)2.5162337662337662E-3)
,((double)4.0584415584415584E-5)
,((double)4.0584415584415584E-5)
,((double)-6.0876623376623375E-4)
,((double)1.5827922077922079E-3)
,((double)1.5827922077922079E-3)
,((double)-6.0876623376623375E-4)
,((double)7.3051948051948055E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)6.5746753246753243E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)8.7662337662337668E-4)

,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)9.8620129870129864E-4)
,((double)2.1915584415584414E-3)
,((double)1.3149350649350649E-3)
,((double)1.497564935064935E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)1.2820616883116883E-2)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.9310064935064934E-3)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.3011363636363637E-3)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.2735389610389614E-3)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)8.5470779220779228E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)5.9172077922077923E-3)

,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)9.8620129870129864E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.9310064935064934E-3)
,((double)1.3149350649350649E-3)
,((double)2.1915584415584414E-3)
,((double)1.497564935064935E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)1.2820616883116883E-2)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)8.5470779220779228E-3)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)-2.3011363636363637E-3)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-4.2735389610389614E-3)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)

,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-1.4204545454545455E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)9.8620129870129864E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-1.0957792207792207E-3)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-6.0876623376623375E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)6.5746753246753243E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)-4.2735389610389614E-3)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)8.5470779220779228E-3)
,((double)1.497564935064935E-3)
,((double)2.1915584415584414E-3)
,((double)1.3149350649350649E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)1.2820616883116883E-2)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-2.3011363636363637E-3)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)5.9172077922077923E-3)

,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)6.4935064935064935E-5)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)9.8620129870129864E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)1.4813311688311687E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)6.5746753246753243E-4)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.5827922077922079E-3)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-1.0957792207792207E-3)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.3011363636363637E-3)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-1.6436688311688311E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)1.497564935064935E-3)
,((double)1.3149350649350649E-3)
,((double)2.1915584415584414E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)1.2820616883116883E-2)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)8.5470779220779228E-3)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-4.2735389610389614E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)

,((double)-6.0876623376623375E-4)
,((double)-1.4204545454545455E-5)
,((double)1.4813311688311687E-4)
,((double)-5.4788961038961036E-4)
,((double)2.1915584415584417E-4)
,((double)-3.1047077922077924E-4)
,((double)7.8530844155844158E-4)
,((double)1.3149350649350649E-3)
,((double)-1.6436688311688311E-3)
,((double)6.5746753246753243E-4)
,((double)-1.4204545454545455E-5)
,((double)4.0584415584415584E-5)
,((double)6.4935064935064935E-5)
,((double)-3.1047077922077924E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)7.8530844155844158E-4)
,((double)1.497564935064935E-3)
,((double)-5.8441558441558442E-4)
,((double)9.8620129870129864E-4)
,((double)1.4813311688311687E-4)
,((double)6.4935064935064935E-5)
,((double)1.5827922077922079E-3)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-5.4788961038961036E-4)
,((double)1.0957792207792207E-3)
,((double)2.1915584415584414E-3)
,((double)-1.6436688311688311E-3)
,((double)-1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.2735389610389614E-3)
,((double)2.1915584415584417E-4)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)-2.3011363636363637E-3)
,((double)-3.1047077922077924E-4)
,((double)-4.7483766233766234E-4)
,((double)-5.4788961038961036E-4)
,((double)6.5746753246753243E-4)
,((double)-3.2873376623376621E-4)
,((double)1.1505681818181819E-3)
,((double)-2.6298701298701297E-3)
,((double)-4.76663961038961E-3)
,((double)1.4793019480519481E-3)
,((double)-4.2735389610389614E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)1.0957792207792207E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)8.5470779220779228E-3)
,((double)1.3149350649350649E-3)
,((double)1.497564935064935E-3)
,((double)2.1915584415584414E-3)
,((double)1.1505681818181819E-3)
,((double)-2.1367694805194807E-3)
,((double)-4.76663961038961E-3)
,((double)7.3965097402597406E-3)
,((double)2.2189529220779221E-2)
,((double)-1.4793019480519481E-3)
,((double)1.2820616883116883E-2)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)5.9172077922077923E-3)

,((double)1.5827922077922079E-3)
,((double)6.4935064935064935E-5)
,((double)1.4813311688311687E-4)
,((double)1.0957792207792207E-3)
,((double)-5.4788961038961036E-4)
,((double)1.8262987012987013E-5)
,((double)-4.7483766233766234E-4)
,((double)-1.6436688311688311E-3)
,((double)2.1915584415584414E-3)
,((double)-1.0957792207792207E-3)
,((double)6.4935064935064935E-5)
,((double)4.0584415584415584E-5)
,((double)-1.4204545454545455E-5)
,((double)7.8530844155844158E-4)
,((double)-4.7483766233766234E-4)
,((double)1.8262987012987013E-5)
,((double)-3.1047077922077924E-4)
,((double)-5.8441558441558442E-4)
,((double)1.497564935064935E-3)
,((double)9.8620129870129864E-4)
,((double)1.4813311688311687E-4)
,((double)-1.4204545454545455E-5)
,((double)-6.0876623376623375E-4)
,((double)7.8530844155844158E-4)
,((double)-3.1047077922077924E-4)
,((double)2.1915584415584417E-4)
,((double)-5.4788961038961036E-4)
,((double)-1.6436688311688311E-3)
,((double)1.3149350649350649E-3)
,((double)6.5746753246753243E-4)
,((double)1.0957792207792207E-3)
,((double)7.8530844155844158E-4)
,((double)7.8530844155844158E-4)
,((double)7.3965097402597406E-3)
,((double)-2.6298701298701297E-3)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.6298701298701297E-3)
,((double)7.3965097402597406E-3)
,((double)8.5470779220779228E-3)
,((double)-5.4788961038961036E-4)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-4.2735389610389614E-3)
,((double)1.8262987012987013E-5)
,((double)1.8262987012987013E-5)
,((double)2.1915584415584417E-4)
,((double)-3.2873376623376621E-4)
,((double)-3.2873376623376621E-4)
,((double)-2.1367694805194807E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.1367694805194807E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.7483766233766234E-4)
,((double)-3.1047077922077924E-4)
,((double)-5.4788961038961036E-4)
,((double)-3.2873376623376621E-4)
,((double)6.5746753246753243E-4)
,((double)1.4793019480519481E-3)
,((double)-4.76663961038961E-3)
,((double)-2.6298701298701297E-3)
,((double)1.1505681818181819E-3)
,((double)-4.2735389610389614E-3)
,((double)-1.6436688311688311E-3)
,((double)-5.8441558441558442E-4)
,((double)-1.6436688311688311E-3)
,((double)-2.6298701298701297E-3)
,((double)1.4793019480519481E-3)
,((double)1.4793019480519481E-3)
,((double)-2.6298701298701297E-3)
,((double)-1.4793019480519481E-3)
,((double)-1.4793019480519481E-3)
,((double)-4.9310064935064934E-3)
,((double)2.1915584415584414E-3)
,((double)1.497564935064935E-3)
,((double)1.3149350649350649E-3)
,((double)7.3965097402597406E-3)
,((double)-4.76663961038961E-3)
,((double)-2.1367694805194807E-3)
,((double)1.1505681818181819E-3)
,((double)-1.4793019480519481E-3)
,((double)2.2189529220779221E-2)
,((double)1.2820616883116883E-2)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)5.9172077922077923E-3)

,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)8.7662337662337668E-4)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)-1.2175324675324675E-5)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)8.7662337662337668E-4)
,((double)-1.2175324675324675E-5)
,((double)-1.2175324675324675E-5)
,((double)7.3051948051948055E-4)
,((double)9.8620129870129864E-4)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)8.7662337662337668E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)5.9172077922077923E-3)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)6.5746753246753243E-4)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)5.9172077922077923E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)-1.0957792207792207E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)5.9172077922077923E-3)
,((double)6.5746753246753243E-4)
,((double)9.8620129870129864E-4)
,((double)-1.0957792207792207E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)8.5470779220779228E-3)
,((double)1.2820616883116883E-2)
,((double)-4.9310064935064934E-3)
,((double)5.9172077922077923E-3)
,((double)-1.0957792207792207E-3)
,((double)9.8620129870129864E-4)
,((double)6.5746753246753243E-4)
,((double)8.5470779220779228E-3)
,((double)-4.2735389610389614E-3)
,((double)-2.3011363636363637E-3)
,((double)-4.2735389610389614E-3)
,((double)-4.9310064935064934E-3)
,((double)1.2820616883116883E-2)
,((double)5.9172077922077923E-3)
,((double)8.7662337662337668E-4)
,((double)8.7662337662337668E-4)
,((double)8.7662337662337668E-4)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)5.9172077922077923E-3)
,((double)0.10650974025974026)

};

static const double tria_L1_L2_L2[108] = {
((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)2.3809523809523812E-3)
,((double)3.9682539682539683E-4)
,((double)-3.1746031746031746E-3)
,((double)-1.5873015873015873E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)2.3809523809523812E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)3.8095238095238099E-2)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)4.7619047619047623E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)3.8095238095238099E-2)
,((double)2.3809523809523812E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)2.3809523809523812E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)-4.7619047619047623E-3)
,((double)3.8095238095238099E-2)
,((double)1.9047619047619049E-2)
,((double)1.2698412698412698E-2)
,((double)-4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)1.9047619047619049E-2)
,((double)3.8095238095238099E-2)
,((double)1.2698412698412698E-2)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)

,((double)2.3809523809523812E-3)
,((double)3.9682539682539683E-4)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)3.9682539682539683E-4)
,((double)2.3809523809523812E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.1904761904761904E-2)
,((double)-1.5873015873015873E-3)
,((double)4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)-1.5873015873015873E-3)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-4.7619047619047623E-3)
,((double)-3.1746031746031746E-3)
,((double)4.7619047619047623E-3)
,((double)1.2698412698412698E-2)
,((double)3.8095238095238099E-2)
,((double)1.9047619047619049E-2)
,((double)-3.1746031746031746E-3)
,((double)-4.7619047619047623E-3)
,((double)4.7619047619047623E-3)
,((double)1.2698412698412698E-2)
,((double)1.9047619047619049E-2)
,((double)3.8095238095238099E-2)

};

static const double tria_L1_L2p_L2p[147] = {
((double)8.4920634920634924E-2)
,((double)2.7777777777777776E-2)
,((double)2.7777777777777776E-2)
,((double)-2.3809523809523808E-2)
,((double)2.3809523809523808E-2)
,((double)-2.3809523809523808E-2)
,((double)5.3571428571428568E-2)
,((double)2.7777777777777776E-2)
,((double)1.6666666666666666E-2)
,((double)8.3333333333333332E-3)
,((double)-1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)2.1428571428571429E-2)
,((double)2.7777777777777776E-2)
,((double)8.3333333333333332E-3)
,((double)1.6666666666666666E-2)
,((double)-1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)2.1428571428571429E-2)
,((double)-2.3809523809523808E-2)
,((double)-1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-2.1428571428571429E-2)
,((double)2.3809523809523808E-2)
,((double)1.1111111111111112E-2)
,((double)1.1111111111111112E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)2.1428571428571429E-2)
,((double)-2.3809523809523808E-2)
,((double)-1.1111111111111112E-2)
,((double)-1.1111111111111112E-2)
,((double)1.2698412698412698E-2)
,((double)-1.2698412698412698E-2)
,((double)1.2698412698412698E-2)
,((double)-2.1428571428571429E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)2.1428571428571429E-2)
,((double)-2.1428571428571429E-2)
,((double)2.1428571428571429E-2)
,((double)-2.1428571428571429E-2)
,((double)4.8214285714285716E-2)
,((double)7.6984126984126988E-2)
,((double)5.0000000000000003E-2)
,((double)3.0555555555555555E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)5.0000000000000003E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)3.3333333333333333E-2)
,((double)2.2222222222222223E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-3.3333333333333333E-2)
,((double)-2.2222222222222223E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)3.214285714285714E-2)
,((double)2.1428571428571429E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

,((double)7.6984126984126988E-2)
,((double)3.0555555555555555E-2)
,((double)5.0000000000000003E-2)
,((double)-4.9206349206349205E-2)
,((double)4.9206349206349205E-2)
,((double)-4.9206349206349205E-2)
,((double)5.3571428571428568E-2)
,((double)3.0555555555555555E-2)
,((double)1.6666666666666666E-2)
,((double)1.6666666666666666E-2)
,((double)-2.2222222222222223E-2)
,((double)2.2222222222222223E-2)
,((double)-2.2222222222222223E-2)
,((double)2.1428571428571429E-2)
,((double)5.0000000000000003E-2)
,((double)1.6666666666666666E-2)
,((double)5.0000000000000003E-2)
,((double)-3.3333333333333333E-2)
,((double)3.3333333333333333E-2)
,((double)-3.3333333333333333E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)4.9206349206349205E-2)
,((double)2.2222222222222223E-2)
,((double)3.3333333333333333E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.214285714285714E-2)
,((double)-4.9206349206349205E-2)
,((double)-2.2222222222222223E-2)
,((double)-3.3333333333333333E-2)
,((double)3.8095238095238099E-2)
,((double)-3.8095238095238099E-2)
,((double)3.8095238095238099E-2)
,((double)-3.214285714285714E-2)
,((double)5.3571428571428568E-2)
,((double)2.1428571428571429E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)3.214285714285714E-2)
,((double)-3.214285714285714E-2)
,((double)4.8214285714285716E-2)

};

static const double tria_L1_L3_L3[300] = {
((double)4.464285714285714E-3)
,((double)3.720238095238095E-4)
,((double)3.720238095238095E-4)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)0.)
,((double)3.720238095238095E-4)
,((double)5.9523809523809529E-4)
,((double)7.4404761904761911E-5)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)7.4404761904761911E-5)
,((double)5.9523809523809529E-4)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)-1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)2.2321428571428571E-4)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)-8.9285714285714283E-4)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)2.6785714285714286E-3)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)4.8214285714285716E-2)
,((double)5.9523809523809529E-4)
,((double)3.720238095238095E-4)
,((double)7.4404761904761911E-5)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)4.464285714285714E-3)
,((double)3.720238095238095E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)0.)
,((double)7.4404761904761911E-5)
,((double)3.720238095238095E-4)
,((double)5.9523809523809529E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)-1.3392857142857143E-3)
,((double)2.2321428571428571E-4)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)1.1160714285714285E-3)
,((double)2.6785714285714286E-3)
,((double)1.5625000000000001E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)1.2053571428571429E-2)
,((double)1.5625000000000001E-3)
,((double)2.6785714285714286E-3)
,((double)1.1160714285714285E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)1.2053571428571429E-2)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)-8.9285714285714283E-4)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)4.8214285714285716E-2)

,((double)5.9523809523809529E-4)
,((double)7.4404761904761911E-5)
,((double)3.720238095238095E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)-8.9285714285714283E-4)
,((double)1.3392857142857143E-3)
,((double)7.4404761904761911E-5)
,((double)5.9523809523809529E-4)
,((double)3.720238095238095E-4)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)-8.9285714285714283E-4)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.2321428571428571E-4)
,((double)1.3392857142857143E-3)
,((double)3.720238095238095E-4)
,((double)3.720238095238095E-4)
,((double)4.464285714285714E-3)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)2.6785714285714286E-3)
,((double)2.6785714285714286E-3)
,((double)-1.3392857142857143E-3)
,((double)0.)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)2.2321428571428571E-4)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-4.4642857142857141E-4)
,((double)2.2321428571428571E-4)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)2.2321428571428571E-4)
,((double)-8.9285714285714283E-4)
,((double)-1.3392857142857143E-3)
,((double)-2.0089285714285712E-3)
,((double)4.0178571428571425E-3)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)0.)
,((double)0.)
,((double)1.5625000000000001E-3)
,((double)1.1160714285714285E-3)
,((double)2.6785714285714286E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)2.4107142857142858E-2)
,((double)1.2053571428571429E-2)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)1.1160714285714285E-3)
,((double)1.5625000000000001E-3)
,((double)2.6785714285714286E-3)
,((double)-2.0089285714285712E-3)
,((double)-2.0089285714285712E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)2.4107142857142858E-2)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)-8.9285714285714283E-4)
,((double)2.2321428571428571E-4)
,((double)-1.3392857142857143E-3)
,((double)4.0178571428571425E-3)
,((double)-2.0089285714285712E-3)
,((double)0.)
,((double)-6.0267857142857146E-3)
,((double)-6.0267857142857146E-3)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)1.3392857142857143E-3)
,((double)1.3392857142857143E-3)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)0.)
,((double)1.2053571428571429E-2)
,((double)1.2053571428571429E-2)
,((double)0.)
,((double)4.8214285714285716E-2)

};

static const double tria_L2_L3_L3[600] = {
  ((double)3.2738095238095239E-3)
  ,((double)2.2321428571428571E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)0.)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.9841269841269841E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-1.9841269841269841E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-1.9841269841269841E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.2738095238095239E-3)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)3.4722222222222222E-5)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.9841269841269841E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.3392857142857144E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)0.)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)0.)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)2.6785714285714287E-4)
  ,((double)0.)
  ,((double)2.6785714285714287E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)-1.9841269841269841E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)2.2321428571428571E-4)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)3.4722222222222222E-5)
  ,((double)-1.9841269841269841E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-4.4642857142857141E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)2.2321428571428571E-4)
  ,((double)3.2738095238095239E-3)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.3392857142857144E-4)
  ,((double)-2.2321428571428571E-4)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.3392857142857144E-4)
  ,((double)0.)
  ,((double)-2.2321428571428571E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.910714285714286E-4)
  ,((double)-4.4642857142857141E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.8035714285714288E-4)
  ,((double)7.5892857142857142E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)7.5892857142857142E-4)
  ,((double)5.8035714285714288E-4)
  ,((double)2.6785714285714286E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-4.4642857142857141E-4)
  ,((double)-4.910714285714286E-4)
  ,((double)-1.3392857142857143E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.6785714285714287E-4)
  ,((double)2.6785714285714287E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)1.1904761904761906E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.5793650793650796E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.142857142857143E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)5.7857142857142857E-2)

  ,((double)3.9682539682539683E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-5)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)2.142857142857143E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)1.1904761904761906E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)2.5793650793650796E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714283E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.7857142857142857E-2)

  ,((double)1.1904761904761906E-3)
  ,((double)3.9682539682539683E-5)
  ,((double)2.5793650793650796E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.7857142857142857E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-5)
  ,((double)3.9682539682539683E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)7.1428571428571429E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)1.25E-3)
  ,((double)1.25E-3)
  ,((double)2.142857142857143E-3)
  ,((double)2.5793650793650796E-4)
  ,((double)3.9682539682539683E-5)
  ,((double)1.1904761904761906E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-8.9285714285714283E-4)
  ,((double)1.7857142857142857E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)1.7857142857142857E-4)
  ,((double)-8.9285714285714283E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)7.1428571428571429E-4)
  ,((double)7.1428571428571429E-4)
  ,((double)0.)
  ,((double)-1.6071428571428571E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)8.0357142857142849E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)1.25E-3)
  ,((double)0.)
  ,((double)8.0357142857142849E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)2.142857142857143E-3)
  ,((double)0.)
  ,((double)6.4285714285714285E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)5.7857142857142857E-2)

};

static const double tria_L2_L2p_L2p[294] = {
  ((double)7.1428571428571426E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-2.7777777777777779E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-8.7301587301587304E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)1.6666666666666666E-2)
  ,((double)0.)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)-5.9523809523809521E-3)
  ,((double)0.)
  ,((double)-5.5555555555555558E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)0.)
  ,((double)-5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)-8.7301587301587304E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-5.9523809523809521E-3)
  ,((double)-5.5555555555555558E-3)
  ,((double)0.)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)1.6666666666666666E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)1.5873015873015873E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-9.6428571428571423E-3)

  ,((double)7.7777777777777779E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)5.7857142857142857E-2)

  ,((double)9.3650793650793651E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)6.9841269841269843E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)2.2222222222222223E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.214285714285714E-2)
  ,((double)4.9206349206349205E-2)
  ,((double)2.2222222222222223E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)6.9841269841269843E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-6.9841269841269843E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)-3.8095238095238099E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-5.7142857142857141E-2)
  ,((double)5.7142857142857141E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)5.7857142857142857E-2)

  ,((double)7.7777777777777779E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)2.3809523809523808E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)4.1269841269841269E-2)
  ,((double)1.1111111111111112E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)6.0714285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)2.8571428571428571E-2)
  ,((double)-2.8571428571428571E-2)
  ,((double)5.7857142857142857E-2)

};

static const double tria_L3_L2p_L2p[490] = {
  ((double)1.1111111111111112E-2)
  ,((double)2.7777777777777779E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)2.142857142857143E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)2.7777777777777779E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)2.1825396825396826E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)2.142857142857143E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)-2.1825396825396826E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-3.1746031746031746E-3)
  ,((double)3.1746031746031746E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)2.142857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)2.142857142857143E-3)
  ,((double)-2.142857142857143E-3)
  ,((double)8.7662337662337668E-4)
  ,((double)3.1746031746031746E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)1.984126984126984E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)3.9682539682539683E-4)
  ,((double)-3.9682539682539683E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)1.984126984126984E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)2.3809523809523812E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)-7.9365079365079365E-4)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.7662337662337668E-4)

  ,((double)3.1746031746031746E-3)
  ,((double)1.984126984126984E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-3.9682539682539683E-4)
  ,((double)3.9682539682539683E-4)
  ,((double)-3.9682539682539683E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.984126984126984E-3)
  ,((double)2.3809523809523812E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)1.3392857142857143E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)5.9523809523809529E-4)
  ,((double)7.1428571428571426E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)1.1904761904761906E-3)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.9682539682539683E-4)
  ,((double)7.9365079365079365E-4)
  ,((double)1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.9682539682539683E-4)
  ,((double)-7.9365079365079365E-4)
  ,((double)-1.1904761904761906E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.3392857142857143E-3)
  ,((double)1.3392857142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.7662337662337668E-4)

  ,((double)2.3214285714285715E-2)
  ,((double)1.7857142857142857E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)5.3571428571428572E-3)
  ,((double)1.7857142857142856E-2)
  ,((double)0.)
  ,((double)-8.9285714285714281E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.7857142857142856E-2)
  ,((double)2.6785714285714284E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.0714285714285714E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)0.)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)-1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)1.0714285714285714E-2)
  ,((double)1.7857142857142857E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)-1.7857142857142857E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)2.5000000000000001E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)2.3214285714285715E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)2.5000000000000001E-2)
  ,((double)2.6785714285714284E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)-2.3214285714285715E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.3214285714285715E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)2.5000000000000001E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)2.3214285714285715E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)0.)
  ,((double)2.5000000000000001E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.3214285714285715E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.3214285714285715E-2)
  ,((double)5.3571428571428572E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-2.3214285714285715E-2)
  ,((double)-5.3571428571428572E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)1.7857142857142856E-2)
  ,((double)-8.9285714285714281E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)0.)
  ,((double)1.7857142857142856E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.0714285714285714E-2)
  ,((double)-1.0714285714285714E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)8.9285714285714281E-3)
  ,((double)1.7857142857142857E-3)
  ,((double)1.0714285714285714E-2)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-8.9285714285714281E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.0714285714285714E-2)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)0.)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.4285714285714285E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)2.3214285714285715E-2)
  ,((double)0.)
  ,((double)1.7857142857142857E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)-1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)0.)
  ,((double)1.7857142857142857E-3)
  ,((double)-1.7857142857142857E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)8.8392857142857145E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.2142857142857142E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)3.2142857142857142E-3)
  ,((double)5.9172077922077923E-3)

  ,((double)0.11785714285714285)
  ,((double)5.3571428571428568E-2)
  ,((double)5.3571428571428568E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)0.10607142857142857)
  ,((double)5.3571428571428568E-2)
  ,((double)3.214285714285714E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)5.3571428571428568E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)6.4285714285714279E-2)
  ,((double)3.214285714285714E-2)
  ,((double)3.214285714285714E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)5.7857142857142857E-2)
  ,((double)-6.4285714285714279E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)-3.214285714285714E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-4.2857142857142858E-2)
  ,((double)4.2857142857142858E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)0.10607142857142857)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)5.7857142857142857E-2)
  ,((double)-5.7857142857142857E-2)
  ,((double)0.10650974025974026)

};

static const double tria_L1dx_L1[9] = {
  ((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

};

static const double tria_L2dx_L2[36] = {
  ((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.10000000000000001)
  ,((double)0.10000000000000001)
  ,((double)-3.3333333333333333E-2)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)0.10000000000000001)
  ,((double)-0.10000000000000001)
  ,((double)0.)
  ,((double)0.)
  ,((double)-0.13333333333333333)
  ,((double)0.13333333333333333)

  ,((double)-3.3333333333333333E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)0.13333333333333333)
  ,((double)0.26666666666666666)
  ,((double)0.26666666666666666)

  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-6.6666666666666666E-2)
  ,((double)-0.13333333333333333)
  ,((double)-0.26666666666666666)
  ,((double)-0.26666666666666666)

};

static const double tria_L3dx_L3[100] = {
  ((double)-3.8095238095238099E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)1.607142857142857E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)6.160714285714286E-2)
  ,((double)6.160714285714286E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-1.607142857142857E-2)

  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-9.6428571428571433E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)-0.14464285714285716)

  ,((double)-3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)0.)
  ,((double)9.6428571428571433E-2)
  ,((double)0.)
  ,((double)-9.6428571428571433E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.14464285714285716)

  ,((double)1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)-2.1428571428571429E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)0.19285714285714287)
  ,((double)2.4107142857142858E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.24107142857142858)

  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)5.8928571428571427E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)0.19285714285714287)
  ,((double)0.19285714285714287)
  ,((double)-7.2321428571428578E-2)
  ,((double)9.6428571428571433E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-5.8928571428571427E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.19285714285714287)
  ,((double)7.2321428571428578E-2)
  ,((double)-9.6428571428571433E-2)

  ,((double)2.6785714285714286E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)2.1428571428571429E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.24107142857142858)

  ,((double)-1.607142857142857E-2)
  ,((double)1.607142857142857E-2)
  ,((double)0.)
  ,((double)0.14464285714285716)
  ,((double)-0.14464285714285716)
  ,((double)-0.24107142857142858)
  ,((double)-9.6428571428571433E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)0.24107142857142858)
  ,((double)0.)

};

static const double tria_L1_L1dx_L1[27] = {
  ((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

};

static const double tria_L2_L2dx_L2[216] = {
  ((double)-3.0952380952380953E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-7.1428571428571426E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)4.3650793650793652E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.8095238095238099E-2)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)3.1746031746031744E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)9.5238095238095247E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-9.5238095238095247E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.0952380952380953E-2)
  ,((double)-3.5714285714285713E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)4.7619047619047623E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.8095238095238099E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.5873015873015873E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.5873015873015873E-3)
  ,((double)-9.5238095238095247E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)

  ,((double)3.5714285714285713E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)4.3650793650793652E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047616E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-4.7619047619047616E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)6.3492063492063489E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.1746031746031744E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-2.5396825396825397E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)

  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063489E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.2698412698412698E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-2.5396825396825397E-2)
  ,((double)-0.10158730158730159)
  ,((double)0.)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)0.15238095238095239)
  ,((double)7.6190476190476197E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-0.15238095238095239)
  ,((double)-7.6190476190476197E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)2.5396825396825397E-2)
  ,((double)0.)
  ,((double)0.10158730158730159)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)7.6190476190476197E-2)
  ,((double)0.15238095238095239)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-7.6190476190476197E-2)
  ,((double)-0.15238095238095239)

};

static const double tria_L3_L3dx_L3[1000] = {
  ((double)-1.7708333333333333E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)2.976190476190476E-3)
  ,((double)1.488095238095238E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)3.3482142857142857E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6785714285714284E-2)
  ,((double)2.2767857142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.0089285714285716E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)2.2098214285714287E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.2767857142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.0133928571428573E-3)
  ,((double)6.4285714285714285E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-6.6964285714285718E-5)
  ,((double)-1.4062499999999999E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.0089285714285712E-3)
  ,((double)8.4375000000000006E-3)
  ,((double)0.)
  ,((double)4.0178571428571428E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-6.6964285714285718E-5)
  ,((double)-2.0089285714285712E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.2098214285714286E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-2.976190476190476E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)1.488095238095238E-3)
  ,((double)1.7708333333333333E-2)
  ,((double)1.488095238095238E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6785714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.2767857142857143E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-6.4285714285714285E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-2.6785714285714284E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)-2.2098214285714287E-2)
  ,((double)1.1049107142857144E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)5.3571428571428572E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.4062499999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-5.3571428571428572E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.4464285714285714E-2)

  ,((double)-1.488095238095238E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.976190476190476E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)1.488095238095238E-3)
  ,((double)2.976190476190476E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-4.2187500000000003E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)0.)
  ,((double)4.2187500000000003E-3)
  ,((double)0.)
  ,((double)-3.8169642857142855E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)3.8169642857142855E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)3.214285714285714E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.214285714285714E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)-6.6964285714285718E-5)
  ,((double)1.4062499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-3.8169642857142855E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-4.0178571428571432E-2)
  ,((double)1.40625E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)3.3482142857142857E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.40625E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)6.8303571428571432E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.0089285714285716E-2)
  ,((double)1.6071428571428571E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)4.2187500000000003E-3)
  ,((double)4.3392857142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-9.0401785714285723E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-2.3504464285714285E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)6.5089285714285711E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.40625E-2)
  ,((double)3.2142857142857142E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-6.8303571428571432E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.40625E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.0312500000000002E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-4.2187500000000003E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-1.6071428571428571E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.8080357142857145E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.3392857142857143E-2)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.3504464285714285E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)-6.5089285714285711E-2)

  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.4866071428571428E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)2.6116071428571429E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)-1.40625E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.2053571428571428E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.700892857142857E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)-2.2098214285714287E-2)
  ,((double)-3.8169642857142855E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.8080357142857145E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-2.3504464285714285E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)0.10848214285714286)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.0089285714285712E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-0.13017857142857142)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-8.6785714285714285E-2)

  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.4866071428571428E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)6.8303571428571432E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)-1.40625E-2)
  ,((double)-3.2142857142857142E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.2053571428571428E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)3.0133928571428573E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.424107142857143E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)0.10848214285714286)
  ,((double)5.424107142857143E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-0.10848214285714286)
  ,((double)-5.424107142857143E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.424107142857143E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.10848214285714286)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.0044642857142858E-2)
  ,((double)-6.8303571428571432E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)3.2142857142857142E-3)
  ,((double)1.40625E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.4866071428571428E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.1049107142857144E-2)
  ,((double)-3.0133928571428573E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)6.4285714285714285E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.2098214285714286E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-5.424107142857143E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)0.10848214285714286)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-0.10848214285714286)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-5.424107142857143E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)0.10848214285714286)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)1.0044642857142858E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.40625E-2)
  ,((double)-4.0178571428571432E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4866071428571428E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.2098214285714287E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-1.4464285714285714E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-4.700892857142857E-2)
  ,((double)0.)
  ,((double)-2.0089285714285712E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)2.3504464285714285E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-0.10848214285714286)
  ,((double)-6.5089285714285711E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)0.)
  ,((double)3.2544642857142855E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.13017857142857142)
  ,((double)8.6785714285714285E-2)

  ,((double)-2.6785714285714286E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.8928571428571428E-2)
  ,((double)1.4732142857142858E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.8928571428571428E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)0.13017857142857142)
  ,((double)8.4375000000000006E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.17357142857142857)
  ,((double)4.8214285714285711E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.424107142857143E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-4.3392857142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)-0.17357142857142857)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)0.)
  ,((double)6.5089285714285711E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)-8.6785714285714285E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)8.6785714285714285E-2)
  ,((double)0.)

};

static const double tria_L1dy_L1[9] = {
  ((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)-0.16666666666666666)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)
  ,((double)0.16666666666666666)

};

static const double tria_L2dy_L2[36] = {
  ((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.10000000000000001)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)-3.3333333333333333E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.10000000000000001)
  ,((double)0.10000000000000001)

  ,((double)3.3333333333333333E-2)
  ,((double)-6.6666666666666666E-2)
  ,((double)3.3333333333333333E-2)
  ,((double)-0.26666666666666666)
  ,((double)-0.26666666666666666)
  ,((double)-0.13333333333333333)

  ,((double)-3.3333333333333333E-2)
  ,((double)6.6666666666666666E-2)
  ,((double)-3.3333333333333333E-2)
  ,((double)0.26666666666666666)
  ,((double)0.26666666666666666)
  ,((double)0.13333333333333333)

  ,((double)0.10000000000000001)
  ,((double)0.)
  ,((double)-0.10000000000000001)
  ,((double)0.13333333333333333)
  ,((double)-0.13333333333333333)
  ,((double)0.)

};

static const double tria_L3dy_L3[100] = {
  ((double)-3.8095238095238099E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-1.1309523809523809E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)3.4821428571428573E-2)
  ,((double)-6.160714285714286E-2)
  ,((double)1.607142857142857E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)

  ,((double)1.1309523809523809E-2)
  ,((double)1.1309523809523809E-2)
  ,((double)3.8095238095238099E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)6.160714285714286E-2)
  ,((double)6.160714285714286E-2)
  ,((double)-3.4821428571428573E-2)
  ,((double)-1.607142857142857E-2)

  ,((double)2.6785714285714286E-3)
  ,((double)2.1428571428571429E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)-0.19285714285714287)
  ,((double)-2.4107142857142858E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)-0.24107142857142858)

  ,((double)-1.3392857142857142E-2)
  ,((double)-5.8928571428571427E-2)
  ,((double)-1.3392857142857142E-2)
  ,((double)7.2321428571428578E-2)
  ,((double)-0.19285714285714287)
  ,((double)-0.19285714285714287)
  ,((double)7.2321428571428578E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-9.6428571428571433E-2)

  ,((double)1.3392857142857142E-2)
  ,((double)5.8928571428571427E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)0.19285714285714287)
  ,((double)0.19285714285714287)
  ,((double)-7.2321428571428578E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)9.6428571428571433E-2)

  ,((double)1.3392857142857142E-2)
  ,((double)-2.1428571428571429E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)-4.8214285714285716E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)0.19285714285714287)
  ,((double)9.6428571428571433E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)0.24107142857142858)

  ,((double)-3.4821428571428573E-2)
  ,((double)0.)
  ,((double)-6.160714285714286E-2)
  ,((double)-4.8214285714285716E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)0.)
  ,((double)9.6428571428571433E-2)
  ,((double)0.14464285714285716)

  ,((double)6.160714285714286E-2)
  ,((double)0.)
  ,((double)3.4821428571428573E-2)
  ,((double)9.6428571428571433E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)4.8214285714285716E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)0.)
  ,((double)-0.14464285714285716)

  ,((double)-1.607142857142857E-2)
  ,((double)0.)
  ,((double)1.607142857142857E-2)
  ,((double)0.24107142857142858)
  ,((double)9.6428571428571433E-2)
  ,((double)-9.6428571428571433E-2)
  ,((double)-0.24107142857142858)
  ,((double)-0.14464285714285716)
  ,((double)0.14464285714285716)
  ,((double)0.)

};

static const double tria_L1_L1dy_L1[27] = {
  ((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)
  ,((double)4.1666666666666664E-2)

  ,((double)-4.1666666666666664E-2)
  ,((double)-4.1666666666666664E-2)
  ,((double)-8.3333333333333329E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.1666666666666664E-2)
  ,((double)4.1666666666666664E-2)
  ,((double)8.3333333333333329E-2)

};

static const double tria_L2_L2dy_L2[216] = {
  ((double)-3.0952380952380953E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)3.5714285714285713E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-4.7619047619047623E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.1428571428571426E-3)
  ,((double)4.3650793650793652E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-9.5238095238095247E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.5873015873015873E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)3.5714285714285713E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)4.3650793650793652E-3)
  ,((double)-7.1428571428571426E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-4.7619047619047616E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047616E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-7.9365079365079361E-3)
  ,((double)0.)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)0.)

  ,((double)3.5714285714285713E-3)
  ,((double)-4.3650793650793652E-3)
  ,((double)7.1428571428571426E-3)
  ,((double)7.9365079365079361E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.5714285714285713E-3)
  ,((double)-3.5714285714285713E-3)
  ,((double)3.0952380952380953E-2)
  ,((double)4.7619047619047623E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.5873015873015873E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-9.5238095238095247E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)1.5873015873015873E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)9.5238095238095247E-3)
  ,((double)-1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)7.9365079365079361E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)-3.1746031746031744E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-3.1746031746031744E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)4.7619047619047623E-3)
  ,((double)-3.8095238095238099E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-0.15238095238095239)
  ,((double)-7.6190476190476197E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)0.15238095238095239)
  ,((double)7.6190476190476197E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)0.10158730158730159)
  ,((double)0.)
  ,((double)2.5396825396825397E-2)

  ,((double)-4.7619047619047623E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)3.8095238095238099E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-7.9365079365079361E-3)
  ,((double)-1.2698412698412698E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063489E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)-7.6190476190476197E-2)
  ,((double)-0.15238095238095239)
  ,((double)-5.0793650793650794E-2)
  ,((double)-1.9047619047619049E-2)
  ,((double)1.9047619047619049E-2)
  ,((double)-1.2698412698412698E-2)
  ,((double)7.6190476190476197E-2)
  ,((double)0.15238095238095239)
  ,((double)5.0793650793650794E-2)
  ,((double)1.2698412698412698E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)-2.5396825396825397E-2)
  ,((double)0.)
  ,((double)-0.10158730158730159)
  ,((double)-2.5396825396825397E-2)

  ,((double)-1.9047619047619049E-2)
  ,((double)7.9365079365079361E-3)
  ,((double)1.2698412698412698E-2)
  ,((double)-3.1746031746031744E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063489E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.2698412698412698E-2)
  ,((double)-7.9365079365079361E-3)
  ,((double)1.9047619047619049E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)3.1746031746031744E-2)
  ,((double)6.3492063492063489E-2)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)6.3492063492063492E-3)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-5.0793650793650794E-2)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)-6.3492063492063492E-3)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)5.0793650793650794E-2)
  ,((double)3.1746031746031744E-2)
  ,((double)0.)
  ,((double)-3.1746031746031744E-2)
  ,((double)2.5396825396825397E-2)
  ,((double)-2.5396825396825397E-2)
  ,((double)0.)

};

static const double tria_L3_L3dy_L3[1000] = {
  ((double)-1.7708333333333333E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.488095238095238E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6785714285714286E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.976190476190476E-3)
  ,((double)8.8541666666666662E-4)
  ,((double)1.488095238095238E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)1.4732142857142858E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)6.6964285714285718E-5)
  ,((double)-4.0178571428571425E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)5.0223214285714289E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-5.0223214285714289E-3)
  ,((double)7.633928571428571E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)-2.0089285714285712E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.4062499999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)8.4375000000000006E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-1.1049107142857144E-2)
  ,((double)6.4285714285714285E-3)
  ,((double)3.0133928571428573E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.6071428571428571E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.6785714285714284E-2)
  ,((double)2.4107142857142856E-3)
  ,((double)2.2767857142857143E-3)
  ,((double)2.2098214285714287E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-7.0312500000000002E-3)
  ,((double)2.0089285714285716E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)5.3571428571428572E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)6.0267857142857146E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-4.8214285714285711E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-1.488095238095238E-3)
  ,((double)-2.976190476190476E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)8.8541666666666662E-4)
  ,((double)2.976190476190476E-3)
  ,((double)1.488095238095238E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)1.2053571428571429E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-3.8169642857142855E-3)
  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.214285714285714E-2)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)-2.4107142857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)3.214285714285714E-2)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)2.4107142857142858E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.6964285714285718E-5)
  ,((double)6.0267857142857146E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)0.)
  ,((double)-1.8080357142857143E-3)
  ,((double)0.)
  ,((double)-2.4107142857142856E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)0.)
  ,((double)4.2187500000000003E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)0.)
  ,((double)1.8080357142857143E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-4.2187500000000003E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)0.)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)0.)

  ,((double)-1.488095238095238E-3)
  ,((double)-8.8541666666666662E-4)
  ,((double)-2.976190476190476E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.488095238095238E-3)
  ,((double)1.488095238095238E-3)
  ,((double)1.7708333333333333E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)1.3392857142857142E-2)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)6.6964285714285718E-5)
  ,((double)1.4062499999999999E-3)
  ,((double)0.)
  ,((double)-7.4330357142857141E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.0089285714285712E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-8.4375000000000006E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-2.3437499999999999E-3)
  ,((double)0.)
  ,((double)2.8124999999999999E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-7.633928571428571E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)-4.8214285714285711E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-2.8124999999999999E-3)
  ,((double)7.4330357142857141E-3)
  ,((double)7.633928571428571E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.6964285714285718E-5)
  ,((double)6.6964285714285718E-5)
  ,((double)5.3571428571428572E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.2098214285714286E-3)
  ,((double)-5.0223214285714289E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.2767857142857143E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6785714285714284E-2)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-2.2098214285714287E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)7.0312500000000002E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.2767857142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)7.0312500000000002E-3)
  ,((double)1.6071428571428571E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-1.4732142857142858E-3)
  ,((double)-5.3571428571428572E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)1.4464285714285714E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.4107142857142856E-3)
  ,((double)-4.0178571428571432E-2)
  ,((double)1.40625E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)-2.0089285714285716E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.3437499999999999E-3)
  ,((double)-1.4062499999999999E-3)
  ,((double)0.)
  ,((double)1.4866071428571428E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)1.2053571428571428E-3)
  ,((double)-4.0178571428571425E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-7.4330357142857141E-3)
  ,((double)-0.10848214285714286)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)2.3504464285714285E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-2.0089285714285712E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.2053571428571428E-3)
  ,((double)-4.700892857142857E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)0.)
  ,((double)2.2098214285714287E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)7.2321428571428578E-2)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.8080357142857145E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)7.2321428571428571E-3)
  ,((double)0.13017857142857142)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)8.6785714285714285E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.40625E-2)
  ,((double)3.2142857142857142E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-2.6116071428571429E-3)
  ,((double)-6.8303571428571432E-3)
  ,((double)1.0044642857142858E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)2.3437499999999999E-3)
  ,((double)0.)
  ,((double)-1.4062499999999999E-3)
  ,((double)1.4866071428571428E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)5.0223214285714289E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)-1.4062499999999999E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-0.10848214285714286)
  ,((double)-5.424107142857143E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)7.633928571428571E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)0.10848214285714286)
  ,((double)5.424107142857143E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)5.424107142857143E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.2053571428571429E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)0.)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.4285714285714285E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)-3.6160714285714286E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.1049107142857144E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-3.0133928571428573E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.10848214285714286)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)

  ,((double)0.)
  ,((double)-2.3437499999999999E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)-1.4866071428571428E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.3482142857142857E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-3.2142857142857142E-3)
  ,((double)-1.40625E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)6.8303571428571432E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)2.0089285714285714E-4)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)5.424107142857143E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-7.4330357142857141E-3)
  ,((double)-2.4107142857142858E-2)
  ,((double)-7.633928571428571E-3)
  ,((double)2.7120535714285715E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-0.10848214285714286)
  ,((double)2.7120535714285715E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-5.424107142857143E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)2.4107142857142858E-2)
  ,((double)7.633928571428571E-3)
  ,((double)-2.7120535714285715E-2)
  ,((double)5.424107142857143E-2)
  ,((double)0.10848214285714286)
  ,((double)-2.7120535714285715E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)5.424107142857143E-2)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.)
  ,((double)5.424107142857143E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)3.0133928571428573E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)1.1049107142857144E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)0.)
  ,((double)-1.2053571428571428E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-6.4285714285714285E-3)
  ,((double)0.)
  ,((double)-5.4241071428571428E-3)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)9.0401785714285723E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)0.)
  ,((double)4.8214285714285711E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-0.10848214285714286)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-2.1696428571428571E-2)

  ,((double)0.)
  ,((double)1.4062499999999999E-3)
  ,((double)-2.3437499999999999E-3)
  ,((double)2.2098214285714286E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)1.4062499999999999E-3)
  ,((double)-1.4866071428571428E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)2.8124999999999999E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.4107142857142856E-3)
  ,((double)2.6116071428571429E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)-2.2098214285714286E-3)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.40625E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)-1.0044642857142858E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-2.0089285714285714E-4)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)2.8124999999999999E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)5.0223214285714289E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)2.7120535714285715E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)-5.0223214285714289E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)-2.7120535714285715E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)7.4330357142857141E-3)
  ,((double)6.0267857142857146E-3)
  ,((double)4.0178571428571425E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)-5.4241071428571428E-3)
  ,((double)0.10848214285714286)
  ,((double)3.6160714285714289E-2)
  ,((double)-2.3504464285714285E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-2.2098214285714287E-2)
  ,((double)5.4241071428571428E-3)
  ,((double)0.)
  ,((double)1.4464285714285714E-2)
  ,((double)-7.2321428571428578E-2)
  ,((double)-1.8080357142857145E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.2053571428571428E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)1.1049107142857144E-2)
  ,((double)-5.4241071428571428E-3)
  ,((double)5.4241071428571428E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)4.700892857142857E-2)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)0.)
  ,((double)-7.2321428571428571E-3)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-3.2544642857142855E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-8.6785714285714285E-2)

  ,((double)6.0267857142857146E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.3482142857142857E-4)
  ,((double)1.0044642857142858E-2)
  ,((double)-6.8303571428571432E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)3.2142857142857142E-3)
  ,((double)1.40625E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)2.6116071428571429E-3)
  ,((double)2.4107142857142856E-3)
  ,((double)1.3392857142857142E-2)
  ,((double)-2.8124999999999999E-3)
  ,((double)-2.2098214285714286E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)4.0178571428571432E-2)
  ,((double)-1.40625E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)-2.0089285714285714E-4)
  ,((double)4.0178571428571425E-3)
  ,((double)2.3504464285714285E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)4.0178571428571428E-4)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)3.8169642857142855E-3)
  ,((double)3.8169642857142855E-3)
  ,((double)2.0089285714285712E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)3.6160714285714289E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-1.6071428571428571E-3)
  ,((double)0.)
  ,((double)-2.0089285714285716E-2)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857145E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)-4.2187500000000003E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)-1.4464285714285714E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)9.0401785714285723E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-4.3392857142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-4.8214285714285711E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)-3.2544642857142855E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-6.5089285714285711E-2)

  ,((double)-1.3392857142857142E-2)
  ,((double)-2.4107142857142856E-3)
  ,((double)-2.6116071428571429E-3)
  ,((double)-2.0089285714285716E-2)
  ,((double)1.0044642857142858E-2)
  ,((double)2.2098214285714286E-3)
  ,((double)2.8124999999999999E-3)
  ,((double)1.40625E-2)
  ,((double)-4.0178571428571432E-2)
  ,((double)-1.2053571428571429E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)3.3482142857142857E-4)
  ,((double)-1.8080357142857143E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)2.6116071428571429E-3)
  ,((double)6.8303571428571432E-3)
  ,((double)-1.0044642857142858E-2)
  ,((double)-1.40625E-2)
  ,((double)-3.2142857142857142E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)-2.0089285714285712E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.8169642857142855E-3)
  ,((double)-3.6160714285714289E-2)
  ,((double)1.2656250000000001E-2)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)1.2656250000000001E-2)
  ,((double)-3.6160714285714289E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)1.0044642857142856E-3)
  ,((double)1.0044642857142856E-3)
  ,((double)4.0178571428571428E-4)
  ,((double)1.8080357142857143E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)7.2321428571428571E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)-1.0044642857142856E-3)
  ,((double)-1.0044642857142856E-3)
  ,((double)-4.0178571428571428E-4)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-1.0848214285714286E-2)
  ,((double)-4.0178571428571425E-3)
  ,((double)2.0089285714285714E-4)
  ,((double)-1.0044642857142856E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-1.8080357142857143E-3)
  ,((double)-2.3504464285714285E-2)
  ,((double)-1.2656250000000001E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-7.0312500000000002E-3)
  ,((double)4.2187500000000003E-3)
  ,((double)7.0312500000000002E-3)
  ,((double)1.8080357142857143E-3)
  ,((double)-9.0401785714285723E-3)
  ,((double)-1.2656250000000001E-2)
  ,((double)1.4464285714285714E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)2.0089285714285716E-2)
  ,((double)0.)
  ,((double)1.6071428571428571E-3)
  ,((double)1.8080357142857145E-2)
  ,((double)-3.6160714285714286E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-7.2321428571428571E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)3.2544642857142855E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)0.)
  ,((double)2.1696428571428571E-2)
  ,((double)6.5089285714285711E-2)

  ,((double)-2.6785714285714286E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.4732142857142858E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)-1.2053571428571428E-3)
  ,((double)9.6428571428571423E-3)
  ,((double)-1.2053571428571429E-2)
  ,((double)2.8928571428571428E-2)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)0.)
  ,((double)1.4732142857142858E-3)
  ,((double)1.4732142857142858E-3)
  ,((double)2.6785714285714286E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)1.2053571428571428E-3)
  ,((double)-9.6428571428571423E-3)
  ,((double)1.2053571428571429E-2)
  ,((double)1.2053571428571429E-2)
  ,((double)-9.6428571428571423E-3)
  ,((double)-2.8928571428571428E-2)
  ,((double)6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-8.4375000000000006E-3)
  ,((double)-6.5089285714285711E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)-0.17357142857142857)
  ,((double)-4.8214285714285711E-3)
  ,((double)0.)
  ,((double)-4.8214285714285711E-3)
  ,((double)2.1696428571428571E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)-5.424107142857143E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)-4.3392857142857143E-2)
  ,((double)4.8214285714285711E-3)
  ,((double)0.)
  ,((double)4.8214285714285711E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)5.424107142857143E-2)
  ,((double)5.424107142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)8.4375000000000006E-3)
  ,((double)0.)
  ,((double)-6.0267857142857146E-3)
  ,((double)-2.1696428571428571E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)4.3392857142857143E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)0.17357142857142857)
  ,((double)-4.8214285714285711E-3)
  ,((double)3.6160714285714286E-3)
  ,((double)-6.0267857142857146E-3)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-1.0848214285714286E-2)
  ,((double)1.0848214285714286E-2)
  ,((double)3.2544642857142855E-2)
  ,((double)0.13017857142857142)
  ,((double)6.0267857142857146E-3)
  ,((double)-3.6160714285714286E-3)
  ,((double)4.8214285714285711E-3)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)1.0848214285714286E-2)
  ,((double)0.)
  ,((double)-3.2544642857142855E-2)
  ,((double)-1.0848214285714286E-2)
  ,((double)-0.13017857142857142)
  ,((double)-1.4464285714285714E-2)
  ,((double)0.)
  ,((double)1.4464285714285714E-2)
  ,((double)8.6785714285714285E-2)
  ,((double)2.1696428571428571E-2)
  ,((double)-2.1696428571428571E-2)
  ,((double)-8.6785714285714285E-2)
  ,((double)-6.5089285714285711E-2)
  ,((double)6.5089285714285711E-2)
  ,((double)0.)

};

#define dim 		((I)2)
#define nfaceinelm 	((I)3)






struct DG_DATA
{
  
  //
  //
  //
  struct matrix_handle m_flux{};

  //
  // submatrix.
  //
  struct matrix_handle m_flux_part_corr{};
  
  //
  // submatrix.
  //
  struct matrix_handle m_flux_part_sol{};

  
  struct matrix_handle uvw_ldofs;

  struct vector_handle m_local_matrices;
  struct matrix_handle m_local_matrix_part_corr;  
  struct matrix_handle m_local_matrix_part_sol;

  struct vector_handle vec_neicorr;
  struct vector_handle vec_neisol;

  struct matrix_handle mat_belm;
  struct vector_handle vec_nrmelm[nfaceinelm];
  struct vector_handle vec_uface[nfaceinelm];

  struct vector_handle m_local_rhs;
  struct vector_handle hsol;

  pR m_udofs;
  pR m_local_matrices_memory;
  pR m_flux_memory{};
  R belm[dim*dim];

  virtual ~DG_DATA()
  {
    if (m_flux_memory)
      {
	free(m_flux_memory);
	m_flux_memory = nullptr;
      }    
    if (m_udofs)
      {
	free(m_udofs);
	m_udofs = nullptr;
      }    
  };

  R nrmelm[dim * nfaceinelm];

  R lcrhs[128];
  DG_DATA(I teta_n_,
	  I trial_n_,
	  I test_n_,
	  I teta_u_n_)
  {
    this->m_flux_memory = (pR)malloc(sizeof(R)*(trial_n_*test_n_+teta_n_*test_n_));
    matrix_handle_def(&this->m_flux_part_corr,test_n_,trial_n_,this->m_flux_memory,test_n_);
    matrix_handle_def(&this->m_flux_part_sol,test_n_,teta_n_,&this->m_flux_memory[trial_n_*test_n_],test_n_);    

    this->m_udofs = (pR)malloc(sizeof(R)*(teta_u_n_*dim));
    matrix_handle_def(&this->uvw_ldofs,teta_u_n_,dim,this->m_udofs,teta_u_n_);
    
    matrix_handle_def(&this->mat_belm,
		      2,2,belm,2);

    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&this->vec_nrmelm[localFaceIndex],dim,&nrmelm[dim * localFaceIndex],1);
      }

    I len  = trial_n_*test_n_+teta_n_*test_n_;
    m_local_matrices_memory = (pR)malloc(sizeof(R)*len);
    vector_handle_def(&this->m_local_matrices,len,m_local_matrices_memory,1);
    matrix_handle_def(&this->m_local_matrix_part_corr,test_n_,trial_n_,  m_local_matrices_memory, test_n_);
    matrix_handle_def(&this->m_local_matrix_part_sol, test_n_,teta_n_,   m_local_matrices_memory + trial_n_ * test_n_, test_n_);

    vector_handle_def(&this->m_local_rhs,trial_n_,lcrhs,1);
  };

};


struct DG_HANDLE
{

  struct matrix_handle m_EVALU;
  struct matrix_handle m_UVWDOFS;
  struct matrix_handle m_BMAT;
  struct vector_handle m_BRHS;

  struct matrix_handle m_mat_tmpbrhs_uelm;
  struct matrix_handle m_brhs_uvw;
  struct matrix_handle m_mat_tmpbrhs_uface[nfaceinelm];
  
};


struct DG
{  
  
public:  
  typedef enum enum_DG
    {
    
      IA_lc=0,
      IA_lc_elm=1,
      IA_lc_face=2,
      IA_lc_elm_a=3,
      IA_lc_u=4,
      IA_lc_elm_u0=5,
      IA_lc_elm_u1=6,
      I_npoints=11,
      I_npoints_element=12,
      I_npoints_boundary=13,
      RA_lc=14,
      I_lc_len=16,
      RA_bmat=17,
      RA_bmatx=18,
      I_bmat_n=19,
      I_bmat_m=20,
      I_bmat_len=21,
      ERR=22,
      RA_EVAL_TETA_U=23,

      I_QFACE_N=28,
      I_QELM_N=29,

      I_TRIAL_FAMILY=30,
      I_TRIAL_DEGREE=31,
      I_TRIAL_NBASIS=32,

    
      I_TEST_FAMILY=33,
      I_TEST_DEGREE=34,
      I_TEST_NBASIS=35,
    
      I_TETA_FAMILY=36,
      I_TETA_DEGREE=37,
      I_TETA_NBASIS=38,
    
      I_TETA_U_FAMILY=39,
      I_TETA_U_DEGREE=40,
      I_TETA_U_NBASIS=41,
    
      I_TETA_A_FAMILY=42,
      I_TETA_A_DEGREE=43,
      I_TETA_A_NBASIS=44,

      I_n=64

    } info_t;

  
public: static DG_DATA * create_data(cst_mkS 		s_teta_a_,
				     cst_mkS 		s_teta_u_,
				     cst_mkS 		s_teta_,
				     cst_mkS 		s_test_,
				     cst_mkS 		s_trial_)
  {
    return new DG_DATA(mkS_n(s_teta_),mkS_n(s_trial_),mkS_n(s_test_),mkS_n(s_teta_u_));
  };
  
public: static void dgadvection_init_with_quadrature(cst_mkS 		s_teta_a_,
						     cst_mkS 		s_teta_u_,
						     cst_mkS 		s_teta_,
						     cst_mkS 		s_test_,
						     cst_mkS 		s_trial_,
						     
						     cst_pI 		qelm_n_,
						     cst_pR 		qelm_p_,
						     cst_pR 		qelm_w_,
						     
						     cst_pI 		qface_n_,
						     cst_pR 		qface_p_,
						     cst_pR 		qface_w_,
						     
						     cst_pI           	iinfo_n_,
						     pI		iinfo_,
						     cst_pI		rinfo_n_,
						     pR 		rinfo_,
						     cst_pI 		rwork_n_,
						     pR		rwork_)
  {
  
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_n_[0] < DG::I_n)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }
  
    iinfo_[DG::I_TETA_A_NBASIS]  = (s_teta_a_) ? mkS_n(s_teta_a_):((I)0);
    iinfo_[DG::I_TETA_U_NBASIS]	 = mkS_n(s_teta_u_);
    iinfo_[DG::I_TETA_NBASIS]	 = mkS_n(s_teta_);
    iinfo_[DG::I_TRIAL_NBASIS]	 = mkS_n(s_trial_);
    iinfo_[DG::I_TEST_NBASIS]	 = mkS_n(s_test_);
    iinfo_[DG::I_QELM_N] 	 = qelm_n_[0];  

    if (ferr_[0]) 
      {
	return;
      }
  
    if (iinfo_[DG::I_TRIAL_NBASIS]!=iinfo_[DG::I_TEST_NBASIS])   
      { 
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_TEST_NBASIS]!=IINFO[DG::I_TRIAL_NBASIS]\n");    
      }

    mkS_st s_teta_a;
    mkS_st s_teta_u;
    mkS_st s_teta;
    mkS_st s_test;
    mkS_st s_trial;

    mkS_copy(&s_teta_a,s_teta_a_);
    mkS_copy(&s_teta_u,s_teta_u_);
    mkS_copy(&s_teta,s_teta_);
    mkS_copy(&s_test,s_test_);
    mkS_copy(&s_trial,s_trial_);

    const I s_teta_a_degree = mkS_k(s_teta_a_);
    const I s_teta_u_degree = mkS_k(s_teta_u_);

    const I trial_n	= iinfo_[DG::I_TRIAL_NBASIS];
    const I test_n	= iinfo_[DG::I_TEST_NBASIS];
    const I teta_n	= iinfo_[DG::I_TETA_NBASIS];
    const I teta_a_n 	= iinfo_[DG::I_TETA_A_NBASIS];
    const I teta_u_n 	= iinfo_[DG::I_TETA_U_NBASIS];

    fprintf(stdout,"  Discontinuous Galerkin Definition\n");
    fprintf(stdout,"  trial_n   " ifmt "\n",trial_n);
    fprintf(stdout,"  test_n    " ifmt "\n",test_n);
    fprintf(stdout,"  teta_n    " ifmt "\n",teta_n);
    fprintf(stdout,"  teta_a_n  " ifmt "\n",teta_a_n);
    fprintf(stdout,"  teta_u_n  " ifmt "\n",teta_u_n);
  
    const I trial_nXtrial_n 	= trial_n * trial_n;

    
    const I npoints_element	= teta_a_n + dim * teta_u_n;
    const I npoints_boundary	= nfaceinelm * qface_n_[0];
    const I npoints 	 	= npoints_element + npoints_boundary;  
    const I bmat_n 		= trial_n * test_n + teta_n * test_n;
    const I bmat_m 		= npoints;
    
    const I bmatx_size 		= bmat_n*(nfaceinelm * qface_n_[0] * nfaceinelm);
    const I evalu_size 		= teta_u_n * (dim*teta_u_n + nfaceinelm * qface_n_[0]);  
    if (rinfo_n_[0] < 2 * npoints + bmat_n*npoints + bmatx_size + evalu_size)
      {
	fprintf(stderr,"*** DGI too small rinfo_ array (" ifmt "<" ifmt ")\n",
		rinfo_n_[0],
		2*npoints+bmat_n*npoints+bmatx_size+evalu_size);
	ferr_[0] = DGERR_MEMORY;
	return;
      }
  
    /* on recopie la quadrature pour les faces */  
    /* on recopie les interpolants */
    
    iinfo_[DG::I_npoints]     		= npoints;
    iinfo_[DG::I_npoints_element]  	= npoints_element;
    iinfo_[DG::I_npoints_boundary] 	= npoints_boundary;
    
    iinfo_[DG::I_bmat_n] 		= bmat_n;
    iinfo_[DG::I_bmat_m] 		= iinfo_[DG::I_npoints];
    iinfo_[DG::I_bmat_len] 		= iinfo_[DG::I_bmat_n]*iinfo_[DG::I_npoints];

    iinfo_[DG::I_lc_len] 		= 2*iinfo_[DG::I_npoints];
  
    /* matrice de flux */  
    /* les points d evaluation en coordonnees locales */
    iinfo_[DG::IA_lc]        	= 0;
    iinfo_[DG::IA_lc_elm]    	= 0;
    iinfo_[DG::IA_lc_elm_a]  	= 0;
    iinfo_[DG::IA_lc_elm_u0] 	= iinfo_[DG::I_TETA_A_NBASIS];
    iinfo_[DG::IA_lc_elm_u1] 	= iinfo_[DG::I_TETA_A_NBASIS]+iinfo_[DG::I_TETA_U_NBASIS];
    iinfo_[DG::IA_lc_face]   	= npoints_element;
    iinfo_[DG::IA_lc_u]      	= iinfo_[DG::I_TETA_A_NBASIS];
  
    iinfo_[DG::RA_bmat] 		= 2*npoints;
    iinfo_[DG::RA_bmatx] 		= iinfo_[DG::RA_bmat]  + iinfo_[DG::I_bmat_len];
    iinfo_[DG::RA_EVAL_TETA_U]    = iinfo_[DG::RA_bmatx] + bmatx_size;

    

    const I npoints_velocity_involved 	= iinfo_[DG::I_npoints]-iinfo_[DG::I_TETA_A_NBASIS];  
    
    pR bmat_x = &rinfo_[iinfo_[DG::RA_bmat]];  
    const I bmat_ld = bmat_n;
    pR teta_a_trial_test 	= &bmat_x[0];

    pR teta_u_nabla_trial_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_trial_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld];
      }
    
    //    pR teta_u_dxtrial_test 	= &bmat_x[teta_a_n*bmat_ld];
    //    pR teta_u_dytrial_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_ld];
    pR teta_u_nabla_teta_test[dim];
    for (I idim =0;idim<dim;++idim)
      {
	teta_u_nabla_teta_test[idim] = &bmat_x[(teta_a_n+idim*teta_u_n)*bmat_ld+trial_nXtrial_n];
      }

    pR teta_a_teta_test  	= &bmat_x[trial_nXtrial_n];
    //    pR teta_u_dxteta_test 	= &bmat_x[teta_a_n*bmat_n+trial_nXtrial_n];
    //    pR teta_u_dyteta_test 	= &bmat_x[(teta_a_n+teta_u_n)*bmat_n+trial_nXtrial_n];
    
    if (teta_a_n>0)
      {
	
	//
	// Form the matrix.
	//
	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_trial),
		mkS_b(&s_test),
		teta_a_trial_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);

#if 0
	I i;for (i=0;i<36;++i) { printf("hh %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }      
#endif

	mkS_kji(mkS_b(&s_teta_a),
		mkS_b(&s_teta),
		mkS_b(&s_test),
		teta_a_teta_test,
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
#if 0
	for (i=0;i<36;++i) { printf("gg %e \n",nsFABS(tria_L2_L2[i]-bmat_x[i]));  }
#endif
      
      }

  
    if (teta_u_n>0)
      {
	
#if DG_CONSERVATIVE_WEAK_FORM
	mkS_kji(mkS_b(&s_teta_u),
		mkS_b(&s_trial),
		mkS_dx(&s_test),
		&bmat_x[teta_a_n*bmat_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
	mkS_kji(mkS_b(&s_teta_u),
		mkS_b(&s_trial),
		mkS_dy(&s_test),
		&bmat_x[(teta_a_n+teta_u_n)*bmat_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
	mkS_kji(mkS_b(&s_teta_u),
		mkS_b(&s_teta),
		mkS_dx(&s_test),
		&bmat_x[teta_a_n*bmat_n+trial_nXtrial_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
	mkS_kji(mkS_b(&s_teta_u),
		mkS_b(&s_teta),
		mkS_dy(&s_test),
		&bmat_x[(teta_a_n+teta_u_n)*bmat_n+trial_nXtrial_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
      
#else

	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_trial,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_trial_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }

	for (I idim =0;idim < dim;++idim)
	  {
	    mkS_kji(mkS_b(&s_teta_u),
		    mkS_derivative(&s_teta,idim),
		    mkS_b(&s_test),
		    teta_u_nabla_teta_test[idim],
		    &bmat_ld,
		    qelm_n_,
		    qelm_w_,
		    qelm_p_,
		    qelm_n_,
		    rwork_n_,
		    rwork_,
		    ferr_);
	  }

	
#if 0
	I i,j;for (j=0;j<6;++j) for (i=0;i<36;++i) { printf("gg %e \n",nsFABS(tria_L2_L2dx_L2[j*36+i]-bmat_x[teta_a_n*bmat_n+j*bmat_n+i]));  }
	mk_err("salut !\n");
#endif
#if 0
	mkS_kji(mkS_b(&s_teta_u),
		mkS_dx(&s_teta),
		mkS_b(&s_test),
		teta_u_dxteta_test,//     &bmat_x[teta_a_n*bmat_n+trial_nXtrial_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);

      
	mkS_kji(mkS_b(&s_teta_u),
		mkS_dy(&s_teta),
		mkS_b(&s_test),
		teta_u_dyteta_test,// &bmat_x[(teta_a_n+teta_u_n)*bmat_n+trial_nXtrial_n],
		&bmat_n,
		qelm_n_,
		qelm_w_,
		qelm_p_,
		qelm_n_,
		rwork_n_,
		rwork_,
		ferr_);
#endif 
#endif
      }  

    /* -------------------------------------------------------------------------------------- */
    /* LOCAL COORDINATES -------------------------------------------------------------------  */
    /* -------------------------------------------------------------------------------------- */            
    { pR lc = &rinfo_[iinfo_[DG::IA_lc]];

      if (iinfo_[DG::I_TETA_A_NBASIS]>0)
	{
	  mkS_lagrange_localspl_tria(&s_teta_a_degree,
				     &lc[iinfo_[DG::IA_lc_elm_a]],
				     &npoints);
	}
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u0]],
					 &npoints);
    
      mkS_lagrange_localspl_tria	(&s_teta_u_degree,
					 &lc[iinfo_[DG::IA_lc_elm_u1]],
					 &npoints);
    
      mkS_bmapping		(qface_n_[0],
				 &lc[iinfo_[DG::IA_lc_face]],
				 &npoints,
				 qface_p_);
    
      /* -------------------------------------------------------------------------------------- */
      /* CALCUL DE L EVALUATION DE TETA_U POUR TOUS LES POINTS D INTEGRATION ELM + 3*NFACES */
      /* -------------------------------------------------------------------------------------- */            
      mkS_basis(mkS_b(&s_teta_u),
		&npoints_velocity_involved,
		&rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		&iinfo_[DG::I_TETA_U_NBASIS],
		&rinfo_[iinfo_[DG::IA_lc_u]],
		&iinfo_[DG::I_npoints],// DG::I_npoints			   			   
		rwork_,
		rwork_n_,
		ferr_);  }
  
#if 0
    {
      I i,j;pR lc = &rinfo_[iinfo_[DG::IA_lc]];
      for (i=0;i<npoints_velocity_involved;++i)
	{
	  printf("allo %e %e " ifmt "\n",lc[iinfo_[DG::IA_lc_face]+i],lc[iinfo_[DG::IA_lc_face]+npoints+i],iinfo_[DG::I_TETA_U_NBASIS]);
	  for (j=0;j<iinfo_[DG::I_TETA_U_NBASIS];++j)
	    {
	      printf("%e\n",rinfo_[iinfo_[DG::RA_EVAL_TETA_U] + iinfo_[DG::I_TETA_U_NBASIS]*i+j]);
	    }
	}
      exit(1);}
#endif


    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD */
    /* -------------------------------------------------------------------------------------- */        

    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_trial),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n],
	     bmat_n,
	     ferr_);
  
    mkS_bwji(qface_n_,
	     qface_p_,
	     qface_w_,
	     mkS_b(&s_teta),
	     mkS_b(&s_test),
	     &bmat_x[(teta_a_n+2*teta_u_n)*bmat_n+trial_nXtrial_n],
	     bmat_n,
	     ferr_);

  
    /* -------------------------------------------------------------------------------------- */
    /* MATRICES PONDEREES SUR LE BORD CONTRAPOSEE AVEC LE BORD */
    /* -------------------------------------------------------------------------------------- */            

    { pR bmatflux = &rinfo_[iinfo_[DG::RA_bmatx]];

      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_trial),
		   mkS_b(&s_test),
		   bmatflux,
		   bmat_n,
		   ferr_);    
    
      mkS_bwji_nei(qface_n_,
		   qface_p_,
		   qface_w_,
		   mkS_b(&s_teta),
		   mkS_b(&s_test),
		   &bmatflux[trial_nXtrial_n],
		   bmat_n,
		   ferr_); }
    return;
#if 0
    I s_teta_a[4]      = {iinfo_[DG::I_TETA_A_FAMILY],0,iinfo_[DG::I_TETA_A_DEGREE],iinfo_[DG::I_TETA_A_NBASIS]};
    I s_teta_u[4]      = {iinfo_[DG::I_TETA_U_FAMILY],0,iinfo_[DG::I_TETA_U_DEGREE],iinfo_[DG::I_TETA_U_NBASIS]};
    I s_teta[4]        = {iinfo_[DG::I_TETA_FAMILY],0,iinfo_[DG::I_TETA_DEGREE],iinfo_[DG::I_TETA_NBASIS]};
    I s_test[4]        = {iinfo_[DG::I_TEST_FAMILY],0,iinfo_[DG::I_TEST_DEGREE],iinfo_[DG::I_TEST_NBASIS]};
    I s_trial[4]       = {iinfo_[DG::I_TRIAL_FAMILY],0,iinfo_[DG::I_TRIAL_DEGREE],iinfo_[DG::I_TRIAL_NBASIS]};
#endif
  };

  static void define(DG_HANDLE* 	handle_,
		     cst_mkS 		s_teta_a_,
		     cst_mkS 		s_teta_u_,
		     cst_mkS 		s_teta_,
		     cst_mkS 		s_test_,
		     cst_mkS 		s_trial_,
		     cst_pI         	iinfo_n_,
		     pI			iinfo_,		  
		     cst_pI		rinfo_n_,
		     pR 		rinfo_,
		     cst_pI 		rwork_n_,
		     pR			rwork_)
  {
    
    pI ferr_ 	= &iinfo_[DG::ERR];  
    ferr_[0] 	= (I)0;  
    if (iinfo_[DG::I_QFACE_N]<1)
      {
	ferr_[0] =
	  DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QFACE_N]<1\n");
	exit(1);
      }
  
    if (iinfo_[DG::I_QELM_N]<1)
      {
	ferr_[0] = DGERR_USER;
	fprintf(stderr,"*** MOK_DGI_INIT ERROR IINFO[DG::I_QELM_N]<1\n");
	exit(1);
      }

    const I qelm_n_	[1] 	= {iinfo_[DG::I_QELM_N]*iinfo_[DG::I_QELM_N]};
    const I qface_n_	[1] 	= {iinfo_[DG::I_QFACE_N]};
    
    R qelm_p_	[64*64*2];
    R qelm_w_	[64*64];
    R qface_p_	[64];
    R qface_w_	[64];
  
    fprintf(stdout,"generate quadrature " ifmt " " ifmt "\n",
	    qelm_n_[0],
	    qface_n_[0]);
  
    {
      const I rwork_n = (rwork_n_[0]>3*iinfo_[DG::I_QELM_N])?rwork_n_[0]-3*iinfo_[DG::I_QELM_N]:(I)0;
      if (rwork_n<1)
	{
	  fprintf(stderr,"*** MOK_DGI_INIT ERROR TOO SMALL RWORK_N\n");
	  ferr_[0] = DGERR_MEMORY;
	  return;
	}
      mkQ_legendre_interval_(qface_n_,
			     qface_p_,
			     &negal1,
			     qface_w_,
			     &negal1,
			     rwork_,
			     rwork_n_,
			     ferr_);
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:first mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      mkQ_legendre_interval_(&iinfo_[DG::I_QELM_N],
			     rwork_,
			     &negal1,
			     &rwork_[2*iinfo_[DG::I_QELM_N]],
			     &negal1,
			     &rwork_[3*iinfo_[DG::I_QELM_N]],
			     rwork_n_,
			     ferr_);    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:second mkQ_legendre_interval_ failed (ferr_ = " ifmt "\n",ferr_[0]);
	  return;
	}
      I qelm_n=0;
      mkQ_collapse_		(&iinfo_[DG::I_QELM_N],
				 rwork_,
				 &negal1,
				 &rwork_[2*iinfo_[DG::I_QELM_N]],
				 &negal1,
				 &qelm_n,
				 qelm_p_,
				 qelm_n_,
				 qelm_w_,
				 &negal1,
				 ferr_);
    
      if (ferr_[0])
	{
	  ferr_[0] = 1;
	  fprintf(stderr,"*** MOK_DGI:mkQ_legendre_interval_ failed\n");
	  return;
	}
    }
    fprintf(stdout,"generate quadrature done\n");



    dgadvection_init_with_quadrature(s_teta_a_,
				     s_teta_u_,
				     s_teta_,
				     s_test_,
				     s_trial_,
		  
				     qelm_n_,
				     qelm_p_,
				     qelm_w_,
		  
				     qface_n_,
				     qface_p_,
				     qface_w_,

				     iinfo_n_,
				     iinfo_,

				     rinfo_n_,
				     rinfo_,

				     rwork_n_,
				     rwork_);


    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	na      = &iinfo_[DG::I_TETA_A_NBASIS];
    const I bmat_n 		= iinfo_[DG::I_bmat_n];  
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  


    
    pR brhs      		= &rwork_[0];
    pR brhs_uelm 		= &rwork_[1];
    I  npts_involving_u 	= bmat_m-na[0];/*nTot-teta_a_n*/    
    pR tmpbrhs        		= &rwork_[bmat_m]; 
    pR tmpbrhs_uelm   		= &tmpbrhs[na[0]]; 
    pR tmpbrhs_u 		= tmpbrhs_uelm;
    const I tmpbrhs_ufaceoff 	= bmat_m;
    
    cst_pR bmat	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    //    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    matrix_handle_def(&handle_->m_BMAT,
		      bmat_n,
		      bmat_m,(pR)bmat,
		      bmat_n);
    

    matrix_handle_def(&handle_->m_mat_tmpbrhs_uelm,
		      nu[0],
		      dim,
		      tmpbrhs_uelm,
		      tmpbrhs_ufaceoff);

    matrix_handle_def(&handle_->m_UVWDOFS,
		      nu[0],
		      dim,
		      brhs_uelm,
		      nu[0]);

    matrix_handle_def(&handle_->m_brhs_uvw,
		      npts_involving_u,
		      dim,
		      tmpbrhs_u,
		      tmpbrhs_ufaceoff);
    
    vector_handle_def(&handle_->m_BRHS,
		      bmat_m,
		      brhs,
		      1);

    matrix_handle_def(&handle_->m_EVALU,
		      nu[0],
		      npts_involving_u,
		      &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]],
		      iinfo_[DG::I_TETA_U_NBASIS]);
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {    
	matrix_handle_def(&handle_->m_mat_tmpbrhs_uface[localFaceIndex],
			  qface_n_[0],
			  dim,
			  &tmpbrhs[na[0]+dim*nu[0]+localFaceIndex*qface_n_[0]],
			  tmpbrhs_ufaceoff);
      }
    
    
  };
  
  static void solve(DG_HANDLE*  handle_,
		    DG_DATA*    data_,
		    cst_pR 	xa_,
		    cst_pR 	xu_,
		    cst_pR	rhs_,
		    cst_pI 	rhsoff_,
		    
		    cst_pI  	cnc_u_,
		    cst_pI  	cncoff_u_,
		    cst_pR 	data_u_,		       
		    cst_pR	data_v_,
		    
		    cst_pR	sol_,
		    cst_pI 	soloff_,
		    
		    pR		corr_,
		    cst_pI 	corroff_,
		    
		    cst_pR 	t_,
		    cst_pI 	nelm_,
		    cst_pR 	coo_,
		    cst_pI 	cooff_,
		    cst_pI 	cnc_,
		    cst_pI 	cncoff_,
		    cst_pI 	adj_,
		    cst_pI 	adjoff_,
		    cst_pI 	vcod_,
		    cst_pI 	noboundary_cod_,
		    
		    cst_pI 	rwork_n_,
		    pR  	rwork_,
		    cst_pI 	iwork_n_,
		    pI	iwork_,
		    
		    cst_pR 	rinfo_,
		    const I  	iinfo_[DG::I_n],
		    R 		rres_[DG_rres_n],
		    I 		ires_[DG_ires_n])
  {  
    
    cst_pI	trial_n = &iinfo_[DG::I_TRIAL_NBASIS];
    cst_pI	test_n  = &iinfo_[DG::I_TEST_NBASIS];
    cst_pI	teta_n  = &iinfo_[DG::I_TETA_NBASIS];
    cst_pI	nu      = &iinfo_[DG::I_TETA_U_NBASIS];
    cst_pI	teta_a_n      = &iinfo_[DG::I_TETA_A_NBASIS];
    cst_pI	qface_n = &iinfo_[DG::I_QFACE_N];
    const I q_nx3 		= qface_n[0]*3;
    const I bmat_n 		= iinfo_[DG::I_bmat_n];
    const I bmat_m 		= iinfo_[DG::I_bmat_m];  
    cst_pR eval_u_ 		= &rinfo_[iinfo_[DG::RA_EVAL_TETA_U]];
    const I eval_uoff_        	= iinfo_[DG::I_TETA_U_NBASIS];
    cst_pR bmat 	 		= &rinfo_[iinfo_[DG::RA_bmat]];
    cst_pR bmatflux 		= &rinfo_[iinfo_[DG::RA_bmatx]];
    cst_pR fpart[3][3];
    { I i,j;
      for (i=0;i<3;++i)
	for (j=0;j<3;++j)
	  fpart[i][j] = &bmatflux[bmat_n*( (i*3+j)*qface_n[0]  ) ];
    }
    
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	vector_handle_def(&data_->vec_uface[localFaceIndex],qface_n[0],&rwork_[iinfo_[DG::IA_lc_face]+localFaceIndex*qface_n[0]],1);
      }
  
    /*_____________________________________________________*/

    static const R
      r1=(R)1.0,
      r0=(R)0.0;
  
    static const I
      n1	= (I)1,
      n2 	= (I)2;
  
    /*_____________________________________________________*/
#if DG_CONSERVATIVE_WEAK_FORM
    const R mxu=xu_[0]*((R)-1.0);
#else
    const R mxu=xu_[0]*((R)1.0);
#endif

    /*_____________________________________________________*/

    I marker_size	= 0;
  
    pI
      lcperm 	= NULL,
      perm 	= NULL,
      graph 	= NULL,
      marker 	= NULL,
      marker_begin= NULL,
      marker_flag = NULL;

    pR
      brhs_a      = NULL,
      uface0 	= NULL,
      uface1 	= NULL,
      uface2 	= NULL;
  
    R mx;  
    I
      id,
      ielm,
      iter_gauss_seidel,
      b1,
      b2,
      pp,
      n,
      //      vcod[3*1],codface[8*1]
      cnc[3*1],neids[8*1],neids_face[8*1];
  
    R
      jacelm[1],
      x,
      longueurs[3],
      jacface[3*1],
      cooelm[6*1],
      lcsol[21];
    pR data_a_ = nullptr;

    //
    // VELOCITY EVALUATION
    //

    //
    // Pointers to bmat and bmatx
    //

#if __mk_debug__
    if (bmat_m!=1 + 2*nu[0]+3*qface_n[0]) 
      {
	fprintf(stderr,"*** DGERR " ifmt " " ifmt "",1 + 2*nu[0]+3*qface_n[0],bmat_m);exit(1);
      }
#endif
    //    I npts_involving_u = bmat_m-1;/*nTot-teta_a_n*/



    //
    //
    //

#if 0
    R 	matflux[21*(21+21)];/*trial_n*(trial_n+teta_n)*/
    pR 		matflux_corr 	= &matflux[0];
    pR 		matflux_residu 	= &matflux[trial_n[0]*trial_n[0]];

    matrix_handle_def(&handle_->m_FC,trial_n[0],trial_n[0],matflux_corr,trial_n[0]);
    matrix_handle_def(&handle_->m_FR,trial_n[0],teta_n[0],matflux_residu,trial_n[0]);
#endif

    // R lcrhs[21*2];
    const I renum_iwork_n = 3*nelm_[0] + 4*(nelm_[0]+1) + trial_n[0];
    ires_[DG_ires_err] = (I)0;
  
    /*_____________________________________________________*/
    if (rwork_n_[0]< bmat_m + 2 * bmat_m)
      {      
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_rw_n] = 3*bmat_m;
	fprintf(stderr,"too small rwork_n_ " ifmt " " ifmt "\n",3*bmat_m,rwork_n_[0]);
	return;
      }


    brhs_a    	= &rwork_[0];  
    uface0 	= &rwork_[iinfo_[DG::IA_lc_face]];
    uface1 	= &rwork_[iinfo_[DG::IA_lc_face]+qface_n[0]];
    uface2 	= &rwork_[iinfo_[DG::IA_lc_face]+2*qface_n[0]];

    pR ufaces[3] = {uface0,uface1,uface2};
  
    const I tmpbrhs_ufaceoff 	= bmat_m;
    pR tmpbrhs        	= &rwork_[bmat_m]; 
    pR tmpbrhs_ufaces 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface0 	= &tmpbrhs[1+2*nu[0]]; 
    pR tmpbrhs_uface1 	= &tmpbrhs[1+2*nu[0]+qface_n[0]];
    pR tmpbrhs_uface2 	= &tmpbrhs[1+2*nu[0]+2*qface_n[0]];
    
    //  pR tmpbrhs_uelm   	= &tmpbrhs[1]; 
    //  pR tmpbrhs_u	= tmpbrhs_uelm;
    /*_____________________________________________________*/
    if (renum_iwork_n>iwork_n_[0])
      {
	ires_[DG_ires_err] = (I)1;
	ires_[DG_ires_required_iw_n] = renum_iwork_n;
	fprintf(stderr,"too small iwork_n_ " ifmt " " ifmt "\n",renum_iwork_n,iwork_n_[0]);
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
	  data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];				
	  x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;						  
	  data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	  x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;  
	  data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	  x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	  jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	  /*____ INFO ELEMENT DONE */ 	
	  /*___________________________________________________________________________________________________________________*/  
	  /* 
	     eval u faces
	  */      
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
	  { I k;
	    for (k=0;k<nu[0];++k)
	      {
		data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
	      } }
  
	  /*___________________________________________________________________________________________________________________*/
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	  Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);


	  
	  Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	  Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	
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
	fprintf(stderr,"*** DGERR RENUM FAILED " ifmt "!=" ifmt "",b2,nelm_[0]);
	goto __state_error;
      }
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    /*___________________________________________________________________________________________________________________*/  
    iter_gauss_seidel = 0;
  __state_next_iter_gauss_seidel:

    fprintf(stdout,"next iter gauss seidel \n");
    if (0)
      {
	dg_print_mesh(mesh,
		      "dgres." ifmt, iter_gauss_seidel);
	dg_print_sol(3 * nelm_[0],
		     corr_,
		     "dgres." ifmt,iter_gauss_seidel);
      }
    mx = (R)0.0;
    /*___________________________________________________________________________________________________________________*/  
    ielm=0;      
  __state_next_elm:
    id 		= perm[ielm+1]-1;
    neids[0]	= adj_[id*adjoff_[0]];
    neids[1]	= adj_[id*adjoff_[0]+1];
    neids[2]	= adj_[id*adjoff_[0]+2];

    cnc[0]        = cnc_[id*cncoff_[0]+0];
    cnc[1]        = cnc_[id*cncoff_[0]+1];
    cnc[2]        = cnc_[id*cncoff_[0]+2];

#if 0    
    vcod[0] 	= vcod_[cnc[0]-1];
    vcod[1] 	= vcod_[cnc[1]-1];
    vcod[2] 	= vcod_[cnc[2]-1];
    codface[0] 	= (neids[0])? noboundary_cod_[0] :  MAX(vcod_[0],vcod_[1]);
    codface[1] 	= (neids[1])? noboundary_cod_[0] :  MAX(vcod_[1],vcod_[2]);
    codface[2] 	= (neids[2])? noboundary_cod_[0] :  MAX(vcod_[2],vcod_[0]);
#endif
    
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
    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
    jacelm[0]  = nsSQRT((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  

    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];
    data_->mat_belm.x[1] = cooelm[0]-cooelm[2];
    data_->mat_belm.x[2] = cooelm[3]-cooelm[4];
    data_->mat_belm.x[3] = cooelm[1]-cooelm[0];    
    x = data_->mat_belm.x[2];
    data_->mat_belm.x[2]=data_->mat_belm.x[1];
    data_->mat_belm.x[1]=x;
    
    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
    /*____ INFO ELEMENT DONE */ 
    /*____ COPY LOCAL RHS */ 
    {I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j] = rhs_[rhsoff_[0]*id+j];}	


    if (data_a_)
      {
	for (I k=0;k<teta_a_n[0];++k)
	  {
	    brhs_a[k] = jacelm[0] * xa_[0] * data_a_[k] ;
	  }
      }
    else
      {
	brhs_a[0] = jacelm[0] * xa_[0];
      }
      
  
    /*___________________________________________________________________________________________________________________*/  
    /* 
       eval u : elm and faces
    */  
    /*___________________________________________________________________________________________________________________*/  

    //
    // Copy the velocity.
    //
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }
    { I k;
      for (k=0;k<nu[0];++k)
	{
	  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[id*cncoff_u_[0]+k]-1];
	} }

    //
    // Wrap the solution.
    //
    vector_handle_def(&data_->hsol,teta_n[0],(pR)&sol_[id*soloff_[0]],1);
  
    //
    // Evaluate u at all the required coordinates.
    //
    handle_->m_brhs_uvw = handle_->m_EVALU.transpose() * data_->uvw_ldofs;
    
    /*____ TRANSFORM ON UELM */
    handle_->m_UVWDOFS = mxu * handle_->m_mat_tmpbrhs_uelm * data_->mat_belm.transpose();

    //
    // On the boundary.
    //
    for (I localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex] = jacface[localFaceIndex] * handle_->m_mat_tmpbrhs_uface[localFaceIndex] * data_->vec_nrmelm[localFaceIndex];
      }
    
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	data_->vec_uface[localFaceIndex].apply([xu_](const R r)
					       {
						 return (r<0.0) ? -r*xu_[0] : ((R)0.0);
					       });
      }
      
    /*____ APPLY FLUX ON UPWIND STREAM FOR CORRECTION AND SOLUTION */ 
    for (I localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
      {
	if (neids[localFaceIndex])
	  {
	    Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[localFaceIndex][neids_face[localFaceIndex]],&bmat_n,ufaces[localFaceIndex],&n1,&r0,data_->m_flux_memory,&n1);

	    
	    vector_handle_def(&data_->vec_neicorr,trial_n[0],  &corr_[(neids[localFaceIndex]-1)*corroff_[0]], 1);
	    vector_handle_def(&data_->vec_neisol,teta_n[0],   (pR)&sol_[(neids[localFaceIndex]-1)*soloff_[0]], 1);

	    
	    data_->m_local_rhs += data_->m_flux_part_corr * data_->vec_neicorr;	 
	    data_->m_local_rhs += data_->m_flux_part_sol  * data_->vec_neisol;	 
	  }
      }
  
#if DG_CONSERVATIVE_WEAK_FORM
    /*____ COMPUTE MAX(U.N,0), APPLY FLUX MATRIX ON DOWNWIND STREAM */ 
    Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface0[j] = (uface0[j]>0.0) ? uface0[j]*(xu_[0]) : (R)0.0;}
    Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface1[j] = (uface1[j]>0.0) ? uface1[j]*(xu_[0]) : (R)0.0;}
    Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);
    {I j;for (j=0;j<qface_n[0];++j)uface2[j] = (uface2[j]>0.0) ? uface2[j]*(xu_[0]) : (R)0.0;}
#endif
  
    /*--- COMPUTE MATRICES  */
    data_->m_local_matrices = handle_->m_BMAT * handle_->m_BRHS;
  
    /*--- COMPUTE LOCAL RESIDUAL  */
    data_->m_local_rhs -= data_->m_local_matrix_part_sol * data_->hsol;
  
    /*--- COMPUTE LOCAL CORRECTION   */
    matrix_handle_gesv(&data_->m_local_matrix_part_corr,&data_->m_local_rhs,lcperm);
  
    /*--- COMPUTE NRMS OF LOCAL CORRECTION */
    { I j; for (j=0;j<test_n[0];++j) {lcsol[j] = corr_[corroff_[0]*id+j];}}
    { I j;for (j=0;j<test_n[0];++j) corr_[corroff_[0]*id+j] = data_->m_local_rhs.x[j];}
    { I j;for (j=0;j<test_n[0];++j) data_->m_local_rhs.x[j]-=lcsol[j];}
    { R xx=(R)0.0; {I j;for (j=0;j<test_n[0];++j) xx+=data_->m_local_rhs.x[j]*data_->m_local_rhs.x[j];}
      if (mx<xx) mx = xx; }		    
  
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
      fprintf(stdout,"*** DGINFO time %.5f iter " ifmt "/" ifmt " tol = %8.15e\n",t_[0],iter_gauss_seidel,nelm_[0],mx);
      if ( (mx>((R)1.0e-24)) AND (iter_gauss_seidel<nelm_[0]) )
	{
	  ++iter_gauss_seidel;
	  goto __state_next_iter_gauss_seidel;
	}
      else
	{	
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

      rres_[DG_rres_areaL1]	= 0.0; 
      rres_[DG_rres_nrmL2]	= 0.0; 
      rres_[DG_rres_nrmLInf]	= 0.0; 
      rres_[DG_rres_jumpL2]	= 0.0; 
      rres_[DG_rres_johnson]	= 0.0;   
      { I jelm;
	for (jelm=0;jelm<nelm_[0];++jelm)
	  { 
	    /**/
	    neids[0]   = adj_[jelm*adjoff_[0]];neids[1]	= adj_[jelm*adjoff_[0]+1];neids[2]	= adj_[jelm*adjoff_[0]+2];

	    cnc[0]        = cnc_[jelm*cncoff_[0]+0];
	    cnc[1]        = cnc_[jelm*cncoff_[0]+1];
	    cnc[2]        = cnc_[jelm*cncoff_[0]+2];	  
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
	    data_->nrmelm[0] = cooelm[4] - cooelm[3];data_->nrmelm[1] = cooelm[0] - cooelm[1];					
	    x = ((R)1.0)/jacface[0];data_->nrmelm[0] *= x;data_->nrmelm[1] *= x;							  
	    data_->nrmelm[2] = cooelm[5] - cooelm[4];data_->nrmelm[3] = cooelm[1] - cooelm[2];
	    x = ((R)1.0)/jacface[1];data_->nrmelm[2] *= x;data_->nrmelm[3] *= x;	  
	    data_->nrmelm[4] = cooelm[3] - cooelm[5];data_->nrmelm[5] = cooelm[2] - cooelm[0] ;
	    x = ((R)1.0)/jacface[2];data_->nrmelm[4] *= x;data_->nrmelm[5] *= x;
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 
	    if (longueurs[0]<longueurs[1]) { x = longueurs[0]; longueurs[0] = longueurs[1]; longueurs[1] = x; } 
	    if (longueurs[1]<longueurs[2]) { x = longueurs[1]; longueurs[1] = longueurs[2]; longueurs[2] = x; } 	  
	    jacelm[0]  = nsSQRT((longueurs[0]+longueurs[1]+longueurs[2])*(longueurs[2]-(longueurs[0]-longueurs[1]))*(longueurs[2]+(longueurs[0]-longueurs[1]) )*(longueurs[0]+(longueurs[1]-longueurs[2]) ))/((R)2.0); 	  
	    data_->mat_belm.x[0] = cooelm[5]-cooelm[3];data_->mat_belm.x[1] = cooelm[0]-cooelm[2];data_->mat_belm.x[2] = cooelm[3]-cooelm[4];data_->mat_belm.x[3] = cooelm[1]-cooelm[0];  	  
	    x = data_->mat_belm.x[2];data_->mat_belm.x[2]=data_->mat_belm.x[1];data_->mat_belm.x[1]=x;	  
	    jacface[0] *=(R)0.5;jacface[1] *=(R)0.5;jacface[2] *=(R)0.5;	  
	    /**/	  

	    Blas_dgemv("N",trial_n,trial_n,jacelm,bmat,trial_n,&corr_[jelm*corroff_[0]],&n1,&r0,data_->m_local_rhs.x,&n1);

	    rres_[DG_rres_nrmL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);

	    { I j;
	      for (j=0;j<trial_n[0];++j)
		{
		  rres_[DG_rres_areaL1] += data_->m_local_rhs.x[j];
		  if (rres_[DG_rres_nrmLInf]<nsFABS(corr_[corroff_[0]*jelm+j]))
		    {
		      rres_[DG_rres_nrmLInf]=nsFABS(corr_[corroff_[0]*jelm+j]);
		    }
		} }

	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	    { I k;
	      for (k=0;k<nu[0];++k)
		{
		  data_->uvw_ldofs.x[data_->uvw_ldofs.ld+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		} }
	  
	    /* u au points de gauss de l element */
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x,&n1,&r0,tmpbrhs_ufaces,&n1);
	    Blas_dgemv("T",nu,&q_nx3,&r1,&eval_u_[eval_uoff_*(iinfo_[DG::IA_lc_face]-1)],&eval_uoff_,data_->uvw_ldofs.x + data_->uvw_ldofs.ld,&n1,&r0,&tmpbrhs_ufaces[tmpbrhs_ufaceoff],&n1);
	  
	    /* u.n sur les faces */
	    Blas_dgemv("N",qface_n,&n2,&jacface[0],tmpbrhs_uface0,&tmpbrhs_ufaceoff,&data_->nrmelm[2*0],&n1,&r0,uface0,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[1],tmpbrhs_uface1,&tmpbrhs_ufaceoff,&data_->nrmelm[2*1],&n1,&r0,uface1,&n1);
	    Blas_dgemv("N",qface_n,&n2,&jacface[2],tmpbrhs_uface2,&tmpbrhs_ufaceoff,&data_->nrmelm[2*2],&n1,&r0,uface2,&n1);	
	  
	    /* mise en a zero */
	    { I j;for (j=0;j<qface_n[0];++j) uface0[j] = (uface0[j]<0.0) ? uface0[j]*(-xu_[0]) : (R)0.0;}       
	    { I j;for (j=0;j<qface_n[0];++j) uface1[j] = (uface1[j]<0.0) ? uface1[j]*(-xu_[0]) : (R)0.0;}	
	    { I j;for (j=0;j<qface_n[0];++j) uface2[j] = (uface2[j]<0.0) ? uface2[j]*(-xu_[0]) : (R)0.0;}
	  	  
	    { I k;
	      for (k=0;k<3;++k)
		{
		  if ( neids[k] )
		    { 
		      /* 
			 on forme la matrice de flux pour l'arete
		      */
		      Blas_dgemv("N",&bmat_n,qface_n,&r1,fpart[k][neids_face[k]],&bmat_n,(k==2)?uface2:((k==1)?uface1:uface0),&n1,&r0,data_->m_flux_memory,&n1);
		      
		      /* 
			 on calcule le produit scalaire
		      */
		      Blas_dgemv("N",trial_n,trial_n,&jacface[k],data_->m_flux_part_corr.x,trial_n,&corr_[(neids[k]-1)*corroff_[0]],&n1,&r1,data_->m_local_rhs.x,&n1);	   
		      rres_[DG_rres_jumpL2]+=Blas_ddot(trial_n,data_->m_local_rhs.x,&n1,&corr_[corroff_[0]*jelm],&n1);/* saut */	
		    } 
		} }
	  
	  } }

      printf("compute nrms\n"
	     "   areaL1  : %e\n"
	     "   nrmL2   : %e\n"
	     "   nrmLInf : %e\n"
	     "   jumpL2  : %e\n"
	     "   jhonson : %e\n",
	     rres_[DG_rres_areaL1],
	     rres_[DG_rres_nrmL2],
	     rres_[DG_rres_nrmLInf],
	     rres_[DG_rres_jumpL2],
	     rres_[DG_rres_johnson]);
    
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

};



int main(int 		argc,
	 const char**	argv)
{
  

  //
  // Set up the monitor
  //
  Monitor_def(0,
	      argv[0],
	      MonitorMode_STD);
  
  mesh = (ns_mesh*)calloc(1,sizeof(ns_mesh));
  STR errmsg;
  Err err;
  ns_mesh_read(mesh,
	       errmsg,
	       &err,
	       argv[1]);



#if 0
  ns_mesh_write_medit(mesh,
		      out);
#endif
  
  
  I 		dg_iinfo[DG::I_n];
  R		dg_rres	[DG_rres_n];
  I		dg_ires	[DG_ires_n];
  I 		dg_rinfo_n;
  I 		dg_iinfo_n = DG::I_n;
  pR		dg_rinfo;
  I 		dg_rwork_n;
  pR		dg_rwork;
  I 		dg_iwork_n;
  pI		dg_iwork;
    
  dg_iinfo[DG::I_QELM_N]	= 10;
  dg_iinfo[DG::I_QFACE_N]  	= 10;  
  dg_rinfo_n			= (I)640000;
  dg_rinfo 			= (pR)malloc(sizeof(R)*dg_rinfo_n);
  dg_rwork_n			= (I)1280000;
  dg_rwork 			= (pR)malloc(sizeof(R)*dg_rwork_n);
  dg_iwork_n			= (I)1280000;
  dg_iwork 			= (pI)malloc(sizeof(I)*dg_iwork_n);

  mkS_st shape_A;
  mkS_st shape_F;
  mkS_st shape_U;
  mkS_definit	(&shape_A,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 0,
		 __emk_discontinuous,
		 &err);
    
  mkS_definit	(&shape_F,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 1,
		 __emk_discontinuous,
		 &err);
    
  mkS_definit	(&shape_U,
		 __eTopology_TRIANGLE,
		 __emkS_FAMILY_lagrange,
		 1,
		 __emk_discontinuous,
		 &err);		    

  DG_HANDLE h;
  DG::define(&h,
	     &shape_A,
	     &shape_U,
	     &shape_F,
	     &shape_F,
	     &shape_F,
	     &dg_iinfo_n,
	     dg_iinfo,
	     &dg_rinfo_n,
	     dg_rinfo,
	     &dg_rwork_n,
	     dg_rwork); 
  
  const I nelm 	= ns_mesh_nelm(mesh);
  const I numNodes 	= ns_mesh_get_numNodes(mesh);
  const I fn = mkS_n(&shape_F);
  const I len = fn * nelm;
  
  double * sol		= (double*)calloc(len,sizeof(R)); 
  double * rhs 		= (double*)calloc(len,sizeof(R)); 
  double * corr 	= (double*)calloc(len,sizeof(R)); 
  double * u 		= (double*)calloc(numNodes,sizeof(R)); 
  double * v 		= (double*)calloc(numNodes,sizeof(R));
  
  for (int i=0;i<numNodes;++i)
    {
      u[i] = 1.0;
    }
    
  const I rhsoff 	= fn;
  const I soloff 	= fn;
  const I corroff 	= fn;
    
  const double xu 		= 1.0;
  const double t 		= 1.0;
  const I noboundary_cod  	= 100;

  const I  cooff = 2;
  const I cncoff = 6;
  const I adjoff = 3;

  cst_pI cnc_u = mesh->cnc;
  const I cncoff_u = 6;
  //    exit(1);


  
  for (I jelm=0;jelm<mesh->nelm;++jelm)
    {
      for (int jadj=0;jadj<3;++jadj)
	{
	  if (mesh->adj[jelm*3+jadj] == 0)
	    {
	      I cncelm[6];
	      ns_mesh_get_cellToNodes(mesh,&jelm,cncelm);
	      const double * p0 = &mesh->coo[cncelm[0]*2+0];
	      const double * p1 = &mesh->coo[cncelm[1]*2+0];
	      const double * p2 = &mesh->coo[cncelm[2]*2+0];
	      const double * p[3] = {p0,p1,p2};
	      const I jedge   = mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
	      const R nx	= mesh->normaledge[2*jedge+0];
	      const R ny 	= mesh->normaledge[2*jedge+1];
	      const R xa      = mesh->jacedge[jedge];	    
	       		
	      //	{I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",s_->coo[cncelm[j]*2+0],s_->coo[cncelm[j]*2+1],s_->cod[cncelm[j]]); } }
	      // sol[jelm*3+jadj] = 1.0;
	      // sol[jelm*3+(jadj+1)%3] = 1.0;
	      //
	      //
	      if (nx * 1 + ny * 0.0<0.0)
		{
		  double iw[6];
		  double ip[6];
		  iw[0] = 1.713244923791704e-1;
		  iw[1] = 3.607615730481386e-1;
		  iw[2] = 4.679139345726911e-1;
		  iw[3] = 4.679139345726911e-1;
		  iw[4] = 3.607615730481386e-1;
		  iw[5] = 1.713244923791704e-1;
		
		  ip[0] = -9.32469514203152e-1 * 0.5 + 0.5;
		  ip[1] = -6.612093864662645e-1 * 0.5 + 0.5;
		  ip[2] = -2.386191860831969e-1 * 0.5 + 0.5;
		  ip[3] = 2.386191860831969e-1 * 0.5 + 0.5;
		  ip[4] = 6.612093864662645e-1 * 0.5 + 0.5;
		  ip[5] = 9.32469514203152e-1 * 0.5 + 0.5;
		
		  double v0 = 0.0,v1=0.0,v2 = 0.0;
		  for (int iq=0;iq<6;++iq)
		    {
		      double r,s;
		      if (jadj==0)
			{
			  r = ip[iq];
			  s = 0.0;
			}
		      else if (jadj==1)
			{
			  r = 1.0 - ip[iq];
			  s = ip[iq];

			}
		      else if (jadj==2)
			{
			  r = 0.0;
			  s = 1.0 - ip[iq];
			}
		      double x = (1.0-r-s) * p[0][0] + r * p[1][0]+ s * p[2][0];
		      double y = (1.0-r-s) * p[0][1] + r * p[1][1]+ s * p[2][1];
		      double psi0 = 1.0-r-s;
		      double psi1 = r;
		      double psi2 = s;
		      //		std::cout << "x= " << x << std::endl;
		      if ( (x < 0.0 ? -x : x)>0.0000001)
			{
			  exit(1);
			}
		      //		std::cout << "y= " << y << std::endl;
		      //		      v0 += psi0 * sin(12.0*y) * iw[iq];
		      //		      v1 += psi1 * sin(12.0*y)* iw[iq];
		      //v2 += psi2 * sin(12.0*y)* iw[iq];
		      v0 += psi0 * 1.0 * iw[iq];
		      v1 += psi1 * 1.0* iw[iq];
		      v2 += psi2 * 1.0* iw[iq];
		    }
		  v0 *= xa;
		  v1 *= xa;
		  v2 *= xa;
		  rhs[3*jelm+0] += v0;
		  rhs[3*jelm+1] += v1;
		  rhs[3*jelm+2] += v2;
		}
	    }
	}
    }
  for (I i=0;i<3;++i)
    {
      std::cout << rhs[i] << std::endl;
    }
  dg_print_mesh(mesh,
		"dg");
  
  dg_print_sol(fn * nelm,
	       rhs,
	       "dg");
  
  double xa = 0.0;

  
  DG_DATA *data  = DG::create_data(&shape_A,
				   &shape_U,
				   &shape_F,
				   &shape_F,
				   &shape_F);
  
  DG::solve(&h,
	    data,
	    &xa,
	    &xu,
	   rhs,
	   &rhsoff,
		      
	   cnc_u,
	   &cncoff_u,
	   u,
	   v,
		      
	   sol,
	   &soloff,
		      
	   corr,
	   &corroff,
		      
	   &t,
		      
	   &nelm,
	   mesh->coo,
	   &cooff,
	   mesh->cnc,
	   &cncoff,
	   mesh->adj,
	   &adjoff,
	   mesh->cod,
	   &noboundary_cod,

	   &dg_rwork_n,
	   dg_rwork,
	   &dg_iwork_n,
	   dg_iwork,
	   dg_rinfo,
	   dg_iinfo,
	   &dg_rres[0],
	   &dg_ires[0]);
  
  dg_print_mesh(mesh,
		"dgres");
  dg_print_sol(fn * nelm,
	       corr,
	       "dgres");
  return 0;
}



