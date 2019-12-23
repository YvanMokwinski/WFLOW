#include "Config.h"
#include "Type.h"

#include "SparseBlock.h"
#include "eSparseIterativeMethod.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include "SparseStokes.h"
#include "SparseIterative.h"
#include "Monitor.h"
#include "Blas.h"
#include <string.h>



#ifndef NDEBUG
#define wrong_param_Sparse(_cond,_i,_s)					\
  if ( (_cond) )							\
    {									\
      fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); \
      exit(1);								\
    } ((void)0)
#endif

I 	Sparse_get_n		(cst_pSparse const self_) 	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_n);
#endif
  return self_->n; 
}

L 	Sparse_is_fortran_index	(cst_pSparse const self_) 	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(!self_,1,Sparse_get_n);
#endif
  return ( self_->format == 1) ? __emnsYES : __emnsNO; 
}


I 	Sparse_get_m		(cst_pSparse const self_) 	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_m);
#endif
  return self_->m; 
}

I 	Sparse_get_nc		(cst_pSparse const self_) 	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_nc);
#endif
  return self_->nc; 
}

cst_pI Sparse_get_i		(cst_pSparse const self_)	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_i);
#endif
  return self_->i; 
}

cst_pI Sparse_get_b		(cst_pSparse const self_) 	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_b);
#endif
  return self_->b; 
}

cst_pR Sparse_get_x		(cst_pSparse const self_)	
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_x);
#endif
  return self_->x; 
}

pI 	Sparse_get_owni	(pSparse const self_)		
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_owni);
#endif
  return self_->own_i; 
}

pI 	Sparse_get_ownb	(pSparse const self_) 		
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_ownb);
#endif
  return self_->own_b; 
}

pR 	Sparse_get_ownx	(pSparse const self_)		
{ 
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_get_ownx);
#endif
  return self_->own_x; 
}

void Sparse_clr(pSparse const self_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_clr);
#endif
  memset(self_,0,sizeof(Sparse));
}

void Sparse_clear(pSparse const self_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_clear);
  wrong_param_Sparse(! self_->own_x,1,Sparse_clear);
#endif
  const I N 	= Sparse_get_nc(self_);
  pR r 		=  Sparse_get_ownx(self_);
  { I i 	= 0;
    for (i=0;i<N;++i)
      {
	r[i]=((R)0.0);
      } }
}


void Sparse_free(pSparse const self_)
{
  if (self_)
    {
      if (self_->own_x)
	free(self_->own_x);
      if (self_->own_i)
	free(self_->own_i);
      if (self_->own_b)
	free(self_->own_b);
      Sparse_clr(self_);
    }
}

pSparse  Sparse_kill(pSparse const self_)
{
  if (self_)
    {
      Sparse_free(self_);      
      free(self_);
    }
  return NULL;
}

void Sparse_free_thread(pSparse const self_,const I iproc_)
{
  if (self_)
    {
      if (self_->own_x)
	free(self_->own_x);
      if (self_->own_i)
	free(self_->own_i);
      if (self_->own_b)
	free(self_->own_b);
      Sparse_clr(self_);
    }
}

void Sparse_fortran_indexation(pSparse const self_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_fortran_indexation);
#endif
  
  { pI b	= Sparse_get_ownb(self_); const I N 	= Sparse_get_n(self_);  { I i; for (i=0;i<=N;++i) *(b++) +=1;} }
  { pI ii 	= Sparse_get_owni(self_); const I NC	= Sparse_get_nc(self_); { I i; for (i=0;i<NC;++i) *(ii++)+=1;} }
  self_->format = 1;

}

void Sparse_c_indexation(pSparse const self_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_c_indexation);
#endif
  { pI b	= Sparse_get_ownb(self_); const I N 	= Sparse_get_n(self_);  { I i; for (i=0;i<=N;++i) *(b++) -=1;} }
  { pI ii 	= Sparse_get_owni(self_); const I NC	= Sparse_get_nc(self_); { I i; for (i=0;i<NC;++i) *(ii++)-=1;} }
  self_->format = 0;
}


void Sparse_spy(cst_pSparse  const 	self_,
		const char* 			filename_,
		...)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_spy);
  wrong_param_Sparse(! filename_,2,Sparse_spy);
#endif
  cst_pI begin 	= Sparse_get_b(self_);
  cst_pI index 	= Sparse_get_i(self_);
  const I N 	= Sparse_get_n(self_);
  FILE * fich 		= fopen(filename_,"w");
  if (self_->format==1)
    {
      { I i;
	for (i=0;i<N;++i)
	  {
	    const I NC = begin[1]-begin[0];
	    { I j;
	      for (j=0;j<NC;++j)
		{
		  fprintf(fich,""ifmt" "ifmt"\n",*(index++)-1,N-i);
		} }
	    ++begin;
	  } }
    }
  else
    {

      { I i;
	for (i=0;i<N;++i)
	  {
	    const I NC = begin[1]-begin[0];
	    { I j;
	      for (j=0;j<NC;++j)
		{
		  fprintf(fich,""ifmt" "ifmt"\n",*(index++),N-i);
		} }
	    ++begin;
	  } }

    }
  fclose(fich);
  return;
}


void Sparse_write(cst_pSparse  const 	self_,
		  const char* 		filename_,...)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_write);
  wrong_param_Sparse(! filename_,2,Sparse_write);
#endif
  cst_pI begin 	= Sparse_get_b(self_);
  cst_pI index 	= Sparse_get_i(self_);
  cst_pR x 	= Sparse_get_x(self_);
  const I N 	= Sparse_get_n(self_);
  FILE * fich	= fopen(filename_,"w");
  if (self_->format==1)
    {
      { I i;
	for (i=0;i<N;++i)
	  {
	    { I j;
	      for (j=begin[i];j<begin[i+1];++j)
		{
		  fprintf(fich,""ifmt" "ifmt" "rfmt"\n",index[j-1],N - i,x[j-1]);
		} }
	  } }
    }
  else
    {
      { I i;
	for (i=0;i<N;++i)
	  {
	    { I j;
	      for (j=begin[i];j<begin[i+1];++j)
		{
		  fprintf(fich,""ifmt" "ifmt" "rfmt"\n",index[j],N - i,x[j]);
		} }
	  }  }
    }
  fclose(fich);
  return;
}

#if 0
static int compar_lexicorder(const void *a_, const void *b_)
{
  cst_pI  a = (cst_pI)a_;
  cst_pI  b = (cst_pI)b_;
  if (b[0]<a[0])
    return 1;
  else if (b[0]>a[0])
    return -1;
  else 
    return 0;
}


static void sort_pour_divector(pI C,const I g,const I  d)
{
  if ( g>=d ) return;
  else
    {
      I i,der;
      I tmpv;      
      tmpv 		=  C[g];
      C[g] 		=  C[( (g+d)/2 )];
      C[( (g+d)/2 )]	= tmpv;
      der=g;
      for(i=g+1; i<=d; i++)
	{
	  if ( C[i] > C[g] )
	    {
	      der++;
	      tmpv 	=  C[der];
	      C[der] 	=  C[i];
	      C[i]	=  tmpv;
	    }
	}/**/
      
      tmpv 	=  C[der];
      C[der] 	=  C[g];
      C[g]	=  tmpv;
      
      sort_pour_divector(C,g,der-1);
      sort_pour_divector(C,der+1,d);

    }
}/*-- --*/
#endif
#define mxmx 32

void mok_iarray_sort1_incr(pI 		C,
			   const I 	g_,
			   const I 	d_)
{
  I pile_g[64];
  I pile_d[64];
  I pile_n=1,g,d;
  pile_g[pile_n-1]=g_;
  pile_d[pile_n-1]=d_;
 __mok__sort1_incr_unstack:
  --pile_n;
  g = pile_g[pile_n];
  d = pile_d[pile_n];
  if ( g<d ) 
    {
       I mid = ( (g+d)/2 );
       I i,der;
       I tmp;
       tmp 		= C[g];
      C[g]		= C[mid];
      C[mid]		= tmp;      
      der=g;
      for(i=g+1; i<=d; i++)
	{
	  if ( C[i] < C[g] )
	    {
	      der++;
	      tmp 	=  C[der];
	      C[der] 	=  C[i];
	      C[i]	=  tmp;
	    }
	}/**/	  
      tmp 		=  C[der];
      C[der] 		=  C[g];
      C[g]		=  tmp;
      pile_g[pile_n]	=  der+1;
      pile_d[pile_n]	=  d;
      ++pile_n;
      pile_g[pile_n]	=g;
      pile_d[pile_n]	=der-1;
      ++pile_n;
    }
  if (pile_n>0)
    goto  __mok__sort1_incr_unstack;
}/*-- --*/





void Sparse_sort(pSparse const self_)
{  
#ifndef NDEBUG
  wrong_param_Sparse(!self_,1,Sparse_sort);
#endif
  { I j;
    const I N 		= Sparse_get_n(self_);
    cst_pI begin 	= Sparse_get_b(self_);
    pI     ownindex 	= Sparse_get_owni(self_);
    for (j=0;j<N;++j,++begin)
      {
#if 0
	qsort(&ownindex[*begin],begin[1]-begin[0],sizeof(I),compar_lexicorder);
#else
	mok_iarray_sort1_incr(&ownindex[*begin],
			      0,
			      begin[1]-begin[0]-1);
#endif
      } }
}



typedef struct
{
  I start;
  I end;
  pSparse self;
  I ithread;
} struct_Sparse_sortParallel;

void* routine_Sparse_sortParallel(void * usrptr_)
{  
  struct_Sparse_sortParallel * pll 	= (struct_Sparse_sortParallel * )usrptr_;
  pSparse const self_ 			= pll->self;
  const I start 			= pll->start;
  const I stop 				= pll->end;
  printf("hello "ifmt" "ifmt"\n",start,stop);
  { I j;
    cst_pI begin 	= Sparse_get_b(self_) + start;
    pI     ownindex 	= Sparse_get_owni(self_);
    for (j=start;j<stop;++j,++begin)
      {
#if 1
	mok_iarray_sort1_incr(&ownindex[*begin],
				       0,
				       begin[1]-begin[0]-1);
#else
	qsort(&ownindex[*begin],begin[1]-begin[0],sizeof(I),compar_lexicorder);
#endif
      } }
  pthread_exit((void*)usrptr_);
  return NULL;
}

void Sparse_sortParallel(pSparse const self_)
{  
#ifndef NDEBUG
  wrong_param_Sparse(!self_,1,Sparse_sort);
#endif
#define nproc 16
  struct_Sparse_sortParallel plls[nproc];
  const I N 		= Sparse_get_n(self_);
  const I n 		= N/nproc;
  I Ns[nproc+1];
  Ns[0] = 0;
  { I i;
    for (i=1;i<nproc;++i)
      {
	Ns[i] = Ns[i-1] + n;
      } }
  Ns[nproc] = N;
  int			status;
  pthread_t 		thread[nproc];
  pthread_attr_t 	attr;
  pthread_mutex_t 	mutex_mkPARALLEL;    
  pthread_mutex_init ( &mutex_mkPARALLEL, NULL );    
  pthread_attr_init ( &attr ); 
  pthread_attr_setdetachstate ( &attr, PTHREAD_CREATE_JOINABLE );        
  { I i;
    for (i=0;i<nproc;++i)
      {
	plls[i].ithread = i;
	plls[i].self  = self_;
	plls[i].start = Ns[i];
	plls[i].end   = Ns[i+1];
	pthread_create(&thread[i],&attr,routine_Sparse_sortParallel,(void*)&plls[i]);  
      } }        
  { I i;
    for (i=0;i<nproc;++i)
      {
	pthread_join(thread[i],(void**)&status);
      } }
  pthread_attr_destroy ( &attr );
  pthread_mutex_destroy (&mutex_mkPARALLEL); 
}



void Sparse_dirichlet(pSparse const 	self_,
		      const I 		n_,
		      cst_pI 		dir_,
		      const I 		dir_dec_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_dirichlet);
  wrong_param_Sparse(! dir_,3,Sparse_dirichlet);
#endif
  { I i; 
    for (i=0;i<n_;++i) 
      { 
	const I k = dir_dec_+dir_[i];
	{ I j;
	  for (j=self_->b[k];j<self_->b[k+1];++j)
	    { 
	      self_->own_x[j-1] = (self_->i[j-1]!=k+1) ? ((R)0.0) : ((R)1.0); 
	    } }
      } }

}

void Sparse_dirichlet21(pSparse const self_,pR Arhs_,const I n_,cst_pI dir_,const I dir_dec_)
{
#ifndef NDEBUG
  wrong_param_Sparse(! self_,1,Sparse_dirichlet21);
  wrong_param_Sparse(! dir_,3,Sparse_dirichlet21);
#endif
  { I i; 
    for (i=0;i<n_;++i) 
      { 
	const I k = dir_dec_+dir_[i]; 
	Arhs_[k] = ((R)0.0);
	{ I j;
	  for (j=self_->b[k];j<self_->b[k+1];++j)
	    {
	      self_->own_x[j-1] = (self_->i[j-1]!=k+1) ? ((R)0.0) : ((R)1.0); 
	    } }
      } } 
}


void Sparse_gemvtr(cst_pR 			a_,
		   cst_pSparse  const 		self_,
		   cst_pR 			rhs_,
		   pR 				sol_)
{

#ifndef NDEBUG
  wrong_param_Sparse(! a_,1,"Sparse_gemvtr");
  wrong_param_Sparse(! self_,2,"Sparse_gemvtr");
  wrong_param_Sparse(! rhs_,3,"Sparse_gemvtr");
  wrong_param_Sparse(! sol_,5,"Sparse_gemvtr");
  wrong_param_Sparse((self_->format!=1)&&(self_->format!=0),2,"Sparse_gemvtr");
#endif
  if (self_->format==1)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {
	    for (j=self_->b[i];j<self_->b[i+1];++j)
	      {
		sol_[self_->i[j-1]-1] +=  a_[0] * self_->x[j-1]*rhs_[i];
	      } 
	  } }
    }
  else if (self_->format==0)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {	    
	    for (j=self_->b[i];j<self_->b[i+1];++j)
	      {
		sol_[self_->i[j]] += a_[0]*self_->x[j]*rhs_[i] ;
	      } 
	  } }
    }
}



void Sparse_atmux(  cst_pSparse  const	self_,
		    pR 			y_,
		    cst_pR 			x_)
{

  { I i;
    for (i=0;i<self_->m;++i)
      {
	y_[i] = ((R)0.0);
      } }
  if (self_->format==1)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {
	    for (j=self_->b[i];j<self_->b[i+1];++j)
	      {
		y_[self_->i[j-1]-1] +=  self_->x[j-1]*x_[i];
	      } 
	  } }
    }
  else if (self_->format==0)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {	    
	    for (j=self_->b[i];j<self_->b[i+1];++j)
	      {
		y_[self_->i[j]] += self_->x[j]*x_[i] ;
	      } 
	  } }
    }
}


void Sparse_asstmux(  cst_pSparse  const	self_,
		      pR 			y_,
		      cst_pR 			x_)
{

  if (self_->format==1)
    {
      { I i;
	for (i=0;i<self_->n;++i)
	  {
	    { I j;
	      for (j=self_->b[i];j<self_->b[i+1];++j)
		{
		  y_[self_->i[j-1]-1] +=  self_->x[j-1]*x_[i];
		} }
	  } }
    }
  else if (self_->format==0)
    {
      { I i;
	for (i=0;i<self_->n;++i)
	  {
	    { I j;
	      for (j=self_->b[i];j<self_->b[i+1];++j)
		{
		  y_[self_->i[j]] += self_->x[j]*x_[i] ;
		}  }
	  } }
    }
}




pSparse 	Sparse_buildReference		(const I 	n_,
						 const I 	m_,
						 const I 	nc_,
						 cst_pI		b_,
						 cst_pI		i_,
						 cst_pR 	x_)
{
  pSparse self  	= (pSparse)calloc(1,sizeof(Sparse));
  self->n 		= n_;
  self->m 		= m_;
  self->nc		= nc_;
  self->b 		= b_;
  self->i 		= i_;
  if (x_)
    {
      self->x 		= x_;
    }
  else
    {
      self->own_x	= (pR)calloc(self->nc,sizeof(R));
      self->x		= self->own_x;
    }
  self->format 		= 1;
  return self;

}



pSparse 	Sparse_clone	(cst_pSparse const self_)
{
  pSparse self  	= (pSparse)calloc(1,sizeof(Sparse));
  self->n 		= self->n;
  self->m 		= self->m;
  self->nc		= self->nc;
  self->b 		= self->b;
  self->i 		= self->i;

  self->own_x		= (pR)calloc(self->nc,sizeof(R));
  self->x		= self->own_x;

  self->format 		= self->format;
  return self;

}



pSparse Sparse_build(const I 	n_,
		     const I 	m_,
		     const I 	nc_,
		     pI 	own_b_,
		     pI 	own_i_,
		     pR 	own_x_)
{
  pSparse self  	= (pSparse)calloc(1,sizeof(Sparse));
  self->n 		= n_;
  self->m 		= m_;
  self->nc		= nc_;
  self->own_b 		= own_b_;
  self->own_i 		= own_i_;
  self->own_x 		= (own_x_)? own_x_ : (pR)calloc(self->nc,sizeof(R));
  self->x		= self->own_x;
  self->i		= self->own_i;
  self->b		= self->own_b;
  self->format 		= 0;
  //  self->format 		= 1;
  return self;
}

pSparse	Sparse_new			(const I 	n_,
					 const I 	m_,
					 const I 	nc_)
{
  pSparse self  	= (pSparse)calloc(1,sizeof(Sparse));
  self->n 		= n_;
  self->m 		= m_;
  self->nc		= nc_;
  self->own_b 		= (pI)calloc(self->n+1,sizeof(I));
  self->own_i 		= (pI)calloc(self->nc,sizeof(I));
  self->own_x 		= (pR)calloc(self->nc,sizeof(R));
  self->x		= self->own_x;
  self->i		= self->own_i;
  self->b		= self->own_b;
  self->format 		= 0;
  return self;
}





void Sparse_amux( cst_pSparse const 		self_,
		  pR 			y_,
		  cst_pR 			x_)
{

  if (self_->format==1)
    {
      { I i,j;
	R s;
	for (i=0;i<self_->n;++i)
	  {
	    for (s=((R)0.0),j=self_->b[i];j<self_->b[i+1];++j)
	      s+=self_->x[j-1]*x_[(self_->i[j-1]-1)];
	    y_[i] = s;
	  } }      
    }
  else if (self_->format==0)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {
	    R s;
	    for (s=(R)0.0,j=self_->b[i];j<self_->b[i+1];++j)
	      {
		s+=self_->x[j]*x_[self_->i[j]];
	      } 
	    y_[i] = s;
	  } }
    }
}

void Sparse_assmux( cst_pSparse const	self_,
		    pR 			y_,
		    cst_pR		x_)
{

  if (self_->format==1)
    {
      { I i,j;
	R s;
	for (i=0;i<self_->n;++i)
	  {
	    for (s=((R)0.0),j=self_->b[i];j<self_->b[i+1];++j)
	      s+=self_->x[j-1]*x_[(self_->i[j-1]-1)];
	    y_[i] += s;
	  } }      
    }
  else if (self_->format==0)
    {
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {
	    R s;
	    for (s=(R)0.0,j=self_->b[i];j<self_->b[i+1];++j)
		s+=self_->x[j]*x_[self_->i[j]];
	    y_[i] += s;
	  } }
    }
}


void       	Sparse_extractDiagonal		(cst_pSparse const 	self_,
						 pR 			values_)
{
  if (1 == self_->format)
    {
      for (I rowIndex = 0; rowIndex < self_->n; ++rowIndex)
	{
	  const I bound = self_->b[rowIndex + 1] - 1;
	  for (I at = self_->b[rowIndex] - 1; at < bound; ++at)
            {
	      if (self_->i[at] == rowIndex + 1)
		{
		  values_[rowIndex] = self_->x[at];
		  break;
               }
            }
         }
    }
  else if (0 == self_->format)
    {
      for (I rowIndex = 0; rowIndex < self_->n; ++rowIndex)
	{
	  const I bound = self_->b[rowIndex + 1];
	  for (I at = self_->b[rowIndex]; at < bound; ++at)
            {
	      if (self_->i[at] == rowIndex)
		{
		  values_[rowIndex] = self_->x[at];
		  break;
               }
            }
         }
    }
}

void Sparse_gemv(cst_pR 			a_,
		 cst_pSparse const 		self_,
		 cst_pR 			rhs_,
		 cst_pR 			b_,
		 pR 				sol_)
{

#ifndef NDEBUG
  wrong_param_Sparse(! a_,1,"Sparse_gemv");
  wrong_param_Sparse(! self_,2,"Sparse_gemv");
  wrong_param_Sparse(! rhs_,3,"Sparse_gemv");
  wrong_param_Sparse(! b_,4,"Sparse_gemv");
  wrong_param_Sparse(! sol_,5,"Sparse_gemv");
  wrong_param_Sparse((self_->format!=1)&&(self_->format!=0),2,"Sparse_gemv");
#endif
  if (self_->format==1)
    {
      if (b_[0]==((R)0.0))
	{
	  { I i,j;
	    R s;
	    for (i=0;i<self_->n;++i)
	      {
		for (s=((R)0.0),j=self_->b[i];j<self_->b[i+1];++j)
		  s+=self_->x[j-1]*rhs_[(self_->i[j-1]-1)];
		sol_[i] = a_[0] * s;
	      } }
	}
      else
	{
	  { I i,j;
	    for (i=0;i<self_->n;++i)
	      {
		R s;
		for (s=(R)0.0,j=self_->b[i];j<self_->b[i+1];++j)
		  {
		    s+=self_->x[j-1]*rhs_[(self_->i[j-1]-1)];
		  } 
		sol_[i] = sol_[i] * b_[0] + a_[0] * s;
	      } }
	}
      
    }
  else if (self_->format==0)
    {
      /*      ns_err("pas ici");*/
      { I i,j;
	for (i=0;i<self_->n;++i)
	  {
	    R s;
	    for (s=(R)0.0,j=self_->b[i];j<self_->b[i+1];++j)
	      {
		s+=self_->x[j]*rhs_[self_->i[j]];
	      } 
	    if (b_[0]==((R)0.0))
	      sol_[i] = a_[0] * s;
	    else
	      sol_[i] = sol_[i] * b_[0] + a_[0] * s;
	  } }
    }
}



void Sparse_ass(pSparse const	self_,
		  const I	i_,
		  const I	j_,
		  const R	a)
  
{

  if (self_->format==1)
    {
      { I i;
	for (i=self_->b[i_];i<self_->b[i_+1];++i)
	  if (self_->i[i-1]==j_+1){self_->own_x[i-1]+=a;break;} 
#ifndef NDEBUG
	if (i>=self_->b[i_+1]) {fprintf(stderr,"Sparse_ass(format=1):erreur non trouve "ifmt" "ifmt"\n",i_,j_);exit(1);}
#endif
      }  
    }
  else if (self_->format==0)
    {
      { I i;
	for (i=self_->b[i_];i<self_->b[i_+1];++i)
	  if (self_->i[i]==j_){self_->own_x[i]+=a;break;} 
#ifndef NDEBUG
	if (i>=self_->b[i_+1]) {fprintf(stderr,"Sparse_ass(format=0):erreur non trouve "ifmt" "ifmt"\n",i_,j_);exit(1);}
#endif
      }  
    }
  else
    {
      fprintf(stderr,"Sparse_ass:erreur format\n");
      exit(1);
    }
}


void Sparse_assmatelm(pSparse const 	self_,
		      cst_pI 		n_,
		      cst_pI 		locnumer_,
		      cst_pR 		matelm_,
		      cst_pI 		matelmoff_,
		      pI 		iwork_,
		      pI 		blank_)
{ 


  I N = 0;  

  { I i;
    for (i=0;i<n_[0];++i)
      {

	if (! blank_[locnumer_[i]])
	  {
	    blank_[locnumer_[i]] 	= ++N;
	    iwork_[N-1] 	  	= locnumer_[i];
	  }

      } }

  { I j;
    for (j=0;j<n_[0];++j)
      {	
	const I jddl 			= locnumer_[j];	
	const I n			= self_->b[jddl+1]-self_->b[jddl];
	const I m 			= self_->b[jddl]-self_->format;
	{ cst_pI it 			= &self_->i[m];
	  pR rt				= &self_->own_x[m];
	  { I i;
	    for (i=0;i<n;++i)
	      {
		I k;
		if ( (k=blank_[it[i]-self_->format]) )
		  {
		    rt[i] += matelm_[j + (k-1)*matelmoff_[0]];
		  }
	      } } }	
      } }  
  { I i;
    for (i=0;i<N;++i) 
      {
	blank_[iwork_[i]] = 0;
      } }
}









