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
#include "ExternPardiso.h"
#include "DG_Jacobian.hpp"
ns_mesh * mesh;

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
  //  const I numNodes = ns_mesh_get_numNodes(s_);
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


void dg_print_mesh(const ns_mesh*s_,cst_pI codelm,const char * name_,...)
{
  //  const I numNodes = ns_mesh_get_numNodes(s_);
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
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


void dg_print_mesh(I numNodes,I nelm,pI cnc,I cncoff,pR coo,pI codnodes,I codnodesoff,pI codelm,I codelmoff,const char * name_,...)
{

  { char ctmp[512];
    { char ctmp2[512];
      va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
    I cncelm[3];
    FILE * fil = fopen(ctmp,"w");
    fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",nelm*3);
    { I i;
      for (i=0;i<nelm;++i)
	{
	  for (I k=0;k<3;++k)
	    {
	      cncelm[k] = cnc[cncoff*i+k];
	    }

	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",coo[cncelm[j]*2+0],coo[cncelm[j]*2+1],codnodes[codnodesoff * cncelm[j]]); } }
	} } 
    
    fprintf(fil,"Triangles\n" ifmt "\n",nelm); 
    { I i;
      for (i=0;i<nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,codelm[i*codelmoff]);
	} } 
    fprintf(fil,"End\n");						
    fclose(fil); }  
}


int comp(const void * a ,const void * b)
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
};



static void coordinates2d(long long int				N_, 
			  double				box_[],
			  double* 				vertex_,
			  long long int                         vertexoff_,
			  unsigned long long int *  		cod2_,
			  long long int				codoff_)
{
  static const unsigned long long m = 1LL<<62;
  static const int BitTab[2] 	= {1,2};
  static const int GeoCod[4]	= {1,2,0,3};
  static const int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
  static const double len 	= 4.611686018427387904e+18;
  double box[4],dbl;
  unsigned long long int IntCrd[2],cod;
  int rot[4];
  box[0] 	= box_[0];
  box[1] 	= box_[1];
  box[2] 	= len / (box_[2] - box_[0]);
  box[3] 	= len / (box_[3] - box_[1]);
  double loc[2];
  { int i;
    for(i=0; i<N_; i++)
      {
	loc[0] = vertex_[vertexoff_*i+0];
	loc[1] = vertex_[vertexoff_*i+1];
	/* Convert double precision coordinates to integers */
	dbl 	= (loc[0] - box[0]) * box[0+2];
	IntCrd[0] = (unsigned long long int)dbl;
	dbl 	= (loc[1] - box[1]) * box[1+2];
	IntCrd[1] = (unsigned long long int)dbl;
	/* Binary hilbert renumbering loop */
	cod = 0;
	rot[0] = GeoCod[0];
	rot[1] = GeoCod[1];
	rot[2] = GeoCod[2];
	rot[3] = GeoCod[3];
	{ int b;
	  for(b=0;b<31;b++)
	    {
	      int GeoWrd = 0;

	      if(IntCrd[0] & m)
		GeoWrd |= BitTab[0];
	      IntCrd[0] = IntCrd[0]<<1;

	      if(IntCrd[1] & m)
		GeoWrd |= BitTab[1];
	      IntCrd[1] = IntCrd[1]<<1;

	      const int NewWrd = rot[ GeoWrd ];
	      cod = cod<<2 | NewWrd;
	      rot[0] = HilCod[ NewWrd ][ rot[0] ];
	      rot[1] = HilCod[ NewWrd ][ rot[1] ];
	      rot[2] = HilCod[ NewWrd ][ rot[2] ];
	      rot[3] = HilCod[ NewWrd ][ rot[3] ];
	    } }
	cod2_[codoff_*i+0] 	= cod;
	cod2_[codoff_*i+1] 	= i;
      } }

}

struct hilbert_reordering_t
{
  unsigned long long int cod;
  unsigned long long int id;
};


int hilbert_comp_t(const void * a_,const void * b_)
{
  const I * a = (const I*)a_;
  const I * b = (const I*)b_;
  if (*a < *b)
    {
      return -1;
    }
  else if (*a > *b)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}



struct topology_elem
{
  I m_partid[1];
  I m_cnc[6];
  I m_adj[3];
  I m_gid[1];
  I m_h[1];
};


void hilbert_reordering(I 		meshnelm_,
			topology_elem * meshelms_,
			pI 		perm_,
			cst_pR 		coo_)
{
  double * vertices = (double *)malloc(sizeof(double) * 2*meshnelm_);
  I * tmp = (I*)malloc(sizeof(I) * 2*meshnelm_);
  for (I ielm=0;ielm<meshnelm_;++ielm)
    {
      //      double x = 0.0;
      //      double y = 0.0;

      I i0 = meshelms_[ielm].m_cnc[0]-1;
      I i1 = meshelms_[ielm].m_cnc[1]-1;
      I i2 = meshelms_[ielm].m_cnc[2]-1;

      cst_pR p0 = &coo_[2*i0];
      cst_pR p1 = &coo_[2*i1];
      cst_pR p2 = &coo_[2*i2];
      
      vertices[2*ielm+0]  = (p0[0]+p1[0]+p2[0])/3.0;
      vertices[2*ielm+1]  = (p0[1]+p1[1]+p2[1])/3.0;
    }
  
  //
  // Compute the centroids of the elements.
  //
  
  double bb[4];
  bb[0] = vertices[0];
  bb[1] = vertices[1];
  bb[2] = vertices[0];
  bb[3] = vertices[1];
  for (I ielm=1;ielm<meshnelm_;++ielm)
    {
      
      bb[0] = std::min(bb[0],vertices[2*ielm+0]);
      bb[1] = std::min(bb[1],vertices[2*ielm+1]);
      bb[2] = std::max(bb[2],vertices[2*ielm+0]);
      bb[3] = std::max(bb[3],vertices[2*ielm+1]);      
    }
  

  bb[0] = 0;
  bb[1] = 0;
  bb[2] = 1;
  bb[3] = 2;

  
  coordinates2d(meshnelm_, 
		bb,
		vertices,
		2,
		(unsigned long long int*)tmp,
		2);
  
  qsort(tmp,meshnelm_,2*sizeof(I),hilbert_comp_t);
  
#if 0
  
  for (I ielm=0;ielm<mesh_->nelm;++ielm)
    {
      // perm_[tmp[2*ielm+1]] = ielm;
      perm_[ielm] = tmp[2*ielm+1];
    }

  FILE * f = fopen("curve.txt","w");
  for (I i=0;i<mesh_->nelm;++i)
    {
      fprintf(f," %e %e\n",vertices[2*perm_[i]+0],vertices[2*perm_[i]+1]);
    }
  fclose(f);
  
#endif  
  for (I ielm=0;ielm<meshnelm_;++ielm)
    {
      // perm_[tmp[2*ielm+1]] = ielm;
      perm_[tmp[2*ielm+1]] = ielm;
    }  
}



I calculate_nghosts(I 			meshnelm,
		    topology_elem * 	meshelms,			    
		    I 			ipart,
		    pI 			marker,
		    pI 			select_n,
		    pI                  select)
{
  select_n[0] = 0;
  I nghost = 0;
  for (I ielm=0;ielm<meshnelm;++ielm)
    {
      if (meshelms[ielm].m_partid[0] == ipart + 1)//partid[ielm * partid_off] 
	{
	  for (I k=0;k<3;++k)
	    {
	      I nei = meshelms[ielm].m_adj[k];
	      if (nei > 0)
		{
		  I neipart = meshelms[nei-1].m_partid[0]; // partid[(nei-1) * partid_off];
		  if (neipart != ipart + 1)
		    {
		      if (0 == marker[nei-1])
			{
			  select[select_n[0]] = nei-1;
			  marker[nei-1] = ++select_n[0];			  
			  ++nghost;
			  //  fprintf(stdout,"nghost "ifmt ", part " ifmt"\n",nghost,ipart);
			}		  
		    }
		}
	    }
	}
    }
  
  for (I i=0;i<select_n[0];++i)
    {
      marker[select[i]] = 0;
    }
  
  return nghost;
}

struct partition_t
{
  
public: I 			m_nnodes;
public: I  			m_nelm;
public: I  			m_nelm_interior_ghost;
public: I  			m_nelm_exterior_ghost;
public: topology_elem * 	m_elms;
public: pR 			m_coo;
public: pI 			m_gidnode;
  
public: void print_mesh(const char * 	name_,...)
  {
    //    const I numNodes = m_nnodes;
    { char ctmp[512];
      { char ctmp2[512];
	va_list args;
      va_start (args,name_);
      vsprintf(ctmp2,name_,args);
      va_end(args);
      sprintf(ctmp,"%s.mesh",ctmp2); }    
      //      I cncelm[3];
      FILE * fil = fopen(ctmp,"w");
      fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n" ifmt "\n",this->m_nelm*3);
      { I i;
      for (i=0;i<this->m_nelm;++i)
	{
	  {I j;for (j=0;j<3;++j){ fprintf(fil,"" rfmt " " rfmt " " ifmt "\n",
					  m_coo[(m_elms[i].m_cnc[j]-1)*2+0],
					  m_coo[(m_elms[i].m_cnc[j]-1)*2+1],((I)0));}}
	} } 
      
      fprintf(fil,"Triangles\n" ifmt "\n",m_nelm); 
      { I i;
      for (i=0;i<m_nelm;++i)
	{
	  fprintf(fil,"" ifmt " " ifmt " " ifmt " " ifmt "\n",i*3+1,i*3+2,i*3+3,m_elms[i].m_partid[0]);
	} } 
      fprintf(fil,"End\n");						
      fclose(fil); }  
    
  };
  
public: void define(I 			ipart,
		    I 			meshnelm,
		    topology_elem * 	meshelms,
		    cst_pI		begin_part,
		    pR 			coo)
  {    
    //    pI marker   	= (pI)calloc(meshnelm,sizeof(I));
    pI select   	= (pI)calloc(meshnelm,sizeof(I));
    // I  select_n 	= 0;
    m_coo = coo;
    //
    // Calculate the size of the ghosts.
    //
#if 0    
    this->m_nelm_exterior_ghost	= calculate_nghosts(meshnelm,
						    meshelms,
						    ipart,
						    marker,
						    &this->m_nelm_exterior_ghost,
						    select);
#endif
    this->m_nelm_exterior_ghost = 0;
    fprintf(stdout,"this->m_nelm_exterior_ghost " ifmt "\n",this->m_nelm_exterior_ghost);
    
    //
    // Allocate.
    //
    I nelm_partition 	= begin_part[ipart+1] - begin_part[ipart];    
    I nelm 		= nelm_partition + this->m_nelm_exterior_ghost;   
    this->m_nelm        = nelm;
    this->m_elms 	= (topology_elem *)malloc(sizeof(topology_elem)*nelm);
    //
    // Copy the partition
    //
    memcpy(this->m_elms,
	   &meshelms[begin_part[ipart]],
	   sizeof(topology_elem)*nelm_partition);


    //
    // Sort the ghost per partitions.
    //
#if 1    
    //
    // Copy the ghost.
    // 
    for (I i=0;i<this->m_nelm_exterior_ghost;++i)
      {
	I ielm = select[i];
	memcpy(&this->m_elms[nelm_partition + i],
	       &meshelms[ielm],
	       sizeof(topology_elem));       
      }
    
#endif
    if (this->m_nelm_exterior_ghost>0)
      {
	qsort(&m_elms[nelm_partition],
	  this->m_nelm_exterior_ghost,
	  sizeof(topology_elem),
	  [](const void * a,const void * b)
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
   //
   //
   // we have renumeroted the ghost.
   //
   //
  };
  
};




int comp_part(const void * a ,const void * b)
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
      //
      // inside the partitions.
      //
      
      
      

      
      return 0;
    }
}


void build_partid(I meshnelm,
		  topology_elem*meshelms,
		   I 		nparts,
		   pI 		partid,
		  I 		partid_off,
		  cst_pR coo,
		  pI begin_part)
 {
   //
   // Space filling curve reordering of the elements.
   //
   pI perm = (pI)malloc(sizeof(I)*meshnelm);
   hilbert_reordering(meshnelm,meshelms,perm,coo);
   
   I nelm_per_partitions = meshnelm / nparts;
   //   pI part_begin = (pI)calloc(nparts+1,sizeof(I));
   begin_part[0] = 0;
   for (I ipart=1;ipart<nparts;++ipart)
     {
       begin_part[ipart] = begin_part[ipart-1]  + nelm_per_partitions;
     }
   begin_part[nparts] = begin_part[nparts-1]  + nelm_per_partitions + (meshnelm % nparts);
   
   for (I ielm=0;ielm<meshnelm;++ielm)
     {
       I jelm = perm[ielm];
       I k = 1;
       for (;k<=nparts;++k)
	 {
	   if (jelm < begin_part[k])
	     {
	       partid[ielm] = k;
	       break;
	     }
	 }
     }
   
   //
   //
   //
   //   topology_elem * elm = (topology_elem *)malloc(sizeof(topology_elem)*meshnelm);
   for (I ielm=0;ielm<meshnelm;++ielm)
     {
//       for (I i=0;i<6;++i)
//	 elm[ielm].m_cnc[i] = mesh->cnc[6*ielm+i];
//       for (I i=0;i<3;++i)
//	 elm[ielm].m_adj[i] = mesh->adj[3*ielm+i];
       meshelms[ielm].m_gid[0] = ielm;
       meshelms[ielm].m_partid[0] = partid[ielm];
       meshelms[ielm].m_h[0] = perm[ielm];
     }
   
   qsort(meshelms,meshnelm,sizeof(topology_elem),
	 [](const void * a,const void * b)
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
	       
	       if (a_[11] < b_[11] )
		 {
		   return -1;
		 }
	       else if (a_[11] > b_[11] )
		 {
		   return 1;
		 }
	       else
		 {
		   return 0;
		 }
	     }
	  
	 });

   free(perm);
   {
     perm = (pI)malloc(sizeof(I)*meshnelm);
     for (I ielm=0;ielm<meshnelm;++ielm)
       {
	 perm[meshelms[ielm].m_gid[0]] = ielm;
       }
     
     for (I ielm=0;ielm<meshnelm;++ielm)
       {
	 //      for (I i=0;i<6;++i)
	 //	mesh->own_cnc[6*ielm+i] = elm[ielm].m_cnc[i];
	 for (I i=0;i<3;++i)	
	   {
	     meshelms[ielm].m_adj[i] = (meshelms[ielm].m_adj[i]==0) ? 0 : perm[ meshelms[ielm].m_adj[i] - 1 ] + 1;
	   }
	 partid[ielm] = meshelms[ielm].m_partid[0];
       }

     
     for (I ielm=0;ielm<meshnelm;++ielm)
       {
	 for (I i=0;i<3;++i)	
	   {
	     I nei = meshelms[ielm].m_adj[i];// = (meshelms[ielm].m_adj[i]==0) ? 0 : perm[ meshelms[ielm].m_adj[i] - 1 ] + 1;
	     if (nei > 0)
	       {
		 bool found = false;
		 for (I j=0;j<3;++j)
		   {
		     if (ielm + 1 == meshelms[nei-1].m_adj[j])
		       {
			 found = true;
			 break;
		       }
		   }
		 if (!found)
		   {
		     fprintf(stderr,"error adjacence ielm = " ifmt "\n",ielm);
		     exit(1);
		   }
	       }
	   }
       }

     
     free(perm);
   }
   
 }





#if 0
void partition(ns_mesh * 	mesh,
	       topology_elem *  elms,
	       pI 		part_begin,
	       I 		ipart,
	       cst_pI 		part_nelm,
	       cst_pI 		partid,
	       I 		partid_off)
{
   pI marker   	= (pI)calloc(mesh->nelm,sizeof(I));
   pI select   	= (pI)calloc(mesh->nelm,sizeof(I));
   I  select_n 	= 0;

   //
   // Calculate the size of the ghosts.
   //
   I nghost 	= calculate_nghosts(mesh,
				    ipart,
				    partid,
				    partid_off,
				    marker,
				    &select_n,
				    select);

   //
   // Allocate.
   //
   I nelm = part_nelm[ipart] + nghost;   
   topology_elem * pelm = (topology_elem *)malloc(sizeof(topology_elem)*nelm);

   //
   // Copy the partition
   //
   memcpy(pelm,&elm[part_begin[ipart]],sizeof(topology_elem)*(part_begin[ipart+1] - part_begin[ipart]));   

   //
   // Copy the ghost.
   // 
   for (I i=0;i<nghost;++i)
     {
       I ielm = select[i];
       memcpy(&pelm[part_nelm[ipart] + i],&elms[ielm],sizeof(topology_elem));       
     }

}
#endif




 
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

  I 	npartitions 	= 5;
  pI 	partid 		= (pI)malloc(sizeof(I)*mesh->nelm);
  I 	partid_off	= 1;


  I meshnelm = mesh->nelm;
  topology_elem * meshelms = (topology_elem *)calloc(mesh->nelm,sizeof(topology_elem));
  for (I ielm=0;ielm<mesh->nelm;++ielm)
    {
      for (I i=0;i<6;++i)
	meshelms[ielm].m_cnc[i] = mesh->cnc[6*ielm+i];
      for (I i=0;i<3;++i)
	meshelms[ielm].m_adj[i] = mesh->adj[3*ielm+i];
    }
  
  //
  // Build partid.
  //
  pI begin_part = (pI)malloc(sizeof(I)*(npartitions+1));
  build_partid(meshnelm,
	       meshelms,
	       npartitions,
	       partid,
	       partid_off,
	       mesh->coo,
	       begin_part);
  
  partition_t * partitions = (partition_t * )calloc(npartitions,sizeof(partition_t));
  for (I ipart=0;ipart < npartitions;++ipart)
    {
      partitions[ipart].define(ipart,
			       meshnelm,
			       meshelms,
			       begin_part,
			       mesh->own_coo);
      
      partitions[ipart].print_mesh("dg." ifmt,ipart);

#if 1
      //      pI 	graph;
      //      I 	graph_off	= 3;
      I 	nbasis 		= 3;
      DG_Jacobian jacobian(partitions[ipart].m_nelm,
			   3,
			   nbasis,
			   nbasis,
			   mesh->adj,
			   3);
#endif
    }
  return 0;
  

}



