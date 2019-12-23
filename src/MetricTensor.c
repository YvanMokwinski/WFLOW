#include "MetricTensor.h"
#include "Blas.h"
void nsRHDLE_def(nsRHDLE_ST*hdle_,const I n_,const I m_,pR x_,const I off_)
{
  hdle_->n=n_;
  hdle_->m=m_;
  hdle_->x=x_;
  hdle_->off=off_;
}

/*#include "libns.h"*/
static const I negal1=(I)1;
#if __ns_debug__

#define wrong_param_MetricTensor(_cond,_i,_s) if ( (_cond) ) { fprintf(stderr,"routine '"#_s"'\nline : '%d'\nfile : '%s'\nwrong parameter %d\n",__LINE__,__FILE__,-(_i)); exit(1);} ((void)0)

cst_pR 	cst_MetricTensor_x	(const MetricTensor*metric_) 	{ wrong_param_MetricTensor(NOT metric_,1,cst_MetricTensor_x);	return metric_->x;}
pR 		MetricTensor_x	(MetricTensor*metric_) 		{ wrong_param_MetricTensor(NOT metric_,1,MetricTensor_x);	return metric_->x;}
I 		MetricTensor_n	(const MetricTensor*metric_) 	{ wrong_param_MetricTensor(NOT metric_,1,MetricTensor_n);	return metric_->n;}
const ns_mesh*	MetricTensor_mesh	(const MetricTensor*metric_) 	{ wrong_param_MetricTensor(NOT metric_,1,MetricTensor_mesh); 	return metric_->mesh;}

#else
#define wrong_param_MetricTensor(_cond,_i,_s) ((void)0)
#endif


void MetricTensor_compute_conformite(const MetricTensor * 	M_,
				 const MetricTensor * 	Q_,
				 R 		mmmv_[4],
				 const char * 		filename_)
{
  const ns_mesh * mesh 	= MetricTensor_mesh(M_);
#if __ns_debug__
  const ns_mesh * mesh_Q= MetricTensor_mesh(Q_);
  if (mesh_Q!=mesh)
    ns_err("MetricTensor_compute_conformite failed");
#endif
  const I numNodes 	= ns_mesh_get_numNodes(mesh);
  R M[__eDim_ALL*__eDim_ALL];
  R Q[__eDim_ALL*__eDim_ALL];  
  mmmv_[0]=(R)0.0;
  mmmv_[1]=(R)1.0e+20;
  mmmv_[2]=(R)0.0;
  mmmv_[3]=(R)0.0;
  char basename[256];
  //  char * basename = NULL;
  FILE*f = NULL;
  if (filename_)
    {
      char ctmp[256];
      ns_basename(filename_,basename);
      ns_mesh_print(mesh,basename);
      sprintf(ctmp,"%s.bb",basename);
      f=fopen(ctmp,"w");
      fprintf(f,"2 1 "ifmt" 2\n",numNodes);
    }

  { I i;
    for (i=(I)0;i<numNodes;++i)
      {     
	M[0] = M_->x[3*i+0];
	M[1] = M_->x[3*i+1];
	M[2] = M_->x[3*i+2];
	Q[0] = Q_->x[3*i+0];
	Q[1] = Q_->x[3*i+1];
	Q[2] = Q_->x[3*i+2];
	const R iM[3]  = {M[2]/(M[0]*M[2]-M[1]*M[1]),-M[1]/(M[0]*M[2]-M[1]*M[1]),M[0]/(M[0]*M[2]-M[1]*M[1])};
	const R iQ[3]  = {Q[2]/(Q[0]*Q[2]-Q[1]*Q[1]),-Q[1]/(Q[0]*Q[2]-Q[1]*Q[1]),Q[0]/(Q[0]*Q[2]-Q[1]*Q[1])};
	const R iMQ[3] = {iM[0]*Q[0]+iM[1]*Q[1],iM[0]*Q[1]+iM[1]*Q[2],iM[1]*Q[1]+iM[2]*Q[2]};
	const R iQM[3] = {iQ[0]*M[0]+iQ[1]*M[1],iQ[0]*M[1]+iQ[1]*M[2],iQ[1]*M[1]+iQ[2]*M[2]};
	const R T[3]   = {((R)2.0)-iMQ[0]-iQM[0],-iMQ[1]-iQM[1],((R)2.0)-iMQ[2]-iQM[2]};
	const R c1     = nsSQRT( T[0]*T[0] + ((R)2.0)*T[1]*T[1] + T[2]*T[2]);
	R c=c1;
	mmmv_[0]+=c;
	if (c<mmmv_[1])
	  mmmv_[1]=c;
	if (c>mmmv_[2])
	  mmmv_[2]=c;
	if (f)
	  {
	    fprintf(f,"%e\n",c1);
	  }
      } }
  if (f)
    {
      fclose(f);
    }
  mmmv_[0]/=((R)numNodes);
  { I i;
    for (i=(I)0;i<numNodes;++i)
      {     
	M[0] = M_->x[3*i+0];
	M[1] = M_->x[3*i+1];
	M[2] = M_->x[3*i+2];
	Q[0] = Q_->x[3*i+0];
	Q[1] = Q_->x[3*i+1];
	Q[2] = Q_->x[3*i+2];
	const R iM[3]  = {M[2]/(M[0]*M[2]-M[1]*M[1]),-M[1]/(M[0]*M[2]-M[1]*M[1]),M[0]/(M[0]*M[2]-M[1]*M[1])};
	const R iQ[3]  = {Q[2]/(Q[0]*Q[2]-Q[1]*Q[1]),-Q[1]/(Q[0]*Q[2]-Q[1]*Q[1]),Q[0]/(Q[0]*Q[2]-Q[1]*Q[1])};
	const R iMQ[3] = {iM[0]*Q[0]+iM[1]*Q[1],iM[0]*Q[1]+iM[1]*Q[2],iM[1]*Q[1]+iM[2]*Q[2]};
	const R iQM[3] = {iQ[0]*M[0]+iQ[1]*M[1],iQ[0]*M[1]+iQ[1]*M[2],iQ[1]*M[1]+iQ[2]*M[2]};
	const R T[3]   = {2.0-iMQ[0]-iQM[0],-iMQ[1]-iQM[1],2.0-iMQ[2]-iQM[2]};
	const R c      = nsSQRT( T[0]*T[0] + 2.0*T[1]*T[1] + T[2]*T[2]);
#if 0
	const R c=nsSQRT((M[0]-Q[0])*(M[0]-Q[0]) + 2.0*(M[1]-Q[1])*(M[1]-Q[1]) + (M[2]-Q[2])*(M[2]-Q[2]));
#endif
	mmmv_[3]+=(mmmv_[0]-c)*(mmmv_[0]-c);
      } }
  mmmv_[3]/=((R)numNodes);
  return ;
}

void MetricTensor_from_size(MetricTensor*	metric_,
			const R 	h_)
{
  const R v = ((R)1.0)/(h_*h_);
  { I i;
    for (i=0;i<metric_->n;++i)
      {
	metric_->x[3*i+0] = v;
	metric_->x[3*i+1] = ((R)0.0);
	metric_->x[3*i+2] = v;
      } }
}


void nsOP_METRIC_INFO_def(nsOP_METRIC_INFO 	op_,
			  cst_pR 		hmin_,
			  cst_pR		hmax_,
			  cst_pR 		itol_)
{
  memset(op_,0,sizeof(nsOP_METRIC_INFO_ST));
  op_->hmin = hmin_;
  op_->hmax = hmax_;
  op_->itol = itol_;
}			  

MetricTensor*ns_mesh_killmetric(ns_mesh*mesh_,MetricTensor*metric_)
{
  if (metric_)
    {
      { I i;
	for (i=0;i<64;++i)
	  {
	    if (mesh_->metrics[i]==metric_)
	      break;
	  } 
	if (i<64)
	  {
	    mesh_->metrics[i] = MetricTensor_kill(mesh_->metrics[i]);
	    return NULL;
	  }
	else
	  {
	    fprintf(stderr,"MetricTensor_killmetric:metric not found (adress=%p)\n",metric_);
	    exit(1);
	    return NULL;
	  } }
    }
  return NULL;
}



pMetricTensor ns_mesh_addmetric(ns_mesh * mesh_)
{
  { I i;
    for (i=0;i<64;++i)
      {
	if (NULL == mesh_->metrics[i])
	  break;
      } 
    if (i<64)
      {
	mesh_->metrics[i] = MetricTensor_new(mesh_);
	return mesh_->metrics[i];
      }
    else
      {
	fprintf(stderr,"MetricTensor_addmetric:maximum metric is reached\n");
	return NULL;
      }  }
}


void MetricTensor_def(pMetricTensor	metric_,
		      const ns_mesh * 	mesh_)
{
  memset(metric_,0,
	 sizeof(MetricTensor_def));  
  metric_->dim			= ns_mesh_dimension(mesh_);
  const I m = ((metric_->dim+1)*metric_->dim)/2;
  metric_->mesh 		= mesh_;
  metric_->n 			= ns_mesh_get_numNodes(mesh_);
#if 0
  printf(" "ifmt" "ifmt"\n",m,metric_->n);
  exit(1);
#endif
  metric_->x = calloc(m*metric_->n,sizeof(R));
  /*  metric_->x 			= ns_memory_rcalloc( m *mesh_->nvertex,0);*/
  nsRHDLE_def(&metric_->data,m,metric_->n,metric_->x,m);
}


MetricTensor *MetricTensor_new(const ns_mesh*mesh_)
{
  MetricTensor* m = malloc(sizeof(MetricTensor));
  MetricTensor_def(m,mesh_);
  return m;
}

void MetricTensor_eig(const R 	T[3],
		      R  	ee_[2],
		      R 	ev_[4])
{ 
  
  R m11 = T[0];
  R m12 = T[1];
  R m22 = T[2];
  pR  lambda1 	= ee_; 
  pR  lambda2	= &ee_[1]; 
  pR  x1	= ev_; 
  pR  y1	= &ev_[1];
  pR  x2	= &ev_[2]; 
  pR  y2	= &ev_[3];
  int NbRotations[1];   
#define MAXITER 50  // Nombre maximum d'iterations
  // la matrice mat est initialisee a une matrice reelle symetrique
  // 2 x 2. A la fin, les elements au-dessus de la diagonale sont detruits.
  R mat[2][2] = { {m11, m12},
		       {m12, m22} };
  
  // Le vecteur lambda contiendra les deux valeurs propres.
  // Il est initialise a la diagonale de mat
  R lambda[2] = { m11, m22 };
  
  // La matrice eigen contiendra les deux vecteurs propres par colonne.
  // Elle est initialisee a l'identite
  R eigen[2][2] = { {1.0, 0.0},
			 {0.0, 1.0} };
  
   // Les vecteurs b et z sont internes a cette methode. Ils accumulent
   // les termes tapq de l'equation 11.1.14....
  R zloc[2] = {0.0, 0.0};
  R bloc[2] = {m11, m22};
  
  // Initialisation du nombre de rotations de Jacobi
  NbRotations[0] = 0;
  
  // Boucle principale
  int    iter;
  R tresh, g, h, s, tt, theta, tau, tmp;
  for ( iter=1; iter<=MAXITER; iter++ )
    {
      // erreur contient la somme des termes hors diagonaux.
      // erreur est le critere de convergence. 
      R erreur = nsFABS( mat[0][1] );
      if ( erreur == 0.0 ) break;

      if ( iter < 4 ) 
      {
         tresh = 0.05 * erreur;
      }
      else
      {
         tresh = 0.0;
      }

      g = ((R)100.0)*nsFABS( mat[0][1] );
      if (    ( iter > 4 )
           && ( nsFABS( lambda[0] ) + g == nsFABS( lambda[0] ) )
           && ( nsFABS( lambda[1] ) + g == nsFABS( lambda[1] ) ) )
      { 
         mat[0][1] = 0.0;
      }
      else if ( nsFABS( mat[0][1] ) > tresh )
      {
         h = lambda[1] - lambda[0];
         if ( nsFABS( h ) + g == nsFABS( h ) )
         {
            tt =  mat[0][1] / h;
         }
         else
         {
            theta = 0.5 * h / mat[0][1];
            tt     = 1.0 / ( nsFABS( theta ) + nsSQRT( 1.0 + theta*theta ) );
            if ( theta < 0.0 ) 
            { 
               tt = -tt;
            }
         }
         tmp  = 1.0 / nsSQRT( 1.0 + tt*tt );
         s    = tt * tmp;
         tau  = s / ( 1.0 + tmp );
         h    = tt * mat[0][1];
         zloc[0]   -= h;
         zloc[1]   += h;
         lambda[0] -= h;
         lambda[1] += h;
         mat[0][1]  = 0.0;
         
         // ROTATE(eigen,k,0,k,1) with k=0;
         g = eigen[0][0];
         h = eigen[0][1];
         eigen[0][0] = g - s*(h + g*tau);
         eigen[0][1] = h + s*(g - h*tau);

         // ROTATE(eigen,k,0,k,1) with k=1;
         g = eigen[1][0];
         h = eigen[1][1];
         eigen[1][0] = g - s*(h + g*tau);
         eigen[1][1] = h + s*(g - h*tau);
         
         NbRotations[0]+=1;
      }
      
      // Update lambda avec la somme dans z et reinitialise z
      bloc[0]  += zloc[0];
      bloc[1]  += zloc[1];
      lambda[0] = bloc[0];
      lambda[1] = bloc[1];
      zloc[0]   = 0.0;
      zloc[1]   = 0.0;
   }

   // Output les resultats.
   lambda1[0] = lambda[0];
   lambda2[0] = lambda[1];
   x1[0] = eigen[0][0];
   y1[0] = eigen[1][0];
   x2[0] = eigen[0][1];
   y2[0] = eigen[1][1];
   
   // Si on est rendu ici, c'est qu'on est sorti de la boucle, soit
   // par le break si on a converge, soit parce qu'on a atteint le nombre
   // maximum d'ierations.
#if 0
   if ( iter <= MAXITER ) return 0;
   else                   return iter;
#endif
}

void sdp_intersect(cst_pR 				A_,
		   cst_pR 				B_,
		   pR 					C_)
{
  R A_ee[4*4],A_ev[4];
  R B_ee[4*4],B_ev[4];

#define __mult(_t1,_t2,_t)			\
  _t[0] = _t1[0] * _t2[0]  + _t1[2] * _t2[1];	\
  _t[1] = _t1[1] * _t2[0]  + _t1[3] * _t2[1];	\
  _t[2] = _t1[0] * _t2[2]  + _t1[2] * _t2[3];	\
  _t[3] = _t1[1] * _t2[2]  + _t1[3] * _t2[3]



	R N1[4],N2[4],N3[4],N4[4],N5[4],M3[4],iN1[4],N3t[4],N5t[4];
	
	MetricTensor_eig(A_,A_ee,A_ev);
	MetricTensor_eig(B_,B_ee,B_ev);	      	      

	
	N1[0] = A_ev[0]*nsSQRT(A_ee[0]);
	N1[2] = A_ev[1]*nsSQRT(A_ee[0]);
	N1[1] = A_ev[2]*nsSQRT(A_ee[1]);
	N1[3] = A_ev[3]*nsSQRT(A_ee[1]);
	
	N2[0] = B_ev[0]*nsSQRT(B_ee[0]);
	N2[2] = B_ev[1]*nsSQRT(B_ee[0]);
	N2[1] = B_ev[2]*nsSQRT(B_ee[1]);
	N2[3] = B_ev[3]*nsSQRT(B_ee[1]);
	
	
	/*
	  a b  d  -b 
	  c d  -c a 
	  
	  ad-bc 
	*/ 
	iN1[0] = N1[3]/(N1[0]*N1[3]-N1[1]*N1[2]);
	iN1[1] = -N1[1]/(N1[0]*N1[3]-N1[1]*N1[2]);
	iN1[2] = -N1[2]/(N1[0]*N1[3]-N1[1]*N1[2]);
	iN1[3] = N1[0]/(N1[0]*N1[3]-N1[1]*N1[2]);
	__mult(N2,iN1,N3);  
#if 0
	printf("N2 %8.5e %8.5e\n   %8.5e %8.5e\n\n",N2[0],N2[2],N2[1],N2[3]);
	printf("N1 %8.5e %8.5e\n   %8.5e %8.5e\n\n",N1[0],N1[2],N1[1],N1[3]);
	printf("T1 %8.5e %8.5e\n   %8.5e %8.5e\n\n",N3[0],N3[2],N3[1],N3[3]);
#endif
	N3t[0]=N3[0];
	N3t[1]=N3[2];
	N3t[2]=N3[1];
	N3t[3]=N3[3];
	
	
	
	__mult(N3t,N3,M3);  
	M3[1] = 0.5 *(M3[1]+M3[2]);
	M3[2]=M3[3];

	R T5[4],T5_ee[4],T5_ev[4];
	T5[0] = M3[0];
	T5[1] = M3[1];
	T5[2] = M3[3];
	MetricTensor_eig(T5,T5_ee,T5_ev);	      	      
	if (T5_ee[0]<1.0)
	  T5_ee[0]=(R)1.0;
	if (T5_ee[1]<1.0)
	  T5_ee[1]=(R)1.0;
	
	N4[0] = T5_ev[0]*nsSQRT(T5_ee[0]);
	N4[2] = T5_ev[1]*nsSQRT(T5_ee[0]);
	N4[1] = T5_ev[2]*nsSQRT(T5_ee[1]);
	N4[3] = T5_ev[3]*nsSQRT(T5_ee[1]);
	
	__mult(N4,N1,N5);  
	N5t[0]=N5[0];
	N5t[1]=N5[2];
	N5t[2]=N5[1];
	N5t[3]=N5[3];
	R C[4];
	__mult(N5t,N5,C);
	C_[0]=C[0];
	C_[1]=(C[1]+C[2])*((R)0.5);
	C_[2]=C[3];

}



void MetricTensor_setidentity(MetricTensor * metric_)
{
#if __ns_debug__
  wrong_param_MetricTensor(NOT metric_,1,MetricTensor_setidentity);
#endif
  I i;
  for (i=0;i<metric_->n;++i) metric_->x[3*i+0]=((R)1.0);
  for (i=0;i<metric_->n;++i) metric_->x[3*i+1]=((R)0.0);
  for (i=0;i<metric_->n;++i) metric_->x[3*i+2]=((R)1.0);
}


Err MetricTensor_print(const MetricTensor * 	metric_,
		       const char * 		name_)
{ 
#if __ns_debug__
  wrong_param_MetricTensor(NOT metric_,1,MetricTensor_print);
  wrong_param_MetricTensor(NOT name_,2,MetricTensor_print);
#endif
  char print_ctmp[512];
  sprintf(print_ctmp,"%s.bb",name_);
  FILE * fil = fopen(print_ctmp,"w");
  if (NOT fil)
    return __eErr_file;
  fprintf(fil,"2 3 "ifmt" 2\n",metric_->n);
  { I i;
    for (i=0;i<metric_->n;++i)
      { 
	fprintf(fil,""rfmt" "rfmt" "rfmt"\n",metric_->x[3*i+0],metric_->x[3*i+1],metric_->x[3*i+2]); 
      } }
  fclose(fil); 
  return __eErr_no;
}




R MetricTensor_density(const ns_mesh * 	mesh_,
			 const MetricTensor * 	metric_)
{
  I cncelm[3];
  static const R untiers = ((R)1.0)/((R)3.0);
  R d_metric = ((R)0.0);
  { I ielm;
    for (ielm=0;ielm<mesh_->nelm;++ielm)
      {
	
	ns_mesh_cncelm(mesh_,
			&ielm,
			cncelm);
	
	d_metric+= mesh_->jacelm[ielm]*untiers*(nsSQRT(nsFABS(metric_->x[3*cncelm[0]+0]*metric_->x[3*cncelm[0]+2]-metric_->x[3*cncelm[0]+1]*metric_->x[3*cncelm[0]+1]))
					     +nsSQRT(nsFABS(metric_->x[3*cncelm[1]+0]*metric_->x[3*cncelm[1]+2]-metric_->x[3*cncelm[1]+1]*metric_->x[3*cncelm[1]+1]))
					     +nsSQRT(nsFABS(metric_->x[3*cncelm[2]+0]*metric_->x[3*cncelm[2]+2]-metric_->x[3*cncelm[2]+1]*metric_->x[3*cncelm[2]+1])) );
      } }
  return d_metric;
}

#if 0
R MetricTensor_densityold(const ns_mesh * 	mesh_,
			 const MetricTensor * 	metric_)
{
  I cncelm[3];
  static const R untiers = ((R)1.0)/((R)3.0);
  R d_metric = ((R)0.0);
  { I i;
    for (i=0;i<mesh_->nelm;++i)
      {
	ns_mesh_cncelm0(mesh_,&i,cncelm);
	d_metric+= mesh_->jacelm[i]*untiers*(nsPOW(nsFABS(metric_->x[3*cncelm[0]+0]*metric_->x[3*cncelm[0]+2]-metric_->x[3*cncelm[0]+1]*metric_->x[3*cncelm[0]+1]),untiers)
				      +nsPOW(nsFABS(metric_->x[3*cncelm[1]+0]*metric_->x[3*cncelm[1]+2]-metric_->x[3*cncelm[1]+1]*metric_->x[3*cncelm[1]+1]),untiers)
				      +nsPOW(nsFABS(metric_->x[3*cncelm[2]+0]*metric_->x[3*cncelm[2]+2]-metric_->x[3*cncelm[2]+1]*metric_->x[3*cncelm[2]+1]),untiers) );
      } }
  return d_metric;
}

#endif


void MetricTensor_copy(MetricTensor * metric_,
		    const MetricTensor*metric2_)
{
#if __ns_debug__
  wrong_param_MetricTensor(NOT metric_,1,MetricTensor_print);
  wrong_param_MetricTensor(NOT metric2_,2,MetricTensor_copy);
#endif
  { const I N = metric_->n*3;
    Blas_dcopy(&N,metric2_->x,&negal1,metric_->x,&negal1); }
}


void MetricTensor_free(MetricTensor*metric_)
{
  metric_->n 	= ((I)0);
  free(metric_->x);
  //  ns_memory_rfree(metric_->x,0);
}

MetricTensor * MetricTensor_kill(MetricTensor*metric_)
{
  if (metric_)
    {
      MetricTensor_free(metric_);
      memset(metric_,0,sizeof(MetricTensor));
      free(metric_);
    }
  return NULL;
}


void ns_vector2metric(MetricTensor * 	metric_,
		      const ns_var * 	vect_,
		      cst_pR 		itol_)
{
  
  R a[2],b[3];
  { I i;
    for (i=0;i<metric_->n;++i)
      {
	a[0] = vect_->x[i]*itol_[0];
	a[1] = vect_->y[i]*itol_[0];
	b[0] = a[0]*a[0];
	b[1] = a[1]*a[0];
	b[2] = a[1]*a[1];
	metric_->x[3*i+0]=b[0];
	metric_->x[3*i+1]=b[1];
	metric_->x[3*i+2]=b[2];
      } }
}


void MetricTensor_hconstraint(MetricTensor * 	metric_,cst_pR vmin_,cst_pR vmax_)
{
  R met2[3],met2_ee[2],met2_ev[4];
  { I i;
    const R vmin = vmin_[0];
    const R vmax = vmax_[0];
    for (i=0;i<metric_->n;++i)
      {
	met2[0] = metric_->x[3*i+0];
	met2[1] = metric_->x[3*i+1];
	met2[2] = metric_->x[3*i+2];	
	MetricTensor_eig(met2,met2_ee,met2_ev);
	met2_ee[0] = nsFABS(met2_ee[0]);
	met2_ee[1] = nsFABS(met2_ee[1]);
	if ( met2_ee[0] < vmin )
	  met2_ee[0] = vmin;
	if ( met2_ee[0] > vmax )
	  met2_ee[0] = vmax;
	if ( met2_ee[1] < vmin )
	  met2_ee[1] = vmin;
	if ( met2_ee[1] > vmax )
	  met2_ee[1] = vmax;	
	/*
	  e0a e1c   a  b
	  e0b e1d   c  d
	  
	  e0aa e1cc  e0ba e1cd
	  e0bb e1dd
	*/
	met2[0] = met2_ee[0] * met2_ev[0]*met2_ev[0] + met2_ee[1] * met2_ev[2]*met2_ev[2];
	met2[1] = met2_ee[0] * met2_ev[0]*met2_ev[1] + met2_ee[1] * met2_ev[2]*met2_ev[3];
	met2[2] = met2_ee[0] * met2_ev[1]*met2_ev[1] + met2_ee[1] * met2_ev[3]*met2_ev[3];		
	metric_->x[3*i+0] = met2[0];
	metric_->x[3*i+1] = met2[1];
	metric_->x[3*i+2] = met2[2];	
      } }  
}





void ns_tensor2metric_intersection(MetricTensor * 		metric_,
				   const ns_var * 	hessian_,
				   cst_pR 		itol_,
				   cst_pR 		vmin_,
				   cst_pR 		vmax_)
{
  const R vmin = vmin_[0];
  const R vmax = vmax_[0];
  R met1[3];
  R met2[3],met2_ee[2],met2_ev[4];
  R met3[3];
  { I i;
    for (i=0;i<metric_->n;++i)
      {
	met2[0] = hessian_->x[i];
	met2[1] = hessian_->y[i];
	met2[2] = hessian_->z[i];
	
	MetricTensor_eig(met2,met2_ee,met2_ev);
	met2_ee[0] = nsFABS(met2_ee[0]);
	met2_ee[1] = nsFABS(met2_ee[1]);
	met2_ee[0] *= itol_[0];
	met2_ee[1] *= itol_[0];	
	if ( met2_ee[0] < vmin )
	  met2_ee[0] = vmin;
	if ( met2_ee[0] > vmax )
	  met2_ee[0] = vmax;
	if ( met2_ee[1] < vmin )
	  met2_ee[1] = vmin;
	if ( met2_ee[1] > vmax )
	  met2_ee[1] = vmax;	
	/*
	  e0a e1c   a  b
	  e0b e1d   c  d
	  
	  e0aa e1cc  e0ba e1cd
	  e0bb e1dd
	*/
	met2[0] = met2_ee[0] * met2_ev[0]*met2_ev[0] + met2_ee[1] * met2_ev[2]*met2_ev[2];
	met2[1] = met2_ee[0] * met2_ev[0]*met2_ev[1] + met2_ee[1] * met2_ev[2]*met2_ev[3];
	met2[2] = met2_ee[0] * met2_ev[1]*met2_ev[1] + met2_ee[1] * met2_ev[3]*met2_ev[3];	
	
	met1[0] = metric_->x[3*i+0];
	met1[1] = metric_->x[3*i+1];
	met1[2] = metric_->x[3*i+2];
	
	sdp_intersect(met1,
		      met2,
		      met3);	
	
	metric_->x[3*i+0] = met3[0];
	metric_->x[3*i+1] = met3[1];
	metric_->x[3*i+2] = met3[2];	

      } }  
}



void ns_tensor2metric(MetricTensor * 	metric_,
		      const ns_var * 	hessian_,
		      cst_pR 		itol_,
		      cst_pR 		vmin_,
		      cst_pR 		vmax_)
{
  const R vmin = vmin_[0];
  const R vmax = vmax_[0];
  R met2[3],met2_ee[2],met2_ev[4];
  { I i;
    for (i=0;i<metric_->n;++i)
      {
	/***/	
	met2[0] = hessian_->x[i];
	met2[1] = hessian_->y[i];
	met2[2] = hessian_->z[i];
	/***/	
	MetricTensor_eig(met2,met2_ee,met2_ev);
	met2_ee[0] = nsFABS(met2_ee[0]);
	met2_ee[1] = nsFABS(met2_ee[1]);
	met2_ee[0] *= itol_[0];
	met2_ee[1] *= itol_[0];	
	if ( met2_ee[0] < vmin )
	  met2_ee[0] = vmin;
	if ( met2_ee[0] > vmax )
	  met2_ee[0] = vmax;
	if ( met2_ee[1] < vmin )
	  met2_ee[1] = vmin;
	if ( met2_ee[1] > vmax )
	  met2_ee[1] = vmax;	
	/*
	  e0a e1c   a  b
	  e0b e1d   c  d
	  
	  e0aa e1cc  e0ba e1cd
	  e0bb e1dd
	*/
	met2[0] = met2_ee[0] * met2_ev[0]*met2_ev[0] + met2_ee[1] * met2_ev[2]*met2_ev[2];
	met2[1] = met2_ee[0] * met2_ev[0]*met2_ev[1] + met2_ee[1] * met2_ev[2]*met2_ev[3];
	met2[2] = met2_ee[0] * met2_ev[1]*met2_ev[1] + met2_ee[1] * met2_ev[3]*met2_ev[3];			
	/***/	
	metric_->x[3*i+0] = met2[0];
	metric_->x[3*i+1] = met2[1];
	metric_->x[3*i+2] = met2[2];	
	/***/	

      } }  
}

void MetricTensor_intersection(const MetricTensor *   metricB_,
			   MetricTensor * 	metric_,
			    cst_pR 		vmin_,
			    cst_pR 		vmax_)
{
  R met1[3];
  R met2[3];
  R met3[3];
  { I i;
    for (i=0;i<metric_->n;++i)
      {
	met2[0] = metricB_->x[3*i+0];
	met2[1] = metricB_->x[3*i+1];
	met2[2] = metricB_->x[3*i+2];
	
	met1[0] = metric_->x[3*i+0];
	met1[1] = metric_->x[3*i+1];
	met1[2] = metric_->x[3*i+2];
	
	sdp_intersect(met1,
		      met2,
		      met3);	
	
	metric_->x[3*i+0] = met3[0];
	metric_->x[3*i+1] = met3[1];
	metric_->x[3*i+2] = met3[2]; 
      } }  
}



void MetricTensor_write_medit(pMetricTensorReadOnly 	self_,
			      const char * 		name_,
			      ...)
{
  char ctmp[512];  
  char ctmp2[512];
  va_list args;
  va_start (args,name_);
  vsprintf(ctmp2,name_,args);
  va_end(args);
  sprintf(ctmp,"%s.sol",ctmp2);
  FILE * fil= fopen(ctmp,"w");
  fprintf(fil,"MeshVersionFormatted\n1\nDimension\n2\nSolAtVertices\n"ifmt"\n1 3\n",self_->n);  

  //  fprintf(fil,"2 3 "ifmt" 2\n",self_->n);  
  { I k;
    for (k=0;k<self_->n;++k)
      {
	fprintf(fil,"" rfmt " " rfmt" " rfmt "\n",
		self_->x[3*k+0],
		self_->x[3*k+1],
		self_->x[3*k+2]);
      } }

  fprintf(fil,"End");
  fclose(fil);
}

void MetricTensor_read_medit(MetricTensor * 		metric_,
			 const char * 		name_,
			 ...)
{
  char ctmp2[512];
  va_list args;
  va_start (args,name_);
  vsprintf(ctmp2,name_,args);
  va_end(args);
  FILE * fil= fopen(ctmp2,"r");
  {  I a1,a2,a3,a4;
    int jj=fscanf(fil,""ifmt" "ifmt" "ifmt" "ifmt"\n",&a1,&a2,&a4,&a3);
    if (jj){}
    if (metric_->x)
      {
	if (a4!=metric_->n)
	  {
	    fprintf(stderr,"file '%s' is corrupted\n",name_);
	    exit(1);
	  }
      }
    else
      {
	metric_->x = malloc(sizeof(R)*3*a4);
	metric_->n = a4;
      }
    { I k;
      for (k=0;k<metric_->n;++k)
	{
	  jj=fscanf(fil,"%lf %lf %lf",
		 &metric_->x[3*k+0],
		 &metric_->x[3*k+1],
		 &metric_->x[3*k+2]);
	} }		
    fclose(fil);
  }
}


void MetricTensor_geometric_statistical(MetricTensor * M_)
{
  const ns_mesh * mesh		= M_->mesh;
  pI 		patchelm 	= malloc(sizeof(I)*mesh->nelm);
  const I nvertex 		= ns_mesh_get_numNodes(mesh);
  pI marker 			= calloc(nvertex,sizeof(I));
  pI select 			= calloc(nvertex,sizeof(I));
  R tens[3];
  R coo_edge[1024];  
  { I i;
    for (i=(I)0;i<nvertex;++i)
      {     
	I patchelm_n;
	cst_pR coo_P = &mesh->coo[2*i];
	{ I j;
	  for (j=mesh->bpatch[i];j<mesh->bpatch[i+1];++j)
	    {
	      patchelm[j-mesh->bpatch[i]] = mesh->patch[j];
	    } }
	patchelm_n = mesh->bpatch[i+1]-mesh->bpatch[i];
	I n=0;
	{ I j;
	  for (j=0;j<patchelm_n;++j)
	    {
	      cst_pI cnc = &mesh->cnc[6*patchelm[j]];
	      { I k;
		for (k=0;k<3;++k)
		  {
		    if ( (cnc[k]!=i+1) AND (NOT marker[cnc[k]-1]) )
		      {
			marker[cnc[k]-1]=++n;
			select[n-1]=cnc[k]-1;
		      }
		  } }
	    } }
	{ I j;
	  for (j=0;j<n;++j)
	    {
	      cst_pR coo_vertex = &mesh->coo[2*select[j]];
	      coo_edge[2*j+0]=coo_vertex[0]-coo_P[0];
	      coo_edge[2*j+1]=coo_vertex[1]-coo_P[1];
	    } }
	tens[0]=((R)0.0);
	tens[1]=((R)0.0);
	tens[2]=((R)0.0);
	{ I j;
	  for (j=0;j<n;++j)
	    {
	    tens[0]+=coo_edge[2*j+0]*coo_edge[2*j+0];
	    tens[1]+=coo_edge[2*j+0]*coo_edge[2*j+1];
	    tens[2]+=coo_edge[2*j+1]*coo_edge[2*j+1];
	    } }

	const R x = ((R)1.0)/((R)n);
	tens[0]*=x;
	tens[1]*=x;
	tens[2]*=x;
	R y=tens[0]*tens[2]-tens[1]*tens[1];
	R a = ((R)1.0)/y;
	R b = tens[0];
	tens[0] = tens[2]*a*((R)0.5);
	tens[2] = b*a*((R)0.5);
	tens[1] *= -a*((R)0.5);
	M_->x[3*i+0]=tens[0];
	M_->x[3*i+1]=tens[1];
	M_->x[3*i+2]=tens[2];
	{ I j;
	  for (j=0;j<n;++j)
	    {
	      marker[select[j]]=(I)0;
	    } }
      } } 	
  free(marker);
  free(select);
  free(patchelm);
  return ;
}
