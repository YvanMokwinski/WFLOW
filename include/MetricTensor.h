#ifndef __header_MetricTensor_h__
#define __header_MetricTensor_h__
#include "ns_config.h"

#include "ns_var.h"
#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
  I 	n;
  I 	off;
  I 	m;
  pR 		x;  
} nsRHDLE_ST;

void nsRHDLE_def(nsRHDLE_ST*hdle_,const I n_,const I m_,pR x_,const I off_);

/**
  \brief Metric information    
*/
typedef struct
{
  cst_pR 	hmin;
  cst_pR 	hmax;
  cst_pR 	itol;
} nsOP_METRIC_INFO_ST,*nsOP_METRIC_INFO;

typedef const nsOP_METRIC_INFO_ST cst_nsOP_METRIC_INFO_ST;
typedef cst_nsOP_METRIC_INFO_ST * cst_nsOP_METRIC_INFO;

void 		nsOP_METRIC_INFO_def	(nsOP_METRIC_INFO 	op_,
					 cst_pR 		hmin_,
					 cst_pR		hmax_,
					 cst_pR 		itol_);

/**
  \brief Metric object
*/
typedef struct
{
  nsRHDLE_ST 	data;
  eDim		dim;
  pR   		x;
  I 		n;
  const ns_mesh*mesh;
} MetricTensor,*pMetricTensor;

typedef const MetricTensor * pMetricTensorReadOnly;

/**
   \brief create a metric from the mesh
   @param mesh_ 	pointer to the mesh
   @return pointer to the created metric
*/
pMetricTensor	ns_mesh_addmetric		(ns_mesh*		mesh_);
/**
   \brief kill metric from mesh
   @param mesh_ 	pointer to the mesh
   @param metric_ 	pointer to the metric
*/
MetricTensor*	ns_mesh_killmetric		(ns_mesh*		mesh_,
						 MetricTensor*		metric_);

/**
   \brief free a  metric
   @param metric_ 	pointer to the metric
*/
void 		MetricTensor_free(MetricTensor*	metric_);

/**
   \brief define a  metric
   @param metric_ 	pointer to the metric
   @param mesh_ 	pointer to the mesh
*/
void 		MetricTensor_def			(pMetricTensor		self_,
							 const ns_mesh *	mesh_);
/**
   \brief kill a  metric
   @param metric_ 	pointer to the metric
   @return NULL
*/
pMetricTensor  	MetricTensor_kill			(pMetricTensor		self_);

/**
   \brief create metric
   @param mesh_ 	pointer to the mesh
   @return pointer to the metric
*/
pMetricTensor 	MetricTensor_new			(const ns_mesh*		mesh_);

/**
   \brief read metric from medit file format
   @param metric_ 	pointer to the struct
   @param name_ 	name file
*/
void 		MetricTensor_read_medit			(pMetricTensor		self_,
							 const char *		name_,
							 ...);

/**
   \brief write metric into medit file format
   @param metric_ 	pointer to the struct
   @param name_ 	name file
*/
void 		MetricTensor_write_medit		(pMetricTensorReadOnly self_,
							 const char * 	name_,
							 ...);

/**
   \brief print metric (medit format)
   @param metric_ 	pointer to the struct
   @param name_ 	name file
*/
Err 	MetricTensor_print			(pMetricTensorReadOnly 	self_,
						 const char * 		name_);


#if __ns_debug__
cst_pR 		cst_MetricTensor_x			(pMetricTensorReadOnly metric_);
pR 		MetricTensor_x				(MetricTensor*metric_);
I 		MetricTensor_n				(pMetricTensorReadOnly metric_);
const ns_mesh*	MetricTensor_mesh			(pMetricTensorReadOnly metric_);
#else
#define 	MetricTensor_mesh(_m) 	(_m)->mesh
#define 	cst_MetricTensor_x(_m) 	(_m)->x
#define 	MetricTensor_x(_m) 	(_m)->x
#define 	MetricTensor_n(_m) 	(_m)->n
#endif

void 		ns_var_vector2metric		(const ns_var * 	nabla_,
						 pMetricTensor		metric_,
						 pMetricTensorReadOnly  	metric_to_respect_,
						 cst_nsOP_METRIC_INFO	metric_info_);

void 		ns_var_tensor2metric		(const ns_var * 	tensor_,
						 MetricTensor*		metric_,
						 pMetricTensorReadOnly 	metric_to_respect_,
						 cst_nsOP_METRIC_INFO	metric_info_);


void 		MetricTensor_geometric_statistical	(pMetricTensor  M_);


void 		MetricTensor_from_size		(pMetricTensor metric_,
						 const R h_);

void *  	MetricTensor_transfert		(pMetricTensor  		v_,
						 pMetricTensorReadOnly  	vold_,
						 void * 		old_localization_,
						 Err *		err_);

void 		MetricTensor_compute_conformite	(pMetricTensorReadOnly  	M_,
						 pMetricTensorReadOnly  	Q_,
						 R 		mmmv_[4],
						 const char * 		filename_);


/*-----------------------------------------------------------------------------*/
/**
   \brief compute eigenvalues and eigenvectors
   @param T tensor
   @param ee_ eigenvalues
   @param ev_ eigenvectors
*/
void MetricTensor_eig(const R	T[3],
		   R  	ee_[2],
		   R 	ev_[4]);

/*-----------------------------------------------------------------------------*/
/**
   \brief compute intersection
   @param A_ metric
   @param B_ metric
   @param C_ intersection
*/
void sdp_intersect(cst_pR 				A_,
		   cst_pR 				B_,
		   pR 					C_);


/*-----------------------------------------------------------------------------*/
/**
   \brief set metric to identity
   @param metric_ 	pointer to the struct
*/
void MetricTensor_setidentity(pMetricTensor  metric_);


/*-----------------------------------------------------------------------------*/
/**
   \brief compute metric density
   @param mesh_ 	mesh
   @param metric_ 	metric
*/
R MetricTensor_density			(const ns_mesh * 	mesh_,
					 pMetricTensorReadOnly 	metric_);


/*-----------------------------------------------------------------------------*/
/**
   \brief copy metric
   @param metric_ pointer to the struct
   @param metric_to_copy_ pointer to the metric to copy
*/
void MetricTensor_copy			(pMetricTensor  		metric_,
					 pMetricTensorReadOnly 	metric_to_copy_);


/*-----------------------------------------------------------------------------*/
/**
   \brief free metric struct
   @param metric_ pointer to the struct
*/


/*-----------------------------------------------------------------------------*/
/**
   \brief transform vector to tensor vect_*vect_'
   @param metric_ 	pointer to the struct
   @param vect_ 	pointer to the vector finite element variable
*/
void ns_vector2metric			(pMetricTensor  	metric_,
					 const ns_var * 	vect_,
					 cst_pR 		itol_);


/*-----------------------------------------------------------------------------*/
/**
   \brief compute absolute values with constraints values'
   @param metric_ 	pointer to the struct
   @param vmin_ 	minimum eigenvalue 
   @param vmax_ 	maximum eigenvalue 
*/
void MetricTensor_hconstraint		(pMetricTensor  		metric_,
					 cst_pR 		vmin_,
					 cst_pR 		vmax_);


/*-----------------------------------------------------------------------------*/
/**
   \brief compute absolute values with constraints practical values and threshold with intersection procedure'
   @param metric_ 	pointer to the struct
   @param tensor_ 	pointer to the tensor
   @param itol_ 	inverse tolerance
   @param vmin_ 	minimum eigenvalue 
   @param vmax_ 	maximum eigenvalue 
*/
void ns_tensor2metric_intersection	(pMetricTensor 		metric_,
					 const ns_var * 	tensor_,
					 cst_pR 		itol_,
					 cst_pR 		vmin_,
					 cst_pR 		vmax_);

/*-----------------------------------------------------------------------------*/
/**
   \brief compute absolute values with constraints practical values and threshold'
   @param metric_ 	pointer to the struct
   @param tensor_ 	pointer to the tensor
   @param itol_ 	inverse tolerance
   @param vmin_ 	minimum eigenvalue 
   @param vmax_ 	maximum eigenvalue 
*/
void ns_tensor2metric			(pMetricTensor  		metric_,
					 const ns_var * 	hessian_,
					 cst_pR 		itol_,
					 cst_pR 		vmin_,
					 cst_pR 		vmax_);

/*-----------------------------------------------------------------------------*/
/**
   \brief compute metric intersection with constraints practical values'
   @param metricB_ 	pointer to the metric
   @param metric_ 	pointer to the metric
   @param vmin_ 	minimum eigenvalue 
   @param vmax_ 	maximum eigenvalue 
*/
void MetricTensor_intersection		(pMetricTensorReadOnly    	metricB_,
					 pMetricTensor  		metric_,
					 cst_pR 			vmin_,
					 cst_pR 			vmax_);

#ifdef __cplusplus
}
#endif

#endif
