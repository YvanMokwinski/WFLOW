#ifndef __header_nsOP_INFO_h__
#define __header_nsOP_INFO_h__
#include "eSmoothingMethod.h"
#include "eOperator.h"

typedef struct
{
  eOperator			op;
  eSmoothingMethod 		smoothing_method;
  I 				qh_degree;
} nsOP_INFO_ST,*nsOP_INFO;

typedef const nsOP_INFO_ST cst_nsOP_INFO_ST;
typedef cst_nsOP_INFO_ST *cst_nsOP_INFO;

void nsOP_INFO_def		(nsOP_INFO 			operator_,
				 cst_eOperator 		op_,
				 cst_eSmoothingMethod 		smoothing_method,
				 const I 			qh_degree_);

static const nsOP_INFO_ST instance_nsOP_INFO_IH[1] = {{__eOperator_ih,__eSmoothingMethod_ERROR,0}};
static const nsOP_INFO_ST instance_nsOP_INFO_PH[1] = {{__eOperator_ph,__eSmoothingMethod_ERROR,0}};

#define nsOP_INFO_ih(_operator) 			nsOP_INFO_def((_operator),__eOperator_ih,__eSmoothingMethod_ERROR,0)
#define nsOP_INFO_ph(_operator) 			nsOP_INFO_def((_operator),__eOperator_ph,__eSmoothingMethod_ERROR,0)
#define nsOP_INFO_qh(_operator,_method,_degree) 	nsOP_INFO_def((_operator),__eOperator_qh,(_method),(_degree))

#endif
