#include "Type.h"

void shape_l2ortho_3_f(R A[],
			const R x,
			const R y)
{
  A[0]=(R)1.41421356237309504880168872421e+0;
  A[1]=(R)6.0e+0*x-2.0e+0;
  A[2]=(R)6.92820323027550917410978536602e+0*y+3.46410161513775458705489268301e+0*x-3.46410161513775458705489268301e+0;
};

void shape_l2ortho_3_fdx(R 	A[],

			  const R 	x,
			  const R 	y)
{
  A[0]=(R)0.0e+0;
  A[1]=(R)6.0e+0;
  A[2]=(R)3.46410161513775458705489268301e+0;
};

void shape_l2ortho_3_fdy(R A[],
			  const R x,
			  const R y)
{
  A[0]=(R)0.0e+0;
  A[1]=(R)0.0e+0;
  A[2]=(R)6.92820323027550917410978536602e+0;
};
void shape_l2ortho_3_fdxx(R A[],
			   const R x,
			   const R y)
{
  A[0]=(R)0.0e+0;
  A[1]=(R)0.0e+0;
  A[2]=(R)0.0e+0;
};
void shape_l2ortho_3_fdxy(R A[],
			   const R x,
			   const R y)
{
  A[0]=(R)0.0e+0;
  A[1]=(R)0.0e+0;
  A[2]=(R)0.0e+0;
};

void shape_l2ortho_3_fdyy(R A[],
			   const R x,
			   const R y)
{
  A[0]=(R)0.0e+0;
  A[1]=(R)0.0e+0;
  A[2]=(R)0.0e+0;
};
