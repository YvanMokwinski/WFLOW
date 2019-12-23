#include <signal.h>
#include <pthread.h>
#include "Cmdline.h"
#include "nsGLOBAL.h"
//#include "nsPBLM.h"
#include "Monitor.h"
#include "ns_sys.h"
#include "MetricTensor.h"
#include <stdio.h>



void MetricTensor_write_pie(pMetricTensor self_,const char * metricName,const char * meshName,FILE * out)
{
  const ns_mesh * mesh = self_->mesh;
  fprintf(out,"CHAMP<double> %s_ValTenseur ()={",metricName);

  const I numNodes = ns_mesh_get_numNodes(mesh);
  { I i;
    for (i=0;i<numNodes;++i)
      {
	fprintf(out,"" rfmt " " rfmt " " rfmt "\n",
		self_->x[3*i+0],
		self_->x[3*i+1],
		self_->x[3*i+2]);
      } }
  
  fprintf(out,"};");
  
  fprintf(out,"CHAMP<int> %s_ConnecP1 ()={",metricName);
  {
    I cncelm[16];
    const I numCells = ns_mesh_nelm(mesh);
    { I cellIndex;
      for (cellIndex=0;cellIndex<numCells;++cellIndex)
	{
	  ns_mesh_cncelm		(mesh,
					 &cellIndex,
					 cncelm);
	  
	  fprintf(out,"" ifmt " " ifmt " " ifmt "\n",
		  1+cncelm[0],
		  1+cncelm[1],
		  1+cncelm[2]);
	} }
  }
  fprintf(out,"};");

  fprintf(out,"SOLUTION SolutionQuiContient%s(%s_Maillage)={\n",metricName,meshName);  
  fprintf(out,"VARIABLE %s(LagrTrian03,%s_ValTenseur%%3,%s_ConnecP1,%s);",metricName,metricName,metricName,meshName);  
  fprintf(out,"};");

}


void ns_mesh_write_pie(ns_mesh * mesh,
		       const char * name,
		       FILE * out)
{

  fprintf(out,"CHAMP<double> %s_coo()=\n{\n",name);
  
  {
    R xc[3];
    const I numNodes = ns_mesh_get_numNodes(mesh);
    { I nodeIndex;
      for (nodeIndex=0;nodeIndex<numNodes;++nodeIndex)
	{
	  ns_mesh_get_coovertex		(mesh,
					 &nodeIndex,
					 xc);
	  fprintf(out,"" rfmt " " rfmt "\n",xc[0],xc[1]);
	} }
  }

  fprintf(out,"};\n");
  fprintf(out,"CHAMP<int> %s_cnc()=\n{\n",name);

  {
    I cncelm[16];
    const I numCells = ns_mesh_nelm(mesh);
    { I cellIndex;
      for (cellIndex=0;cellIndex<numCells;++cellIndex)
	{
	  ns_mesh_cncelm(mesh,
			 &cellIndex,
			 cncelm);
	  
	  fprintf(out,"" ifmt " " ifmt " " ifmt "\n",
		  1+cncelm[0],
		  1+cncelm[1],
		  1+cncelm[2]);	  
	} }
  }


  fprintf(out,"};\n");
  fprintf(out,"CHAMP<int> %s_cod()=\n{\n",name);

  {
    const I numNodes = ns_mesh_get_numNodes(mesh);
    { I nodeIndex;
      for (nodeIndex=0;nodeIndex<numNodes;++nodeIndex)
	{
	  fprintf(out,"" ifmt "\n",mesh->own_cod[nodeIndex]);
	} }
  }


  // MOK_GEOM_NAME
  fprintf(out,"};\n");
  fprintf(out,"Maillage %s_Maillage( LaGeometrie ,  %s_coo%%2,%s_cod%%1 ) =\n",
	  name,
	  name,
	  name);
  fprintf(out," { ZONE %s (LagrTrian03,%s_coo%%2,%s_cnc,100,%s_cod%%1); };\n",
	  name,
	  name,
	  name,
	  name);
}


#if 0
static void ns_main_signal(int sigid) 
{

  fflush(stdout);
  fflush(stderr);
  fprintf(stderr,"\n");
  Monitor_errmsg(0,"ns_main_signal:Unexpected error"); 
  switch(sigid) {
    case SIGABRT:
      Monitor_errmsg(0,"ns_main_signal:<Abnormal stop>");  
      exit(1);
      break;
    case SIGFPE:
      Monitor_errmsg(0,"ns_main_signal:<Floating-point exception>"); 
      exit(1);
      break;
    case SIGILL:
      Monitor_errmsg(0,"ns_main_signal:<Illegal instruction>"); exit(1);break;
    case SIGSEGV:
      Monitor_errmsg(0,"ns_main_signal:<Segmentation fault>");  exit(1);break;
    case SIGTERM:
      {
	Monitor_errmsg(0,"<Program killed>2");  
	exit(1);
	break;
      }
    case SIGINT:
      {
	Monitor_errmsg(0,"<Program killed>");  
	exit(1);
	break;
      }
  default:
    {
      Monitor_errmsg(0,"  ??? unrecognized error");  
      exit(1);
      break;
    }
  }

}

static void ns_main_atexit()
{
#ifndef NDEBUG
  fprintf(stderr,"ns_main_atexit ...\n");
#endif

#ifndef NDEBUG
  fprintf(stderr,"ns_main_atexit done.\n");
#endif
}
#endif
void  Global_checkError(const char * 	filename_,
			const int 	fileline_,
			STR 		errmsg_,
			Err*		err_)
{
  if (err_[0])
    {
      Monitor_errmsg(0,"ERROR DETECTED");
      Monitor_errmsg(0,"ERROR CODE      : '%s' (err=" efmt ")",eErr_get_string(err_[0]),err_[0]);
      Monitor_errmsg(0,"ERROR MESSAGE   : '%s'",errmsg_);
      Monitor_errmsg(0,"ERROR FROM FILE : '%s'",filename_);
      Monitor_errmsg(0,"ERROR AT LINE   : %d",fileline_);
      exit(err_[0]);
    }
}






static int index(int i,int j)
{
  return j*2+i;
};


static void Mul(cst_pR A,cst_pR B,pR C)
{
  C[index(0,0)] = A[index(0,0)] * B[index(0,0)] + A[index(0,1)] * B[index(1,0)];
  C[index(0,1)] = A[index(0,0)] * B[index(0,1)] + A[index(0,1)] * B[index(1,1)]; 

  C[index(1,0)] = A[index(1,0)] * B[index(0,0)] + A[index(1,1)] * B[index(1,0)];
  C[index(1,1)] = A[index(1,0)] * B[index(0,1)] + A[index(1,1)] * B[index(1,1)];

};

#define CHECK_ERROR(_errmsg,_err) Global_checkError(__FILE__,__LINE__,_errmsg,_err)
  static void ComputeMetric(cst_pR xyz_,pR tensor)
  {
    const R x = xyz_[0];
    const R y = xyz_[1];

    static const R ox = 1.0;
    static const R oy = 0.5;
    static const R eps = 0.125/8.0;
    
    const R radius = sqrt((x-ox)*(x-ox)+(y-oy)*(y-oy));

    const R h = 0.005;

    R P[4];
    R PT[4];
    R D[4];
    R C[4];
    R C1[4];

    R normal[2];
#if 1
    if ( (radius <= 0.25 + 8.0*eps) && (radius >= 0.25 - 8.0*eps) )
      {
#endif
	normal[0] = (x-ox)/radius;
	normal[1] = (y-oy)/radius;

	P[0] = normal[0];
	P[1] = normal[1];
	P[2] = -P[1];
	P[3] = P[0];

	PT[0] = P[0];
	PT[1] = P[2];
	PT[2] = P[1];
	PT[3] = P[3];
	
	D[0] = exp(  -512.0 * (radius - 0.25) * (radius - 0.25) ) / h / h ;
	D[1] = 0.0;
	D[2] = 0.0;
	D[3] = D[0] / 32.0 / 32.0;// 1.0 / (h*32.0) /(h*32.0);

	Mul(D,PT,C);
	Mul(P,C,C1);
	
	//
	//  P0 P2  D0 0   P0 P1
	//  P1 P3  0  D3  P2 P3
	//
	//
	//  P0 P2  D0*P0 D0*P1
	//  P1 P3  D3*P2 D3*P3
	//
	
	tensor[0] = C1[0];
	tensor[1] = (C1[1] + C1[2])*0.5;
	tensor[2] = C1[3];
#if 1
      }
    else
      {
	tensor[0] = 64.0;
	tensor[1] = 0.0;
	tensor[2] = 64.0;	
      }
#endif
    
  };




static pMetricTensor Read(const char * name,
			  const char * out,
			  ns_mesh ** mesh_)
  {
  
    //
    // READ MESH
    //
    ns_mesh *mesh = (ns_mesh*)calloc(1,sizeof(ns_mesh));
    STR errmsg;
    Err err;
      ns_mesh_read(mesh,
		   errmsg,
		   &err,
		   name);
    
#if 0
      CHECK_ERROR(errmsg,
		  &err);
#endif
      pMetricTensor metric = ns_mesh_addmetric(mesh);
      
      MetricTensor_geometric_statistical(metric);

      {
	R xc[3];
	R tensor[3];
	const I numNodes = ns_mesh_get_numNodes(mesh);
	I nodeIndex;
	for (nodeIndex=0;nodeIndex<numNodes;++nodeIndex)
	  {
	    ns_mesh_get_coovertex		(mesh,
						 &nodeIndex,
						 xc);

	    ComputeMetric(xc,tensor);
	    
	    metric->x[3*nodeIndex+0] = tensor[0];
	    metric->x[3*nodeIndex+1] = tensor[1];
	    metric->x[3*nodeIndex+2] = tensor[2];
	  }
      }

      
      ns_mesh_write_medit(mesh,
			  out);

      MetricTensor_write_medit(metric,
			       out);

      mesh_[0] = mesh;
      return metric;
  };

  

//  atexit(ns_main_atexit);
//  signal(SIGABRT,ns_main_signal);
//  signal(SIGFPE,ns_main_signal);
//  signal(SIGILL,ns_main_signal);
//  signal(SIGSEGV,ns_main_signal);
//  signal(SIGTERM,ns_main_signal);
//  signal(SIGINT,ns_main_signal);

#include "Shape.hpp"
int main(int 		argc,
	 const char**	argv)
{


  

  
  //
  // Set up the monitor
  //
  Monitor_def(0,
	      argv[0],
	      MonitorMode_STD);


  char out[512];
  const char * filename = argv[1];
  pMetricTensor metricTensor;
  for (int iter=1;iter <= 10;++iter)
  {
    sprintf(out,"allo%d",iter);
    ns_mesh * mesh;
    metricTensor = Read(filename,out,&mesh);

#if 1
    
    { FILE * fich = fopen("mesh2adapt.pie","w");
      ns_mesh_write_pie(mesh,
			"mesh2adapt",
			fich);      
      fclose(fich); }

    
    { FILE * fich = fopen("meshmetric.pie","w");
      ns_mesh_write_pie(mesh,
			"meshmetric",
			fich);
      
      MetricTensor_write_pie(metricTensor,
			     "myMetric",
			     "meshmetric",
			     fich);
      fclose(fich); }

    { FILE * fich = fopen("tout.pie","w");
      fprintf(fich,"#include \"GEOMETRIE.pie\"\n");
      fprintf(fich,"#include \"mesh2adapt.pie\"\n");
      fprintf(fich,"#include \"meshmetric.pie\"\n");
      fclose(fich); }
    
    {
      char cmd[512];   
      sprintf(cmd,"oort.exe -N -T myMetric -o res.pie tout.pie"); 
      int i = system(cmd);
      i = system("pie2medit.exe res.pie -o zzz");
      if (i!=0)
	{

	}
    }


    


#else
    {
      char cmd[512];   
      sprintf(cmd,"../ThirdParty/mmg-master/build/bin/mmg2d_O3 -hmin 0.000000001  -in %s.mesh -sol %s.sol -out zzz.mesh",
	      out,
	      out);
      
      system(cmd);
    }
#endif

    metricTensor = MetricTensor_kill(metricTensor);
    ns_mesh_free(mesh);

    
    filename = "zzz.mesh";
    
		
  }  
    return 0;


#if 0
  //
  // Call main
  // 
  {
    STR errmsg;
    Err err;
    ns_main(argc,
	    argv,
	    errmsg,
	    &err);
    
    return err;
  }
#endif

  return 0;
}
