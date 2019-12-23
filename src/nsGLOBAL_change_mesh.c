#include "nsGLOBAL.h"
#include "Monitor.h"

void nsGLOBAL_change_mesh(nsGLOBAL_ST*const 	self_,
			  const char * 		basename_)
{   
  char ctmp[512];	  	  	  
  int err_system;

  switch(self_->meshAdaptivityMethod)
    {
    case __eMeshAdaptivity_STATIC2:
      {
	/*on subdivise le maillage*/
#if 0
	ns_mesh_print(self_->metric->mesh,"remeshinput");
	nsMETRIC_print(self_->metric,"remeshinput");
#endif
	sprintf(ctmp,"mkmake_remesh remeshinput.mesh -o remesh.mesh --remesh-l2 --remesh-log mkmake_remesh.log");	  
	err_system = system(ctmp);
	if ( err_system )
	  {
	    Monitor_errmsg(self_->iproc,"pre_integrator failed, err = %d\n",err_system);
	  }  
	break;
      }
    case __eMeshAdaptivity_STATIC3:
      {
	/*on subdivise le maillage*/
#if 0
	ns_mesh_print(self_->metric->mesh,"remeshinput");
	nsMETRIC_print(self_->metric,"remeshinput");
#endif
	sprintf(ctmp,"mkmake_remesh remeshinput.mesh -o remesh.mesh --remesh-l3 --remesh-log mkmake_remesh.log");	  
	err_system = system(ctmp);
	if ( err_system )
	  {
	    Monitor_errmsg(self_->iproc,"pre_integrator failed, err = %d\n",err_system);
	  }	  
	break;
      }
    case __eMeshAdaptivity_ISOTROPIC:
      {
#if 0
	ns_mesh_print(self_->metric->mesh,"remeshinput");
	nsMETRIC_print(self_->metric,"remeshinput");	  
#endif
	sprintf(ctmp,"mkmake_remesh remeshinput.bb -o remesh.mesh --remesh-geom %s_geometry.pie --remesh-log mkmake_remesh.log",basename_);	  
	err_system = system(ctmp);
	if ( err_system )
	  {
	    Monitor_errmsg(self_->iproc,"mkmake_remesh failed, err = %d, dommand was '%s'\n",err_system,ctmp);
	    exit(1);
	  }	  
	break;
      }
    case __eMeshAdaptivity_ANISOTROPIC:
      {
#if 0
	ns_mesh_print(self_->metric->mesh,"remeshinput");
	nsMETRIC_print(self_->metric,"remeshinput");	  
#endif
	/*--remesh-smooth-borouchaki 2.0*/
	sprintf(ctmp,"mkmake_remesh remeshinput.bb -o remesh.mesh  --remesh-geom %s_geometry.pie --remesh-log mkmake_remesh.log",basename_);	  
	err_system = system(ctmp);
	if ( err_system )
	  {
	    Monitor_errmsg(self_->iproc,"mkmake_remesh failed, err = %d, dommand was '%s'\n",err_system,ctmp);
	    exit(1);
	  }	  
	break;
      }
    case __eMeshAdaptivity_ALL:
    case __eMeshAdaptivity_ERROR:
      {
	Monitor_errmsg(self_->iproc,"nsGLOBAL_change_mesh:wrong mesh adaptivity method");
	break;
      }
    }

  //
  // Free the mesh
  //
  ns_mesh_free		((ns_mesh*)self_->mesh_usrptr);  

  //
  // Free
  //
  Global_free		(self_); 

  Err err;
  STR errmsg;

  //
  // Read the new mesh
  //
  ns_mesh_read		((ns_mesh*)self_->mesh_usrptr,
			 errmsg,
			 &err,
			 "remesh.mesh");

  //
  // Initialize
  //

#if 0
  Global_initialize	(self_,
			 (ns_mesh*)self_->mesh_usrptr,
			 errmsg,&err);
#endif

  
}



