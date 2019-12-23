#include "nsGLOBAL.h"
#include "ns_constantes.h"

void nsGLOBAL_transient_print(const nsGLOBAL_ST*const G_)
{	  
  pTimeReadOnly const gTimeInfo	= GlobalReadOnly_get_Time(&G_[0]);
  ns_var_print		(&G_[0].var_p,"p",
			 &gTimeInfo->itime_interval,
			 &gTimeInfo->itime);
  ns_var_print		(&G_[0].var_u,"u",
			 &gTimeInfo->itime_interval,
			 &gTimeInfo->itime);

#if 0
  const nsPARAMS_ST*const Params = &G_[0].params;
#endif

#if 0  
  if (  G_[0].params.linfo[__ens_linfo_vnsdg])
    {
#if 0
      { mkSTRING aa;
	sprintf(aa,"integrator.dg.%.5"nsFORMAT_INTEGER".bin",tim->itime);
	ns_var_write_bin(&G_[0].var_dg,aa); }
#endif
#if 0
      print_dscalar1	("fdisc",x,_nf);
#endif
      if (0)
	{

#if 0       
	  { mkSTRING aa;
	    sprintf(aa,"integrator.dg.%.5"nsFORMAT_INTEGER".bb",tim->itime);
	    FILE * fich = fopen(aa,"w");
	    fprintf(fich,"2 1 "ifmt" 1\n",G_[0].pointeur_maillage->nelm);
	    {I i;for(i=0;i<G_[0].pointeur_maillage->nelm;++i)fprintf(fich,"%e\n",G_[0].var_dg.x[i*G_[0].nddlelmdg]*sqrt(2.0));}
	    fclose(fich);
#endif
	  }
      
	}
#if 0
      mkSTRING aa;
      sprintf(aa,"integrator.dg.%.5"nsFORMAT_INTEGER"",tim->itime);
      ns_var_write_medit_with_mesh(&G_[0].var_dg,aa); 
      
#endif
      
#if 0
      if (pele_pele==0)
	{
	  pele_pele=1;
	  ns_print_dgspace	("dg.mesh",
				 _shape_dg);
	  
	}
#endif
      
#if 0
      nsGLOBAL_print_dgspace	(&G_[0],
				 "visudg.mesh",
				 _shape_dg);
      nsGLOBAL_dg_print_sol	(&G_[0],
				 "visudg",		 
				 _shape_dg,
				 G_[0].var_dg.x,
				 &G_[0].dg_n);      
#endif
      
    }
#if 0
  if  ( (Params->linfo[__ens_linfo_axisymetric_x]) OR (Params->linfo[__ens_linfo_axisymetric_y]) )
    {
      const I axi_axe 	= (Params->linfo[__ens_linfo_axisymetric_x])?0:1;
      I * axi_perm  	= imalloc(G_[0].pointeur_maillage->nvertex);
      I * axi_rperm 	= imalloc(G_[0].pointeur_maillage->nvertex);
      const I axi_nrot 	= 40;
      I axi_nvertex		= 0;
      I noaxi_nvertex	= 0;
      print_mesh_extrude_axi_permutation	(G_[0].pointeur_maillage,
						 &axi_axe,
						 axi_perm,
						 axi_rperm,
						 &axi_nvertex,
						 &noaxi_nvertex);	  
      char ctmp[512];
      sprintf(ctmp,"integrator.faxi.%.5"nsFORMAT_INTEGER"",tim->itime);
      print_scalar_extrude_axi		(G_[0].pointeur_maillage,
					 ctmp,
					 &axi_nrot,
					 axi_rperm,			   			   
					 &axi_nvertex,
					 &noaxi_nvertex,
					 G_[0].var_f.x);
      ifree(axi_perm);
      ifree(axi_rperm); 
    }	
#endif
#endif
}
