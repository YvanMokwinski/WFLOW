#if 0
#include "DgTransport.h"
#include "Err.h"
int DgTransport_test(int 		argc,
		     const char ** 	argv)
{
  /*_____________________________HELP______________________________*/
  char ofilename	[512];
  Err err = __eErr_no;
  L info 		= __emnsYES;
  if (info)
    {
      printf("mokmake_dg ... \n");
    }

  
  mkS_st 	teta_a;
  mkS_st 	teta_u;
  mkS_st 	teta;
  mkS_st 	test;
  mkS_st 	trial;
  const I 	test_degree = 2;
  mkS_definit(&test  ,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,test_degree,__emk_discontinuous,&err);
  mkS_definit(&trial ,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,test_degree,__emk_discontinuous,&err);
  mkS_definit(&teta  ,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,test_degree,__emk_discontinuous,&err);
  mkS_definit(&teta_u,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,test_degree,__emk_discontinuous,&err);
  mkS_definit(&teta_a,__eTopology_TRIANGLE,__emkS_FAMILY_lagrange,test_degree,__emk_discontinuous,&err);
  I dg_iinfo[DG_I_n];
  dg_iinfo[DG_I_QELM_N]	  = 10;
  dg_iinfo[DG_I_QFACE_N]  = 10;  
  /*
    dg_iinfo[DG_I_TRIAL_FAMILY]   = DGS_L2ORTHONORMAL;
    dg_iinfo[DG_I_TRIAL_NBASIS]   = fn;
    dg_iinfo[DG_I_TRIAL_DEGREE]   = fk;
    dg_iinfo[DG_I_TEST_FAMILY]    = DGS_L2ORTHONORMAL;
    dg_iinfo[DG_I_TEST_NBASIS]    = fn;
    dg_iinfo[DG_I_TEST_DEGREE]    = fk;
    dg_iinfo[DG_I_TETA_FAMILY]    = DGS_L2ORTHONORMAL;
    dg_iinfo[DG_I_TETA_NBASIS]    = fn;
    dg_iinfo[DG_I_TETA_DEGREE]    = fk;
    dg_iinfo[DG_I_TETA_A_FAMILY]  = DGS_CANONIC;
    dg_iinfo[DG_I_TETA_A_DEGREE]  = 0;
    dg_iinfo[DG_I_TETA_A_NBASIS]  = 1;     
    dg_iinfo[DG_I_TETA_U_FAMILY]  = DGS_LAGRANGE;
    dg_iinfo[DG_I_TETA_U_NBASIS]  = fn;
    dg_iinfo[DG_I_TETA_U_DEGREE]  = fk;
    dg_iinfo[DG_ERR] = 0
  */
  R 	rinfo[1024];
  I 	rinfo_n = 1024;
  I 	iinfo_n = DG_I_n;
  R rwork[128000];
  I rwork_n = 128000;
  
#if 0
void  DgTransport_init	(mkS 		s_teta_a_,
			 mkS 		s_teta_u_,
			 mkS 		s_teta_,
			 mkS 		s_test_,
			 mkS 		s_trial_,
			 I		iinfo_[DG_I_n],
			 cst_pI		rinfo_n_,
			 pR 		rinfo_,
			 cst_pI 	rwork_n_,
			 pR		rwork_);
#endif
 dgadvection(&teta_a,
		   &teta_u,
		   &teta,
		   &test,
		   &trial,
		   &iinfo_n,
		   dg_iinfo,
		   &rinfo_n,
		   rinfo,
		   &rwork_n,
		   rwork);
  //  zone = mkZ_kill(zone);
  //  free(obasename);
  return err;  
}
#endif
