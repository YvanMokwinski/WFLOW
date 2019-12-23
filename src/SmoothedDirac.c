#include "Type.h"
#include "SmoothedDirac.h"
#include "SmoothedHeaviside.h"
#include <math.h>

void SmoothedDirac_def(pSmoothedDirac 	const 	self_,
		       cst_eDirac  		dirac_,
		       cst_pR 	const 	eps_,
		       pErr const err_) 
{
  err_[0] = __eErr_no;
  memset(self_,
	 0,
	 sizeof(SmoothedDirac));
  switch(dirac_)
    {
    case __eDirac_m0p2:
    case __eDirac_m0p3:
    case __eDirac_m0p4:
    case __eDirac_m0p5:
    case __eDirac_m0p6:
    case __eDirac_m2p0:
    case __eDirac_m2p1:
    case __eDirac_m2p2:
    case __eDirac_m2p3:
    case __eDirac_m4p0:
    case __eDirac_m4p1:
    case __eDirac_m4p2:
    case __eDirac_m4p3:
    case __eDirac_sin:
      {
	self_->dirac 	= dirac_;
	self_->eps	= eps_;
	break;
      }

    case __eDirac_ERROR:
    case __eDirac_ALL:
      {
	err_[0] = __eErr_switch;
	break;
      }
    }
}



					  
#if 0
cst_pS __ensDIRAC_rinfo_strings[__ensDIRAC_rinfo_n] 
= { "__ensDIRAC_rinfo_error",
    "eps",
    "epsmin",
    "epsmax"};

cst_pS __ensDIRAC_iinfo_strings[__ensDIRAC_iinfo_n] 
= {"__ensDIRAC_info_error",
   "lev"};


cst_pS __ensHEAVISIDE_rinfo_strings[__ensHEAVISIDE_rinfo_n] 
= { "__ensHEAVISIDE_rinfo_error",
    "eps",
    "epsmin",
    "epsmax"};

cst_pS __ensHEAVISIDE_iinfo_strings[__ensHEAVISIDE_iinfo_n] 
= {"__ensHEAVISIDE_info_error",
   "lev"};
#endif
/*__remk2(ensDIRAC,"dirac kind");*/





#define  r238  			((R)238.0)
#define  r325  			((R)325.0)
#define  r345  			((R)345.0)
#define  r256  			((R)256.0)
#define  r25   			((R)25.0)
#define  r35   			((R)35.0)
#define  r63   			((R)63.0)
#define  r231  			((R)231.0)
#define  r843  			((R)843.0)
#define  r1617 			((R)1617.0)
#define  r4830 			((R)4830.0)
#define  r938  			((R)938.0)
#define  r1212 			((R)1212.0)
#define  r378  			((R)378.0)
#define  r4165 			((R)4165.0)
#define  r7651 			((R)7651.0)
#define  r7938 			((R)7938.0)
#define  r1218 			((R)1218.0)
#define  r1024 			((R)1024.0)
#define  r0    			((R)0.0)
#define  r1    			((R)1.0)
#define  r2    			((R)2.0)
#define  r3    			((R)3.0)
#define  r5    			((R)5.0)
#define  r8    			((R)8.0)
#define  r16   			((R)16.0)
#define  r20   			((R)20.0)
#define  r29   			((R)29.0)
#define  r4    			((R)4.0)
#define  r13   			((R)13.0)
#define  r32   			((R)32.0)
#define  r42   			((R)42.0)
#define  r64   			((R)64.0)
#define  r79   			((R)79.0)
#define  r123  			((R)123.0)
#define  r128  			((R)128.0)
#define  r135  			((R)135.0)
#define  r161  			((R)161.0)
#define  r189  			((R)189.0)
#define  r385  			((R)385.0)
#define  r1220  		((R)1220.0)
#define  r1653  		((R)1653.0)
#define  r269  			((R)269.0)
#define  r666  			((R)666.0)
#define  r495  			((R)495.0)
#define  r8031  		((R)8031.0)
#define  r5005  		((R)5005.0)
#define  r21100 		((R)21100.0)
#define  r22565 		((R)22565.0)
#define  r25480 		((R)25480.0)
#define  r67795 		((R)67795.0)
#define  r3410  		((R)3410.0)
#define  r12210 		((R)12210.0)
#define  r2048 			((R)2048.0)
#define  r49140			((R)49140.0)
#define  r12285 		((R)12285.0)
#define  r990 			((R)990.0)
#define  r2230 			((R)2230.0)
#define  r1540 			((R)1540.0)
#define  r15015 		((R)15015.0)
#define  r162 			((R)162.0)
#define  r2203 			((R)2203.0)				  
#define  r1_2   		((R)1.0/2.0) 
#define  r1_8   		((R)1.0/8.0)
#define  r1_16   		((R)1.0/16.0)
#define  r1_32  		((R)1.0/32.0)
#define  r1_64  		((R)1.0/64.0)
#define  r1_128 		((R)1.0/128.0)
#define  r1_256 		((R)1.0/256.0)
#define  r1_512 		((R)1.0/512.0)
#define  r1_1024 		((R)1.0/1024.0)
#define  r1_2048 		((R)1.0/2048.0)
#define  r1_4096 		((R)1.0/4096.0)
#define  mr1_2   		((R)-0.5)
#define  r21     		((R)21.0)
#define  r45     		((R)45.0)
#define  mr45    		((R)-45.0)
#define  r50     		((R)50.0)
#define  r9      		((R)9.0)
#define  r6      		((R)6.0)
#define  r7      		((R)7.0)
#define  r10     		((R)10.0)
#define  r15     		((R)15.0)
#define  r105    		((R)105.0)
#define  r175    		((R)175.0)
#define  r147    		((R)147.0)





static const R table_bn[__eHeaviside_ALL]={0,0,3,4,5,6,7,3,4,5,6,5,6,7,8};
static const R table_b[8*__eHeaviside_ALL] = {
  r0,r0,r0,r0,r0,r0,r0,r0,/*err*/
  r0,r0,r0,r0,r0,r0,r0,r0,/*sin*/
  r1_16*r8,-r9*r1_16,r3*r1_16,r0,r0,r0,r0,r0,
  r1_32*r16,-r29*r1_32,r20*r1_32,-r5*r1_32,r0,r0,r0,r0,
  r1_256*r128,-r325*r1_256,r345*r1_256,-r175*r1_256,r35*r1_256,r0,r0,r0,
  r1_512*r256,-r1_512*r843,r1_512*r1218,-r1_512*r938,r1_512*r378,-r1_512*r63,r0,r0,
  r1_2048*r1024,-r1_2048*r4165,r1_2048*r7651,-r1_2048*r7938,r1_2048*r4830,-r1_2048*r1617,r1_2048*r231,r0,
  
  r1_8*r4,-r1_8*r5,r1_8*r5,r0,r0,r0,r0,r0,
  r1_32*r16,r1_32*r13,-r1_32*r42,r1_32*r21,r0,r0,r0,r0,
  r1_64*r32,r1_64*r9,-r1_64*r123,r1_64*r135,-r1_64*r45,r0,r0,r0,
  r1_512*r256,-r1_512*r79,-r1_512*r1220,r1_512*r2230,-r1_512*r1540,r1_512*r385,r0,r0,
  
  r1_128*r64,r1_128*r161,-r1_128*r161,-r1_128*r189,r1_128*r189,r0,r0,r0,
  r1_256*r128, r1_256*r269, -r1_256*r666,-r1_256*r162,r1_256*r990,-r1_256*r495,r0,r0,
  r1_2048*r1024, r1_2048*r1653,  -r1_2048*r8031, r1_2048*r3410,  r1_2048*r12210,-r1_2048*r15015,r1_2048*r5005,r0,
  r1_4096* r2048,  r1_4096*r2203,-r1_4096*r21100 ,r1_4096*r22565 ,r1_4096*r25480, -r1_4096*r67795,r1_4096*r49140 ,-r1_4096*r12285
  };

static   const I  	table_b_degree[__eHeaviside_ALL]= {(I)0,(I)1,
							   (I)3,(I)4,(I)5,(I)6,(I)7,
							   (I)1,(I)2,(I)3,(I)4,
							   (I)1,(I)2,(I)3,(I)4};

static const I  table_a_degree[8]={(I)0,(I)1,(I)2,(I)3,(I)4,(I)5,(I)6,(I)7};
static const R 	table_a[8*8] = { 
  r1,r0,r0,r0,r0,r0,r0,r0,
  r1,r1,r0,r0,r0,r0,r0,r0,
  r1,r2,r1,r0,r0,r0,r0,r0,
  r1,r3,r3,r1,r0,r0,r0,r0,  
  r1,r4,r6,r4,r1,r0,r0,r0,
  r1,r5,r10,r10,r5,r1,r0,r0,
  r1,r6,r15,r20,r15,r6,r1,r0,
  r1,r7,r21,r35,r35,r21,r7,r1  };    

#define horner_b(_kind,_x,_t) { const I _n = table_bn[_kind]; cst_pR _y = &table_b[_kind*8]; { I _j; for (_x=_y[_n-1]*_t,_j=_n-2;_j>0;--_j) {_x+=_y[_j];_x*=_t;} _x+=_y[0]; } }
#define horner_a(_kind,_x,_t) { const I _n = table_a_degree[table_b_degree[_kind]]; cst_pR _y = &table_a[table_b_degree[_kind]*8]; { I _j; for (_x=_y[_n]*_t,_j=_n-1;_j>0;--_j) {_x+=_y[_j];_x*=_t;} _x+=_y[0]; } }

#define hornerdx_b(_kind,_x,_t) { const I _n = table_bn[_kind]; cst_pR _y = &table_b[_kind*8]; { I _j; for (_x=((R)(_n-1))*_y[_n-1]*_t,_j=_n-2;_j>1;--_j) {_x+=((R)_j)*_y[_j];_x*=_t;} _x+=_y[1]; } }

#define hornerdx_a(_kind,_x,_t) { const I _n = table_a_degree[table_b_degree[_kind]]; cst_pR _y = &table_a[table_b_degree[_kind]*8]; { I _j; for (_x=((R)_n)*_y[_n]*_t,_j=_n-1;_j>1;--_j) {_x+=((R)_j)*_y[_j];_x*=_t;} _x+=_y[1]; } }







void SmoothedDirac_eval	(pSmoothedDiracReadOnly	const	self_, 
			 const I 			rn_,
			 pR 			const 	r_,
			 pErr 			const 	err_)
{
  const R eps 	= self_->eps[0];
  const R meps 	= -eps;
  const R ieps 	= ((R)1.0)/eps;
  cst_eDirac kdirac = self_->dirac;
  switch(kdirac)
    {
    case __eDirac_m0p2:
    case __eDirac_m0p3:
    case __eDirac_m0p4:
    case __eDirac_m0p5:
    case __eDirac_m0p6:
    case __eDirac_m2p0:
    case __eDirac_m2p1:
    case __eDirac_m2p2:
    case __eDirac_m2p3:
    case __eDirac_m4p0:
    case __eDirac_m4p1:
    case __eDirac_m4p2:
    case __eDirac_m4p3:
      {
	{ I _jj;
	  for (_jj=0;_jj<rn_;++_jj)
	    {
	      const R t = r_[_jj];
	      if (t<=meps)
		{
		  r_[_jj]=(R)0.0;
		}
	      else if (t>=eps)
		{
		  r_[_jj]=(R)0.0;
		}
	      else
		{
		  R a,b,c,d;
		  const R t2=t*ieps;
		  hornerdx_a(kdirac,d,t2);
		  horner_b(kdirac,c,t2);
		  horner_a(kdirac,a,t2);
		  hornerdx_b(kdirac,b,t2);
		  r_[_jj] = (a*b+c*d)*ieps;
		}
	    } }
	break;
      }
    case __eDirac_sin:
      {
	const R pidiveps = MNS_PI*ieps;
	{ I _jj;
	  for (_jj=0;_jj<rn_;++_jj)
	    {
	      const R t = r_[_jj];
	      if (t<=meps)
		r_[_jj] = (R)0.0;
	      else if (t>=eps)
		r_[_jj] = (R)0.0;
	      else
		{
		  r_[_jj] = ieps*((R)0.5)*(((R)1.0)+nsCOS(pidiveps*t));
		}
	    } }
	break;
      }
    case __eDirac_ERROR:
    case __eDirac_ALL:
      {
	err_[0] = __eErr_switch;
	break;
      }
    }
}




#undef  r238
#undef  r325
#undef  r345
#undef  r256
#undef  r25 
#undef  r35 
#undef  r63 
#undef  r231
#undef  r843
#undef  r1617
#undef  r4830
#undef  r938 
#undef  r1212
#undef  r378 
#undef  r4165
#undef  r7651
#undef  r7938
#undef  r1218
#undef  r1024
#undef  r0   
#undef  r1   
#undef  r2   
#undef  r3   
#undef  r5   
#undef  r8   
#undef  r16  
#undef  r20  
#undef  r29  
#undef  r4   
#undef  r13  
#undef  r32  
#undef  r42  
#undef  r64  
#undef  r79  
#undef  r123 
#undef  r128 
#undef  r135 
#undef  r161 
#undef  r189 
#undef  r385 
#undef  r1220
#undef  r1653
#undef  r269 
#undef  r666 
#undef  r495 
#undef  r8031
#undef  r5005
#undef  r21100
#undef  r22565
#undef  r25480
#undef  r67795
#undef  r3410 
#undef  r12210
#undef  r2048 
#undef  r49140
#undef  r12285
#undef  r990 
#undef  r2230 
#undef  r1540 
#undef  r15015 
#undef  r162 
#undef  r2203 
#undef  r1_2  
#undef  r1_8  
#undef  r1_16 
#undef  r1_32  
#undef  r1_64  
#undef  r1_128 
#undef  r1_256 
#undef  r1_512 
#undef  r1_1024 
#undef  r1_2048 
#undef  r1_4096 
#undef  mr1_2   
#undef  r21     
#undef  r45     
#undef  mr45    
#undef  r50     
#undef  r9      
#undef  r6      
#undef  r7      
#undef  r10     
#undef  r15     
#undef  r105    
#undef  r175
#undef  r147
