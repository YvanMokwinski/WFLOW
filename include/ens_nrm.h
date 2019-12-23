#ifndef __HEADER_MOK_ENUM_NRM_H__
#define __HEADER_MOK_ENUM_NRM_H__
#if __mok_pure_enum__

enum __emk_nrm 	{ __emk_nrm_error=0,
		  __emk_nrm_L1,
		  __emk_nrm_L2,
		  __emk_nrm_LInf,
		  __emk_nrm_H1,
		  __emk_nrm_H2,	
		  __emk_nrm_L1_jump,
		  __emk_nrm_L2_jump,
		  __emk_nrm_LInf_jump,		 
		  __emk_nrm_n};

typedef enum __emk_nrm emk_nrm;

#else

#define  __emk_nrm_error 	0
#define	 __emk_nrm_L1 		1
#define  __emk_nrm_L2 		2
#define	 __emk_nrm_LInf 	3
#define	 __emk_nrm_H1 		4
#define  __emk_nrm_H2 		5
#define	 __emk_nrm_L1_jump 	6
#define	 __emk_nrm_L2_jump 	7
#define  __emk_nrm_LInf_jump 	8	 
#define	 __emk_nrm_n 		9
typedef I  emk_nrm;

#endif

typedef const emk_nrm cst_emk_nrm;

const char * 		__emk_nrm_string		(cst_emk_nrm    i_);
emk_nrm 		__emk_nrm_enum			(const char* name_);
const char * 		__emk_nrm_string_latex		(cst_emk_nrm i);
void 			__emk_nrm_strings_copy		(const char * names_[__emk_nrm_n]);
emk_nrm 		__emk_nrm_next			(cst_emk_nrm i);

#endif
