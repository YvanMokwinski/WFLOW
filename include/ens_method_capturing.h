#ifndef __HEADER_emkCAPTURING_H__
#define __HEADER_emkCAPTURING_H__

#if __mok_pure_enum__


enum __ens_method_capturing { __ens_method_capturing_error = 0,
			      __ens_method_capturing_linear,
			      __ens_method_capturing_c2,
			      __ens_method_capturing_n };

typedef enum __ens_method_capturing ens_method_capturing;

#else

#define __ens_method_capturing_error  	0
#define __ens_method_capturing_linear 	1
#define __ens_method_capturing_c2 	2
#define __ens_method_capturing_n  	3

typedef I ens_method_capturing;

#endif

typedef const ens_method_capturing cst_ens_method_capturing;
#if 0
enum __emk_capturing_iinfo { __emk_capturing_iinfo_error=0,
			     __emk_capturing_iinfo_lev,
			     __emk_capturing_iinfo_n };

enum __emk_capturing_rinfo { __emk_capturing_rinfo_error=0,
			     __emk_capturing_rinfo_iso,
			     __emk_capturing_rinfo_tol,
			     __emk_capturing_rinfo_epsmin,
			     __emk_capturing_rinfo_epsmax,
			     __emk_capturing_rinfo_area,
			     __emk_capturing_rinfo_n };
ens_method_capturing  		__ens_method_capturing_mpsint2enum	(cst_mpsint e_);
#endif

ens_method_capturing 		__ens_method_capturing_enum		(const char * name);
void 				f__ens_method_capturing_enum_		(const char * name,I*e_);
const char * 			__ens_method_capturing_string		(cst_ens_method_capturing flag_);
void 				__ens_method_capturing_strings_copy	(const char * names_[__ens_method_capturing_n]);

#endif

#ifdef __define__ens_method_capturing_strings__
/*#######################################################################################*/

#if 0
const char * __ens_method_capturing_strings[__ens_method_capturing_n] 	= {"__ens_method_capturing_error",
									   "c0",
									   "c2"};
enum __emk_capturing_iinfo 	__emk_capturing_iinfo_enum		(const char * name);
enum __emk_capturing_rinfo 	__emk_capturing_rinfo_enum		(const char * name);
const char * 			__emk_capturing_iinfo_string		(const enum __emk_capturing_iinfo flag_);
const char * 			__emk_capturing_rinfo_string		(const enum __emk_capturing_rinfo flag_);
const char * __emk_capturing_rinfo_strings[__emk_capturing_rinfo_n] 	= {"__emk_capturing_rinfo_error",
									   "capturing_tol",
									   "capturing_iso",
									   "capturing_epsmin",
									   "capturing_epsmax",
									   "capturing_area"};

const char * __emk_capturing_iinfo_strings[__emk_capturing_iinfo_n] 	= {"__emk_capturing_iinfo_error",
									   "capturing_lev"};
/*#######################################################################################*/
#endif
#endif
