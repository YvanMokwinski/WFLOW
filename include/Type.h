#ifndef __header_Type_h__
#define __header_Type_h__

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <float.h>
#if 1
#ifdef __cplusplus

#include <cstdlib>				
#include <math.h>				
#include <cstdio>				
#include <cstring>				
#include <fstream>				
#include <iostream>				
#include <vector>				
#include <iomanip>
#include <memory>
using namespace std; 

#define MNS_THROW(_msg,_err)  { std::cerr << "// ERROR FILE: " << __FILE__ << std::endl <<  "//       LINE: "  << __LINE__ << std::endl << "//       MSG : " << (_msg) << std::endl;throw(_err); } (void)0


#endif
#endif
#define __MK_WITH_INTEL_LAPACK__ 1

#ifndef AND
#define AND &&
#else
#error AND is already defined
#endif

#ifndef OR
#define OR ||
#else
#error OR is already defined
#endif

#ifndef NOT
#define NOT !
#else
#error NOT is already defined
#endif


#ifndef  MAX
#define  MAX(a,b)           ((a) > (b) ? (a) : (b))
#else
#error MAX is already defined
#endif

#ifndef  MIN
#define  MIN(a,b)           ((a) < (b) ? (a) : (b))
#else
#error MIN already defined
#endif


/**
   \brief Maximum length of a character array
*/
#define __mSTR_maxlen 256
/**
   \brief Define a type for a character array
*/
typedef char STR[__mSTR_maxlen];

#ifndef __MNS_ILP__
#ifndef __MNS_LP__
#warning  __MNS_ILP__ and __MNS_LP__ undefined, set __MNS_LP__ by default
#define __MNS_ILP__
#endif
#endif

#ifndef __MNS_DLP__
#ifndef __MNS_DP__
#ifndef __MNS_SP__
#warning  __MNS_DLP__ and __MNS_DP__ undefined, set __MNS_DP__ by default
#define __MNS_DP__
#endif
#endif
#endif

#ifdef __MNS_ILP__
/**
   \brief Define a type for a long long integer
*/
typedef long long int I;
typedef unsigned long long int uI;
/**
   \brief Define the format to convert a long long integer
*/
#ifndef iformat
#define iformat "Ld"
#endif
#ifndef uiformat
#define uiformat "Lu"
#endif
/**
   \brief Define the (full) format to convert a long long integer
*/

#define nsFORMAT_INTEGER 	"Ld"
#define nsFORMAT_UINTEGER 	"Lu"
#define nsPRINT_SIZET 		"%Zu"

#ifndef ifmt
#define ifmt "%Ld"
#endif
#ifndef uifmt
#define uifmt "%Lu"
#endif
#endif

#ifdef __MNS_LP__
/**
   \brief Define a type for a long integer
*/
typedef long int I;
typedef long int uI;
/**
   \brief Define the format to convert a long  integer
*/
#ifndef iformat
#define iformat "ld"
#endif
#ifndef uiformat
#define uiformat "lu"
#endif
/**
   \brief Define the (full) format to convert a  long integer
*/
#define nsFORMAT_INTEGER 	"ld"
#define nsFORMAT_UINTEGER 	"lu"
#define nsPRINT_SIZET 		"%zu"
#ifndef ifmt
#define ifmt "%ld"
#endif
#ifndef uifmt
#define uifmt "%lu"
#endif
#endif


#ifdef __MNS_DLP__
/**
   \brief Define a type for a long double
*/
typedef long double R;

#define RMAX LDBL_MAX
#define RMIN LDBL_MIN

#define header(_name) q##_name##_

/**
   \brief Define the (full) format to convert a  long double
*/
#ifndef rfmt
#define rfmt "%8.17Le"
#endif
#define nsFORMAT_REAL "Le"


#define nsFLOOR 	floorl
#define nsROUND 	roundl
#define nsGAMMA 	tgammal
#define nsSQRT 		sqrtl
#define nsEXP 		expl
#define nsLOG 		logl
#define nsSIN 		sinl
#define nsCOS 		cosl
#define nsTAN 		tanl
#define nsSINH 		sinhl
#define nsCOSH 		coshl
#define nsTANH 		tanhl
#define nsASINH		asinhl
#define nsACOSH		acoshl
#define nsATANH		atanhl
#define nsASIN 		asinl
#define nsACOS 		acosl
#define nsATAN 		atanl
#define nsFABS 		fabsl
#define nsPOW 		powl


/**
   \brief Define the format to convert a  long double
*/
#ifndef rformat
#define rformat "Le"
#endif
/**
   \brief Overwrite the name of fabs function for long double
*/
#define nsFABS 	fabsl
/**
   \brief Overwrite the name of sqrt function for long double
*/
#define nsSQRT 	sqrtl
#endif


#ifdef __MNS_DP__
/**
   \brief Define a type for a double
*/
typedef double		R;

#define RMAX DBL_MAX
#define RMIN DBL_MIN

#define header(_name) d##_name##_
/**
   \brief Define the format to convert a double
*/
#ifndef rfmt
#define rfmt "%8.15e"
#endif
#define nsFORMAT_REAL "e"


#define nsGAMMA 	tgamma
#define nsFLOOR 	floor
#define nsROUND 	round
#define nsSQRT 		sqrt
#define nsEXP 		exp
#define nsLOG 		log
#define nsASINH		asinh
#define nsACOSH		acosh
#define nsATANH		atanh
#define nsSIN 		sin
#define nsCOS 		cos
#define nsTAN 		tan
#define nsSINH 		sinh
#define nsCOSH 		cosh
#define nsTANH 		tanh
#define nsASIN 		asin
#define nsACOS 		acos
#define nsATAN 		atan
#define nsFABS 		fabs
#define nsPOW 		pow

/**
   \brief Define the (full) format to convert a double
*/
#ifndef rformat
#define rformat "le"
#endif
/**
   \brief Overwrite the name of fabs function for double
*/
#define nsFABS 	fabs
/**
   \brief Overwrite the name of sqrt function for double
*/
#define nsSQRT 	sqrt
#endif


#ifdef __MNS_SP__
/**
   \brief Define a type for a float
*/
typedef float		R;

#define RMAX FLT_MAX
#define RMIN FLT_MIN

#define header(_name) s##_name##_
/**
   \brief Define the (full) format to convert a float
*/
#ifndef rfmt
#define rfmt "%8.5e"
#endif
#define nsFORMAT_REAL "e"
/**
   \brief Define the (full) format to convert a float
*/
#ifndef rformat
#define rformat "e"
#endif
/**
   \brief Overwrite the name of fabs function for float
*/
#define nsFABS 	fabsf
/**
   \brief Overwrite the name of sqrt function for float
*/
#define nsSQRT 	sqrtf
#endif



#ifndef RESTRICT 


/**
   \brief Macro for RESTRICT directive (c99)
*/
#define RESTRICT 

#ifdef __cplusplus
#undef  RESTRICT
#define RESTRICT __restrict__
#else

#ifdef __STDC_VERSION__

#if __STDC_VERSION__ >= 199901L
#undef  RESTRICT
#define RESTRICT restrict

#else
#warning MNS DO NOT APPLY RESTRICT DIRECTIVE
#endif

#else

#warning MNS DO NOT APPLY RESTRICT DIRECTIVE, __STDC_VERSION__ not defined

#endif
#endif


#endif

/**
   \brief Define a type for I pointers
*/
typedef 	I* RESTRICT  	pI;
/**
   \brief Define a type for uI pointers
*/
typedef 	uI* RESTRICT  	puI;
/**
   \brief Define a type for R pointers
*/
typedef 	R* RESTRICT  	pR;
/**
   \brief Define a type for constant I pointers
*/
typedef  const 	I* RESTRICT 	cst_pI;
/**
   \brief Define a type for constant uI pointers
*/
typedef  const 	uI* RESTRICT 	cst_puI;
/**
   \brief Define a type for constant R pointers
*/
typedef  const 	R* RESTRICT 	cst_pR;



/**
   \brief Define basic logical values 
*/
enum __emnsLOGIC {__emnsNO=0,
		  __emnsYES};
/**
   \brief Define a type for logic values
*/
typedef enum __emnsLOGIC L;
/**
   \brief Define a type for logic value pointer
*/
typedef L * RESTRICT pL;
/**
   \brief Define a type for constant logic value pointer
*/
typedef const L * RESTRICT cst_pL;


/**
   \brief Define a type for a character pointer
*/
typedef char * 	RESTRICT pS;
/**
   \brief Define a type for a constant character pointer
*/
typedef const char * RESTRICT cst_pS;


#ifndef NDEBUG
/**
   \brief Define a usefull macro for debug
*/
#define DebugVerif(_cond) if (!(_cond)) { fprintf(stderr,"condition failed line %d file %s\n",__LINE__,__FILE__);   }
#endif

#define NotYetImplemented fprintf(stderr,"piece of code not yet implemented line %d file %s\n",__LINE__,__FILE__);exit(1)


#define MNS_PI      		((R)3.141592653589793238512808959406e+00)
#define MNS_HALFPI  		((R)1.570796326794896619256404479703e+00)
#define MNS_TWICEPI 		((R)6.283185307179586477025617918812e+00)
#define MNS_PIINVERSE 		((R)3.183098861837906715381747713156e-01)
#define MNS_TWICEPIINVERSE 	((R)1.591549430918953357690873856578e-01)
#define MNS_INVLOG2             ((R)1.442695040888963407387651782798e+00)
#define MNS_SQRT2		((R)1.41421356237309504880168872420970)
#define MNS_HALFSQRT2 		((R)0.707106781186547524400844362104849)

typedef void * nsOBJ;
typedef int nsFLAG;
#define efmt "%d"



#define mStringVariable(_target,_source)				\
  { va_list _mSTR_GATHER_args;						\
    va_start (_mSTR_GATHER_args,_source);				\
    vsprintf(_target,_source,_mSTR_GATHER_args);			\
    va_end(_mSTR_GATHER_args); } (void)0


#endif
