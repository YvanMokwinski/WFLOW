
/*! \mainpage MNS (Multiphase Navier Stokes), a C/C++/Fortran library
 * 
 * \section intro_sec Introduction
 *
 * This library is a part of my Ph.D work. 
 * It is able to perform numerical simulation of multiphase flow using finite element methods.
 *
 * \section install_sec Installation
 *
 * \subsection step0 Step 0: type 'make' to get help
 * \subsection step1 Step 1: make sure libConfig/Config.h is properly configured
 * \subsection step2 Step 2: type 'make dep' to re-generate dependencies
 * \subsection step3 Step 3: type 'make all' to compile both optimized and debug version
 *  
 * etc...
 */
#define YES 1
#define NO  0

/**
   \brief Macro to specify the path of the installed version of MNS
*/
#define MNS_PATH       			"/home/yvan/MNS"
/**
   \brief Macro to specify the version 
*/
#define MNS_VERSION    			"0.1"
/**
   \brief Macro to specify the date
*/
#define MNS_DATE       			"2013"
/**
   \brief Macro to specify the company
*/
#define MNS_COMPANY    			"Ecole Polytechnique de Montreal, Canada"
/**
   \brief Macro to force enum types to be replaced by integer types
*/
#define MNS_FORCE_FORTRAN_COMPATIBILITY NO
/**
   \brief Macro to specify the use of SSE INSTRUCTION (must have an Intel compiler)
*/
#define MNS_SSE_INSTRUCTION 		NO
/**
   \brief Macro to specify the use of the Intel MKL
*/
#define MNS_WITH_INTEL_LAPACK 		YES
/**
   \brief Macro to specify the maximum number of threads
*/
#define MNS_MAX_THREAD 			8



#define __MNS_WITH_SPARSEKIT__ 		YES
#define __MNS_WITH_PARDISO__ 		YES

/**
   \brief Macro to specify how to type enum values
   @note set to YES for optimized version (NDEBUG)
*/
#ifndef NDEBUG
#define MNS_DEFINE_ENUM_AS_INTEGER MNS_FORCE_FORTRAN_COMPATIBILITY
#else
#define MNS_DEFINE_ENUM_AS_INTEGER YES
#endif

