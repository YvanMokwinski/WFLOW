#ifndef __header_Cmdline_h__
#define __header_Cmdline_h__

#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define __mCmdline_MAXARG 128

  /**
     \brief Struct to define a Cmdline
  */
  typedef struct
  {
    int 		argc;			/*!< number of arguments in command line */
    STR 		argv[__mCmdline_MAXARG];/*!< arguments in command line */
    int 		ref_argc;		/*!< reference number of argument in command line */
    cst_pS *  		ref_argv;		/*!< reference arguments in command line */
  } Cmdline,* RESTRICT pCmdline;
  
  typedef const Cmdline * RESTRICT cst_pCmdline;

  /**
     \brief Define the command line
     @param self_ self object
     @param argc_ nb arguments
     @param argv_ array of arguments
   */
  void 		Cmdline_def		(pCmdline 	self_,
					 const int	argc_,
					 cst_pS* 	argv_);

  /**
     \brief Display the command line
     @param self_ self object
   */
  void 		Cmdline_disp		(cst_pCmdline 	self_);

  /**
     \brief Free the Cmdline object
     @param self_ self object
  */
  void 		Cmdline_free		(pCmdline 	self_);

  /**
     \brief Get a logical option
     @param self_ self object
     @param opt_ name of the option
     @return __emnsNO if not found, __emnsYES otherwise
  */
  L 		Cmdline_get_logical	(pCmdline  	self_,					 
					 const STR 	opt_);

  /**
     \brief Get a string option
     @param self_ self object
     @param opt_ name of the option
     @param v_ value of the option
     @return __emnsNO if not found, __emnsYES otherwise
  */
  L 		Cmdline_get_string	(pCmdline  	self_,
					 const STR 	opt_,
					 STR 	       	v_);
  
  /**
     \brief Get a real option
     @param self_ self object
     @param opt_ name of the option
     @param v_ value of the option
     @return __emnsNO if not found, __emnsYES otherwise
  */
  L  		Cmdline_get_real	(pCmdline  	self_,
					 const STR 	opt_,
					 pR 		v_);
  
  /**
     \brief Get an integer option
     @param self_ self object
     @param opt_ name of the option
     @param v_ value of the option
     @return __emnsNO if not found, __emnsYES otherwise
  */
  L  		Cmdline_get_integer	(pCmdline  	self_,
					 const STR 	opt_,
					 pI 		v_);
  
  /**
     \brief Get the ith argument of the command line
     @param self_ self object
     @param ith_  ith argument
     @return the ith argument 
  */
  cst_pS 	Cmdline_get_arg		(pCmdline  	self_,
					 const int 	ith_);
  
  /**
     \brief Check invalid argument
     @param self_ self object
     @return 1 if an invalid argument is found, 0 otherwise
  */
  Err 		Cmdline_check_invalid	(cst_pCmdline  self_);

  /**
     \brief Is empty ?
     @param self_ self object
     @return __emnsNO if not empty, __emnsYES otherwise
  */
  L		Cmdline_isempty		(cst_pCmdline  self_);

  
#ifdef __cplusplus
}
#endif


#endif
