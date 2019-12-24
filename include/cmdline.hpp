#pragma once 
#include <string.h>  
#include <stdio.h>  
//!
//!   @brief Class to define a command line.
//!
class cmdline
{
  
public:
    
  //!
  //! @brief Define the command line
  //! @param self_ self object
  //! @param argc_ nb arguments
  //! @param argv_ array of arguments
  //! 
  cmdline(int	argc_,
	  char** 	argv_);
  
  //!  
  //! @brief Free the cmdline object
  //! @param self_ self object
  //!  
  ~cmdline();

    
  //!  
  //! @brief Display the command line
  //! @param self_ self object
  //!   
  void disp() const;

  //!  
  //! @brief Get a logical option
  //! @param self_ self object
  //! @param opt_ name of the option
  //! @return false if not found, true otherwise
  //!  
  bool	get_logical(const char 	*opt_);

  //!  
  //! @brief Get a string option
  //! @param self_ self object
  //! @param opt_ name of the option
  //! @param v_ value of the option
  //! @return false if not found, true otherwise
  //!  
  bool get_string	(const char* 	opt_,
			 char*   	v_);
    
  //!  
  //! @brief Get a real option
  //! @param self_ self object
  //! @param opt_ name of the option
  //! @param v_ value of the option
  //! @return false if not found, true otherwise
  //!  
  template<typename real_t>
  bool  		get_real	(const char* 	opt_,
					 real_t*	v_);
    
  //!  
  //! @brief Get an integer option
  //! @param self_ self object
  //! @param opt_ name of the option
  //! @param v_ value of the option
  //! @return false if not found, true otherwise
  //!
  template<typename int_t>
  bool  		get_integer	(const char* 	opt_,
					 int_t*		v_);
    
  //!  
  //! @brief Get the ith argument of the command line
  //! @param self_ self object
  //! @param ith_  ith argument
  //! @return the ith argument 
  //!  
  const char* 		get_arg		(const int 	ith_) const;
    
  //!  
  //! @brief Check invalid argument
  //! @param self_ self object
  //! @return 1 if an invalid argument is found, 0 otherwise
  //!  
  bool 			check_invalid() const;
    
  //!  
  //! @brief Is empty ?
  //! @param self_ self object
  //! @return false if not empty, true otherwise
  //!  
  bool			isempty() const;

private:
  
  static constexpr int MAXLEN = 256;
  static constexpr int MAXARG = 128;

  int 		argc{};			//!< number of arguments in command line 
  char 		argv[MAXARG][MAXLEN]; //!< arguments in command line 
  int 		ref_argc{};		//!< reference number of argument in command line 
  char **	ref_argv{};		//!< reference arguments in command line
  
  int search(int 		argc_,
	     char argv_[cmdline::MAXARG][cmdline::MAXLEN] ,
	     const char* 	opt_,
	     const int 	nchoice_ = 0,
	     const char** opt_choice_ = NULL);
  
};
