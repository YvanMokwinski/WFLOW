#ifndef __header_Err_h__
#define __header_Err_h__
#include "Type.h"

typedef enum __eErr {  __eErr_no=0,/*!< value for no error */
		       __eErr_user,/*!< value for an user error */
		       __eErr_file,/*!< value for a file error */
		       __eErr_switch,/*!< value for a switch error */
		       __eErr_blank,/*!< value for a blank error */
		       __eErr_iwork,/*!< value for an iwork error */
		       __eErr_intern,/*!< value for an internal error */
		       __eErr_memory,/*!< value for a memory error */
		       __eErr_ALL /*!< total number of values for errors */ } Err;

typedef Err * RESTRICT pErr;



#ifdef __cplusplus

/**
   \brief Overload of increment operator for C++
*/
#if 0
Err& operator ++(Err& d_);
#endif
/**
   \brief Overload of output stream operator for C++
*/
std::ostream& operator<<(std::ostream&s_,const Err d_);

#endif

#ifdef __cplusplus
extern "C"
{
#endif

/**
   \brief Get string related to an error
   @param err_ error
   @return related string
*/
cst_pS eErr_get_string	(const Err 	err_);


/**
   \brief Get string related to an error
   @param err_ error
   @param file_ file
   @param line_ line number
*/
void   eErr_check	(const Err 	err_,
			 cst_pS 	file_,
			 const I 	line_);


#define  __mMNS_ERROR_CHECK(_err) eErr_check( (_err) , __FILE__ , __LINE__ )

#ifdef __cplusplus
}
#endif


#endif
