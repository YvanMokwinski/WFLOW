#ifndef __header_Monitor_h__
#define __header_Monitor_h__

#ifdef __cplusplus
extern "C"
{
#endif
  
  /**
     \brief Define mode for monitoring   
     @note an example should be a=MonitorMode_STD | MonitorMode_LOG
  */
  enum __eMonitorMode { 
    MonitorMode_OFF		= 0,/*!< default logical value*/
    MonitorMode_STD		= 1,/*!< mode to display messages on standard output*/
    MonitorMode_LOG		= 2,/*!< mode to display messages on log output*/
    MonitorMode_EXITIFERROR	= 4 /*!< mode to exit if an error message is received */
  };

  /**
     \brief Type MonitorMode enumeration
  */
  typedef enum __eMonitorMode MonitorMode;


  /**
     \brief Shutdown the monitor
     @note it will close all log files.
  */
  void 	Monitor_shutdown(void);


  /**
     \brief Define a monitor
     @param ithread_  thread number
     @param progname_ name of the program
     @param mode_ mode of the monitor
  */
  void 	Monitor_def	(const int 		ithread_,
			 const char *  		progname_,
			 const unsigned char 	mode_);


  /**
     \brief Send a message to the monitor
     @param ithread_  thread number
     @param msg_ message to send
     @note additional arguments are similar to printf
  */
  void 	Monitor_msg	(const int 		ithread_,
			 const char *  		msg_,
			 ...);

  /**
     \brief Send a warning to the monitor
     @param ithread_  thread number
     @param msg_ message to send
     @note additional arguments are similar to printf
  */
  void 	Monitor_warn	(const int		ithread_,
			 const char *  		msg_,
			 ...);


  /**
     \brief Send an error message to the monitor
     @param ithread_  thread number
     @param msg_ message to send
     @note additional arguments are similar to printf
  */
  void 	Monitor_errmsg		(const int 		ithread_,
				 const char * 		msg_,
				 ...);
  
  void 	Monitor_createDirectory	(const char *	name_,
				 ...);
  

#ifdef __cplusplus
}
#endif

#endif
