#ifndef __header_ens_method_transfert_h__
#define __header_ens_method_transfert_h__

enum __ens_method_transfert {__ens_method_transfert_error = 0,
			     __ens_method_transfert_static,			/** \brief Methode par simple evaluation */
			     __ens_method_transfert_L2,			/** \brief Methode par projection L2 globale */
			     __ens_method_transfert_annreinitialization,	/** \brief Methode par reinitialisation  */
			     __ens_method_transfert_n};


enum __ens_method_transfert 	__ens_method_transfert_enum		(const char *  name);
const char * 			__ens_method_transfert_string		(const enum __ens_method_transfert flag_);
void 				__ens_method_transfert_strings_copy	(const char * names_[__ens_method_transfert_n]);


#endif
