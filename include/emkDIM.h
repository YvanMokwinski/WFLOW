#ifndef __HEADER_EMK_DIM_H__
#define __HEADER_EMK_DIM_H__

#if __mok_pure_enum__

enum __emkDIM 	{   __emkDIM_error = 0,   
		    __emkDIM_1,
		    __emkDIM_2,
		    __emkDIM_3,
		    __emkDIM_n };

typedef enum __emkDIM emkDIM;

#else

#define __emkDIM_error 0
#define __emkDIM_1 1
#define __emkDIM_2 2
#define __emkDIM_3 3
#define __emkDIM_n 4

typedef nsINT emkDIM;

#endif

typedef const emkDIM cst_emkDIM;

emk_err emkDIM_fwrite	(FILE*fich_,cst_emkDIM *dim_);
emk_err emkDIM_fread	(FILE*fich_,emkDIM * dim_);

#endif
