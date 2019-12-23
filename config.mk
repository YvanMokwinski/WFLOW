SHMMAX=`cat /proc/sys/kernel/shmmax`
SHMMAX=1240000

PLATFORM=Linux
# 1/CHOOSE YOUR FORTRAN COMPILER
F77=g77
# 2/CHOOSE FFLAGS
FFLAGS= -i8   -pipe -m64 $(OPTIMIZE) -fno-second-underscore  -funroll-loops -funroll-all-loops
# 3/CHOOSE YOUR C COMPILER
CC=gcc
# 4/CHOOSE CFLAG
CFLAGS=-std=c99  -mtune=native  -pipe -m64 $(OPTIMIZE)   -Wall -funroll-loops -funroll-all-loops 
# 5/GIVE BLAS HEADERS
HEADERBLAS= -Iblas/ -D__MNS_WITH_MKL__=0
# 6/GIVE BLAS LIBRARY
LIBBLAS=-Lblas/ -lblas


F77=ifort
CC=icc
CFLAGS=  -pipe -m64  -std=c99  -Wall  -funroll-loops -funroll-all-loops   $(OPTIMIZE)
FFLAGS=-i8  -override-limits -funroll-loops -funroll-all-loops -m64 $(OPTIMIZE)
HEADERBLAS= -I/home/p901195/ifort/mkl/include -D__MNS_WITH_MKL__=1
LIBBLAS= -L/home/p901195/ifort/lib/intel64  -liomp5 -L/home/p901195/ifort/mkl/include/../lib/em64t  -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_intel_thread -lmkl_core    -lpthread #-lmkl_solver_ilp64
LD=$(F77) -nofor-main







F77=gfortran -fdefault-integer-8
CC=gcc
CPP=g++
CFLAGS= -std=c99   -mtune=native -pipe -m64 $(OPTIMIZE)   -Wall  -funroll-loops -funroll-all-loops   #-std=gnu99      -ansi
CPPFLAGS=    -mtune=native -pipe -m64 $(OPTIMIZE)   -Wall  -funroll-loops -funroll-all-loops   #-std=gnu99      -ansi
FFLAGS=    -pipe -m64 $(OPTIMIZE) -fno-second-underscore  -funroll-loops -funroll-all-loops
HEADERBLAS= -I/home/p901195/ifort/mkl/include -D__MNS_WITH_MKL__=1
HEADERBLAS= #-Iblas/ -D__MNS_WITH_MKL__=0
LIBBLAS=-Lblas/ -lblas
LD=$(F77)




# 1/CHOOSE YOUR FORTRAN COMPILER
F77=ifort
# 2/CHOOSE FFLAGS
FFLAGS= -i8  -override-limits -funroll-loops -funroll-all-loops -m64 $(OPTIMIZE)
# 3/CHOOSE YOUR C COMPILER
CC=icc
CPP=icpc
# 4/CHOOSE CFLAG
CFLAGS=   -pipe -m64  -std=c99    -Wall  -funroll-loops -funroll-all-loops   $(OPTIMIZE)
CPPFLAGS=  -pipe -m64 -Wall  -funroll-loops -funroll-all-loops   $(OPTIMIZE)
# 5/GIVE BLAS HEADERS
HEADERBLAS= -I/home/p901195/ifort/mkl/include -D__MNS_WITH_MKL__=1
# 6/GIVE BLAS LIBRARY
LIBBLAS= -L/home/p901195/ifort/lib/intel64  -liomp5 -L/home/p901195/ifort/mkl/include/../lib/em64t  -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_intel_thread -lmkl_core    -lpthread #-lmkl_solver_ilp64
LD=$(F77) -nofor-main
LDFLAGS=   -funroll-loops -funroll-all-loops   $(OPTIMIZE)





##############################################################################################################
# 1/CHOOSE YOUR FORTRAN COMPILER
F77=gfortran
# 2/CHOOSE FFLAGS
FFLAGS= -fdefault-integer-8 #-i8  -override-limits -funroll-loops -funroll-all-loops -m64 $(OPTIMIZE)
# 3/CHOOSE YOUR C COMPILER
CC=gcc
CPP=g++
# 4/CHOOSE CFLAG 
CFLAGS=    -pipe -m64 -std=c99      -Wall  -funroll-loops -funroll-all-loops   $(OPTIMIZE) #-Werror
CPPFLAGS= -std=c++11 -pipe -m64              -Wall  -funroll-loops -funroll-all-loops   $(OPTIMIZE) #-Werror
# 5/GIVE BLAS HEADERS
HEADERBLAS= -I/home/p901195/ifort/mkl/include -D__MNS_WITH_MKL__=1
# 6/GIVE BLAS LIBRARY
LIBBLAS= -L/home/p901195/ifort/lib/intel64  -liomp5 -L/home/p901195/ifort/mkl/include/../lib/em64t  -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core   -lpthread
LD=g++

#$(F77) #-nofor-main  -lmkl_solver_ilp64
LDFLAGS=   -funroll-loops -funroll-all-loops   $(OPTIMIZE)
##############################################################################################################






#-D__MNS_SHMMAX=$(SHMMAX)
#-D__MNS_SHMMAX=$(SHMMAX)
CFLAGS+= -D__MNS_DP__ -D__MNS_ILP__ -D__MNS_WITH_METIS__
CPPFLAGS+=  -D__MNS_DP__ -D__MNS_ILP__ -D__MNS_WITH_METIS__ 
