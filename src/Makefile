#
# A simple makefile to build the program"
# written by Yvan Mokwinski, 13/08/2013"
#
VPATH=$(DISTRIB)/src

#HEADERMETIS=-I$(DISTRIB)/../Kernel/libExtern/metis/Lib/
#LIBMETIS=  -L$(DISTRIB)/../Kernel/libExtern/metis/ -lmetis
#LIBSPARSEKIT2=  -L$(DISTRIB)/../../../Kernel/libExtern/sparsekit2/ -lskit

#-I$(DISTRIB)/../include/
#-I$(DISTRIB)/extern -I$(DISTRIB)/extern/WLS/include
HEADERS=-I$(DISTRIB)/../  -I$(DISTRIB)/../WLS/include  -I$(DISTRIB)/include/  $(HEADERMETIS) $(HEADERBLAS) 

LIBS=$(LIBBLAS) #-L$(DISTRIB)/../lib -lMns_$(CC)$(DEBUG) 

include $(DISTRIB)/config.mk

TRASH=	DgTransport_routine.o \
	DgTransport_test.o \
	DgTransport_init.o \
	DgTransport_correction.o \
	DgTransport_compute.o \


#	SparseIterative.o Cmdline.o \

OBJ=	cmdline.o mkQ.o \
	DgTransport_shape.o \
	LinsysVectors.o \
	Sparse.o \
	SparseFactorization.o \
	SparseBlock.o \
	SparseStokes.o  \
	ExternSparsekit.o \
	ExternPardiso.o \
	Medit.o \
	eKindEquation.o \
	eLinearSolver.o \
	eTransientMethod.o \
	eTensionMethod.o \
	eFileSuffix.o \
	Points.o Matrix.o TransientScheme.o  extern_medit.o  \
	\
	SmoothedHeaviside.o \
	SmoothedDirac.o \
	\
	Monitor.o \
	Time.o \
	\
	ens_rinfo.o \
	ens_linfo.o \
	ens_sinfo.o \
	ens_iinfo.o \
	ensBASIS.o \
	ensBASIS_EDGE_ORTHO.o \
	ensBASIS_TRIA_ORTHO.o \
	ensBASIS_TRIA_LAGR.o \
	\
	eVariable.o \
	\
	\
	Workelm_init.o \
	\
	Err.o \
	\
	Fem.o \
	\
	GlobalReadOnly_get_Workelm.o \
	GlobalReadOnly_get_Time.o \
	GlobalReadOnly_get_Fem.o \
	GlobalReadOnly_get_CsfInfo.o \
	\
	\
	Global_get_Workelm.o \
	Global_get_Time.o \
	Global_get_Fem.o \
	Global_get_CsfInfo.o \
	\
	Global_initialize.o \
	Global_free.o \
	\
	nsGLOBAL_transient_print.o \
	nsGLOBAL_change_mesh.o \
	nsGLOBAL_print.o \
	\
	nsSPACE.o \
	nsSPACE_PRODUCT.o \
	MetricTensor.o \
	ns_mesh.o \
	ns_var.o  \
	ns_basename.o \
	Parameters.o \
	Global_main.o \
	mkmake_ns.o \
	libns_exact.o \
	ensBASIS_EDGE_LAGR.o \


#	ensBASIS_TETRA_ORTHO.o \
#	ensBASIS_TETRA_LAGR.o \

# REMOVED FILES 
# 	Global_precompute.o \
#	Global_newton.o \
#	Global_run.o \
#	nsPBLM_canal.o \
#	nsPBLM.o \
	Linsys.o \
#	nsGLOBAL_init.o \
# GlobalReadOnly_get_Variables.o \
#	Global_get_Variables.o
#	VariablesReadOnly_get.o \
#	\
#	Variables_clear.o \
#	Variables_def.o \
#	Variables_get.o \
#	Variables_set.o \
#



LIB=libMnsStokes.a

all: MnsStokes.exe MnsDg.exe FixedPointMesh.exe TestFiniteElementShape.exe TestDgTransport.exe mkpartition.exe TestDgTransport.version2.exe TestDgTransport.version3.exe

MnsStokes.exe:MnsStokes.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/MnsStokes.exe $(DISTRIB)/bin/MnsStokes$(DEBUG)

FixedPointMesh.exe:FixedPointMesh.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/FixedPointMesh.exe $(DISTRIB)/bin/FixedPointMesh$(DEBUG)

TestFiniteElementShape.exe:TestFiniteElementShape.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/TestFiniteElementShape.exe $(DISTRIB)/bin/TestFiniteElementShape$(DEBUG)


MnsDg.exe:MnsDg.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/MnsDg.exe $(DISTRIB)/bin/MnsDg$(DEBUG)

TestDgTransport.exe:TestDgTransport.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/TestDgTransport.exe $(DISTRIB)/bin/TestDgTransport$(DEBUG)

TestDgTransport.version2.exe:TestDgTransport.version2.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/TestDgTransport.version2.exe $(DISTRIB)/bin/TestDgTransport.version2$(DEBUG)

TestDgTransport.version3.exe:TestDgTransport.version3.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/TestDgTransport.version3.exe $(DISTRIB)/bin/TestDgTransport.version3$(DEBUG)

mkpartition.exe:mkpartition.o $(OBJ) libMnsStokes.a
	$(LD) $(CFLAGS) $(FLAGS) $< $(OBJ) -o $@ $(LIBS) $(LIBMETIS) $(LIBSPARSEKIT2) 
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/mkpartition.exe $(DISTRIB)/bin/mkpartition$(DEBUG)




%.o:%.c
	$(CC) $(CFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

%.o:%.cpp
	$(CC) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

%.o:%.f
	$(F77) $(FFLAGS) $(FLAGS) -c $<  -o $@

libMnsStokes.a:$(OBJ)
	\rm -f $@
	ar crsu $@ $(OBJ)
	ranlib $@
	echo "'$(PLATFORM)_$(CC)'"
	ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/$@  $(DISTRIB)/lib/libMnsStokes_$(CC)$(DEBUG).a
