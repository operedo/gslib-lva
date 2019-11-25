
#FC=gfortran

#FFLAGS= -fpp -DDEBUG -mkl -openmp -O3 -mtune=native -march=native -opt-report=3 -vec-report=3 
#CCFLAGS = -g -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-vectorizer-verbose=2 
#FFLAGS= -g -fpp -mkl -openmp -O3 -mtune=native -march=native -opt-report=3 -vec-report=3 

#FFLAGS= -g -traceback -fpp -mkl -openmp -O3 -mtune=native -march=native 
FFLAGS= -fpp -mkl -openmp -O3 -mtune=native -march=native 
#CCFLAGS = -dumpmachine -g -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-vectorizer-verbose=2 
#CCFLAGS = -g -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize  
CCFLAGS = -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize  

#LDFLAGS = -lstdc++ 
LDFLAGS = -lstdc++ -lgomp
#LDFLAGS= -lstdc++ -L/Soft/gcc/6.2.0/lib64/ -lgomp
#LDFLAGS = -lstdc++ /Soft/intelcc/11.1-0.72/lib/intel64/libiomp5.so
#LINKS=/Soft/intelf/composer_xe_2013.3.163/compiler/lib/ia32/libiomp5.a
#LINKS=/Soft/intelf/composer_xe_2013.3.163/compiler/lib/intel64/libiomp5.a
LINKS=/Soft/intelf/composer_xe_2013.3.163/compiler/lib/intel64/libiomp5.a $(LDFLAGS)
#LINKS=$(LDFLAGS)

#FFLAGS=  -O3 -xSSE4.2
#FC=ifort
FC=/Soft/intelf/composer_xe_2013.3.163/bin/intel64/ifort
#FC=/Soft/intelf/composer_xe_2013.3.163/bin/ia32/ifort
#FC=/Soft/gcc/6.2.0/bin/gfortran
#F90=quicksort.f90  glb.f90  grid_info.f90 class_laf.f90 graph_subs.f90 sqdist.f90 GAMV_LVA.f90 Scaling\ -\ ISOMAP_landmark.f90
CXX=/Soft/gcc/6.2.0/bin/g++
#CXX=g++
#CXX=/Soft/intelcc/11.1-0.72/bin/intel64/icpc

KIND_OBJS=  
F90_OBJS= quicksort.o gauinv.o random.o glb.o  grid_info.o aniso_d.o graph_subs.o kdtree3.o cova_1D.o sortem.o ktsol.o exhaustive_srch.o sgs_lva.o Scaling-ISOMAP_landmark.o
F90_LEVELS_OBJS= quicksort.o gauinv.o random.o glb.o  grid_info.o aniso_d.o graph_subs.o kdtree3.o cova_1D.o sortem.o ktsol.o exhaustive_srch.o sgs_lva_levels.o Scaling-ISOMAP_landmark.o
F90_LEVELS_SPEC_OBJS= quicksort.o gauinv.o random.o glb.o  grid_info.o aniso_d.o graph_subs.o kdtree3.o cova_1D.o sortem.o ktsol.o exhaustive_srch.o sgs_lva_levels_speculative.o Scaling-ISOMAP_landmark.o
F90_OPTLAND_OBJS= quicksort.o gauinv.o random.o glb.o  grid_info.o aniso_d.o graph_subs.o kdtree3.o cova_1D.o sortem.o ktsol.o exhaustive_srch.o sgs_lva_optlandmarks.o Scaling-ISOMAP_landmark.o
#FOR=chknam.for  getindx.for
FOR_OBJS=chknam.o  getindx.o
#LAPACK=LAPACK/eig.f90
LAPACK_OBJS=LAPACK/eig.o LAPACK/solve.o
#LAPACK_OBJS=LAPACK/eig.o LAPACK/solve.o \
#	LAPACK/BLAS/daxpy.o \
#	LAPACK/BLAS/dcopy.o \
#	LAPACK/BLAS/ddot.o \
#	LAPACK/BLAS/dgemm.o \
#	LAPACK/BLAS/dgemv.o \
#	LAPACK/BLAS/dger.o \
#	LAPACK/BLAS/dnrm2.o \
#	LAPACK/BLAS/dscal.o \
#	LAPACK/BLAS/dswap.o \
#	LAPACK/BLAS/dsymv.o \
#	LAPACK/BLAS/dsyr2.o \
#	LAPACK/BLAS/dsyr2k.o \
#	LAPACK/BLAS/dtrmm.o \
#	LAPACK/BLAS/dtrmv.o \
#	LAPACK/BLAS/dtrsm.o \
#	LAPACK/BLAS/idamax.o \
#	LAPACK/BLAS/lsame.o \
#	LAPACK/BLAS/xerbla.o \
#	LAPACK/DOUBLE/dgesv.o \
#	LAPACK/DOUBLE/dgetf2.o \
#	LAPACK/DOUBLE/dgetrf.o \
#	LAPACK/DOUBLE/dgetrs.o \
#	LAPACK/DOUBLE/dlae2.o \
#	LAPACK/DOUBLE/dlaev2.o \
#	LAPACK/DOUBLE/dlanst.o \
#	LAPACK/DOUBLE/dlansy.o \
#	LAPACK/DOUBLE/dlapy2.o \
#	LAPACK/DOUBLE/dlarfb.o \
#	LAPACK/DOUBLE/dlarf.o \
#	LAPACK/DOUBLE/dlarfg.o \
#	LAPACK/DOUBLE/dlarft.o \
#	LAPACK/DOUBLE/dlartg.o \
#	LAPACK/DOUBLE/dlascl.o \
#	LAPACK/DOUBLE/dlaset.o \
#	LAPACK/DOUBLE/dlasr.o \
#	LAPACK/DOUBLE/dlasrt.o \
#	LAPACK/DOUBLE/dlassq.o \
#	LAPACK/DOUBLE/dlaswp.o \
#	LAPACK/DOUBLE/dlatrd.o \
#	LAPACK/DOUBLE/dorg2l.o \
#	LAPACK/DOUBLE/dorg2r.o \
#	LAPACK/DOUBLE/dorgql.o \
#	LAPACK/DOUBLE/dorgqr.o \
#	LAPACK/DOUBLE/dorgtr.o \
#	LAPACK/DOUBLE/dsteqr.o \
#	LAPACK/DOUBLE/dsterf.o \
#	LAPACK/DOUBLE/dsyev.o \
#	LAPACK/DOUBLE/dsytd2.o \
#	LAPACK/DOUBLE/dsytrd.o \
#	LAPACK/DOUBLE/ieeeck.o \
#	LAPACK/DOUBLE/ilaenv.o \
#	LAPACK/DOUBLE/iparmq.o \
#	LAPACK/INSTALL/dlamch.o 
##	LAPACK/invert_routines/sgetf2.o \
##	LAPACK/invert_routines/sgetrf.o \
##	LAPACK/invert_routines/sgetri.o \
##	LAPACK/invert_routines/slaswp.o \
##	LAPACK/invert_routines/strti2.o \
##	LAPACK/invert_routines/strtri.o

OBJS=   eig.o solve.o
#OBJS=   eig.o solve.o \
#	dlamch.o \
#	dgesv.o \
#	dgetf2.o \
#	dgetrf.o \
#	dgetrs.o \
#	dlae2.o \
#	dlaev2.o \
#	dlanst.o \
#	dlansy.o \
#	dlapy2.o \
#	dlarfb.o \
#	dlarf.o \
#	dlarfg.o \
#	dlarft.o \
#	dlartg.o \
#	dlascl.o \
#	dlaset.o \
#	dlasr.o \
#	dlasrt.o \
#	dlassq.o \
#	dlaswp.o \
#	dlatrd.o \
#	dorg2l.o \
#	dorg2r.o \
#	dorgql.o \
#	dorgqr.o \
#	dorgtr.o \
#	dsteqr.o \
#	dsterf.o \
#	dsyev.o \
#	dsytd2.o \
#	dsytrd.o \
#	ieeeck.o \
#	ilaenv.o \
#	iparmq.o \
#	daxpy.o \
#	dcopy.o \
#	ddot.o \
#	dgemm.o \
#	dgemv.o \
#	dger.o \
#	dnrm2.o \
#	dscal.o \
#	dswap.o \
#	dsymv.o \
#	dsyr2.o \
#	dsyr2k.o \
#	dtrmm.o \
#	dtrmv.o \
#	dtrsm.o \
#	idamax.o \
#	lsame.o \
#	xerbla.o 
##	sgetf2.o \
##	sgetrf.o \
##	sgetri.o \
##	slaswp.o \
##	strti2.o \
##	strtri.o 

CPP_OBJS= boostdijkstra.o boostdijkstra_capi.o Boostdijkstra.o

ALL_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_OBJS)
ALL_LEVELS_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_LEVELS_OBJS)
ALL_LEVELS_SPEC_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_LEVELS_SPEC_OBJS)
ALL_OPTLAND_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_OPTLAND_OBJS)
OUT=SGS_LVA
OUT_LEVELS=SGS_LVA_levels_opt01
OUT_LEVELS_SPEC=SGS_LVA_levels_speculative
OUT_OPTLAND=SGS_LVA_optland

#Scaling-ISOMAP_landmark.o : boostdijkstra_mod.o

default: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_OBJS) -o $(OUT)
	cp $(OUT) ../bin/
#	cp $(OUT) ../test/swiss-roll/

levels: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_LEVELS_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_LEVELS_OBJS) -o $(OUT_LEVELS)
	cp $(OUT_LEVELS) ../bin/
	cp $(OUT_LEVELS) ../test/swiss-roll/

levels-spec: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_LEVELS_SPEC_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_LEVELS_SPEC_OBJS) -o $(OUT_LEVELS_SPEC)
	cp $(OUT_LEVELS_SPEC) ../bin/
	cp $(OUT_LEVELS_SPEC) ../test/swiss-roll/


optland: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_OPTLAND_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_OPTLAND_OBJS) -o $(OUT_OPTLAND)
	cp $(OUT_OPTLAND) ../bin/

clean:
	rm *.mod *.o $(OUT) ../bin/$(OUT)

.SUFFIXES: .o .F90 .f90 .for .f .fpp .cpp

.F90.o : ; $(FC) -c $(FFLAGS) $*.F90

.f90.o : ; $(FC) -c $(FFLAGS) $*.f90

.for.o : ; $(FC) -c $(FFLAGS) $*.for

.f.o : ; $(FC) -c $(FFLAGS) $*.f

.fpp.o : ; $(FC) -c $(FFLAGS) $*.fpp

.cpp.o : ; ${CXX} -c ${CCFLAGS} $*.cpp


