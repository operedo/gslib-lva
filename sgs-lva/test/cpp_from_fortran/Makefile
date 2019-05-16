FC = /Soft/intelf/composer_xe_2013.3.163/bin/intel64/ifort 
#FC = gfortran
#CXX = /Soft/intelcc/11.1-0.72/bin/intel64/icpc 
CXX = g++

UNAME := $(shell uname -s)

FCFLAGS = -g -cpp -Wall -Wextra
CCFLAGS = -g -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-vectorizer-verbose=1 
#CCFLAGS = -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp  -mcmodel=medium -Wall -pedantic -O3 -ftree-vectorize -ftree-vectorizer-verbose=1
#CCFLAGS = -fopenmp  -O3 -ftree-vectorize -ftree-vectorizer-verbose=1
#ifeq ($(UNAME_S),Darwin)
LDFLAGS = -lstdc++ -lgomp
#else
#	LDFLAGS = -lc++
#endif

all: test.x
test.o : foo_mod.o

%.x : %.o foo_mod.o foo_capi.o Foo.o
	${FC} $^ -o $@ ${LDFLAGS}

%.o : %.f90
	${FC} ${FCFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $^ -o $@

.PHONY : clean

clean :
	${RM} -rf *.o *.mod test.x
