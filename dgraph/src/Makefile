CC	= g++
#SRC	= Boost_dijkstra.cpp
#OUT	= Boost_dijkstra
SRC	= dgraph.cpp
OUT	= dgraph
SRCMULTI = dgraph_multi.cpp
OUTMULTI = dgraph_multi
#SRC	= test.cpp
#OUT	= test
#SRC     = distributed_shortest_paths_test.cpp
#OUT     = distributed_shortest_paths_test
#BOOST_LIBS = -L/home/apps/lib/boost/1.57
#BOOST_HEADERS = -I/home/operedo/boost_1_57_0
BOOST_LIBS = -L/gpfs/apps/MN3/BOOST/1.56.0/OPENMPI/lib
BOOST_HEADERS = -I/gpfs/apps/MN3/BOOST/1.56.0/OPENMPI/include
FLAG= -mcmodel=medium -Wall -pedantic -O3 -mtune=native -march=native -m64 -funroll-all-loops -ftree-vectorize -ftree-vectorizer-verbose=1 -DNDEBUG

MPICC   = mpic++
ifeq ($(shell uname), Darwin)
LIBS    = -lboost_mpi-mt -lboost_serialization-mt -lboost_graph_parallel-mt -lboost_system-mt -lboost_filesystem-mt -lpthread
else
LIBS    = -lboost_mpi -lboost_serialization -lboost_graph_parallel -lboost_system -lboost_filesystem -lpthread
endif

all: mpi mpimulti

mpi:
	$(MPICC) $(BOOST_HEADERS) $(FLAG) $(SRC) -o $(OUT) $(BOOST_LIBS) $(LIBS)
	cp $(OUT) ../bin

mpimulti:
	$(MPICC) $(BOOST_HEADERS) $(FLAG) $(SRCMULTI) -o $(OUTMULTI) $(BOOST_LIBS) $(LIBS)
	cp $(OUTMULTI) ../bin

clean:
	rm $(OUT) $(OUTMULTI) ../bin/$(OUT) ../bin/$(OUTMULTI)
