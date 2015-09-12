CC    = g++
FLAGS = -mcmodel=medium -Wall -pedantic -O3 -ftree-vectorize -ftree-vectorizer-verbose=1  

all: seq par

seq:
	$(CC) $(FLAGS) Boost_dijkstra.cpp -o Boost_dijkstra

par:
	$(CC) $(FLAGS) -fopenmp Boost_dijkstra.cpp -o Boost_dijkstra_openmp

seq-gprof:
	$(CC) $(FLAGS) -pg Boost_dijkstra.cpp -o Boost_dijkstra_pg

clean:
	rm Boost_dijkstra Boost_dijkstra_openmp 
