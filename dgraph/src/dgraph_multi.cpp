// Copyright (C) 2004-2006 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine

// Example usage of dijkstra_shortest_paths algorithm

// Enable PBGL interfaces to BGL algorithms
#include <boost/graph/use_mpi.hpp>

// Communication via MPI
#include <boost/graph/distributed/mpi_process_group.hpp>

// Dijkstra's single-source shortest paths algorithm
//#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/distributed/dijkstra_shortest_paths.hpp>

// Distributed adjacency list
#include <boost/graph/distributed/adjacency_list.hpp>

// METIS Input
#include <boost/graph/metis.hpp>

// Graphviz Output
#include <boost/graph/distributed/graphviz.hpp>

// Standard Library includes
#include <fstream>
#include <string>
#include <iostream> 
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>     
#include <stdexcept>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <sys/time.h>





#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace std; 
using namespace boost;
using boost::graph::distributed::mpi_process_group;
//using boost::graph::parallel::synchronize;

/* An undirected, weighted graph with distance values stored on the
   vertices. */
typedef adjacency_list<
//			vecS, 
//			distributedS<mpi_process_group, vecS>, 
//			undirectedS,
//                       /*Vertex properties=*/property<vertex_distance_t, float>,
//                       /*Edge properties=*/property<edge_weight_t, float>
			listS, 
			distributedS<mpi_process_group, vecS>, 
			directedS,
                       /*Vertex properties=*/no_property,//property<vertex_distance_t, float>,
                       /*Edge properties=*/property<edge_weight_t, float> 
			>
  			graph_t;

typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
typedef std::pair<int, int> Edge;

//typedef property_map < graph_t, vertex_index_t >::type IndexMap;
//typedef iterator_property_map < float*, IndexMap, float, float& > DistanceMap;

// Define a distance map
typedef property_map < graph_t, vertex_index_t>::const_type IndexMap;  

//typedef local_property_map<mpi_process_group,IndexMap,vertex_index_t> distribIndexMap;
typedef iterator_property_map < std::vector<float>::iterator, IndexMap> DistanceMap;
//typedef iterator_property_map < float*, IndexMap, float, float& > DistanceMap;


int main(int argc, char* argv[])
{
	struct timeval start_time, end_time;
	time_t diff_time;

	int runId=atoi(argv[1]);
	int numRuns=atoi(argv[2]);

  	boost::mpi::environment env(argc,argv);
  	boost::mpi::communicator world;


  	long int num_nodes;
  	long int num_edges;
  	char gfile[20] ; //graph data file
  	char nfile[20] ; //nodes to cal dist for data file
  	char  buffer[128] ;
  	std::istringstream instream ;

  	//string gfile_str = "gird 1.6m nodes.out";
  	//string gfile_str = "girdsmall.out";
  	std::string gfile_str = "grid.out";
  	strcpy(gfile,gfile_str.c_str());

  	std::ifstream infile( gfile, ios::in );
	//if (process_id(process_group(g)) == 0){
	//if (world.rank() == 0){
  		infile.getline(buffer, 128);
  		instream.clear() ;
  		instream.str(buffer) ;
  		instream >> num_nodes >> num_edges;
	//}
	//boost::mpi::broadcast(world, num_nodes ,0);
	//boost::mpi::broadcast(world, num_edges ,0);
	//world.barrier();


  	Edge* edge_array;   // Pointer to edge, initialize to nothing.
  	edge_array = new Edge[num_edges];  // Allocate n ints and save ptr in a.
	//float weighttmp;
  	float* weights;   // Pointer to edge, initialize to nothing.
  	weights = new float[num_edges];  // Allocate n ints and save ptr in a.
  	int num_arcs = num_edges;
  	int n1,n2;

	//if (process_id(process_group(g)) == 0){
	if (world.rank() == 0){
  		std::cerr << "READING graph into C++ program from file grid.out";  std::cerr << std::endl;std::cerr << std::endl;  std::cerr << std::endl;
  		std::cerr << "Found ";std::cerr << num_nodes; std::cerr << " nodes, and "; std::cerr << num_edges; std::cerr << " edges";std::cerr  << std::endl; 
	}
	//}

  	for (int i=0; i<num_edges; ++i) 
  	{ 
      		//get a line from the file
      		if (num_edges/2 == i || (num_edges+1)/2 == i )
      		{
			//if (process_id(process_group(g)) == 0)
			if (world.rank() == 0)
          			std::cerr << "half way done reading the graph ... please wait a bit longer ...";  std::cerr << std::endl; 
      		}
		//if (process_id(process_group(g)) == 0){
		//if (world.rank() == 0){
      			infile.getline(buffer, 128);
      			instream.clear() ;
      			instream.str(buffer) ;
      			instream >> n1 >> n2 >> weights[i];
      			//instream >> n1 >> n2 >> weighttmp;
		//}
		//boost::mpi::broadcast(world, n1 ,0);
		//boost::mpi::broadcast(world, n2 ,0);
		//boost::mpi::broadcast(world, weighttmp ,0);
		//world.barrier();

		//weights[i]=weighttmp;
      		edge_array[i]=Edge(n1-1,n2-1);
  	}


	//if (process_id(process_group(g)) == 0){
	if (world.rank() == 0){
  		std::cerr << "graph read in"; std::cerr << std::endl;  std::cerr << std::endl; 
  		std::cerr << "setting up dijkstra call for your graph"; std::cerr << std::endl; 
	}

  	graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);



  	// Get vertex 0 in the graph
  	//graph_traits<Graph>::vertex_descriptor start = vertex(0, g);
  	//graph_traits<Graph>::vertex_descriptor start = vertex(1, g);

  	//if (process_id(process_group(g)) == 0)
	if (world.rank() == 0)
    		std::cerr << "Graph allocated." << std::endl;


  	int nodes2cal,cur_node;
  	std::string nfile_str = "nodes2cal.out";
  	strcpy(nfile, nfile_str.c_str());
  	std::ifstream infile2( nfile, ios::in ); 
  	//if (process_id(process_group(g)) == 0){
	if (world.rank() == 0){
  		infile2.getline(buffer, 128);
  		instream.clear() ;
  		instream.str(buffer) ;
  		instream >> nodes2cal;
	}
	boost::mpi::broadcast(world, nodes2cal ,0);
	world.barrier();

	if (world.rank() == 0){
  		std::cerr << "working out distances to ";std::cerr << nodes2cal; std::cerr << " nodes"; std::cerr << std::endl; 
  		std::cerr << "NODE:"; std::cerr << std::endl; 
	}

  	//now need to be able to write out the distances

  	ofstream outfile;
	//if (process_id(process_group(g)) == 0){
	//if (world.rank() == 0){
  	//	//outfile.open ("dist_cpp.out");
  	//	//outfile.close();
	//}

	//typedef property_map < graph_t, vertex_index_t >::type IndexMap;
	//typedef iterator_property_map < float*, IndexMap, float, float& > DistanceMap;

	std::vector<vertex_descriptor> p(num_vertices(g));
	//std::vector<int> d(num_vertices(g));
	std::vector<float> d(num_vertices(g));


	int num_vertices_per_process=num_vertices(g);
	int sum;

	std::cerr << "Process " << process_id(process_group(g)) << ": " << num_vertices_per_process << " nodes" << std::endl;

	boost::mpi::all_reduce(world, num_vertices_per_process, sum, std::plus<int>());
	world.barrier();

	//std::cout << "Sum of vertices=" << sum << std::endl;

	//for (int i = 0; i < sum; ++i) {
	//	//std::cout << 0 << "\t" << i << "\t" << get(distance, vertex(i, g)) << std::endl;  // return value unused
	//	std::cout << 1 << "\t" << (i+1) << "\t" << d[i] << std::endl;  // return value unused
	//}
	////synchronize(distance);

	DistanceMap distance;

	diff_time=0;
	double elapsed;

	int landmarkId[nodes2cal];
	if (world.rank() == 0){
		for (int i=0; i<nodes2cal; ++i) 
		{
			std::cerr << i + 1 << flush; 
			std::cerr << std::endl; 
		
			infile2.getline(buffer, 128);
			instream.clear() ;
			instream.str(buffer) ;
			instream >> cur_node;
			landmarkId[i]=cur_node;
		}
	}
	boost::mpi::broadcast(world, landmarkId, nodes2cal ,0);
	world.barrier();
	

	int blockSize=(nodes2cal-numRuns+1)/numRuns;
	int iIni=runId*blockSize;
	int iFin=(runId+1)*blockSize;
	if(runId==numRuns-1){
		iFin=nodes2cal;
	}

	//for (int i=0; i<nodes2cal; ++i) 
	for (int i=iIni; i<iFin; ++i) 
	{
		//if (process_id(process_group(g)) == 0){
		//if (world.rank() == 0){
		//	std::cerr << i + 1 << flush; 
		//	std::cerr << std::endl; 
		//
		//	infile2.getline(buffer, 128);
		//	instream.clear() ;
		//	instream.str(buffer) ;
		//	instream >> cur_node;
		//}
		//boost::mpi::broadcast(world, cur_node ,0);
		cur_node=landmarkId[i];
		//world.barrier();
		vertex_descriptor s = vertex(cur_node-1, g);  //node number output from fortran program will always be +1 (C++ starts at 0 not 1)
		world.barrier();
	    	std::cerr << "Process " << process_id(process_group(g)) << ": " << "starting Dijkstra for node " << cur_node << std::endl;


		//synchronize(process_group(g));

		//start_time = time(NULL);
		gettimeofday(&start_time, NULL);
	  	dijkstra_shortest_paths(g, s,  
	       		predecessor_map(
	        	 	make_iterator_property_map(p.begin(), get(vertex_index, g))).
	       		distance_map(
	        		make_iterator_property_map(d.begin(), get(vertex_index, g))).
	  		weight_map(get(edge_weight, g)).
	  		vertex_index_map(get(vertex_index, g)).
			distance_compare(std::less<float>()).
			distance_combine(closed_plus<float>()).
			distance_inf((std::numeric_limits<float>::max)()).
			distance_zero(0).
			visitor(default_dijkstra_visitor())
		);
		synchronize(process_group(g));

	    	std::cerr << "Process " << process_id(process_group(g)) << ": " <<  "Dijkstra completed." << std::endl;
		//end_time=time(NULL);
		gettimeofday(&end_time, NULL);
		elapsed = ((end_time.tv_sec - start_time.tv_sec) * 1000) + (end_time.tv_usec / 1000 - start_time.tv_usec / 1000);
		std::cerr << "time dijkstra("<< world.rank() <<"):  " << (elapsed/1000.0) <<  " s" << std::endl ;
	
	  	//if (process_id(process_group(g)) == 0){
		//synchronize(distance);
            
            
		int outflag=0;
  		//ofstream outfile;
		
		//start_time=time(NULL);
		gettimeofday(&start_time, NULL);


		char filename[60];		

		for(int procId=0;procId<world.size();procId++){
			world.barrier();
			if (world.rank() == procId){

				sprintf(filename,"dist_cpp.out_%05d_%05d_%05d" ,i,procId,runId);
				// Descomentar esto para guardar resultados en disco
				// Comentar esto para hacer pruebas de escalabilidad de Boost, sin escritura
				
				outfile.open (filename);
				for (int j=0; j<num_vertices_per_process; ++j) 
				{
					//std::cout << d[j] <<  std::endl ; 
					outfile << d[j] <<  std::endl ; 
				} 
				outfile.close();
 			} 
			world.barrier();
			usleep(100);
		}
		gettimeofday(&end_time, NULL);
		elapsed = ((end_time.tv_sec - start_time.tv_sec) * 1000) + (end_time.tv_usec / 1000 - start_time.tv_usec / 1000);
		std::cerr << "time writing("<< world.rank() <<"):  " << (elapsed/1000.0) <<  " s" << std::endl ;

		world.barrier();
	
	}

	//if (process_id(process_group(g)) == 0){
	if (world.rank() == 0){
  		std::cerr << "done dijkstra" << std::endl; 
  		//outfile.close();
		char filenameEnd[60];		
		sprintf(filenameEnd,"lock.%05d" ,runId);
		outfile.open (filenameEnd);
		outfile.close();
	}


  	delete(edge_array) ;
  	delete(weights) ;

  	return 0;
}

