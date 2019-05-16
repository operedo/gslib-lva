#include "Foo.hpp"

#include <stdexcept>
#include <iostream> 
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each

#include <boost/graph/adjacency_list.hpp> 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>     
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std; 
using namespace boost;



//Foo::Foo(int _a, int _b): a(_a), b(_b){
//    cout << "C++ side, constructor" << endl;
//}
//
//Foo::~Foo(){
//    cout << "C++ side, destructor" << endl;
//}
//
//int Foo::bar(int c) const{
//    return a + c;
//}
//
//double Foo::baz(double d) const{
//    return d + b;
//}

void foo_speaker(string s){
//    Foo f(4, 2);
//    cout << s << " Foo(4, 2).bar(3) is: " <<  f.bar(3) << endl;
     cout << s << endl;
}

void dijkstra2_cpp(int *n, int* a) 
{ 
    cout << "Begin dijkstra2 " << *n << " " << a[0] << " " << a[1] << endl;
//     start the timer
}



void dijkstra(string s) 
{ 
    cout << "Begin dijkstra" << endl;
//     start the timer
    time_t start_time, end_time;
    time_t diff_time;


    //typedef adjacency_list < listS, vecS, directedS, no_property, property < edge_weight_t, float > > graph_t;
    typedef adjacency_list < listS, vecS, undirectedS, no_property, property < edge_weight_t, float > > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
    typedef std::pair<int, int> Edge;

    long int num_nodes;
    long int num_edges;
    char gfile[20] ; //graph data file
    char nfile[20] ; //nodes to cal dist for data file
    char  buffer[128] ;
    istringstream instream ;

    //string gfile_str = "gird 1.6m nodes.out";
    //string gfile_str = "girdsmall.out";
    //string gfile_str = "grid.out";
    string gfile_str = "grid.out";
    strcpy(gfile,gfile_str.c_str());

    ifstream infile( gfile, ios::in );
    infile.getline(buffer, 128);
    instream.clear() ;
    instream.str(buffer) ;
    instream >> num_nodes >> num_edges;

    Edge* edge_array;   // Pointer to edge, initialize to nothing.
    edge_array = new Edge[num_edges];  // Allocate n ints and save ptr in a.
    float* weights;   // Pointer to edge, initialize to nothing.
    weights = new float[num_edges];  // Allocate n ints and save ptr in a.
    int num_arcs = num_edges;
    int n1,n2;

    cout << "READING graph into C++ program from file grid.out";  cout << endl;cout << endl;  cout << endl;
    cout << "Found ";cout << num_nodes; cout << " nodes, and "; cout << num_edges; cout << " edges";cout  << endl; 

    for (int i=0; i<num_edges; ++i) 
    { 
        //get a line from the file
        //if (num_edges/2 == i || (num_edges+1)/2 == i )
        //{
        //    cout << "half way done reading the graph ... please wait a bit longer ...";  cout << endl; 
        //}
        infile.getline(buffer, 128);
        instream.clear() ;
        instream.str(buffer) ;
        instream >> n1 >> n2 >> weights[i];

        edge_array[i]=Edge(n1-1,n2-1);
    }

//     _strdate( dateStr1);
//     _strtime( timeStr1 );

    cout << "graph read in"; cout << endl;  cout << endl; 
    cout << "setting up dijkstra call for your graph"; cout << endl; 

    graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);

    cout << "graph_t g created"; cout << endl; 

    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);

    cout << "property_map weightmap created"; cout << endl; 

    std::vector<vertex_descriptor> p(num_vertices(g));

    cout << "vector<vertex_descriptor> p created"; cout << endl; 

    std::vector<float> d(num_vertices(g));

    cout << "vector<float> d created"; cout << endl; 
    
    property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);

    cout << "property_map indexmap created"; cout << endl; 

//need to find out what distances are required?
//in file 'node2cal.out'
//first line is number of nodes to cal
//all other lines contain the index of the node to calculate
//remember in c++ the ind = indFORTRAN-1

    int nodes2cal,cur_node;
    string nfile_str = "nodes2cal.out";
    strcpy(nfile, nfile_str.c_str());
    ifstream infile2( nfile, ios::in ); 
    infile2.getline(buffer, 128);
    instream.clear() ;
    instream.str(buffer) ;
    instream >> nodes2cal;
    cout << "working out distances to ";cout << nodes2cal; cout << " nodes"; cout << endl; 
    cout << "NODE:"; cout << endl; 


#ifndef _OPENMP
    //now need to be able to write out the distances
    ofstream outfile;
    outfile.open ("dist_cpp.out");
    diff_time=0;
    for (int i=0; i<nodes2cal; ++i) 
    //for (int i=0; i<4; ++i) 
    {
//        cout << i + 1 << flush; 
//        cout << endl; 

        infile2.getline(buffer, 128);
        instream.clear() ;
        instream.str(buffer) ;
        instream >> cur_node;

        vertex_descriptor s = vertex(cur_node-1, g);  //node number output from fortran program will always be +1 (C++ starts at 0 not 1)
         start_time = time(NULL);

//        dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap, 
//                          std::less<float>(), closed_plus<float>(), 
//                          (std::numeric_limits<float>::max)(), 0,
//                          default_dijkstra_visitor());


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



        //write these distances out
	// Descomentar para escribir resultado
	// Comentar para medir speedup y escalabilidad
	
        for (int i=0; i<num_nodes; ++i) 
        {
            outfile << d[i] <<  endl ; 
        }  
	

         end_time=time(NULL);
         diff_time=diff_time + end_time-start_time;
    }

    cout << "time for one path:  " << float(diff_time)/float(nodes2cal) <<  " s" << endl ;
    cout << "done dijkstra"; 
    outfile.close();

#else
    //now need to be able to write out the distances

    int num_threads=1;
    int this_thread=0;
#pragma omp parallel
{
    num_threads = omp_get_num_threads();
}

    //cout << num_threads << endl; 

    int cur_node_array[nodes2cal];
    string stringout = "dist_cpp.out"; 

    for (int i=0; i<nodes2cal; ++i) 
    {
        infile2.getline(buffer, 128);
        instream.clear() ;
        instream.str(buffer) ;
        instream >> cur_node_array[i];
    }

    diff_time=0;
    start_time = time(NULL);

#pragma omp parallel default(none) firstprivate(this_thread,stringout,cur_node,p,d) \
                                   shared(num_nodes,num_threads,nodes2cal,cout,cur_node_array,g)
{
    this_thread = omp_get_thread_num();
    char num2str[21];

//    if(num_threads>1000){
//        cout << "Too many threads (max=1000). Exit" << endl;
//        return 1;
//    }

    if(this_thread<10)
    	sprintf(num2str,"_00%d",this_thread);
    else if(this_thread<100)
    	sprintf(num2str,"_0%d",this_thread);
    else if(this_thread<1000)
    	sprintf(num2str,"_%d",this_thread);

    string stringoutLocal = stringout + string(num2str);

    //cout << this_thread << ":" << stringoutLocal << endl;

    ofstream outfile;
    outfile.open(stringoutLocal.c_str());

#pragma omp for nowait
    for (int i=0; i<nodes2cal; ++i) 
    {

//#pragma omp critical 
//{
//        cout << this_thread << " : "<< i + 1 << flush; 
//        cout << endl; 
//}

        cur_node=cur_node_array[i];

        vertex_descriptor s = vertex(cur_node-1, g);  //node number output from fortran program will always be +1 (C++ starts at 0 not 1)

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



        //write these distances out
        // Descomentar para escribir resultado
        // Comentar para medir speedup y escalabilidad
        
        for (int j=0; j<num_nodes; ++j) 
        {
            outfile << d[j] <<  endl ; 
        }  
        

    }


    outfile.close();

}

    end_time=time(NULL);
    diff_time=diff_time + end_time-start_time;
    cout << "time for one path:  " << float(diff_time)/float(nodes2cal) <<  " s" << endl ;
    cout << "done dijkstra"; 


#endif

    delete(edge_array) ;
    delete(weights) ;

//    return 0; 
} 

