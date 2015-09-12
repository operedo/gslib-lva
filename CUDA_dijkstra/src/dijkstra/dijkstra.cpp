/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

/* Template project which demonstrates the basics on how to setup a project 
* example application.
* Host code.
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <multithreading.h>
#include <sstream>

// includes, project
//#include <cutil_inline.h>
//#include "shrUtils.h"

// includes, kernels
#include "dijkstra_kernel.h"

///
//  Macro Options
//
//#define CITY_DATA

///
//  Some test data
//      http://en.literateprograms.org/Dijkstra%27s_algorithm_%28Scala%29

///
//  Generate a random graph
//
void generateRandomGraph(GraphData *graph, int numVertices, int neighborsPerVertex)
{
    graph->vertexCount = numVertices;
    graph->vertexArray = (int*) malloc(graph->vertexCount * sizeof(int));
    graph->edgeCount = numVertices * neighborsPerVertex;
    graph->edgeArray = (int*)malloc(graph->edgeCount * sizeof(int));
    graph->weightArray = (float*)malloc(graph->edgeCount * sizeof(float));

    for(int i = 0; i < graph->vertexCount; i++)
    {
        graph->vertexArray[i] = i * neighborsPerVertex;
    }

    for(int i = 0; i < graph->edgeCount; i++)
    {
        graph->edgeArray[i] = (rand() % graph->vertexCount);
        graph->weightArray[i] = (float)(rand() % 1000) / 1000.0f;
    }
}

//
// Generate an LVA graph from file
//
void generateLVAGraph(GraphData *graph, int numVertices, int neighborsPerVertex)
{

	FILE *fp;
	char buffer[20];
	fp = fopen("grid.out", "r");

	char numVerticesBuff[20], edgeCountBuff[20]; 
	int numVerticesAux, edgeCountAux; 
	int orig, dest;
	float cost; 

	graph->vertexCount = -1;
	graph->edgeCount = -1;


	int currVertex=0, pastVertex=0, accumEdges=0;
	int i;

	if(fp){
		//printf("Reading file grid.out...\n");
		if(fscanf(fp, "%d %d", &numVerticesAux, &edgeCountAux)!=EOF){
			//sscanf(buffer,"%[^ ] %[^ \n]\n",numVerticesBuff, edgeCountBuff);

			//printf("%d %d\n",numVerticesAux,edgeCountAux);
			graph->vertexCount = numVerticesAux;
			graph->edgeCount = edgeCountAux;

			if(numVerticesAux>0 && edgeCountAux>0){
				graph->vertexArray = (int*) malloc(graph->vertexCount * sizeof(int));
				graph->edgeArray = (int*)malloc(graph->edgeCount * sizeof(int));
				graph->weightArray = (float*)malloc(graph->edgeCount * sizeof(float));
			}
			else{
				printf("grid.out must have numVertices>0 and edgeCount>0.\n");
				graph->vertexCount = -1;
				graph->edgeCount = -1;
				return;
			}
		}
		else{
			printf("grid.out must have as first line 'numVertices edgeCount'.\n");
			return;
		}

		i=0;

                graph->vertexArray[0]=0;
		//printf("v[%d]=%d\n",i,graph->vertexArray[i]);

		while(fscanf(fp, "%d %d %f", &orig, &dest, &cost)!=EOF){
			//printf("%d %d %f\n",orig,dest,cost);
			

			if(orig>numVerticesAux || orig<1 || dest>numVerticesAux || dest<1){
				printf("grid.out must have all vertices between 1 and %d.\n",numVerticesAux);
				graph->vertexCount = -1;
				graph->edgeCount = -1;
				return;
			}
			edgeCountAux--;
			if(edgeCountAux<0){
				printf("grid.out must have exactly %d edges.\n",edgeCountAux);
				graph->vertexCount = -1;
				graph->edgeCount = -1;
				return;
			}

			++accumEdges;
			currVertex=orig-1;
			//if(currVertex==0){
			//	pastVertex = currVertex;
			//}
			//if(pastVertex!=currVertex && i>0){
			//	graph->vertexArray[currVertex] = 0;
			//}

			graph->edgeArray[i] = dest-1;
			graph->weightArray[i] = cost;
				
			//graph->vertexArray[currVertex] = graph->vertexArray[currVertex] + 1 + (currVertex>0 ? graph->vertexArray[currVertex-1] : 0);
			
			if(pastVertex!=currVertex && i>0){
				//graph->vertexArray[pastVertex] = accumEdges-1;
				graph->vertexArray[currVertex] = accumEdges-1 + graph->vertexArray[pastVertex];
				//printf("v[%d]=%d\n",currVertex,graph->vertexArray[currVertex]);
				accumEdges=1;
			}

			pastVertex=currVertex;

			i++;
		}
		if(edgeCountAux!=0){
			printf("grid.out must have exactly %d edges.\n",edgeCountAux);
			graph->vertexCount = -1;
			graph->edgeCount = -1;
			return;
		}
		//graph->vertexArray[pastVertex] = accumEdges;
		//printf("v[%d]=%d\n",pastVertex,graph->vertexArray[pastVertex]);

	}
	else{
		printf("grid.out not exists!.\n");
		return;
	}

/*    
    for(int i = 0; i < graph->vertexCount; i++)
    {
        printf("v[%d]=%d\n",i,graph->vertexArray[i]);
    }

    for(int i = 0; i < graph->edgeCount; i++)
    {
        printf("e[%d]=%d\n",i,graph->edgeArray[i]);
        printf("w[%d]=%f\n",i,graph->weightArray[i]);
    }
    */

	fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) 
{
//@TEMP - why do I need this for link to work on Linux?
//cutWaitForThreads(NULL,0);
//@TEMP
	// use command-line specified CUDA device, otherwise use device with highest Gflops/s
	//if( cutCheckCmdLineFlag(argc, (const char**)argv, "device") )
	//	cutilDeviceInit(argc, argv);
	//else
	//	cudaSetDevice( cutGetMaxGflopsDeviceId() );

    bool doGPU;
    bool doRef;
    bool doMultiGPU;
    int generateVerts;
    int generateEdgesPerVert;
    int numSources;
    int sourceId;

    //parseCommandLineArgs(argc, (const char**)argv, doGPU, doMultiGPU, doRef, &numSources, &generateVerts, &generateEdgesPerVert);

    doGPU = 1;
    doMultiGPU = 0;
    doRef = 0;

    // Allocate memory for arrays
    GraphData graph;

    //generateRandomGraph(&graph, generateVerts, generateEdgesPerVert);
    generateLVAGraph(&graph, generateVerts, generateEdgesPerVert);

    //printf("Vertex Count: %d\n", graph.vertexCount);
    //printf("Edge Count: %d\n", graph.edgeCount);

    std::vector<int> sourceVertices;

    //for(int source = 0; source < numSources; source++)
    //{
    //    sourceVertices.push_back(source % graph.vertexCount);
    //}


	FILE *fp;
	fp = fopen("nodes2cal.out", "r");

	if(fp){
		//printf("Reading file nodes2cal.out...\n");
		if(fscanf(fp, "%d", &numSources)!=EOF){
			if(numSources<1){
				printf("nodes2cal.out must have numSources>0.\n");
				return 1;
			}
		}
		else{
			printf("nodes2cal.out must have as first line 'numSources'.\n");
			return 1;
		}
		while(fscanf(fp, "%d", &sourceId)!=EOF){
			sourceVertices.push_back(sourceId-1);
		}

	}
	else{
		printf("nodes2cal.out not exists!.\n");
		return 1;
	}



    int *sourceVertArray = (int*) malloc(sizeof(int) * sourceVertices.size());
    std::copy(sourceVertices.begin(), sourceVertices.end(), sourceVertArray);

    float *results = (float*) malloc(sizeof(float) * sourceVertices.size() * graph.vertexCount);


    unsigned int gpuTimer = 0;
    //cutilCheckError(cutCreateTimer(&gpuTimer));
    //cutilCheckError(cutStartTimer(gpuTimer));

    // Run Dijkstra's algorithm
    if ( doGPU )
    {
        runDijkstra(&graph, sourceVertArray, results, sourceVertices.size() );
    }

    //cutilCheckError(cutStopTimer(gpuTimer));


    unsigned int multiGPUTimer = 0;
    //cutilCheckError(cutCreateTimer(&multiGPUTimer));
    //cutilCheckError(cutStartTimer(multiGPUTimer));

    if ( doMultiGPU )
    {
        runDijkstraMultiGPU(&graph, sourceVertArray, results, sourceVertices.size() );
    }

    //cutilCheckError(cutStopTimer(multiGPUTimer));

    unsigned int refTimer = 0;
    //cutilCheckError(cutCreateTimer(&refTimer));
    //cutilCheckError(cutStartTimer(refTimer));

    if ( doRef )
    {
        runDijkstraRef(&graph, sourceVertArray, results, sourceVertices.size() );
    }

    //cutilCheckError(cutStopTimer(refTimer));

    for (unsigned int i = 0; i < sourceVertices.size(); i++)
    {
        for (int j = 0; j < graph.vertexCount; j++)
        {
            //if (i != j)
            //{
                //printf("%d --> %d: %f\n", sourceVertArray[i], j, results[i * graph.vertexCount + j] );
                printf("%f\n", results[i * graph.vertexCount + j] );
            //}
        }
    }

    free(sourceVertArray);
    free(results);

    //cudaThreadExit();

    return 0;
}

