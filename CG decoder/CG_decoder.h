#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <climits>
#include <chrono>
#include <random>
#include "blossom5/PMexpand.cpp"
#include "blossom5/PMinterface.cpp"
#include "blossom5/PMmain.cpp"
#include "blossom5/PMduals.cpp"
#include "blossom5/PMshrink.cpp"
#include "blossom5/PMinit.cpp"
#include "blossom5/MinCost/MinCost.cpp"
#include "min_heap.h"

using namespace std;

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// Structure to save the coordinates of a block
typedef struct block_node {

	int t0, t1, vertex;

} block_node;

// Structure to save the coordinates of an anyon block
typedef struct syndrome_node {

	int row, column, t0, t1, vertex;

} syndrome_node;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function for a biased coin
bool random_error (double error) {

	if ( ((double)rand() / ((double)(RAND_MAX) + (double)(1))) < error)
		return true;
	else
		return false;
}
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////

// Function that outputs the probability of successful decoding for a given set of inputs
double CG_decoder (int L, double physical_error, double measu_error, double synchronicity, int time, int precision) {

	int success = 0;

	int i, j, t, m;
	int x0, y0, t0, x1, y1, t1, v;
	double p;
	bool syndrome0, syndrome1, flag;

	int n_blocks, n_syn;
	int current, block_size;
	bool error_horizontal, error_vertical;

	bool lattice_H[L][L][time];
	bool lattice_V[L][L][time];

	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			lattice_H[i][j][0] = false;
			lattice_V[i][j][0] = false;
		}
	}

	for (int n = 0; n < precision; n++) {

///////////////////////////////////////// CREATE LATTICE

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {
				for (t = 1; t < time; t++) {
					lattice_H[i][j][t] = random_error(physical_error) ^ lattice_H[i][j][t-1];
					lattice_V[i][j][t] = random_error(physical_error) ^ lattice_V[i][j][t-1];
				}
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET SYNDROME

		n_blocks = 0;

		vector<syndrome_node> syndrome;
		vector<vector<vector<block_node>>> block (L, vector<vector<block_node>>(L));
	    
		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {

				syndrome0 = false;
				t0 = 0;

				for (t = 1; t < time - 1; t++) {
					if (random_error(synchronicity)) {	// Check if stabilizer measurement is successful

						block[i][j].push_back(block_node({t0,t,n_blocks}));	// All blocks are numbered by their order of appearance = n_blocks

						syndrome1 = lattice_H[i][j][t] ^ lattice_H[i][(j-1+L)%L][t] ^ lattice_V[i][j][t] ^ lattice_V[(i-1+L)%L][j][t] ^ random_error(measu_error);
						if (syndrome0 != syndrome1)	// Different consecutive measurements lead to an anyon block
							syndrome.push_back(syndrome_node({i,j,t0,t,n_blocks}));
					
						t0 = t;
						syndrome0 = syndrome1;
						n_blocks++;
					}
				}

				block[i][j].push_back(block_node({t0,time-1,n_blocks}));	// All blocks are numbered by their order of appearance = n_blocks

				// Last measurement round with perfect stabilizer measurement
				syndrome1 = lattice_H[i][j][time-1] ^ lattice_H[i][(j-1+L)%L][time-1] ^ lattice_V[i][j][time-1] ^ lattice_V[(i-1+L)%L][j][time-1];
				if (syndrome0 != syndrome1) 	// Different consecutive measurements lead to an anyon block
					syndrome.push_back(syndrome_node({i,j,t0,time-1,n_blocks}));
			
				n_blocks++;
			}
	   	}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// BUILD CONTRACTED GRAPH

		vector<unordered_set<int>> block_graph(n_blocks);
		vector<vector<double>> weight_graph(n_blocks, vector<double>(n_blocks));

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {

				t0 = block[i][j][0].t0;
				t1 = block[i][j][0].t1;
				v  = block[i][j][0].vertex;


				for (m = 0; block[(i+1+L)%L][j][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block
	
				// Connect all adjacent blocks overlapping in time with the current block
				while (m < block[(i+1+L)%L][j].size() and block[(i+1+L)%L][j][m].t0 < t1) {
					block_graph[v].insert(block[(i+1+L)%L][j][m].vertex);
					block_graph[block[(i+1+L)%L][j][m].vertex].insert(v);

					p = 0.5*(1 - pow(1 - 2*physical_error, MIN(block[(i+1+L)%L][j][m].t1,t1) - MAX(block[(i+1+L)%L][j][m].t0,t0)));	// Contracted edge's new weight
					weight_graph[v][block[(i+1+L)%L][j][m].vertex] = log((1-p)/p);
					weight_graph[block[(i+1+L)%L][j][m].vertex][v] = log((1-p)/p);

					m++;
				}


				for (m = 0; block[i][(j+1+L)%L][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block
	
				// Connect all adjacent blocks overlapping in time with the current block
				while (m < block[i][(j+1+L)%L].size() and block[i][(j+1+L)%L][m].t0 < t1) {
					block_graph[v].insert(block[i][(j+1+L)%L][m].vertex);
					block_graph[block[i][(j+1+L)%L][m].vertex].insert(v);

					p = 0.5*(1 - pow(1 - 2*physical_error, MIN(block[i][(j+1+L)%L][m].t1,t1) - MAX(block[i][(j+1+L)%L][m].t0,t0)));	// Contracted edge's new weight
					weight_graph[v][block[i][(j+1+L)%L][m].vertex] = log((1-p)/p);
					weight_graph[block[i][(j+1+L)%L][m].vertex][v] = log((1-p)/p);

					m++;
				}
			
				block_size = block[i][j].size();
				for (t = 1; t < block_size; t++) {	// The first time step was done separately above because in the next ones we shall also connect the blocks vertically in time

					t0 = block[i][j][t].t0;
					t1 = block[i][j][t].t1;
					v  = block[i][j][t].vertex;


					for (m = 0; block[(i+1+L)%L][j][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block

					// Connect all adjacent blocks overlapping in time with the current block
					while (m < block[(i+1+L)%L][j].size() and block[(i+1+L)%L][j][m].t0 < t1) {
						block_graph[v].insert(block[(i+1+L)%L][j][m].vertex);
						block_graph[block[(i+1+L)%L][j][m].vertex].insert(v);

						p = 0.5*(1 - pow(1 - 2*physical_error, MIN(block[(i+1+L)%L][j][m].t1,t1) - MAX(block[(i+1+L)%L][j][m].t0,t0)));	// Contracted edge's new weight
						weight_graph[v][block[(i+1+L)%L][j][m].vertex] = log((1-p)/p);
						weight_graph[block[(i+1+L)%L][j][m].vertex][v] = log((1-p)/p);

						m++;
					}


					for (m = 0; block[i][(j+1+L)%L][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block
		
					// Connect all adjacent blocks overlapping in time with the current block
					while (m < block[i][(j+1+L)%L].size() and block[i][(j+1+L)%L][m].t0 < t1) {
						block_graph[v].insert(block[i][(j+1+L)%L][m].vertex);
						block_graph[block[i][(j+1+L)%L][m].vertex].insert(v);

						p = 0.5*(1 - pow(1 - 2*physical_error, MIN(block[i][(j+1+L)%L][m].t1,t1) - MAX(block[i][(j+1+L)%L][m].t0,t0)));	// Contracted edge's new weight
						weight_graph[v][block[i][(j+1+L)%L][m].vertex] = log((1-p)/p);
						weight_graph[block[i][(j+1+L)%L][m].vertex][v] = log((1-p)/p);

						m++;
					}

					// Connect consective blocks in time
					block_graph[v].insert(block[i][j][t-1].vertex);
					block_graph[block[i][j][t-1].vertex].insert(v);

					weight_graph[v][block[i][j][t-1].vertex] = log((1-measu_error)/measu_error);
					weight_graph[block[i][j][t-1].vertex][v] = log((1-measu_error)/measu_error);
				}
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET MATCHING GRAPH AND PERFECT MATCHING

		vector<double> distance(n_blocks);

		struct MinHeap* minHeap = createMinHeap(n_blocks);	// Create min heap
		for (j = 0; j < n_blocks; j++)
			minHeap->array[j] = (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode)); 

		n_syn = syndrome.size();
	    	PerfectMatching *pm = new PerfectMatching(n_syn, 10*n_syn);

	    	for (i = 0; i < n_syn; i++) {

			for (j = 0; j < n_blocks; j++) {	// Dijkstra's algorithm initialisation
				distance[j] = INT_MAX;
				minHeap->array[j]->v = j;
				minHeap->array[j]->dist = INT_MAX; 
				minHeap->pos[j] = j;
			}

			distance[syndrome[i].vertex] = 0.;
			decreaseKey(minHeap, syndrome[i].vertex, 0.); 	// Set current anyon at the root of min heap 
			minHeap->size = n_blocks; 

			while (minHeap->size != 0) { 		// Dijkstra's algorithm
				current = extractMin(minHeap); 

				for (auto child: block_graph[current]) {	// Get all neighbours of current
					if (minHeap->pos[child] < minHeap->size) {	// Check if child should be visited

						if (distance[current] + weight_graph[current][child] < distance[child]) { 
							distance[child] = distance[current] + weight_graph[current][child];
							decreaseKey(minHeap, child, distance[child]); 
						}
					}
				}
			}

	   		for (j = i+1; j < n_syn; j++)
		   		pm->AddEdge(i, j, int(round(1000*distance[syndrome[j].vertex])));	// Add the distance between anyons i and j
	    	}

		for (j = 0; j < n_blocks; j++)
			free(minHeap->array[j]);
		free(minHeap->array);
		free(minHeap->pos);
		free(minHeap);

///////////////////////////////////////////////////////////
/////////////////////////////////////////// CORRECT SYNDROME

	    	pm->Solve();	// Use Edmonds' matching algorithm

	    	for (i = 0; i < n_syn; i++) {
			j = pm->GetMatch(i);	// Get the anyon matched with i

			if (i < j) {
		    		x0 = syndrome[i].row;
		    		y0 = syndrome[i].column;
	      	    		x1 = syndrome[j].row;
		    		y1 = syndrome[j].column;

				// Flip the lattice edges between anyons i and j
		    		if ((x1-x0+L) % L < (x0-x1+L) % L) {
					for (m = 0; m < (x1-x0+L) % L; m++)
			    			lattice_V[(x1-m-1+L)%L][y1][time-1] ^= true;
		    		}
		    		else {
					for (m = 0; m < (x0-x1+L) % L; m++)
			    			lattice_V[(x1+m+L)%L][y1][time-1] ^= true;
		    		}

		    		if ((y1-y0+L) % L < (y0-y1+L) % L) {
					for (m = 0; m < (y1-y0+L) % L; m++)
			    			lattice_H[x0][(y1-m-1+L)%L][time-1] ^= true;
		    		}
		    		else {
					for (m = 0; m < (y0-y1+L) % L; m++)
			    			lattice_H[x0][(y1+m+L)%L][time-1] ^= true;
		    		}
			}
	    	}
		delete pm;

///////////////////////////////////////////////////////////
/////////////////////////////////////////// VERIFY LOGICAL ERROR

		flag = true;
		error_horizontal = false, error_vertical = false;

		for (i = 0; i < L; i++) {

			// Check for the parity of the number of errors along any column or row
			for (j = 0; j < L; j++) {
		    		error_horizontal ^= lattice_H[j][i][time-1];
		    		error_vertical ^= lattice_V[i][j][time-1];
			}
			// Check if there is a logical error
			if (error_horizontal or error_vertical) {
		    		flag = false;
				break;
			}
		}
	    	if (flag)
			success++;
	}

	return double(success)/precision;
}
