#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <queue>
#include <stack>
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

	double t0, t1;
	int vertex;

} block_node;

// Structure to save the coordinates of an anyon block
typedef struct syndrome_node {

	int row, column;
	double t0, t1;
	int vertex;

} syndrome_node;

// Structure for the shortest paths in Yen's algorithm. It stores the path, its length (number of edges), and its distance (sum of edge weights)
typedef struct MinHeapNodePaths { 

	int *path;
	int length;
	double dist;

} MinHeapNodePaths;

// Function to check which path has less distance
struct Comp {

	bool operator()(MinHeapNodePaths* s1, MinHeapNodePaths* s2) {
		return s1->dist < s2->dist;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to check if two paths arr1[] and arr2[] with lengths length1 and length2, respectively, are equal or not
bool areEqual (int arr1[], int arr2[], int length1, int length2)  { 

    	if (length1 != length2) 
        	return false; 
  
    	for (int i = 0; i < length1; i++) 
        	if (arr1[i] != arr2[i]) 
            		return false; 

    	return true; 
}

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

// Function that outputs the number of anyon pairs considered (average), and the ratios P1/P0 and P2/P0 between the second-shortest and shortest paths, and third-shortest and shortest paths
tuple<long int, double, double> P2_and_P1_by_P0_continuous (int L, double physical_error, double measu_error, int time, int precision) {

	long int average = 0;
	double P1_by_P0 = 0., P2_by_P0 = 0.;

	int n_shortest_paths = 3;
	int i, j, k, t, l, m, n;
	int x0, y0, x1, y1, v, v1, v2;
	double t0, t1;
	double p, cost;
	bool syndrome0, syndrome1, aux_bool, flag;

	int n_syn, n_blocks, block_size;
	int aux_count, count, current;
	double aux_time;

	stack<int> erasedEdges, mystack;
	unordered_set<int> erasedNodes;
	struct MinHeapNodePaths *testPath;
	int spurNode;

	// List of measurement times
	priority_queue<double, vector<double>, greater<double>> measu_list;
	
	// List of error times for horizontal and vertical lattices. "lattice_H_time2" will be a copy of "lattice_H_time1". The same is true for "lattice_V_time2" and "lattice_V_time1"
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_H_time1 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_H_time2 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_V_time1 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_V_time2 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));

	// Poisson distributions
	random_device rd;
	default_random_engine generator (rd());
	uniform_real_distribution<double> unif_distribution (0.0001,time-0.0001);
	poisson_distribution<int> p_error (-time*0.5*log(1 - 2*physical_error));
	poisson_distribution<int> p_measu (time);


	for (int loop = 0; loop < precision; loop++) {

///////////////////////////////////////// CREATE LATTICE

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {

				aux_count = p_error(generator);	// Sample number of errors from Poisson distribution

				for (t = 0; t < aux_count; t++) {	// Obtain the list of times of the physical errors
					aux_time = unif_distribution(generator);
					lattice_H_time1[i][j].push(aux_time);
					lattice_H_time2[i][j].push(aux_time);
				}


				aux_count = p_error(generator);	// Sample number of errors from Poisson distribution

				for (t = 0; t < aux_count; t++) {	// Obtain the list of times of the physical errors
					aux_time = unif_distribution(generator);
					lattice_V_time1[i][j].push(aux_time);
					lattice_V_time2[i][j].push(aux_time);
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

				aux_count = p_measu(generator);	// Sample number of errors from Poisson distribution
				for (t = 0; t < aux_count; t++)
					measu_list.push(unif_distribution(generator));

				syndrome0 = false;
				aux_bool = false;
				t0 = 0.0;

				while (!measu_list.empty()) {

					t1 = measu_list.top();
					measu_list.pop();

					block[i][j].push_back(block_node({t0,t1,n_blocks})); // All blocks are numbered by their order of appearance equal to n_blocks
		
					// Checks the parity of the number of physical errors prior to the measurement time
					while (!lattice_H_time1[i][j].empty() and lattice_H_time1[i][j].top() <= t1) {
						aux_bool = !aux_bool;
						lattice_H_time1[i][j].pop();
					}

					while (!lattice_V_time1[i][j].empty() and lattice_V_time1[i][j].top() <= t1) {
						aux_bool = !aux_bool;
						lattice_V_time1[i][j].pop();
					}

					while (!lattice_H_time2[i][(j-1+L)%L].empty() and lattice_H_time2[i][(j-1+L)%L].top() <= t1) {
						aux_bool = !aux_bool;
						lattice_H_time2[i][(j-1+L)%L].pop();
					}

					while (!lattice_V_time2[(i-1+L)%L][j].empty() and lattice_V_time2[(i-1+L)%L][j].top() <= t1) {
						aux_bool = !aux_bool;
						lattice_V_time2[(i-1+L)%L][j].pop();
					}

					syndrome1 = aux_bool ^ random_error(physical_error);	// Add measurement error
					if (syndrome0 != syndrome1)	// Different consecutive measurements lead to an anyon block
						syndrome.push_back(syndrome_node({i,j,t0,t1,n_blocks}));
				
					t0 = t1;
					syndrome0 = syndrome1;
					n_blocks++;
				}

				block[i][j].push_back(block_node({t0,double(time),n_blocks}));	// All blocks are numbered by their order of appearance equal to n_blocks

				// Last measurement round with perfect stabilizer measurement
				while (!lattice_H_time1[i][j].empty() and lattice_H_time1[i][j].top() <= time) {
					aux_bool = !aux_bool;
					lattice_H_time1[i][j].pop();
				}

				while (!lattice_V_time1[i][j].empty() and lattice_V_time1[i][j].top() <= time) {
					aux_bool = !aux_bool;
					lattice_V_time1[i][j].pop();
				}

				while (!lattice_H_time2[i][(j-1+L)%L].empty() and lattice_H_time2[i][(j-1+L)%L].top() <= time) {
					aux_bool = !aux_bool;
					lattice_H_time2[i][(j-1+L)%L].pop();
				}

				while (!lattice_V_time2[(i-1+L)%L][j].empty() and lattice_V_time2[(i-1+L)%L][j].top() <= time) {
					aux_bool = !aux_bool;
					lattice_V_time2[(i-1+L)%L][j].pop();
				}

				if (syndrome0 != aux_bool)	// Different consecutive measurements lead to an anyon block
					syndrome.push_back(syndrome_node({i,j,t0,double(time),n_blocks}));

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

					weight_graph[v][block[i][j][t-1].vertex] = log((1-physical_error)/physical_error);
					weight_graph[block[i][j][t-1].vertex][v] = log((1-physical_error)/physical_error);
				}
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET PERFECT MATCHING AND COMPUTE RATIOS

		n_syn = syndrome.size();
		
		vector<double> distance(n_blocks);
		vector<vector<double>> distance_matrix(n_syn, vector<double>(n_syn));
		
		vector<vector<int>> prev(n_syn, vector<int>(n_blocks));

		struct MinHeap* minHeap = createMinHeap(n_blocks);	// Create min heap
		for (j = 0; j < n_blocks; j++)
			minHeap->array[j] = (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode)); 


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
							prev[i][child] = current;	// Set current as the previous vertex of child in the shortest path
							decreaseKey(minHeap, child, distance[child]); 
						}
					}
				}
			}

	   		for (j = i+1; j < n_syn; j++) {
		  		pm->AddEdge(i, j, int(round(1000*distance[syndrome[j].vertex])));
		  		distance_matrix[i][j] = distance[syndrome[j].vertex];
		  		distance_matrix[j][i] = distance[syndrome[j].vertex];
		  	}		
	    	}

		vector<int> prev_aux(n_blocks);
		vector<MinHeapNodePaths*> path_list;
		vector<MinHeapNodePaths*> minHeapPaths;

		pm->Solve();	// Use Edmonds' matching algorithm
		for (i = 0; i < n_syn; i++) {
		
			j = pm->GetMatch(i);	// Get the anyon matched with i
			if (i < j) {
				
				count = 1;
				m = syndrome[j].vertex;
				while (m != syndrome[i].vertex) {	// Obtain the path between anyons i and j in 'mystack' (in reverse order) and its length in 'count'
					mystack.push(m);
					m = prev[i][m];
					count++;
				}

				// Input the path from 'mystack' into the tentative shortest path 'testPath'
				testPath = new struct MinHeapNodePaths;
				testPath->path = new int[count];
				testPath->path[0] = syndrome[i].vertex;
				for (l = 1; l < count; l++) {
					testPath->path[l] = mystack.top();
					mystack.pop();
				}

				testPath->length = count;	// Set 'testPath' length and distance
				testPath->dist = distance_matrix[i][j];

				path_list.push_back(testPath);	// Input 'testPath' into a list of tentative shortest paths

				for (k = 1; k < n_shortest_paths; k++) {

					for (m = 0; m < path_list[k-1]->length - 1; m++) {	// Remove the links that are part of the previous shortest paths which share the same root path

						cost = 0.;
						spurNode = path_list[k-1]->path[m];	// Spur node is retrieved from the previous k-shortest path, k-1

						testPath = new struct MinHeapNodePaths;
						testPath->path = (int*) malloc((m+1) * sizeof(int));

						for (l = 0; l < m; l++) {
							testPath->path[l] = path_list[k-1]->path[l];	// The sequence of nodes from the source to the spur node of the previous k-shortest path
							cost += weight_graph[path_list[k-1]->path[l]][path_list[k-1]->path[l+1]];
							erasedNodes.insert(path_list[k-1]->path[l]);	// The nodes from 'testPath' are erased from the contracted graphs
						}
						testPath->path[m] = spurNode;

						for (l = 0; l < k; l++) {	// For each path in the tentative path list
							if (path_list[l]->length >= m+1) {	

								flag = true;
								for (n = 0; n <= m; n++) {
									if (testPath->path[n] != path_list[l]->path[n]) {	// Check if both paths are equal
										flag = false;
										break;
									}
								}

								if (flag) {	// If paths are equal, delete edge(m, m+1) from contracted syndrome graph. Erased nodes are stored in 'erasedEdges'
									erasedEdges.push(path_list[l]->path[m]);
									erasedEdges.push(path_list[l]->path[m+1]);
									block_graph[path_list[l]->path[m]].erase(path_list[l]->path[m+1]);
									block_graph[path_list[l]->path[m+1]].erase(path_list[l]->path[m]);
								}
							}
						}

						// Calculate the spur path from the spur node to the sink (anyon j)
						for (l = 0; l < n_blocks; l++) {	// Dijkstra's algorithm initialisation
							distance[l] = INT_MAX;
							prev_aux[l] = -1;
							minHeap->array[l]->v = l;
							minHeap->array[l]->dist = INT_MAX;
							minHeap->pos[l] = l;
						}

						distance[spurNode] = 0.;
						prev_aux[spurNode] = spurNode;
						decreaseKey(minHeap, spurNode, 0.);
						minHeap->size = n_blocks; 

						while (minHeap->size != 0) { 		// Dijkstra's algorithm for Yen's algorithm
							current = extractMin(minHeap);

							for (auto child: block_graph[current]) {
								if (minHeap->pos[child] < minHeap->size and erasedNodes.find(child) == erasedNodes.end()) {	// Check if child was not erased

									if (distance[current] + weight_graph[current][child] < distance[child]) {
										distance[child] = distance[current] + weight_graph[current][child];
										prev_aux[child] = current;
										decreaseKey(minHeap, child, distance[child]);
									}
								}
							}
						}


						while (!erasedEdges.empty()) {	// Return erased edges to contracted syndrome graph
							v1 = erasedEdges.top();
							erasedEdges.pop();
							v2 = erasedEdges.top();
							erasedEdges.pop();
							block_graph[v1].insert(v2);
							block_graph[v2].insert(v1);
						}
						unordered_set<int>().swap(erasedNodes);

						count = 0;
						l = syndrome[j].vertex;
						while (prev_aux[l] != -1 and l != spurNode) {		// Obtain path from sink (anyon j) to spurNode
							mystack.push(l);
							l = prev_aux[l];
							count++;
						}

						if (prev_aux[l] != -1) {

							// Entire path is made up of the root path and spur path
							testPath->path = (int*) realloc(testPath->path, (m+1+count) * sizeof(int));

							for (l = 0; l < count; l++) {		// Add spur path to root path
								testPath->path[m+1+l] = mystack.top();
								mystack.pop();
							}

							testPath->length = m+1+count; 
							testPath->dist = cost + distance[syndrome[j].vertex];
							minHeapPaths.push_back(testPath);	// Add the potential k-shortest path to the heap
						}
						else {
							delete [] testPath->path;
							delete testPath;
							while (!mystack.empty())
								mystack.pop();
						}
					}

					if (minHeapPaths.size() == 0) 	// This handles the case of there being no spur paths, or no spur paths left
						break;

					while (1) {
					
						make_heap(minHeapPaths.begin(), minHeapPaths.end(), Comp());	
						sort_heap(minHeapPaths.begin(), minHeapPaths.end(), Comp());	// Order tentative paths inside min heap according to their distance

						flag = true;
						pop_heap(minHeapPaths.begin(), minHeapPaths.end());
						testPath = minHeapPaths.back();	// Pick path with smaller distance from min heap
						minHeapPaths.pop_back();

						// Check if path obtained from min heap is already in the shortest-paths list
						for (l = 0; l < k; l++) {

							if (areEqual(path_list[l]->path, testPath->path, path_list[l]->length, testPath->length)) {	
								flag = false;
								break;
							}
						}

						if (flag) {	// If the path is not in the list, add it
							path_list.push_back(testPath);
							break;
						}
					}
				}
				
				P1_by_P0 += exp(path_list[0]->dist - path_list[1]->dist);
				P2_by_P0 += exp(path_list[0]->dist - path_list[2]->dist);
				average++;
				

				for (l = 0; l < path_list.size(); l++) {
					delete [] path_list[l]->path;
					delete path_list[l];
				} 
				vector<MinHeapNodePaths*>().swap(path_list);
				
				for (l = 0; l < minHeapPaths.size(); l++) {
					delete [] minHeapPaths[l]->path;
					delete minHeapPaths[l];
				} 
				vector<MinHeapNodePaths*>().swap(minHeapPaths);
			}		
		}

		for (j = 0; j < n_blocks; j++)
			delete minHeap->array[j];
		delete [] minHeap->array;
		delete [] minHeap->pos;
		delete minHeap;
	}

	return make_tuple(average, P1_by_P0, P2_by_P0);
}
