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
double CG_decoder_first_second_order_continuous (int L, double physical_error, double measu_error, int time, int precision, double tau) {

	int success = 0;

	int i, j, t, m;
	int x0, y0, x1, y1, v;
	double t0, t1;
	bool syndrome0, syndrome1, flag, aux_bool;
	double beta = log((1-physical_error)/physical_error);
	double exp_beta = physical_error/(1-physical_error);
	double p;

	int n_syn, n_blocks, block_size;
	int aux_count, current;
	double aux_time;
	bool error_horizontal, error_vertical;

	// Lattice state after all measurement rounds
	bool lattice_H[L][L];
	bool lattice_V[L][L];

	// List of measurement times
	priority_queue<double, vector<double>, greater<double>> measu_list;
	
	// List of error times for horizontal and vertical lattices. "lattice_H_time2" will be a copy of "lattice_H_time1". The same is true for "lattice_V_time2" and "lattice_V_time1"
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_H_time1 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_H_time2 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_V_time1 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));
	vector<vector< priority_queue<double, vector<double>, greater<double>> >> lattice_V_time2 (L, vector< priority_queue<double, vector<double>, greater<double>> >(L));

	queue<int> queue;

	// Poisson distributions
	random_device rd;
	default_random_engine generator (rd());
	uniform_real_distribution<double> unif_distribution (0.0001,time-0.0001);
	poisson_distribution<int> p_error (-time*0.5*log(1 - 2*physical_error));
	poisson_distribution<int> p_measu (time);


	for (int n = 0; n < precision; n++) {

///////////////////////////////////////// CREATE LATTICE

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {

				aux_count = p_error(generator);	// Sample number of errors from Poisson distribution
				lattice_H[i][j] = bool(aux_count % 2);	// Lattice state at the end

				for (t = 0; t < aux_count; t++) {	// Obtain the list of times of the physical errors
					aux_time = unif_distribution(generator);
					lattice_H_time1[i][j].push(aux_time);
					lattice_H_time2[i][j].push(aux_time);
				}


				aux_count = p_error(generator);	// Sample number of errors from Poisson distribution
				lattice_V[i][j] = bool(aux_count % 2);	// Lattice state at the end

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

					block[i][j].push_back(block_node({t0,t1,n_blocks}));	// All blocks are numbered by their order of appearance equal to n_blocks
		
					// Check the parity of the number of physical errors prior to the measurement time
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
		vector<vector<double>> degeneracy_graph(n_blocks, vector<double>(n_blocks));

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

					p = MIN(block[(i+1+L)%L][j][m].t1,t1) - MAX(block[(i+1+L)%L][j][m].t0,t0);	// Adjacent-blocks time overlap
					degeneracy_graph[v][block[(i+1+L)%L][j][m].vertex] = p;
					degeneracy_graph[block[(i+1+L)%L][j][m].vertex][v] = p;

					m++;
				}


				for (m = 0; block[i][(j+1+L)%L][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block

				// Connect all adjacent blocks overlapping in time with the current block
				while (m < block[i][(j+1+L)%L].size() and block[i][(j+1+L)%L][m].t0 < t1) {
					block_graph[v].insert(block[i][(j+1+L)%L][m].vertex);
					block_graph[block[i][(j+1+L)%L][m].vertex].insert(v);

					p = MIN(block[i][(j+1+L)%L][m].t1,t1) - MAX(block[i][(j+1+L)%L][m].t0,t0);	// Adjacent-blocks time overlap
					degeneracy_graph[v][block[i][(j+1+L)%L][m].vertex] = p;
					degeneracy_graph[block[i][(j+1+L)%L][m].vertex][v] = p;

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

						p = MIN(block[(i+1+L)%L][j][m].t1,t1) - MAX(block[(i+1+L)%L][j][m].t0,t0);	// Adjacent-blocks time overlap
						degeneracy_graph[v][block[(i+1+L)%L][j][m].vertex] = p;
						degeneracy_graph[block[(i+1+L)%L][j][m].vertex][v] = p;

						m++;
					}


					for (m = 0; block[i][(j+1+L)%L][m].t1 <= t0; m++);	// Find the first adjacent block overlapping in time with the current block
	
					// Connect all adjacent blocks overlapping in time with the current block
					while (m < block[i][(j+1+L)%L].size() and block[i][(j+1+L)%L][m].t0 < t1) {
						block_graph[v].insert(block[i][(j+1+L)%L][m].vertex);
						block_graph[block[i][(j+1+L)%L][m].vertex].insert(v);

						p = MIN(block[i][(j+1+L)%L][m].t1,t1) - MAX(block[i][(j+1+L)%L][m].t0,t0);	// Adjacent-blocks time overlap
						degeneracy_graph[v][block[i][(j+1+L)%L][m].vertex] = p;
						degeneracy_graph[block[i][(j+1+L)%L][m].vertex][v] = p;

						m++;
					}

					// Consecutive blocks have time overlap equal to 1
					block_graph[v].insert(block[i][j][t-1].vertex);
					block_graph[block[i][j][t-1].vertex].insert(v);

					degeneracy_graph[v][block[i][j][t-1].vertex] = 1;
					degeneracy_graph[block[i][j][t-1].vertex][v] = 1;
				}
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET MATCHING GRAPH AND PERFECT MATCHING

		vector<int> distance(n_blocks);
		vector<long double> degeneracy_first(n_blocks);	// First-order degeneracy parameter
		vector<long double> degeneracy_second(n_blocks);	// Second-order degeneracy parameter

		n_syn = syndrome.size();
	    	PerfectMatching *pm = new PerfectMatching(n_syn, 10*n_syn);

	    	for (i = 0; i < n_syn; i++) {

			for (j = 0; j < n_blocks; j++) {	// Dijkstra's algorithm initialisation
				distance[j] = -1;
				degeneracy_second[j] = 0;
			}

			distance[syndrome[i].vertex] = 0;
			degeneracy_first[syndrome[i].vertex] = 1;
			queue.push(syndrome[i].vertex);	// Set current anyon as the first vertex to be visited

			while (!queue.empty()) { 		// Simplified Dijkstra's algorithm
				current = queue.front();
				queue.pop(); 

				for (auto child: block_graph[current]) {	// Get all neighbours of current
				
					if (distance[child] == -1) {	// Check if child has been visited
						distance[child] = distance[current] + 1;
						degeneracy_first[child] = degeneracy_first[current] * degeneracy_graph[child][current];
						queue.push(child);	// Add child to the queue of vertices to be enquired
					}
					else if (distance[child] == distance[current] + 1)
						degeneracy_first[child] += degeneracy_first[current] * degeneracy_graph[child][current];

					else if (distance[child] == distance[current])
						degeneracy_second[child] += degeneracy_first[current] * degeneracy_graph[child][current];

					else if (distance[child] == distance[current] - 1)
						degeneracy_second[current] += degeneracy_second[child] * degeneracy_graph[child][current];
				}
			}

	   		for (j = i+1; j < n_syn; j++)
		  		pm->AddEdge(i, j, int(round(1000*(beta*distance[syndrome[j].vertex] - tau*log(degeneracy_first[syndrome[j].vertex] + exp_beta*degeneracy_second[syndrome[j].vertex])))));	// Add the distance between anyons i and j
	    	}

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
			    			lattice_V[(x1-m-1+L)%L][y1] ^= true;
		    		}
		    		else {
					for (m = 0; m < (x0-x1+L) % L; m++)
			    			lattice_V[(x1+m+L)%L][y1] ^= true;
		    		}

		    		if ((y1-y0+L) % L < (y0-y1+L) % L) {
					for (m = 0; m < (y1-y0+L) % L; m++)
			    			lattice_H[x0][(y1-m-1+L)%L] ^= true;
		    		}
		    		else {
					for (m = 0; m < (y0-y1+L) % L; m++)
			    			lattice_H[x0][(y1+m+L)%L] ^= true;
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
		    		error_horizontal ^= lattice_H[j][i];
		    		error_vertical ^= lattice_V[i][j];
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
