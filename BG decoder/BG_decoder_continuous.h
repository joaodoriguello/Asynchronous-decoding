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

// Structure to save the coordinates of a block
typedef struct block_node {

	int row, column;
	double t0, t1;

} block_node;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function for a biased coin
bool random_error (double error) {

	if ( ((double)rand() / ((double)(RAND_MAX) + (double)(1))) < error)
		return true;
	else
		return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function that outputs the probability of successful decoding for a given set of inputs
double ApSG_decoder_continuous (int L, double physical_error, double measu_error, int time, int precision) {

	int success = 0;

	int i, j, t, m;
	int x0, y0, x1, y1, v;
	double t0, t1;
	int dx, dy;
	double dt;
	bool syndrome0, syndrome1, flag, aux_bool;

	int n_syn, aux_count;
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

		vector<block_node> syndrome;

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
						syndrome.push_back(block_node({i,j,t0,t1}));
				
					t0 = t1;
					syndrome0 = syndrome1;
				}

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
					syndrome.push_back(block_node({i,j,t0,double(time)}));
			}
	   	}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET MATCHING GRAPH AND PERFECT MATCHING

		n_syn = syndrome.size();

	    	PerfectMatching *pm = new PerfectMatching(n_syn, 10*n_syn);
	    	for (i = 0; i < n_syn; i++) {

			x0 = syndrome[i].row;
			y0 = syndrome[i].column;

	   		for (j = i+1; j < n_syn; j++) {

		    		x1 = syndrome[j].row;
		    		y1 = syndrome[j].column;

				dx = MIN((x1-x0+L)%L,(x0-x1+L)%L);	// X-distance between anyons
				dy = MIN((y1-y0+L)%L,(y0-y1+L)%L);	// Y-distance between anyons

				// Time distance between anyons
				if (syndrome[i].t0 >= syndrome[j].t1) 
					dt = syndrome[i].t0 - syndrome[j].t1;
				else if (syndrome[j].t0 >= syndrome[i].t1) 
					dt = syndrome[j].t0 - syndrome[i].t1;
				else
					dt = 0.;

				pm->AddEdge(i, j, int(round(1000*(dx + dy + dt))));	// Add the distance between anyons i and j
			}
	    	}
	    	pm->Solve();	// Use Edmonds' matching algorithm

///////////////////////////////////////////////////////////
/////////////////////////////////////////// CORRECT SYNDROME

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
