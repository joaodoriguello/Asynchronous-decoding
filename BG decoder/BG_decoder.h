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

	int row, column, t0, t1;

} block_node;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function for a biased coin
bool random_error (double error) {

	if ( ((double)rand() / ((double)(RAND_MAX) + (double)(1))) < error)
		return true;
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function that outputs the probability of successful decoding for a given set of inputs
double ApSG_decoder (int L, double physical_error, double measu_error, double synchronicity, int time, int precision) {

	int success = 0;

	int i, j, t, m;
	int x0, y0, t0, x1, y1, t1;
	int dx, dy, dt;
	bool syndrome0, syndrome1, flag;

	int n_syn;
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

		vector<block_node> syndrome;

		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {

				syndrome0 = false;
				t0 = 0;

				for (t = 1; t < time - 1; t++) {
					if (random_error(synchronicity)) {	// Check if stabilizer measurement is successful

						syndrome1 = lattice_H[i][j][t] ^ lattice_H[i][(j-1+L)%L][t] ^ lattice_V[i][j][t] ^ lattice_V[(i-1+L)%L][j][t] ^ random_error(measu_error);
						if (syndrome0 != syndrome1)	// Different consecutive measurements lead to an anyon block
							syndrome.push_back(block_node({i,j,t0,t}));

						t0 = t;
						syndrome0 = syndrome1;
					}
				}

				// Last measurement round with perfect stabilizer measurement
				syndrome1 = lattice_H[i][j][time-1] ^ lattice_H[i][(j-1+L)%L][time-1] ^ lattice_V[i][j][time-1] ^ lattice_V[(i-1+L)%L][j][time-1];
				if (syndrome0 != syndrome1)
					syndrome.push_back(block_node({i,j,t0,time-1}));
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET MATCHING GRAPH

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
					dt = syndrome[i].t0 - syndrome[j].t1 + 1;
				else if (syndrome[j].t0 >= syndrome[i].t1) 
					dt = syndrome[j].t0 - syndrome[i].t1 + 1;
				else
					dt = 0;
		
		   		pm->AddEdge(i, j, int(10*(dx + dy + synchronicity*dt)));	// Add the distance between anyons i and j
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
		    		if ((x1-x0+L)%L < (x0-x1+L)%L) {
					for (m = 0; m < (x1-x0+L)%L; m++)
			    			lattice_V[(x1-m-1+L)%L][y1][time-1] ^= true;
		    		}
		    		else {
					for (m = 0; m < (x0-x1+L)%L; m++)
			    			lattice_V[(x1+m+L)%L][y1][time-1] ^= true;
		    		}

		    		if ((y1-y0+L)%L < (y0-y1+L)%L) {
					for (m = 0; m < (y1-y0+L)%L; m++)
			    			lattice_H[x0][(y1-m-1+L)%L][time-1] ^= true;
		    		}
		    		else {
					for (m = 0; m < (y0-y1+L)%L; m++)
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
