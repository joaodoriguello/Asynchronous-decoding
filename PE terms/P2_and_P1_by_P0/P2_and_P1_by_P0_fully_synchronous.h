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


typedef struct block_node {

	int row, column, t;

} block_node;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool random_error (double error) {

	if ( ((double)rand() / ((double)(RAND_MAX) + (double)(1))) < error)
		return true;
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
tuple<int, double, double> P2_and_P1_by_P0_fully_synchronous (int L, double physical_error, double measu_error, int time, int precision) {

	unsigned long long int degeneracy, average = 0;
	long double P1_by_P0 = 0., P2_by_P0 = 0.;

	int i, j, t, m;
	int x0, y0, x1, y1;
	int dx, dy, dt;
	int L2 = L/2;
	bool syndrome0, syndrome1;

	int n_syn;
	double p = physical_error*physical_error/((1-physical_error)*(1-physical_error));

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

				for (t = 1; t < time - 1; t++) {

					syndrome1 = lattice_H[i][j][t] ^ lattice_H[i][(j-1+L)%L][t] ^ lattice_V[i][j][t] ^ lattice_V[(i-1+L)%L][j][t] ^ random_error(measu_error);
					if (syndrome0 != syndrome1)
						syndrome.push_back(block_node({i,j,t}));

					syndrome0 = syndrome1;
				}

				syndrome1 = lattice_H[i][j][time-1] ^ lattice_H[i][(j-1+L)%L][time-1] ^ lattice_V[i][j][time-1] ^ lattice_V[(i-1+L)%L][j][time-1];
				if (syndrome0 != syndrome1)
					syndrome.push_back(block_node({i,j,time-1}));
			}
		}

///////////////////////////////////////////////////////////
/////////////////////////////////////////// GET MATCHING GRAPH

		n_syn = syndrome.size();

		for (i = 0; i < n_syn; i++) {

			x0 = syndrome[i].row;
			y0 = syndrome[i].column;

			for (j = i+1; j < n_syn; j++) {

		    		x1 = syndrome[j].row;
		    		y1 = syndrome[j].column;

				dx = MIN((x1-x0+L)%L, (x0-x1+L)%L);
				dy = MIN((y1-y0+L)%L, (y0-y1+L)%L);
				dt = abs(syndrome[i].t - syndrome[j].t);

				degeneracy = (1 + (dx == L2)) * (1 + (dy == L2)) * tgamma(dx + dy + dt + 1)/(tgamma(dx + 1) * tgamma(dy + 1) * tgamma(dt + 1));

				if (degeneracy > 1)
					P1_by_P0 += 1.;
				else
					P1_by_P0 += p;

				if (degeneracy > 2)
					P2_by_P0 += 1.;
				else
					P2_by_P0 += p;

				average++;		
			}
		}
	}
		
	return make_tuple(average, P1_by_P0, P2_by_P0);

}
