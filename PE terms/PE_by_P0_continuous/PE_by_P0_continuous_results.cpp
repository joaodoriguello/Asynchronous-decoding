#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "PE_by_P0_continuous.h"

int main (void) {

	srand(time(NULL));

    	int precision = 5;
    	int L = 14;
	int n_shortest_paths = 11;
	int time_factor = 2;
	double physical_error = 0.01695;

	vector<double> PE_by_P0;
	int average;

	std::tie(average, PE_by_P0) = PE_by_P0_continuous(L, physical_error, physical_error, L*time_factor, precision, n_shortest_paths);

	for (int i = 1; i < n_shortest_paths; i++) 
		printf("P%d %lf\n", i, PE_by_P0[i]/average);

    	return 0;
}

