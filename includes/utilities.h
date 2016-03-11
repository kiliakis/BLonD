/*
 * utilities.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "configuration.h"

#ifndef INCLUDES_UTILITIES_H_
#define INCLUDES_UTILITIES_H_

#define dprintf(...)    fprintf(stdout, __VA_ARGS__)     // Debug printf

ftype *linspace(const ftype start, const ftype end, const int n) {
	ftype * a = new ftype[n];
	ftype step = (end - start) / (n - 1);
	ftype value = start;
	for (int i = 0; i < n; ++i) {
		a[i] = value;
		value += step;
	}
	return a;
}

void dump(ftype* a, int n, const char* s) {
	dprintf("\n");
	for (int i = 0; i < n; ++i) {
		dprintf("%s[%d] = %.8lf\n", s, i, a[i]);
	}
	dprintf("\n");

}

void dump(int* a, int n, const char* s) {
	dprintf("\n");
	for (int i = 0; i < n; ++i) {
		dprintf("%s[%d] = %d\n", s, i, a[i]);
	}
	dprintf("\n");

}

ftype mean(const ftype data[], const int n) {
	ftype m = 0;
	for (int i = 0; i < n; ++i) {
		m += data[i];
	}
	return m / n;
}

ftype standard_deviation(const ftype data[], const int n, const ftype mean) {
	ftype sum_deviation = 0.0;
	int i;
	for (i = 0; i < n; ++i)
		sum_deviation += (data[i] - mean) * (data[i] - mean);
	return sqrt(sum_deviation / n);
}

#endif /* INCLUDES_UTILITIES_H_ */
