#include <iostream>

#include "gmres.h"
#include "least_squares.h"
#include "utilities.h"

using namespace std;

void SolveGMRES(double** m, int rows, int columns, double* v, double* x) {
	double** K = new double*[rows];

	for (int i = 0; i < rows; i++)
		K[i] = new double[columns];
	
	for (int i = 0; i < columns; i++)
		K[i][0] = v[i];

	double* y = new double[columns];
	
	double* difference = new double[columns];

	int k = 1;
	
	while (k < 256) {
		for (int i = 0; i < rows; i++) {
			K[i][k] = 0;

			for (int j = 0; j < columns; j++)
				K[i][k] += m[i][j] * K[j][k - 1]; // A * Ki
		}

		k++;

		SolveLeastSquares(K, rows, k - 1, v, y, 1);

		Multiply(K, rows, k - 1, y, rows, x);

		Multiply(m, rows, columns, x, columns, difference);
		Subtract(difference, v, columns, difference);

		if (EuclideanNorm(difference, columns) < EPS_GMRES)
			break;
	}

	delete[] difference;
	   
	delete[] y;

	for (int i = 0; i < rows; i++)
		delete[] K[i];

	delete[] K;
}