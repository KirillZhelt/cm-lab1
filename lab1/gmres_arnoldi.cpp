#include "gmres_arnoldi.h"
#include "least_squares.h"
#include "utilities.h"

void SolveGMRESArnoldi(double** m, int rows, int columns, double* v, double* x) {
	double** Q = new double*[rows];

	for (int i = 0; i < rows; i++)
		Q[i] = new double[columns];

	double* y = new double[columns];

	double* difference = new double[columns];

	int k = 0;

	while (true) {
		// TODO: generate qk

		k++;

		SolveLeastSquares(Q, rows, k, v, y);

		Multiply(Q, rows, k, y, rows, x);

		Multiply(m, rows, columns, x, columns, difference);
		Subtract(difference, v, columns, difference);

		if (EuclideanNorm(difference, columns) < EPS_GMRES_ARNOLDI * 1000)
			break;
	}

	delete[] difference;

	delete[] y;

	for (int i = 0; i < rows; i++)
		delete[] Q[i];

	delete[] Q;
}