#include "least_squares.h"

void SolveLeastSquares(double** m, int rows, int columns, double* v, double* x) {
	double** mtm = new double*[columns];

	for (int i = 0; i < columns; i++)
		mtm[i] = new double[rows];

	

	for (int i = 0; i < columns; i++)
		delete[] mtm[i];

	delete[] mtm;
}