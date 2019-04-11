#include "utilities.h"
#include "least_squares.h"

#include <iostream>

using namespace std;

void SolveLeastSquares(double** m, int rows, int columns, double* v, double* x) {
	double** mtm = new double*[columns];

	for (int i = 0; i < columns; i++)
		mtm[i] = new double[columns];

	for (int i = 0; i < columns; i++) {
		for (int j = 0; j < columns; j++) {
			mtm[i][j] = 0;

			for (int k = 0; k < rows; k++)
				mtm[i][j] += m[k][j] * m[k][i];
		}
	}

	

	for (int i = 0; i < columns; i++)
		delete[] mtm[i];

	delete[] mtm;
}