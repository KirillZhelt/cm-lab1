#include "utilities.h"
#include "least_squares.h"
#include "gauss.h"

#include <iostream>

using namespace std;

void SolveLeastSquares(double** m, int rows, int columns, double* v, double* x, int start_column) {
	double** mtm = new double*[columns];

	for (int i = 0; i < columns; i++)
		mtm[i] = new double[columns];

	for (int i = 0; i < columns; i++) {
		for (int j = start_column; j < columns + start_column; j++) {
			mtm[i][j - start_column] = 0;

			for (int k = 0; k < rows; k++)
				mtm[i][j - start_column] += m[k][j] * m[k][i];
		}
	}

	double* b = new double[columns];

	for (int i = 0; i < columns; i++) {
		b[i] = 0;

		for (int k = 0; k < rows; k++)
			b[i] += v[k] * m[k][i];
	}

	Gauss(mtm, columns, columns, b, x);

	delete[] b;

	for (int i = 0; i < columns; i++)
		delete[] mtm[i];

	delete[] mtm;
}