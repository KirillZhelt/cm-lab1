#include <cmath>
#include <iostream>

#include "utilities.h"
#include "cholesky.h"

using namespace std;

void BuildCholesky(double** m, int rows, int columns, double** lt, double* d) {
	CopyMatrix(m, lt, rows, columns);
	
	for (int i = 0; i < rows; i++)
		d[i] = 1;
	
	for (int k = 0; k < rows; k++) {
		if (lt[k][k] < 0) {
			d[k] = -1;
			lt[k][k] = -lt[k][k];
		}
		  
		double sqrt_divider = sqrt(lt[k][k]);

		for (int j = k; j < columns; j++)
			lt[k][j] /= sqrt_divider;

		for (int i = k + 1; i < rows; i++) {
			auto l = lt[i][k] / lt[k][k];

			for (int j = i; j < columns; j++)
				lt[j][i] = lt[i][j] = lt[i][j] - l * lt[k][j];
		}
	}
}

void SolveCholesky(double** lt, double* d, int rows, int columns, double* v, double* x) {
	double* b = new double[rows];

	CopyVector(v, b, rows);

	double* y = new double[columns];
	
	for (int i = 0; i < columns; i++) {
		double sum = 0;

		for (int j = 0; j < i; j++)
			sum += lt[j][i] * y[j];

		y[i] = (b[i] - sum) / lt[i][i];
	}
	
	double* z = new double[columns];
	for (int i = 0; i < columns; i++)
		z[i] = d[i] * y[i];

	for (int i = columns - 1; i >= 0; i--) {
		double sum = 0;

		for (int j = i + 1; j < columns; j++)
			sum += lt[i][j] * x[j];

		x[i] = (z[i] - sum) / lt[i][i];
	}

	delete[] z;
	delete[] y;
	delete[] b;
}