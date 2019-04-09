#include <cmath>
#include <iostream>

#include "utilities.h"
#include "gauss.h"

using namespace std;

void Gauss(double** m, int rows, int columns, double* v, double* x) {
	double** a = new double*[rows];

	for (int i = 0; i < rows; i++) 
		a[i] = new double[columns];
	
	CopyMatrix(m, a, rows, columns);

	double* b = new double[rows];
	CopyVector(v, b, rows);

	int* x_positions = new int[columns];
	for (int i = 0; i < columns; i++)
		x_positions[i] = i;

	for (int k = 0; k < rows - 1; k++) {
		// выбор главного элемента по матрице
		double max = a[k][k];
		int max_row = k;
		int max_column = k;
		for (int i = k; i < rows; i++) {
			for (int j = k; j < columns; j++) {
				if (fabs(a[i][j]) > fabs(max)) {
					max_row = i;
					max_column = j;

					max = a[i][j];
				}
			}
		}

		swap(b[max_row], b[k]);
		SwapRows(a, max_row, k);

		swap(x_positions[max_column], x_positions[k]);
		SwapColumns(a, rows, max_column, k);

		for (int i = k + 1; i < rows; i++) {
			auto l = a[i][k] / a[k][k];

			for (int j = k; j < columns; j++) 
				a[i][j] -= l * a[k][j];

			b[i] -= l * b[k];
		}
	}

	for (int i = rows - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < columns; j++)
			sum += a[i][j] * x[j];

		x[i] = (b[i] - sum) / a[i][i];
	}
	
	Rearrange(x, x_positions, columns);
	
	delete[] x_positions;

	delete[] b;

	for (int i = 0; i < rows; i++)
		delete[] a[i];

	delete[] a;
}
