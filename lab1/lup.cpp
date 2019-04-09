#include <cmath>
#include <algorithm>

#include "utilities.h"
#include "lup.h"

using namespace std;

void BuildLUP(double** m, int rows, int columns, double** lu, int* p) {
	CopyMatrix(m, lu, rows, columns);
	
	for (int i = 0; i < rows; i++)
		p[i] = i;
	
	for (int k = 0; k < rows - 1; k++) {
		// выбор главного элемента по матрице
		double max = lu[k][k];
		int max_row = k;
		for (int j = k + 1; j < rows; j++) {
			if (fabs(lu[j][k]) > fabs(max)) {
				max_row = j;
				max = lu[j][k];
			}
		}

		swap(p[k], max_row);
		SwapRows(lu, max_row, k);

		for (int i = k + 1; i < rows; i++) {
			lu[i][k] /= lu[k][k];

			for (int j = k + 1; j < columns; j++)
				lu[i][j] -= lu[i][k] * lu[k][j];
		}
	}
}

void SolveLUP(double** lu, int* p, int rows, int columns, double* v, double* x) {
	double* b = new double[rows];

	CopyVector(v, b, rows);
	Rearrange(b, p, rows);

	double* y = new double[columns];

	for (int i = 0; i < columns; i++) {
		double sum = 0;

		for (int j = 0; j < i; j++)
			sum += lu[i][j] * y[j];

		y[i] = b[i] - sum;
	}

	for (int i = columns - 1; i >= 0; i--) {
		double sum = 0;

		for (int j = i + 1; j < columns; j++)
			sum += lu[i][j] * x[j];
		
		x[i] = (y[i] - sum) / lu[i][i];
	}

	delete[] y;
	delete[] b;
}