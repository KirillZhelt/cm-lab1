
#include "utilities.h"
#include "relaxation.h"

void SolveRelaxation(double** m, int rows, int columns, double* v, double w, double* x) {
	CopyVector(v, x, rows);

	double* x_prev = new double[rows];
	double* difference = new double[rows];

	while (true) {
		CopyVector(x, x_prev, rows);

		for (int i = 0; i < columns; i++) {
			double sum = 0;

			for (int j = 0; j < columns; j++) {
				if (j != i)
					sum += m[i][j] * x[j];
			}

			x[i] = (1 - w) * x[i] + (w * (v[i] - sum)) / m[i][i];
		}

		Subtract(x, x_prev, rows, difference);
		if (EuclideanNorm(difference, rows) < EPS)
			break;
	}

	delete[] difference;
	delete[] x_prev;
}