#include "utilities.h"
#include "qr.h"

#include <iostream>

using namespace std;

int Sign(double d) {
	return (d >= 0) - (d < 0);
}

void Transpose(double** m, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = i + 1; j < columns; j++)
			swap(m[i][j], m[j][i]);
	}
}

void BuildQR(double** m, int rows, int columns, double** qr, double* diag_r) {
	// copy transpose
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			qr[i][j] = m[j][i];
	}

	for (int i = 0; i < columns; i++) {
		/* 
			1) итерируемся по каждому столбцу
				1a) находим new_a = -signa1 * ||a|| (diag_r[i] = new_a[i][i])
				1b) находим w = (a - new_a) / ||a - new_a||
				1c) для каждого следующего столбца
					1c-a)  new_a[i] = a[i] - 2 * (a[i], w) * w
		*/
		diag_r[i] = Sign(qr[i][i]) * EuclideanNorm(&qr[i][i], rows - i);

		qr[i][i] -= diag_r[i];
		double norm = EuclideanNorm(&qr[i][i], rows - i);
		qr[i][i] /= norm; // wi


	}

	Transpose(qr, rows, columns); // ? maybe return transposed qr
}

void SolveQR(double** qr, double* diag_r, int rows, int columns, double* v, double* x) {
}
