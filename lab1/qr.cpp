#include "utilities.h"
#include "qr.h"

#include <iostream>

using namespace std;

int Sign(double d) {
	return (d < 0) ? -1 : 1;
}

void BuildQR(double** m, int rows, int columns, double** qr, double* diag_r) {
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
		//diag_r[i] =

	}
}

void SolveQR(double** qr, double* diag_r, int rows, int columns, double* v, double* x) {

}
