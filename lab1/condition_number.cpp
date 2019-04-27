#include <cmath>
#include <algorithm>
#include <iostream>

#include "utilities.h"
#include "condition_number.h"
#include "fill.h"

using namespace std;

double CountLNorm(double** m, int rows, int columns) {
	double result = 0;

	double row_sum = 0;
	for (int i = 0; i < rows; i++) {
		row_sum = 0;

		for (int j = 0; j < columns; j++)
			row_sum += fabs(m[i][j]);

		result = max(result, row_sum);
	}

	return result;
}

void GaussJordanInverse(double** m, int rows, int columns, double** inverse) {
	double** a = new double*[rows];

	for (int i = 0; i < rows; i++)
		a[i] = new double[columns];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			a[i][j] = m[i][j];

		inverse[i][i] = 1;
	}

	for (int k = 0; k < rows - 1; k++) {
		/* выбор главного элемента по столбцу

		   диагональный элемент в нашей матрице равен сумме модулей всех остальных элементов строки
		   диагональный элемент будет равен нулю <=> каждый элемент строки равен нулю
		   такое практически невозможно, следовательно деления на ноль не должно произойти

		   запустил раз 10 матрицу 256х256 не было ни одного свапа при выборе главного элемента по столбцу
		*/

		double max = a[k][k];
		int max_row = k;
		for (int j = k + 1; j < rows; j++) {
			if (fabs(a[j][k]) > fabs(max)) {
				max_row = j;
				max = a[j][k];
			 }
		}
		
		SwapRows(a, max_row, k);

		for (int i = k + 1; i < rows; i++) {
			auto l = a[i][k] / a[k][k];

			for (int j = k; j < columns; j++) {
				a[i][j] -= l * a[k][j];
				inverse[i][j] -= l * inverse[k][j];
			}
		}
	}

	for (int i = rows - 1; i >= 0; i--) {
		for (int j = 0; j < columns; j++) {
			inverse[i][j] /= a[i][i];

			for (int k = i - 1; k >= 0; k--)
				inverse[k][j] -= inverse[i][j];
		}
	}

	for (int i = 0; i < rows; i++)
		delete[] a[i];

	delete[] a;
}

double CountConditionNumber(double** m, int rows, int columns) {
	double** inverse = new double*[rows];

	for (int i = 0; i < rows; i++)
		inverse[i] = new double[columns] {0};

	GaussJordanInverse(m, rows, columns, inverse);

	auto cond_number = CountLNorm(m, rows, columns) * CountLNorm(inverse, rows, columns);

	for (int i = 0; i < rows; i++)
		delete[] inverse[i];

	delete[] inverse;

	return cond_number;
}