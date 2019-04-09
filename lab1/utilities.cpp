#include <iostream>

#include "utilities.h"

using namespace std;

void SwapRows(double** m, int k, int j) {
	if (k == j)
		return;

	swap(m[k], m[j]);
}

void SwapColumns(double** m, int rows, int k, int j) {
	if (k == j)
		return;

	for (int i = 0; i < rows; i++)
		swap(m[i][k], m[i][j]);
}

void Subtract(double* v1, double* v2, int length, double* result) {
	for (int i = 0; i < length; i++)
		result[i] = v1[i] - v2[i];
}

void Rearrange(double* v, int* p, int length) {
	double* tmp = new double[length];

	CopyVector(v, tmp, length);
	for (int i = 0; i < length; i++)
		v[p[i]] = tmp[i];

	delete[] tmp;
}

double MaxNorm(double* v, int length) {
	double max = v[0];

	for (int i = 1; i < length; i++) {
		if (fabs(v[i]) > max)
			max = fabs(v[i]);
	}

	return max;
}

void PrintMatrix(double** m, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			cout << m[i][j] << ' ';

		cout << '\n';
	}
}

void PrintVector(double* v, int length) {
	for (int i = 0; i < length; i++)
		cout << v[i] << ' ';

	cout << "\n";
}

void CopyMatrix(double** src, double** dst, int rows, int columns) {
	for (int i = 0; i < rows; i++)
		memcpy(dst[i], src[i], columns * sizeof(double));
}

void CopyVector(double* src, double* dst, int length) {
	memcpy(dst, src, length * sizeof(double));
}