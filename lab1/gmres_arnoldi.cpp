#include <iostream>

#include "gmres_arnoldi.h"
#include "least_squares.h"
#include "utilities.h"

using namespace std;


void SolveGMRESArnoldi(double** m, int rows, int columns, double* v, double* x) {
	double** Q = new double*[rows];

	for (int i = 0; i < rows; i++)
		Q[i] = new double[columns];

	double** H = new double*[rows];

	for (int i = 0; i < rows; i++)
		H[i] = new double[columns] {};

	double* y = new double[columns];

	double* z = new double[rows];

	double* difference = new double[columns];

	double* hq = new double[rows];

	double* d = new double[rows] {};

	d[0] = EuclideanNorm(v, rows);

	int k = 1;
	for (int i = 0; i < rows; i++)
		Q[i][0] = v[i] / d[0]; // q1

	while (k < 256) {
		// generate new q and h
		for (int i = 0; i < rows; i++) {
			z[i] = 0;

			for (int j = 0; j < columns; j++)
				z[i] += m[i][j] * Q[j][k - 1];
		}

		for (int i = 0; i < k; i++) {
			H[i][k - 1] = 0;

			for (int j = 0; j < rows; j++)
				H[i][k - 1] += z[j] * Q[j][i];

			for (int j = 0; j < rows; j++)
				hq[j] = H[i][k - 1] * Q[j][i];

			Subtract(z, hq, rows, z);
		}

		H[k][k - 1] = EuclideanNorm(z, rows);
		if (H[k][k - 1] < EPS_GMRES_ARNOLDI)
			break;
		
		for (int j = 0; j < rows; j++)
			Q[j][k] = z[j] / H[k][k - 1];

		SolveLeastSquares(H, k + 1, k, d, y);

		Multiply(H, k + 1, k, y, k + 1, difference);

		Subtract(difference, d, k + 1, difference);

		if (EuclideanNorm(difference, k + 1) < EPS_GMRES_ARNOLDI)
			break;
		
		k++;
	}

	Multiply(Q, rows, k, y, rows, x);

	delete[] d;

	delete[] hq;

	delete[] difference;

	delete[] z;

	delete[] y;

	for (int i = 0; i < rows; i++)
		delete[] H[i];

	delete[] H;

	for (int i = 0; i < rows; i++)
		delete[] Q[i];

	delete[] Q;
}


/*
void SolveGMRESArnoldi(double** m, int rows, int columns, double* v, double* x) {
	double** AQ = new double*[rows];

	for (int i = 0; i < rows; i++)
		AQ[i] = new double[columns];

	double** Q = new double*[rows];

	for (int i = 0; i < rows; i++)
		Q[i] = new double[columns];

	double** H = new double*[rows];

	for (int i = 0; i < rows; i++)
		H[i] = new double[columns];

	double* y = new double[columns];
	
	double* z = new double[rows];

	double* difference = new double[columns];

	double* hq = new double[rows];

	double b_norm = EuclideanNorm(v, rows);

	int k = 1;
	for (int i = 0; i < rows; i++)
		Q[i][0] = v[i] / b_norm; // q1

	while (k < 256) {
		for (int i = 0; i < rows; i++) {
			AQ[i][k - 1] = 0;

			for (int j = 0; j < columns; j++)
				AQ[i][k - 1] += m[i][j] * Q[j][k - 1];
		}
		
		SolveLeastSquares(AQ, rows, k, v, y);

		Multiply(Q, rows, k, y, rows, x);

		Multiply(m, rows, columns, x, columns, difference);
		Subtract(difference, v, columns, difference);

		if (EuclideanNorm(difference, columns) < EPS_GMRES_ARNOLDI)
			break;

		// generate qk
		for (int i = 0; i < rows; i++) {
			z[i] = 0;

			for (int j = 0; j < columns; j++)
				z[i] += m[i][j] * Q[j][k - 1];
		}

		for (int i = 0; i < k; i++) {
			H[i][k - 1] = 0;

			for (int j = 0; j < rows; j++)
				H[i][k - 1] += z[j] * Q[j][i];

			for (int j = 0; j < rows; j++)
				hq[j] = H[i][k - 1] * Q[j][i];

			Subtract(z, hq, rows, z);
		}
		
		H[k][k - 1] = EuclideanNorm(z, rows);
		if (H[k][k - 1] == 0)
			break;

		for (int j = 0; j < rows; j++)
			Q[j][k] = z[j] / H[k][k - 1];

		k++;
	}

	delete[] hq;

	delete[] difference;

	delete[] z;
	
	delete[] y;

	for (int i = 0; i < rows; i++)
		delete[] H[i];

	delete[] H;

	for (int i = 0; i < rows; i++)
		delete[] Q[i];

	delete[] Q;

	for (int i = 0; i < rows; i++)
		delete[] AQ[i];

	delete[] AQ;
}
*/