#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>

#include "utilities.h"

#include "fill.h" // task1
#include "condition_number.h" // task2
#include "gauss.h" // task3
#include "lup.h" // task4
#include "cholesky.h" // task5
#include "relaxation.h" // task6
#include "qr.h" // task7
#include "least_squares.h" // task8
#include "gmres.h" // task9
#include "gmres_arnoldi.h" // task10

using namespace std;

const int N = 4;

const int ROWS = 256;
const int COLUMNS = 256;

int main() {
	// 1 TASK (Fill)
	cout << fixed << setprecision(13);

	double** A = new double*[ROWS];

	for (int i = 0; i < ROWS; i++)
		A[i] = new double[COLUMNS];

	Fill(A, ROWS, COLUMNS, N);
	
	double* y = new double[COLUMNS];

	Fill(y, COLUMNS, N);

	double* b = new double[ROWS];

	Multiply(A, ROWS, COLUMNS, y, COLUMNS, b);

	// 2 TASK (Condition number)
	cout << "Condition number: " << CountConditionNumber(A, ROWS, COLUMNS) << endl;

	// 3 TASK (Gauss)
	double* x_gauss = new double[COLUMNS];
	Gauss(A, ROWS, COLUMNS, b, x_gauss);

	/*
	cout << "Vector x_gauss: " << endl;
	PrintVector(x_gauss, COLUMNS);
	cout << endl << endl;

	cout << "Vector y: " << endl;
	PrintVector(y, COLUMNS);
	*/

	double* difference_gauss = new double[COLUMNS];

	Subtract(y, x_gauss, COLUMNS, difference_gauss);
	cout << "gauss max norm: " << MaxNorm(difference_gauss, COLUMNS) << endl;

	// 4 TASK (LUP)
	double** lu = new double*[ROWS];

	for (int i = 0; i < ROWS; i++)
		lu[i] = new double[COLUMNS];

	int* p = new int[ROWS];

	BuildLUP(A, ROWS, COLUMNS, lu, p);

	double* x_lup = new double[COLUMNS];

	SolveLUP(lu, p, ROWS, COLUMNS, b, x_lup);

	double* difference_lup = new double[COLUMNS];
	
	Subtract(y, x_lup, COLUMNS, difference_lup);
	cout << "lup max norm: " << MaxNorm(difference_lup, COLUMNS) << endl;

	// 5 TASK (CHOLESKY)
	double** lt = new double*[ROWS];
	double* d = new double[ROWS];

	for (int i = 0; i < ROWS; i++)
		lt[i] = new double[COLUMNS];

	BuildCholesky(A, ROWS, COLUMNS, lt, d);

	double* x_cholesky = new double[COLUMNS];

	SolveCholesky(lt, d, ROWS, COLUMNS, b, x_cholesky);

	double* difference_cholesky = new double[COLUMNS];

	Subtract(y, x_cholesky, COLUMNS, difference_cholesky);
	cout << "cholesky max norm: " << MaxNorm(difference_cholesky, COLUMNS) << endl;

	// TASK 6 (Relaxation)
	double* x_relaxation = new double[COLUMNS];
	SolveRelaxation(A, ROWS, COLUMNS, b, (N + 1) / 6.0, x_relaxation);

	double* difference_relaxation = new double[COLUMNS];

	Subtract(y, x_relaxation, COLUMNS, difference_relaxation);
	cout << "relaxation max norm: " << MaxNorm(difference_relaxation, COLUMNS) << endl;
	
	// TASK 7 (QR)
	double** qr = new double*[ROWS];
	double* diag_r = new double[ROWS];

	for (int i = 0; i < ROWS; i++)
		qr[i] = new double[COLUMNS];

	BuildQR(A, ROWS, COLUMNS, qr, diag_r);

	double* x_qr = new double[COLUMNS];

	SolveQR(qr, diag_r, ROWS, COLUMNS, b, x_qr);

	double* difference_qr = new double[COLUMNS];

	Subtract(y, x_qr, COLUMNS, difference_qr);
	cout << "qr max norm: " << MaxNorm(difference_qr, COLUMNS) << endl;

	
	// TASK 8 (Least Squares)
	double* x_least_squares = new double[20 * N];
	SolveLeastSquares(A, ROWS, 20 * N, b, x_least_squares);

	double* Ax = new double[ROWS];
	Multiply(A, ROWS, 20 * N, x_least_squares, 20 * N, Ax);

	double* discrepancy = new double[ROWS];
	Subtract(Ax, b, ROWS, discrepancy);
	 
	cout << "euclidean norm for least squares method: " << EuclideanNorm(discrepancy, ROWS) << endl;

	// TASK 9 (GMRES)
	double* x_gmres = new double[COLUMNS];
	SolveGMRES(A, ROWS, COLUMNS, b, x_gmres);

	double* difference_gmres = new double[COLUMNS];

	Subtract(y, x_gmres, COLUMNS, difference_gmres);
	cout << "gmres max norm: " << MaxNorm(difference_gmres, COLUMNS) << endl;

	// TASK 10 (GMRES Arnoldi)
	double* x_gmres_arnoldi = new double[COLUMNS];
	SolveGMRESArnoldi(A, ROWS, COLUMNS, b, x_gmres_arnoldi);

	double* difference_gmres_arnoldi = new double[COLUMNS];

	Subtract(y, x_gmres_arnoldi, COLUMNS, difference_gmres_arnoldi);
	cout << "gmres arnoldi max norm: " << MaxNorm(difference_gmres_arnoldi, COLUMNS) << endl;

	delete[] difference_gmres_arnoldi;
	delete[] x_gmres_arnoldi;

	delete[] difference_gmres;
	delete[] x_gmres;

	delete[] discrepancy;
	delete[] Ax;

	delete[] difference_qr;
	delete[] x_qr;

	delete[] diag_r;

	for (int i = 0; i < ROWS; i++)
		delete[] qr[i];

	delete[] qr;

	delete[] x_least_squares;

	delete[] difference_relaxation;

	delete[] x_relaxation;

	delete[] difference_cholesky;

	delete[] x_cholesky;

	delete[] d;

	for (int i = 0; i < ROWS; i++)
		delete[] lt[i];

	delete[] lt;

	delete[] difference_lup;

	delete[] x_lup;

	delete[] p;

	for (int i = 0; i < ROWS; i++)
		delete[] lu[i];

	delete[] lu;

	delete[] difference_gauss;

	delete[] x_gauss;

	delete[] b;
	delete[] y;

	for (int i = 0; i < ROWS; i++)
		delete[] A[i];

	delete[] A;

	system("pause");

	return 0;
}