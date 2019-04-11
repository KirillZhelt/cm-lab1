#pragma once

void SwapRows(double** m, int k, int j);

void SwapColumns(double** m, int rows, int k, int j);

void Subtract(double* v1, double* v2, int length, double* result);

void Rearrange(double* v, int* p, int length);

double MaxNorm(double* v, int length);

double EuclideanNorm(double* v, int length);

void PrintMatrix(double** m, int rows, int columns);

void PrintVector(double* v, int length);

void CopyMatrix(double** src, double** dst, int rows, int columns);

void CopyVector(double* src, double* dst, int length);

void Multiply(double** m, int rows, int columns,
	double* v, int length, double* b);