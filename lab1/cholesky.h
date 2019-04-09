#pragma once

void BuildCholesky(double** m, int rows, int columns, double** lt, double* d);

void SolveCholesky(double** lt, double* d, int rows, int columns, double* v, double* x);