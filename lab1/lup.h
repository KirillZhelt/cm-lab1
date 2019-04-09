#pragma once

void BuildLUP(double** m, int rows, int columns, double** lu, int* p);

void SolveLUP(double** lu, int* p, int rows, int columns, double* v, double* x);