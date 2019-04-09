#pragma once

#include <cmath>


double CountLNorm(double** m, int rows, int columns);

void GaussJordanInverse(double** m, int rows, int columns, double** inverse);

// число обусловленности, в реале если матрица плохо обусловлена, то мы можем получить не то, 
// если считаем плохо обратную
double CountConditionNumber(double** m, int rows, int columns);