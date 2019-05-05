#pragma once

const double EPS_RELAXATION = 0.0000000000001;

void SolveRelaxation(double** m, int rows, int columns, double* v, double w, double* x);