#pragma once

const double EPS_GMRES = 0.001;

void SolveGMRES(double** m, int rows, int columns, double* v, double* x);