#pragma once

const double EPS_GMRES_ARNOLDI = 0.000000001;

void SolveGMRESArnoldi(double** m, int rows, int columns, double* v, double* x);