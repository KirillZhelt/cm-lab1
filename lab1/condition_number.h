#pragma once

#include <cmath>


double CountLNorm(double** m, int rows, int columns);

void GaussJordanInverse(double** m, int rows, int columns, double** inverse);

// ����� ���������������, � ����� ���� ������� ����� �����������, �� �� ����� �������� �� ��, 
// ���� ������� ����� ��������
double CountConditionNumber(double** m, int rows, int columns);