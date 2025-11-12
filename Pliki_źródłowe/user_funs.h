#pragma once

#include"ode_solver.h"

const double g = 9.80665;
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix x, matrix ud1, matrix ud2);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix lab1dY(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff3T(matrix x, matrix ud1, matrix ud2);
matrix lab3dY(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix ud2);
