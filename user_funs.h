#pragma once

#include"ode_solver.h"

const double g = 9.80665;
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix lab1(matrix);
matrix ff3T(matrix x, matrix ud1, matrix ud2);