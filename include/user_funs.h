#pragma once

#include"ode_solver.h"

const double g = 9.80665;
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix x, matrix ud1, matrix ud2);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix lab1dY(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2T(matrix x, matrix ud1, matrix ud2);
matrix lab2dY(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix ff3T(matrix x, matrix ud1, matrix ud2);
matrix ff3T_outside(matrix x, matrix ud1, matrix ud2);
matrix ff3T_inside(matrix x, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1);
matrix ff4T(matrix x, matrix ud1, matrix ud2);
matrix gf4T(matrix x, matrix ud1, matrix ud2);
matrix hf4T(matrix x, matrix ud1, matrix ud2);