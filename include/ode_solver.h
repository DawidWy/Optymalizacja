//Ten plik nie powinien by� edytowany

#pragma once
#include <utility>
#include <functional>

#include"matrix.h"
#include"user_funs.h"

std::pair<matrix,matrix> solve_ode( std::function<matrix(double, matrix, matrix, matrix)>, double, double, double, matrix, matrix = NAN, matrix = NAN); // throw (string);
