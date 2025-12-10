//Ten plik nie powinien by� edytowany

#pragma once
#include <functional>

#include"solution.h"

solution MC(std::function<matrix(matrix,matrix,matrix)> ff, int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);

std::pair<double,double> expansion(std::function<matrix(matrix,matrix,matrix)> ff, double x0, double d, double alpha, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution fib(std::function<matrix(matrix,matrix,matrix)> ff, double a, double b, double epsilon, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution lag(std::function<matrix(matrix,matrix,matrix)> ff, double a, double b, double epsilon, double gamma, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution HJ(std::function<matrix(matrix,matrix,matrix)> ff, matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution HJ_trial(std::function<matrix(matrix,matrix,matrix)> ff, solution XB, double s, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Rosen(std::function<matrix(matrix,matrix,matrix)> ff, const matrix &x0, const matrix &s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution pen(std::function<matrix(matrix,matrix,matrix)> ff, matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution sym_NM(std::function<matrix(matrix,matrix,matrix)> ff, matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution SD(std::function<matrix(matrix,matrix,matrix)> ff, matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN, bool h_golden); // throw (string);
solution CG(std::function<matrix(matrix,matrix,matrix)> ff, matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN, bool h_golden); // throw (string);
solution Newton(std::function<matrix(matrix,matrix,matrix)> ff, matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN, bool h_golden); // throw (string);
solution golden(std::function<matrix(matrix,matrix,matrix)> ff, double a, double b, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution Powell(std::function<matrix(matrix,matrix,matrix)> ff, matrix x0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution EA(std::function<matrix(matrix,matrix,matrix)> ff, int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
