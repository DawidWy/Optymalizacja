#include "lagrange_optimizer.h"
#include <cmath>
#include <limits>
#include <map>
#include <algorithm>
#include <stdexcept>

LagrangeOptimizer::LagrangeOptimizer(std::function<double(double)> obj_func,
                                     double eps_, int max_evals_, double tol_A_)
    : f(std::move(obj_func)), eps(eps_), max_evals(max_evals_), tolA(tol_A_),
      last_x(0.0), last_f(0.0), have_last(false)
{
    if (!f) throw std::invalid_argument("Brak funkcji celu (nullptr).");
    if (eps <= 0) eps = 1e-6;
    if (max_evals <= 3) max_evals = 20;
    if (tolA <= 0) tolA = 1e-12;
}

void LagrangeOptimizer::setEpsilon(double e) { if (e>0) eps = e; }
void LagrangeOptimizer::setMaxEvals(int m) { if (m>3) max_evals = m; }
void LagrangeOptimizer::setTolA(double t) { if (t>0) tolA = t; }

double LagrangeOptimizer::eval_cached(double x) {
    if (have_last && x == last_x) return last_f;
    double v = f(x);
    last_x = x;
    last_f = v;
    have_last = true;
    return v;
}

void LagrangeOptimizer::sort3(double &x1, double &x2, double &x3,
                              double &f1, double &f2, double &f3)
{

    if (x1 > x2) { std::swap(x1,x2); std::swap(f1,f2); }
    if (x2 > x3) { std::swap(x2,x3); std::swap(f2,f3); }
    if (x1 > x2) { std::swap(x1,x2); std::swap(f1,f2); }
}

std::pair<double,bool> LagrangeOptimizer::lagrange_quad_argmin(double x1, double x2, double x3,
                                                               double f1, double f2, double f3)
{
    // metoda różnic dzielonych
    double denom12 = x2 - x1;
    double denom23 = x3 - x2;
    double denom13 = x3 - x1;

    if (denom12 == 0.0 || denom23 == 0.0 || denom13 == 0.0)
        return {0.0, false};

    double D1 = (f2 - f1) / denom12;
    double D2 = (f3 - f2) / denom23;
    double A = (D2 - D1) / denom13;

    if (std::fabs(A) < tolA) {
        return {0.0, false};
    }

    double B = D1 - A*(x1 + x2);

    double xstar = -B / (2.0 * A);
    return {xstar, true};
}

std::pair<double,double> LagrangeOptimizer::golden_section(double a, double b, int maxCalls) {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double resphi = 2.0 - phi; // = 1/phi^2
    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    double f1 = eval_cached(x1);
    double f2 = eval_cached(x2);
    int calls = 2;
    while (std::fabs(b - a) > eps && calls < maxCalls) {
        if (f1 < f2) {
            b = x2;
            x2 = x1; f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = eval_cached(x1);
        } else {
            a = x1;
            x1 = x2; f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = eval_cached(x2);
        }
        calls++;
    }
    if (f1 < f2) return {x1, f1};
    else return {x2, f2};
}

std::pair<double,double> LagrangeOptimizer::optimize(double a, double b, double c) {
    if (!(a < b && b < c)) throw std::invalid_argument("Wymagane: a < b < c.");

    double x1 = a, x2 = b, x3 = c;
    double f1 = eval_cached(x1), f2 = eval_cached(x2), f3 = eval_cached(x3);
    int calls = 3;

    while (true) {
        sort3(x1,x2,x3,f1,f2,f3);

        if (std::fabs(x3 - x1) < eps) {
            if (f1 <= f2 && f1 <= f3) return {x1, f1};
            if (f2 <= f1 && f2 <= f3) return {x2, f2};
            return {x3, f3};
        }
        if (calls >= max_evals) {
            if (f1 <= f2 && f1 <= f3) return {x1, f1};
            if (f2 <= f1 && f2 <= f3) return {x2, f2};
            return {x3, f3};
        }

        auto [xstar, ok] = lagrange_quad_argmin(x1,x2,x3,f1,f2,f3);

        if (!ok || !(xstar > x1 && xstar < x3)) {
            int remain = max_evals - calls;
            if (remain <= 3) {
                if (f1 <= f2 && f1 <= f3) return {x1,f1};
                if (f2 <= f1 && f2 <= f3) return {x2,f2};
                return {x3,f3};
            }
            auto [xg, fg] = golden_section(x1, x3, remain);
            calls = max_evals;
            return {xg, fg};
        }

        if (xstar <= x1) xstar = x1 + 1e-16*(x3-x1);
        if (xstar >= x3) xstar = x3 - 1e-16*(x3-x1);

        double fstar = eval_cached(xstar);
        calls++;


        if (fstar < f2) {
            if (xstar > x2) {
                x1 = x2; f1 = f2;
                x2 = xstar; f2 = fstar;
            } else {
                x3 = x2; f3 = f2;
                x2 = xstar; f2 = fstar;
            }
        } else {
            if (f1 > f3) {
                x1 = xstar; f1 = fstar;
            } else {
                x3 = xstar; f3 = fstar;
            }
        }

        if (std::fabs(x2 - xstar) < eps) {
            if (f1 <= f2 && f1 <= f3) return {x1,f1};
            if (f2 <= f1 && f2 <= f3) return {x2,f2};
            return {x3,f3};
        }
    }
}
