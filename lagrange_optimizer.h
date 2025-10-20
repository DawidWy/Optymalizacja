#ifndef LAGRANGE_OPTIMIZER_H
#define LAGRANGE_OPTIMIZER_H

#include <functional>
#include <tuple>

class LagrangeOptimizer {
public:
    LagrangeOptimizer(std::function<double(double)> obj_func,
                      double eps = 1e-6,
                      int max_evals = 200,
                      double tol_A = 1e-12);

    std::pair<double,double> optimize(double a, double b, double c);

    void setEpsilon(double e);
    void setMaxEvals(int m);
    void setTolA(double t);

private:
    std::function<double(double)> f;
    double eps;
    int max_evals;
    double tolA;

    double eval_cached(double x);
    void sort3(double &x1, double &x2, double &x3, double &f1, double &f2, double &f3);

    std::pair<double,bool> lagrange_quad_argmin(double x1, double x2, double x3,
                                                double f1, double f2, double f3);

    std::pair<double,double> golden_section(double a, double b, int maxCalls);

    double last_x;
    double last_f;
    bool have_last;
};

#endif
