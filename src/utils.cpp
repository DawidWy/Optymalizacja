#include "utils.h"
#include <random>

matrix rand_range(matrix lower_inc, matrix upper_exc){
    static std::default_random_engine gen;
    int n = get_len(lower_inc);
    int m = get_len(upper_exc);
    if (n != m)
        throw string("matrixrand_range(matrix,matrix):\nwymiary macierzy musza byc takie same");
    matrix result(n, 1);
    for (int i = 0; i < n; ++i) {
        std::uniform_real_distribution<double> dist(lower_inc(i), upper_exc(i));
        result(i) = dist(gen);
    }
    
    return result;

}