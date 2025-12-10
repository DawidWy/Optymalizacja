#include "matrix.h"
#include <cstdlib>
#include <ostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#include"user_funs.h"
#include <algorithm>

#define M_PI 3.14159265358979323846

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	auto Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y.first);									// d�ugo�� rozwi�zania
	double teta_max = Y.second(0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y.second(i, 0))
			teta_max = Y.second(i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2){
	matrix y;
	y = -cos(0.1 * x(0)) * exp(-pow(0.1*x(0)-2*M_PI, 2)) + 0.002 * pow(0.1 * x(0), 2);
	return y;
}

matrix lab1dY(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3, 1);							// Wektor pochodnych
	double g = 9.81, PA = 2.0, PB = 1.0;		// Parametry
	double a = 0.98, b = 0.63;					// współczynniki lepkości i zwężenia strumienia
	double TA_in = 95.0, TB_in_temp = 20.0;		// temperatury dopływu
	double Fin_B = 10.0 / 1000.0;				// dopływ zewnętrzny do zbiornika B [m^3/s]
	double DB = 36.5665 * 1e-4;					// pole powierzchni odpływu ze zbiornika B [m^2]
	double DA = ud2(0);							// pole powierzchni przekroju DA [m^2]
	double VA = Y(0), VB = Y(1), TB = Y(2);		// obecna objętość i temperatura
	double hA = VA / PA, hB = VB / PB;			// wysokości słupa wody w zbiornikach
	double Fout_A;								// odpływ ze zbiornika A
	double Fout_B;								// odpływ ze zbiornika B
	if (hA > 0){
		Fout_A = a * b * DA * sqrt(2 * g * hA);
	}
	else {
		Fout_A = 0;
	};

	if (hB > 0){ 
		Fout_B = a * b * DB * sqrt(2 * g * hB);
	}
	else {
		Fout_B = 0;
	};

	dY(0) = -Fout_A;							// zmiana objętości wody w zbiorniku A
	dY(1) = Fout_A + Fin_B - Fout_B;			// zmiana objętości wody w zbiorniku B
	dY(2) = (abs(VB) < 1e-9) ? 0 : (Fout_A * (TA_in - TB) + Fin_B * (TB_in_temp - TB)) / VB;	// zmiana temperatury
	return dY;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, 1);								//Warunki początkowe
	matrix MT = matrix(1, new double[1]{ m2d(x)*1e-4}); // Pole przekroju DA w m^2
	Y0(0) = 5.0;											// VA = 5 m^3
	Y0(1) = 1.0;											// VB = 1 m^3  
	Y0(2) = 20.0;											// TB = 20 C
	auto Y = solve_ode(lab1dY, 0, 1, 2000, Y0, ud1, MT);	// rozwiązanie równania różniczkowego
	int n = get_len(Y.first);									// długość rozwiązania
	double T_max = 0;
										
	for (int i = 0; i < n; ++i) {							// szukamy maksymalnej temperatury w zbiorniku B
		if (Y.second(i, 2) > T_max) {
			T_max = Y.second(i, 2);
		}
	}

	y = abs(T_max - 50.0);									// wartość funkcji celu (minimalizujemy różnicę od 50°C)
	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);
	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;
	return matrix(1, 1, result);
}

matrix lab2dY(double t, matrix Y, matrix ud1, matrix ud2)
{
    double alpha = Y(0);
    double omega = Y(1);
    double k1 = ud1(0);
    double k2 = ud1(1);
    const double mr = 1.0;
    const double mc = 5.0;
    const double l  = 2.0;
    const double b  = 0.25;
    const double alpha_ref  = 3.141592653589793;
    const double omega_ref  = 0.0;
    const double I = ((1.0 / 3.0) * mr + mc) * l * l;
    double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
    matrix dY(2, 1);
    dY(0) = omega;					//da/dt
    dY(1) = (M - b * omega) / I;    //d^2a/dt^2
    return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0(2, 1);
    Y0(0) = 0.0; // a0
    Y0(1) = 0.0; // da/dt |0

    double t0   = 0.0;
    double dt   = 0.1;
    double tend = 100.0;
    auto Y = solve_ode(lab2dY, t0, dt, tend, Y0, x, NAN);
    const double mr = 1.0;
    const double mc = 5.0;
    const double l  = 2.0;
    const double b  = 0.25;
    const double alpha_ref  = M_PI;
    const double omega_ref  = 0.0;
    const double I = (1.0 / 3.0) * mr * l * l + mc * l * l;
    double Q = 0.0;

    for (int i = 0; i < tend; ++i) {
        double alpha = Y.second(0);
        double omega = Y.second(1);
        double k1 = x(0);
        double k2 = x(1);
        double M  = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
        double e_alpha = alpha_ref - alpha;
        double e_omega = omega_ref - omega;

        Q += (10.0 * e_alpha * e_alpha + e_omega * e_omega + M * M) * dt;
    }
    y = Q;
    return y;
}

matrix ff3T_outside(matrix x, matrix ud1, matrix ud2)
{
    double x0 = x(0); 
    double x1 = x(1);
    double a = ud1(0);
    double c = ud2(0);

    double term = sqrt(pow(x0 / M_PI, 2) + pow(x1 / M_PI, 2));
    
    double denominator = M_PI * term;

    double y = 0;
    if (abs(denominator) < 1e-10) y = 1.0;
    else y = sin(M_PI * term) / denominator;

    double g1 = -x0 + 1.0;
    double g2 = -x1 + 1.0;
    double g3 = sqrt(pow(x0, 2) + pow(x1, 2)) - a;

    double penalty = 0.0;

    if (g1 > 0) penalty += c * pow(g1, 2);
    if (g2 > 0) penalty += c * pow(g2, 2);
    if (g3 > 0) penalty += c * pow(g3, 2);

    return matrix(y + penalty);
}

matrix ff3T_inside(matrix x, matrix ud1, matrix ud2)
{
    double x0 = x(0);
    double x1 = x(1);
    double a = ud1(0);
    double c = ud2(0);

    double g1 = -x0 + 1.0;
    double g2 = -x1 + 1.0;
    double g3 = sqrt(pow(x0, 2) + pow(x1, 2)) - a;

    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return matrix(1e20);
    }

    double term = sqrt(pow(x0 / M_PI, 2) + pow(x1 / M_PI, 2));
    double denominator = M_PI * term;
    
    double y = 0;
    if (abs(denominator) < 1e-10) y = 1.0;
    else y = sin(M_PI * term) / denominator;

    double barrier = 0.0;
    barrier += -1.0 / g1;
    barrier += -1.0 / g2;
    barrier += -1.0 / g3;

    return matrix(y + c * barrier);
}

matrix lab3dY(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(4,1);
    const double m = 0.6; // masa w kilogramach
    const double r = 0.12; //promień w metrach

    const double vx = Y(2); //początkowa prędkość w kierunku x
    const double vy = Y(3); //początkowa prędkość w kierunku y
    const double omega = m2d(ud1); //początkowa rotacja w rad/s

    // Konstansy użyte w obliczeniach
    const double C = 0.47;
    const double S = M_PI * r * r;
    const double rho = 1.2;

    const double Dx = C*rho*S*vx*abs(vx)/2;
    const double Dy = C*rho*S*vy*abs(vy)/2;

    const double Fmx = rho*vy*omega*M_PI*r*r*r;
    const double Fmy = rho*vx*omega*M_PI*r*r*r;
    dY(0) = vx; //dx/dt
    dY(1) = vy; //dy/dt
    dY(2) =-(Dx+Fmx)/m; //dvx/dt
    dY(3) =-(Dy+Fmy)/m-g; //dvy/dt

    return dY;
}

matrix ff3R(matrix x, matrix ud1) {
    
    //x to vx początkowe, ud1 to omega
    matrix Y(4,1);
    Y(0) = 0; //wstępne przesunięcie w metrach
    Y(1) = 100; //wstępna wysokość w metrach
    Y(2) = m2d(x);
    Y(3) = 0;

    pair<matrix,matrix> sol = solve_ode(lab3dY, 0, 0.01, 7, Y, ud1);
    
    int i;
    double x_end;
    double t_end;
    double penalty = 0;

    double y_prev;
    double y_curr;
    double x_prev;
    double x_curr;
    double t_prev;
    double t_curr;
    double alpha;

    int n = get_len(sol.second[1]);

    //Szukamy momentu zderzenia z ziemią
    if (sol.second(n-1,1)<0) {
        for (i=0;sol.second(i,1)>0; i++);

        //interpolacja liniowa
        y_prev = sol.second(i-1,1);
        y_curr = sol.second(i,1);
        x_prev = sol.second(i-1,0);
        x_curr = sol.second(i,0);
        t_prev = sol.first(i-1,0);
        t_curr = sol.first(i,0);

        alpha = -y_prev / (y_curr - y_prev);
        x_end = x_prev + alpha * (x_curr - x_prev);
        t_end = t_prev + alpha * (t_curr - t_prev);
    }
    else { //jeśli go nie ma, zwracamy najniższy punkt
        x_end = sol.second(n-1,0);
        t_end = sol.first(n-1,0);
    }

    // Jeśli piłka nigdy nie minęła y=50, coś poszło nie tak, duże penalty
    if(sol.second(n-1,1)>50) {
        penalty += 5000;
        return -x_end + penalty;
    }

    //Szukamy momentu przecięcia y=50
    for (i=0;sol.second(i,1)>50; i++);
    y_prev = sol.second(i-1,1);
    y_curr = sol.second(i,1);
    x_prev = sol.second(i-1,0);
    x_curr = sol.second(i,0);
    alpha = (50 - y_prev) / (y_curr - y_prev);
    double x_50 = x_prev + alpha * (x_curr - x_prev);
    
    if(x_50 < 3) {
        penalty += 1000 * (3 - x_50)*(3 - x_50);
    }

    if(x_50 > 7) {
        penalty += 1000 * (7 - x_50)*(7 - x_50);
    }
    
    return -x_end + penalty;
}
matrix ff4T(matrix x, matrix ud1, matrix ud2){
    double x1 = x(0);
	double x2 = x(1);
    matrix result = (1.0/6.0) * pow(x1, 6) - 1.05 * pow(x1, 4) + 2.0 * pow(x1, 2) + pow(x2, 2) + x1 * x2;
	return result;
}

double sigmoid(double z) {
    return 1.0 / (1.0 + exp(-z));
}

matrix hf4R(matrix theta, matrix X, matrix Y) {
    try {
        auto x_size = get_size(X);
        int n = x_size.first;   // liczba cech (3)
        int m = x_size.second;  // liczba przykładów (100)
        
        // Obliczanie hipotezy dla każdego przykładu
        matrix h(1, m);
        double J = 0.0;
        
        for (int i = 0; i < m; i++) {
            //theta^T * x^i
            double z = 0.0;
            for (int j = 0; j < n; j++) {
                z += theta(j) * X(j, i);
            }
            
            //h_theta(x^i) = sigmoid(z)
            double h_i = sigmoid(z);
            
            //Składowa funkcji kosztu
            double y_i = Y(0, i);
            
            //Sprawdzanie log(0)
            double epsilon = 1e-15;
            if (h_i < epsilon) h_i = epsilon;
            if (h_i > 1 - epsilon) h_i = 1 - epsilon;
            
            J += y_i * log(h_i) + (1 - y_i) * log(1 - h_i);
        }
        
        //Średni koszt
        J = -J / m;
        
        matrix cost(1, 1);
        cost(0) = J;
        
        return cost;
    }
    catch (string ex_info) {
        throw ("matrix hf4R(...):\n" + ex_info);
    }
}

// Funkcja gradientu
matrix gf4R(matrix theta, matrix X, matrix Y) {
    try {
        auto x_size = get_size(X);
        int n = x_size.first;   // liczba cech (3)
        int m = x_size.second;  // liczba przykładów (100)
        
        // Inicjalizacja gradientu
        matrix grad(n, 1);
        
        // Dla każdego parametru theta_j
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            
            // Dla każdego przykładu treningowego
            for (int i = 0; i < m; i++) {
                //h_theta(x^i)
                double z = 0.0;
                for (int k = 0; k < n; k++) {
                    z += theta(k) * X(k, i);
                }
                double h_i = sigmoid(z);
                
                //(h_theta(x^i) - y^i) * x_j^i
                sum += (h_i - Y(0, i)) * X(j, i);
            }
            
            //średnia
            grad(j) = sum / m;
        }
        
        return grad;
    }
    catch (string ex_info) {
        throw ("matrix gf4R(...):\n" + ex_info);
    }
}