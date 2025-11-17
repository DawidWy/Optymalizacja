#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#include"user_funs.h"
#include <algorithm>


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

// matrix ff1R(matrix x, matrix ud1, matrix ud2){
// 	matrix Y0(3, 1), y;
// 	Y0(0) = 5.0;
// 	Y0(1) = 1.0;
// 	Y0(2) = 20.0;
// 	double t0 = 0.0;
// 	double tend = 2000.0;
// 	double dt = 1.0;
// 	ud2 = m2d(x);
// 	matrix* S = solve_ode(lab1dY, t0, dt, tend, Y0, ud1, ud2);
// 	int n = get_len(S[0]);									
// 	double T_max = 0;
// 	for (int i = 0; i < n; ++i)									
// 		if (S[1](i, 2) > T_max)
// 			T_max = S[1](i, 2);
// 	y = abs(T_max - 50.0);									
// 	S[0].~matrix();
// 	S[1].~matrix();
// 	return y;
// }

matrix lab1dY(double t, matrix Y, matrix ud1, matrix ud2)
{
	double Fin = 0.01;
	double Pa = 2; 				// Pole podstawy zbiornika A
	double Va0 = 5; 			// Objętość wody w temperaturze Ta0
	double Ta0 = 95; 			// Temperatura wody w C w zbiorniku B
	double Pb = 1; 				// Pole podstawy zbiornika B
	double Vb0 = 1; 			// Objętość wody w temperaturze Tb0
	double Tb0 = 20; 			// Temperatura wody w C w zbiorniku B
	double TinB = 20; 			// Temperatura wlewającej się wody do zbiornika B
	double FinB = 10; 			// Prędkość wlewania się wody do zbiornika B w l/s
	double DB = 0.00365665; 	// Pole przekroju otworu z którego wylewa się woda ze zbiornika B
	double a = 0.98; 			// Współczynnik lepkości cieczy
	double b = 0.63; 			// Współczynnik zawężenia strumienia cieczy
	double t0 = 0;
	double tend = 2000;
	double dt = 1;
	double Tmax = 50; 			// Maksymalna porządana temperatura w zbiorniku
	matrix dY(Y);
	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);
	double Faout = VA > 0 ? a * b * m2d(ud1) * sqrt((2 * g * Va0) / Pa) : 0;
	double Fbout = VB > 0 ? a * b * m2d(ud1) * sqrt(2 * g* Vb0 / Pb) : 0;
	dY(0) = -Faout;
	dY(1) = Faout - Fbout + Fin;
	dY(2) = Faout/VB * (Ta0 - TB) + Fin/VB * (TinB - TB);
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


matrix lab1dY(matrix x, matrix ud1, double a, double b, double Va, double Pa, double Db, double Pb, double Fin, double Tinb, double Ta0){
	matrix y;
	double Faout = a * b * m2d(ud1) * sqrt((2 * g * Va) / Pa);
	y(0) = -1 * Faout;
	y(1) = Faout - a * b * Db * sqrt((2 * g * y(0))/ Pb) + Fin;
	y(2) = Fin/y(0) * (Tinb - y(1)) + Faout/y(0) * (Ta0 - y(1));
	return y;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);
	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;
	return matrix(1, 1, result);
}

matrix lab3dY(double t, matrix Y, matrix ud1, matrix ud2)
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
    const double I = (1.0 / 3.0) * mr * l * l + mc * l * l;
    double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - b * omega) / I;
    return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    double t0   = 0.0;
    double dt   = 0.01;
    double tend = 10.0;
    auto Y = solve_ode(lab3dY, t0, dt, tend, Y0, x, NAN);

    int n = get_len(Y.first);
    const double mr = 1.0;
    const double mc = 5.0;
    const double l  = 2.0;
    const double b  = 0.25;
    const double alpha_ref  = 3.141592653589793;
    const double omega_ref  = 0.0;
    const double I = (1.0 / 3.0) * mr * l * l + mc * l * l;
    double Q = 0.0;

    for (int i = 0; i < n; ++i) {
        double alpha = Y.second(i, 0);
        double omega = Y.second(i, 1);
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