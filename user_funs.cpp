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
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
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

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
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
		double Fout_A = a * b * DA * sqrt(2 * g * hA);
	}
	else {
		double Fout_A = 0;
	};

	if (hB > 0){ 
		double Fout_B = a * b * DB * sqrt(2 * g * hB);
	}
	else {
		double Fout_B = 0;
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
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, MT);	// rozwiązanie równania różniczkowego
	int n = get_len(Y[0]);									// długość rozwiązania
	double T_max = 0;
										
	for (int i = 0; i < n; ++i) {							// szukamy maksymalnej temperatury w zbiorniku B
		if (Y[1](i, 2) > T_max) {
			T_max = Y[1](i, 2);
		}
	}

	y = abs(T_max - 50.0);									// wartość funkcji celu (minimalizujemy różnicę od 50°C)
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
	return y;
}


// matrix lab1dY(matrix x, matrix ud1, double a, double b, double Va, double Pa, double Db, double Pb, double Fin, double Tinb, double Ta0){
// 	matrix y;
// 	double Faout = a * b * m2d(ud1) * sqrt((2 * g * Va) / Pa);
// 	y(0) = -1 * Faout;
// 	y(1) = Faout - a * b * Db * sqrt((2 * g * y(0))/ Pb) + Fin;
// 	y(2) = Fin/y(0) * (Tinb - y(1)) + Faout/y(0) * (Ta0 - y(1));
// }

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);
	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;
	return matrix(1, 1, result);
}