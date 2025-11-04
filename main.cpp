/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include "matrix.h"
#include"opt_alg.h"
#include "solution.h"
#include <cmath>
#include <cstdlib>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1()
{
	std::ofstream Sout("symulacja_lab1.csv");
	double* res = new double[2] {0,0};
	double x0 = -45, d = 5, alpha = 1.1, epsilon = 0.0001, gamma = 0.0001, minimum = 62.74818;
	int Nmax = 10000;
	solution wynik1;
	solution wynik2;
	for(int i=0;i<100;i++){
		res = expansion(ff1T, x0, d, alpha, Nmax);
		x0 = x0 + 1;
		Sout << x0 << "," << res[0] << "," << res[1] << "," << solution::f_calls << ",";
		solution::clear_calls();
		wynik1 = fib(ff1T, res[0], res[1], epsilon);
		Sout << wynik1.x(0,0) << "," << wynik1.y(0,0) << "," << solution::f_calls;
		solution::clear_calls();
		if (abs(wynik1.x(0,0)-minimum)<0.001) {
			Sout << ",globalne,";
		}
		else {
			Sout << ",lokalne,";
		}
		wynik2 = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
		Sout << wynik2.x(0,0) << "," << wynik2.y(0,0) << "," << solution::f_calls;
		solution::clear_calls();
		if (abs(wynik1.x(0,0)-minimum)<0.001) {
			Sout << ",globalne\n";
		}
		else {
			Sout << ",lokalne\n";
		}
	}

	cout << ff1R(2,matrix(),matrix());

}

void lab2(){
	srand(time(NULL));
	std::ofstream Sout("symulacja_lab2.csv");
	matrix X;
	double step = 0.01, alpha = 0.8, beta = 0.1, epsilon = 0.0001;
	double a, b;
	int Nmax = 1000;
	for (int i = 0; i < 100; i++)
	{
	a = ((rand() % 200) / 100.0) - 1;
	b = ((rand() % 200) / 100.0) - 1;
	alpha = 0.8;
	X = matrix(2, new double[2] {a, b});
	solution hooke = HJ(ff3T, X, step, alpha, epsilon, Nmax);
	Sout << "x" << a << ";" << "x" << b << ";" << "x" <<
	hooke.x(0) << ";" << "x" << hooke.x(1) << ";" << "x" << hooke.y <<
	";" << "x" << solution::f_calls << ";\n";
	cout << hooke;
	}
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
