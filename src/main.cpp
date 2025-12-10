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
#include "opt_alg.h"
#include "solution.h"
#include "csv.h"
#include "utils.h"
#include "ode_solver.h"
#include <cmath>
#include <cstdlib>
#include <random>

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
		lab4();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl
			 << endl;
	}
	return 0;
}

void lab0()
{
	// Funkcja testowa
	double epsilon = 1e-2;			  // dok�adno��
	int Nmax = 10000;				  // maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5), // dolne oraz g�rne ograniczenie
		a(2, 1);					  // dok�adne rozwi�zanie optymalne
	solution opt;					  // rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); // wywo�anie procedury optymalizacji
	cout << opt << endl
		 << endl;			 // wypisanie wyniku
	solution::clear_calls(); // wyzerowanie licznik�w

	// Wahadlo
	Nmax = 1000;										// dok�adno��
	epsilon = 1e-2;										// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;										// dolne oraz g�rne ograniczenie
	double teta_opt = 1;								// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt); // wywo�anie procedury optymalizacji
	cout << opt << endl
		 << endl;			 // wypisanie wyniku
	solution::clear_calls(); // wyzerowanie licznik�w

	// Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),							 // Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2]{m2d(opt.x), 0.5});	 // MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	auto Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT); // rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");				 // definiujemy strumie� do pliku .csv
	Sout << hcat(Y.first, Y.second);							 // zapisyjemy wyniki w pliku
	Sout.close();										 // zamykamy strumie�
}

void lab1()
{

	CSVStream Sout("symulacja_lab1.csv", ',', {"x0", "x_L", "x_H", "fib_wynik.x", "fib_wynik.y", "fib_wynik.f_calls", "fib_wynik.flag", "lag_wynik.x", "lag_wynik.y", "lag_wynik.f_calls", "lag_wynik.flag"});
	// Sout << fixed;
	// cout << fixed;
	// Problem teoretyczny
	std::pair<double,double> res = {0,0};
	double d = 1, alpha = 1.1, epsilon = 0.0001, gamma = 0.0001, minimum = 62.74818;
	int Nmax = 10000;
	solution wynik1;
	solution wynik2;
	for (int i = 1; i <= 100; i++)
	{
		res = expansion(ff1R, i, d, alpha, Nmax);
		Sout << i << res.first << res.second;
		solution::clear_calls();
		wynik1 = fib(ff1R, res.first, res.second, epsilon);
		Sout << wynik1.x(0, 0) << wynik1.y(0, 0) << solution::f_calls;
		if (abs(wynik1.x(0, 0) - minimum) < 0.001)
		{
			Sout << "globalne";
		}
		else
		{
			Sout << "lokalne";
		}
		solution::clear_calls();
		wynik2 = lag(ff1R, res.first, res.second, epsilon, gamma, Nmax);
		Sout << wynik2.x(0, 0) << wynik2.y(0, 0) << solution::f_calls;
		if (abs(wynik1.x(0, 0) - minimum) < 0.001)
		{
			Sout << "globalne";
		}
		else
		{
			Sout << "lokalne";
		}
		Sout.newline();
		// cout << x0 << "," << res.first << "," << res.second << "," << solution::f_calls << "\n";
		// Sout << x0 << "," << res.first << "," << res.second << "," << solution::f_calls << "\n";
		// cout <<"Przedzial <"<< res.first << " " << res.second << ">, wywaloania " << solution::f_calls << "\n";
		// wynik = fib(ff1T, res.first, res.second, epsilon);
		// cout<<"Wynik fib : "<<wynik<<"\n";
		cout << wynik2.y(0) << endl;
	}

	// Problem rzeczywsity
	double Pa = 2;			// Pole podstawy zbiornika A
	double Va0 = 5;			// Objętość wody w temperaturze Ta0
	double Ta0 = 95;		// Temperatura wody w C w zbiorniku B
	double Pb = 1;			// Pole podstawy zbiornika B
	double Vb0 = 1;			// Objętość wody w temperaturze Tb0
	double Tb0 = 20;		// Temperatura wody w C w zbiorniku B
	double TinB = 20;		// Temperatura wlewającej się wody do zbiornika B
	double FinB = 10;		// Prędkość wlewania się wody do zbiornika B w l/s
	double DB = 0.00365665; // Pole przekroju otworu z którego wylewa się woda ze zbiornika B
	double a = 0.98;		// Współczynnik lepkości cieczy
	double b = 0.63;		// Współczynnik zawężenia strumienia cieczy
	double t0 = 0;
	double tend = 2000;
	double dt = 1;
	double Tmax = 50; // Maksymalna porządana temperatura w zbiorniku
}

void lab2()
{
	srand(time(NULL));
	CSVStream Sout("symulacja_lab2.csv", ',', {"i","a","b","x0","x1","y","calls","czy_global"});
	matrix X;
	double step = 0.01, alpha = 0.8, beta = 0.1, epsilon = 0.0001;
	double a, b;
	int Nmax = 1000;
	for (int i = 0; i < 100; i++)
	{
		a = ((rand() % 200) / 100.0) - 1;
		b = ((rand() % 200) / 100.0) - 1;
		alpha = 0.9;
		X = matrix(2, new double[2]{a, b});
		solution hooke = HJ(ff2T, X, step, alpha, epsilon, Nmax);
		Sout << i << a  << b  << hooke.x(0)  << hooke.x(1)  << hooke.y << solution::f_calls;
		if (abs(hooke.x(0) - 0.0) < 0.001 && abs(hooke.x(1) - 0.0) < 0.001)
		{
			Sout << "Tak";
		}
		else
		{
			Sout << "Nie";
		}
	}

	matrix x0(2, 1);
	x0 = rand_range({0,0},{20,20});
	step = 1.0;
	alpha = 0.5;
	epsilon = 1e-1;
	Nmax = 1000;
	solution::clear_calls();
	solution opt = HJ(ff2R, x0, step, alpha, epsilon, Nmax);
	cout << "Problem rzeczywisty (ramie):" << endl;
	cout << opt << endl;
	matrix k_opt = opt.x;
	matrix Y0(2, 1);
	Y0(0) = 0.0;
	Y0(1) = 0.0;

	double t0 = 0.0;
	double dt = 0.01;
	double tend = 10.0;
	auto Y = solve_ode(lab2dY, t0, dt, tend, Y0, k_opt, NAN);
	int n = get_len(Y.first);

	CSVStream Sram("symulacja_lab2_ramie.csv",',',{"t","alpha","omega","M"});

	const double mr = 1.0;
	const double mc = 5.0;
	const double l = 2.0;
	const double b_r = 0.25;
	const double alpha_ref = 3.141592653589793;
	const double omega_ref = 0.0;
	const double I = (1.0 / 3.0) * mr * l * l + mc * l * l;

	for (int i = 0; i < n; ++i)
	{
		double t = Y.first(i, 0);
		double alpha_t = Y.second(i, 0);
		double omega_t = Y.second(i, 1);

		double k1 = k_opt(0);
		double k2 = k_opt(1);
		double M = k1 * (alpha_ref - alpha_t) + k2 * (omega_ref - omega_t);

		Sram << t << alpha_t << omega_t  << M;
	}

}

void lab3() {
	double epsilon = 1E-3;
	int Nmax = 10000;
	double c_inside = 100;
	double dc_inside = 0.2;
	double c_outside = 1.0;
	double dc_outside = 1.5;
	std::ofstream Sout("symulacja_lab3.csv");
	Sout << "x0_1;x0_2;x1_out;x2_out;norm_out;y_out;f_calls_out;x1_in;x2_in;norm_"
	"in;y_in;f_calls_in\n";
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist;
	std::stringstream test_ss;
	solution test_sol;
	matrix a = matrix(4.0);
	matrix test_x0{};
	for (int i = 0; i < 3; ++i) {
		if (i == 0) {
			a = matrix(4.0);
			x0_dist = std::uniform_real_distribution<>(1, 2);
		} else if (i == 1) {
			a = matrix(4.4934);
			x0_dist = std::uniform_real_distribution<>(1, sqrt(4.4934));
		} else {
			a = matrix(5.0);
			x0_dist = std::uniform_real_distribution<>(1, sqrt(5));
		}
		for (int j = 0; j < 100; ++j) {
			test_x0 = matrix(2, new double[2]{x0_dist(gen), x0_dist(gen)});
			test_ss << test_x0(0) << ";" << test_x0(1) << ";";
			test_sol = pen(ff3T_outside, test_x0, c_outside, dc_outside, epsilon, Nmax, a);
			test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";"
			<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";"
			<< test_sol.y(0) << ";" << solution::f_calls << ";";
			solution::clear_calls();
			test_sol =
			pen(ff3T_inside, test_x0, c_inside, dc_inside, epsilon, Nmax, a);
			test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";"
			<< sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";"
			<< test_sol.y(0) << ";" << solution::f_calls << "\n";
			solution::clear_calls();
		}
	}
	Sout << test_ss.str();
	Sout.close();
}

void load_data(const string& x_filename, const string& y_filename, 
               matrix& X, matrix& Y) {
    // Wczytywanie danych X (3 cechy x 100 przykładów)
    ifstream x_file(x_filename);
    if (!x_file.is_open()) {
        throw string("Nie można otworzyć pliku " + x_filename);
    }
    
    int n = 3;    // liczba cech + bias term
    int m = 100;  // liczba przykładów
    
    X = matrix(n, m);
    
    // Wczytujemy 3 wiersze z pliku XData.txt
    for (int i = 0; i < n; i++) {
        string line;
        getline(x_file, line);  // wczytaj cały wiersz
        
        // Podziel linię na wartości oddzielone średnikami
        stringstream ss(line);
        string value;
        
        for (int j = 0; j < m; j++) {
            getline(ss, value, ';');
            // Usuń spacje z początku i końca
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            
            // Konwertuj na double i zapisz
            X(i, j) = stod(value);
        }
    }
    x_file.close();
    
    // Wczytywanie danych Y (etykiety 0/1)
    ifstream y_file(y_filename);
    if (!y_file.is_open()) {
        throw string("Nie można otworzyć pliku " + y_filename);
    }
    
    Y = matrix(1, m);
    
    // Wczytujemy 1 wiersz z pliku YData.txt
    string line;
    getline(y_file, line);
    
    // Podziel linię na wartości oddzielone średnikami
    stringstream ss(line);
    string value;
    
    for (int j = 0; j < m; j++) {
        getline(ss, value, ';');
        // Usuń spacje z początku i końca
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        // Konwertuj na double i zapisz
        Y(0, j) = stod(value);
    }
    y_file.close();
}

void lab4()
{
	double epsilon = 1e-6;
	double h0 = 0;
	int Nmax = 2147483647;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(-2,2);
	std::ofstream Sout("symulacja_lab4.csv");
	solution grad_result;
	std::stringstream result;
	matrix ud1 = NAN;
	matrix ud2 = NAN;
	Sout<<"x1(0);x2(0);x1;x2;y;f_calls;g_calls;minimum;x1;x2;y;f_calls;g_calls;minimum;x1;x2;y;f_calls;g_calls;H_calls;minimum;\n";
	for (int i = 0; i < 2; i++) {
		if (i == 0) h0 = 0.05;
		else if (i == 1) h0 = 0.25;
		else h0 = 0;
		for (int j=0;j<=100;j++){
			matrix x0(2, 1);
			x0(0) = x0_dist(gen);
            x0(1) = x0_dist(gen);
			grad_result = SD(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
			bool is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1)) < 0.5 && grad_result.y(0) < 0.1);
			result<<x0(0)<<";"<<x0(1)<<";"<<grad_result.x(0)<<";"<<grad_result.x(1)<<";"<<grad_result.y<<";"<<grad_result.f_calls<<grad_result.g_calls;
			if (is_global) result << ";Tak;";
			else result << ";Nie;";
			cout<<i<<"   "<<j<<"   "<<"\n";
			solution::clear_calls();
			grad_result = CG(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
			is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1)) < 0.5 && grad_result.y(0) < 0.1);
			result<<x0(0)<<";"<<x0(1)<<";"<<grad_result.x(0)<<";"<<grad_result.x(1)<<";"<<grad_result.y<<";"<<solution::f_calls<<solution::g_calls;
			if (is_global) result << ";Tak;";
			else result << ";Nie;";
			cout<<i<<"   "<<j<<"   "<<"\n";
			solution::clear_calls();
			grad_result = Newton(ff4T, gf4T, hf4T, x0, h0, epsilon, Nmax, ud1, ud2);
			is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1)) < 0.5 && grad_result.y(0) < 0.1);
			result<<x0(0)<<";"<<x0(1)<<";"<<grad_result.x(0)<<";"<<grad_result.x(1)<<";"<<grad_result.y<<";"<<solution::f_calls<<solution::g_calls<<";"<<solution::H_calls;
			if (is_global) result << ";Tak\n";
			else result << ";Nie\n";
			cout<<i<<"   "<<j<<"   "<<"\n";
			solution::clear_calls();
		}
	}
	Sout<<result.str();
	Sout.close();
}

void lab5()
{
}

void lab6()
{
}
