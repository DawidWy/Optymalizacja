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
void lab4_csv();
void lab5();
void lab6();

int main()
{   

	try
	{
        lab4();
		//lab4_csv();
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

	CSVStream Sout("symulacja_lab1.csv", {"x0", "x_L", "x_H", "fib_wynik.x", "fib_wynik.y", "fib_wynik.f_calls", "fib_wynik.flag", "lag_wynik.x", "lag_wynik.y", "lag_wynik.f_calls", "lag_wynik.flag"}, ',');
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
	CSVStream Sout("symulacja_lab2.csv", {"i","a","b","x0","x1","y","calls","czy_global"},',');
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

	CSVStream Sram("symulacja_lab2_ramie.csv",{"t","alpha","omega","M"},',');

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

matrix calculate_probabilities(matrix theta, matrix X) {
    try {
        auto x_size = get_size(X);
        int n = x_size.first;
        int m = x_size.second;
        
        // Prawdopodobieństwa dla każdego przykładu
        matrix probabilities(1, m);
        
        for (int i = 0; i < m; i++) {
            // Oblicz theta^T * x^i
            double z = 0.0;
            for (int j = 0; j < n; j++) {
                z += theta(j) * X(j, i);
            }
            
            // Oblicz h_theta(x^i) = sigmoid(z) - prawdopodobieństwo przyjęcia
            probabilities(0, i) = 1.0 / (1.0 + exp(-z));
        }
        
        return probabilities;
    }
    catch (string ex_info) {
        throw ("matrix calculate_probabilities(...):\n" + ex_info);
    }
}

matrix calculate_predictions(matrix theta, matrix X, double threshold = 0.5) {
    try {
        matrix probabilities = calculate_probabilities(theta, X);
        auto size = get_size(probabilities);
        int m = size.second;
        
        matrix predictions(1, m);
        
        for (int i = 0; i < m; i++) {
            if (probabilities(0, i) >= threshold) {
                predictions(0, i) = 1;
            } else {
                predictions(0, i) = 0;
            }
        }
        
        return predictions;
    }
    catch (string ex_info) {
        throw ("matrix calculate_predictions(...):\n" + ex_info);
    }
}


double calculate_accuracy_percentage(matrix theta, matrix X, matrix Y) {
    try {
        auto x_size = get_size(X);
        int m = x_size.second;  // liczba przykładów
        
        // Oblicz predykcje
        matrix predictions = calculate_predictions(theta, X);
        
        // Policz poprawne klasyfikacje
        int correct = 0;
        for (int i = 0; i < m; i++) {
            if (predictions(0, i) == Y(0, i)) {
                correct++;
            }
        }
        
        // Oblicz procent
        double accuracy_percentage = (static_cast<double>(correct) / m) * 100.0;
        
        return accuracy_percentage;
    }
    catch (string ex_info) {
        throw ("matrix calculate_accuracy_percentage(...):\n" + ex_info);
    }
}


void lab4()
{
	double epsilon = 1e-6;
	double h0 = 0;
	int Nmax = 2147483647;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(-2,2);
	solution grad_result;
	//std::stringstream result;
	matrix ud1 = NAN;
	matrix ud2 = NAN;
	bool h_golden = false;
        const std::vector<std::string> headersT = {
            "x1(0)",      "x2(0)",      "SD_x1",      "SD_x2",      "SD_y",
            "SD_f_calls", "SD_g_calls", "SD_minimum", "CG_x1",      "CG_x2",
            "CG_y",       "CG_f_calls", "CG_g_calls", "CG_minimum", "N_x1",
            "N_x2",       "N_y",        "N_f_calls",  "N_g_calls",  "N_H_calls",
            "N_minimum"};
        CSVStream Sout("symulacja_lab4T.csv", headersT);
	for (int i = 0; i < 3; i++) {
		if (i == 0) h0 = 0.05;
		else if (i == 1) h0 = 0.25;
		else h_golden = true;
		for (int j=0;j<100;j++){
			matrix x0(2, 1);
			x0(0) = x0_dist(gen);
            x0(1) = x0_dist(gen);
			grad_result = SD(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2, h_golden);
			bool is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1))
			 < 0.5 && grad_result.y(0) < 0.1);
			Sout << x0(0) << x0(1) << grad_result.x(0) << grad_result.x(1)
			<< grad_result.y(0) << solution::f_calls << solution::g_calls << (is_global ? "Tak" : "Nie");
			//cout<<i<<"   "<<j<<"   "<<"\n";

			solution::clear_calls();
			grad_result = CG(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2, h_golden);
			is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1))
			 < 0.5 && grad_result.y(0) < 0.1);
			Sout << grad_result.x(0) << grad_result.x(1)
			<< grad_result.y(0) << solution::f_calls << solution::g_calls << (is_global ? "Tak" : "Nie");
			//cout<<i<<"   "<<j<<"   "<<"\n";

			solution::clear_calls();
			grad_result = Newton(ff4T, gf4T, hf4T, x0, h0, epsilon, Nmax, ud1, ud2, h_golden);
			is_global = (abs(grad_result.x(0)) < 0.5 && abs(grad_result.x(1))
			 < 0.5 && grad_result.y(0) < 0.1);
			Sout << grad_result.x(0) << grad_result.x(1)
			<< grad_result.y(0) << solution::f_calls << solution::g_calls << solution::H_calls << (is_global ? "Tak" : "Nie");
			//cout<<i<<"   "<<j<<"   "<<"\n";
			solution::clear_calls();
		}}
	// Sout.close();

	// matrix X, Y;
    // load_data("XData.txt", "YData.txt", X, Y);
	// CSVStream Rout("symulacjaR_lab4.csv");
	
	// epsilon = 1e-6;

	// //Parametry theta
	// matrix theta0(3, 1);
    // theta0(0) = 0;
    // theta0(1) = 0;
    // theta0(2) = 0;

	// // Początkowy koszt
    // matrix initial_cost = hf4R(theta0, X, Y);
	
	// solution::g_calls=0;
	// solution res = CG(hf4R, gf4R, theta0, 0, 1e-14, 10000, X, Y);
	// Rout << res.x(0) << res.x(1) << res.x(2) << m2d(hf4R(res.x, X, Y)) << calculate_accuracy_percentage(res.x, X, Y) << solution::g_calls << "\n";

	// solution::g_calls=0;
	// res = CG(hf4R, gf4R, theta0, 0.001, 1e-14, 10000, X, Y);
	// Rout << res.x(0) << res.x(1) << res.x(2) << m2d(hf4R(res.x, X, Y)) << calculate_accuracy_percentage(res.x, X, Y) << solution::g_calls << "\n";

	// solution::g_calls=0;
	// res = CG(hf4R, gf4R, theta0, 0.0001, 1e-14, 10000, X, Y);
	// Rout << res.x(0) << res.x(1) << res.x(2) << m2d(hf4R(res.x, X, Y)) << calculate_accuracy_percentage(res.x, X, Y) << solution::g_calls << "\n";

	
}

using traj_t = std::vector<matrix>;

traj_t SD_traj(std::function<matrix(matrix, matrix, matrix)> ff,
               matrix (*gf)(matrix, matrix, matrix),
               matrix x0, double h0, double epsilon, int Nmax,
               matrix ud1, matrix ud2, bool h_golden)
{
    traj_t traj;
    try
    {
        solution Xopt;
        solution::clear_calls();
        matrix x = x0;
        int iter = 0;
        double fy_old = NAN;
		traj.push_back(x);
        while (true)
        {
            matrix grad = gf(x, ud1, ud2);
            solution::g_calls++;
            
            // Warunek stopu - norma gradientu
            double grad_norm = norm(grad);
            if (grad_norm < epsilon) {
                break;
            }
            
            // Kierunek najszybszego spadku
            matrix d = -grad;
            
            // Określenie długości kroku
            double h;
            if (h_golden || h0 == 0) {
                // Użyj metody złotego podziału do znalezienia optymalnego kroku
                h = find_step_length(x, d, ff, ud1, ud2, epsilon, Nmax);
            } else {
                // Użyj stałego kroku
                h = h0;
            }
            
            // Zapamiętaj poprzedni punkt i wartość funkcji
            matrix x_old = x;
            
            // Oblicz nową wartość funkcji w starym punkcie
            matrix f_old_mat = ff(x, ud1, ud2);
            double f_old = m2d(f_old_mat);
            solution::f_calls++;
            
            // Wykonaj krok
            x = x + h * d;
			traj.push_back(x);
            
            // Oblicz nową wartość funkcji
            matrix f_new_mat = ff(x, ud1, ud2);
            double f_new = m2d(f_new_mat);
            solution::f_calls++;
            
            // Sprawdź zmianę w x
            double x_change = norm(x - x_old);
            
            // Sprawdź zmianę wartości funkcji
            double f_change = abs(f_new - f_old);
            
            iter++;
            
            // Warunki stopu
            if (x_change < epsilon || f_change < epsilon) {
                break;
            }
            
            // Sprawdzenie maksymalnej liczby wywołań
            if (solution::g_calls > Nmax || solution::f_calls > Nmax) {
                throw std::string("Przekroczono maksymalną liczbę wywołań w metodzie najszybszego spadku.");
            }
            
            // Maksymalna liczba iteracji
            if (iter > Nmax) {
                throw std::string("Przekroczono maksymalną liczbę iteracji w metodzie najszybszego spadku.");
            }
        }
        
        Xopt.x = x;
        Xopt.y = ff(x, ud1, ud2);
        solution::f_calls++;
        
        return traj;
    }
	catch (string ex_info)
    {
        throw("solution SD(...):\n" + ex_info);
    }
}

traj_t CG_traj(std::function<matrix(matrix, matrix, matrix)> ff,
               matrix (*gf)(matrix, matrix, matrix),
               matrix x0, double h0, double epsilon, int Nmax,
               matrix ud1, matrix ud2, bool h_golden)
{
    traj_t traj;
    try {
        solution::clear_calls();
        solution Xopt;
        Xopt.x = x0;
        Xopt.fit_fun(ff, ud1, ud2);
        Xopt.grad(gf, ud1, ud2);

        int n = get_len(x0);
        int iter = 0;
        matrix d = -Xopt.g;

        traj.push_back(Xopt.x);

        while (true) {
            matrix x_old = Xopt.x;
            matrix g_old = Xopt.g;
            double f_old = Xopt.y(0);

            double alpha;
            if (h0 == 0 || h_golden) {
                alpha = find_step_length(x_old, d, ff, ud1, ud2, 1e-6, 10000);
                if (alpha < 1e-10) alpha = 1e-4;
            } else {
                alpha = h0;
            }

            Xopt.x = x_old + alpha * d;

            Xopt.grad(gf, ud1, ud2);
            Xopt.fit_fun(ff, ud1, ud2);

            iter++;
            traj.push_back(Xopt.x);

            double grad_norm = norm(Xopt.g);
            if (grad_norm < epsilon) break;
            if (abs(Xopt.y(0) - f_old) < epsilon) break;
            if (iter >= Nmax) break;
            if (solution::g_calls >= Nmax) break;

            double g_old_dot_g_old = 0.0, g_new_dot_g_new = 0.0;
            for (int i = 0; i < n; ++i) {
                g_old_dot_g_old += g_old(i) * g_old(i);
                g_new_dot_g_new += Xopt.g(i) * Xopt.g(i);
            }
            double beta = 0.0;
            if (g_old_dot_g_old > 0) beta = g_new_dot_g_new / g_old_dot_g_old;

            d = -Xopt.g + beta * d;

            double directional_derivative = 0.0;
            for (int i = 0; i < n; ++i) directional_derivative += Xopt.g(i) * d(i);
            if (directional_derivative >= 0) d = -Xopt.g;
            if (iter % n == 0) d = -Xopt.g;
        }

        return traj;
    }
    catch (std::string ex) {
        std::cerr << "CG_traj exception: " << ex << std::endl;
        return traj;
    }
}

traj_t Newton_traj(std::function<matrix(matrix, matrix, matrix)> ff,
                   matrix (*gf)(matrix, matrix, matrix),
                   matrix (*Hf)(matrix, matrix, matrix),
                   matrix x0, double h0, double epsilon, int Nmax,
                   matrix ud1, matrix ud2, bool h_golden)
{
    traj_t traj;
    try
    {
        solution Xopt;
        solution::clear_calls();
        matrix x = x0;
        int iter = 0;
        traj.push_back(x);
        while (true)
        {
            // Oblicz gradient i hesjan
            matrix g = gf(x, ud1, ud2);
            Xopt.g_calls++;
            matrix H = Hf(x, ud1, ud2);
            Xopt.H_calls++;
            
            // Sprawdź warunek stopu - norma gradientu
            double grad_norm = norm(g);
            if (grad_norm < epsilon) {
                break;
            }
            
            // Oblicz kierunek Newtona: d = -inv(H) * g
            matrix d;
            try {
                // Spróbuj obliczyć kierunek Newtona
                d = -inv(H) * g;
                
                // Sprawdź czy kierunek jest kierunkiem spadku
                // Jeśli g^T * d >= 0, to kierunek nie jest kierunkiem spadku
                matrix gT_d = trans(g) * d;
                if (m2d(gT_d) >= 0) {
                    // Użyj kierunku gradientu zamiast Newtona
                    d = -g;
                }
            } catch (string) {
                // Jeśli obliczenie odwrotności się nie powiodło, użyj gradientu
                d = -g;
            }
            
            // Określenie długości kroku
            double h;
            if (h_golden || h0 == 0) {
                // Użyj metody złotego podziału do znalezienia optymalnego kroku
                h = find_step_length(x, d, ff, ud1, ud2, epsilon, Nmax);
            } else {
                // Użyj stałego kroku (zwykle 1 dla metody Newtona)
                h = h0;
            }
            
            // Zapamiętaj poprzedni punkt
            matrix x_old = x;
            
            // Oblicz wartość funkcji w starym punkcie
            matrix f_old_mat = ff(x, ud1, ud2);
            solution::f_calls++;
            double f_old = m2d(f_old_mat);
            
            // Wykonaj krok
            x = x + h * d;
			traj.push_back(x);
            
            // Oblicz nową wartość funkcji
            matrix f_new_mat = ff(x, ud1, ud2);
            solution::f_calls++;
            double f_new = m2d(f_new_mat);
            
            // Sprawdź zmianę w x
            double x_change = norm(x - x_old);
            
            // Sprawdź zmianę wartości funkcji
            double f_change = abs(f_new - f_old);
            
            iter++;
            
            // Warunki stopu
            if (x_change < epsilon || f_change < epsilon) {
                break;
            }
            
            // Sprawdzenie maksymalnej liczby wywołań
            if (solution::g_calls > Nmax || solution::H_calls > Nmax || solution::f_calls > Nmax) {
                throw std::string("Przekroczono maksymalną liczbę wywołań w metodzie Newtona.");
            }
            
            // Maksymalna liczba iteracji
            if (iter > Nmax) {
                throw std::string("Przekroczono maksymalną liczbę iteracji w metodzie Newtona.");
            }
        }
        
        Xopt.x = x;
        Xopt.y = ff(x, ud1, ud2);
        solution::f_calls++;
        
        return traj;
    }
    catch (std::string ex) {
        std::cerr << "Newton_traj exception: " << ex << std::endl;
        return traj;
    }
}

void save_trajs_to_csv(const std::string &filename,
                       const traj_t &sd05, const traj_t &sd25, const traj_t &sdg,
                       const traj_t &cg05, const traj_t &cg25, const traj_t &cgg,
                       const traj_t &n05, const traj_t &n25, const traj_t &ng)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Nie mozna otworzyc pliku: " << filename << std::endl;
        return;
    }

    out <<
        "SD_x1_005;SD_x2_005;"
        "SD_x1_025;SD_x2_025;"
        "SD_x1_var;SD_x2_var;"
        "CG_x1_005;CG_x2_005;"
        "CG_x1_025;CG_x2_025;"
        "CG_x1_var;CG_x2_var;"
        "Newton_x1_005;Newton_x2_005;"
        "Newton_x1_025;Newton_x2_025;"
        "Newton_x1_var;Newton_x2_var\n";

    size_t max_len = 0;
    auto update_max = [&](const traj_t &t){ if (t.size() > max_len) max_len = t.size(); };
    update_max(sd05); update_max(sd25); update_max(sdg);
    update_max(cg05); update_max(cg25); update_max(cgg);
    update_max(n05); update_max(n25); update_max(ng);

    for (size_t k = 0; k < max_len; ++k) {
        std::stringstream row;

        auto write_pair = [&](const traj_t &t){
            if (k < t.size()) {
                const matrix &m = t[k];
                row << m(0) << ";" << m(1) << ";";
            } else {
                row << ";" << ";";
            }
        };

        // SD
        write_pair(sd05);
        write_pair(sd25);
        write_pair(sdg);
        // CG
        write_pair(cg05);
        write_pair(cg25);
        write_pair(cgg);
        // Newton
        write_pair(n05);
        write_pair(n25);
        write_pair(ng);

        std::string s = row.str();
        out << s << "\n";
    }

    out.close();
    std::cout << "Zapisano plik CSV: " << filename << std::endl;
}

void lab4_csv()
{
    double epsilon = 1e-6;
    int Nmax = 2147483647;
    matrix ud1 = NAN, ud2 = NAN;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-2.0, 2.0);

    matrix x0(2,1);
    x0(0) = dist(gen);
    x0(1) = dist(gen);

    std::cout << "Start x0: (" << x0(0) << ", " << x0(1) << ")\n";

    traj_t sd05 = SD_traj(ff4T, gf4T, x0, 0.05, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t sd25 = SD_traj(ff4T, gf4T, x0, 0.25, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t sdg  = SD_traj(ff4T, gf4T, x0, 0.0,  epsilon, Nmax, ud1, ud2, true);

    solution::clear_calls();
    traj_t cg05 = CG_traj(ff4T, gf4T, x0, 0.05, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t cg25 = CG_traj(ff4T, gf4T, x0, 0.25, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t cgg  = CG_traj(ff4T, gf4T, x0, 0.0,  epsilon, Nmax, ud1, ud2, true);

    solution::clear_calls();
    traj_t n05 = Newton_traj(ff4T, gf4T, hf4T, x0, 0.05, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t n25 = Newton_traj(ff4T, gf4T, hf4T, x0, 0.25, epsilon, Nmax, ud1, ud2, false);
    solution::clear_calls();
    traj_t ng  = Newton_traj(ff4T, gf4T, hf4T, x0, 0.0,  epsilon, Nmax, ud1, ud2, true);

    save_trajs_to_csv("lab4_wykresy.csv", sd05, sd25, sdg, cg05, cg25, cgg, n05, n25, ng);

}

void lab5()
{
}

void lab6()
{
}
