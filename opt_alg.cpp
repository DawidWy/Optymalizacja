#include"opt_alg.h"
#include <algorithm>
#include <stdexcept>
#include <system_error>
#include<vector>
#include<utility>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		vector <double> x = {x0};
		double i = 0;
		x[1] = x0 + d;
		if(lab1(x[1]) == lab1(x[0]))
		{
			p[0] = x[0];
			p[1] = x[1];
			return p;
		}
		if(lab1(x[1])>lab1(x[0]))
		{
			d = -d;
			x[1] = x0+d;
			if(lab1(x[1])>=lab1(x[0])){
				p[0] = x[1];
				p[1] = x[0] -d;
				return p;
			}
		}
		while (lab1(x[i]) <= lab1(x[i+1]))
		{
			if(i>Nmax) throw "Przekroczono limit wywołań";
			i++;
			x[i+1] = x[0] + pow(alpha, i)*d;
		}
		if(d>0)
		{
			p[0] = x[i-1];
			p[1] = x[i+1];
		}
		p[0] = x[i+1];
		p[1] = x[i-1];
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
		// Generacja ciągu Fibbonacciego
        vector<int> fibs = {1, 1};
        int k = 2;
        double L = b - a;
        
        while (fibs[k-1] < L / epsilon) {
            fibs.push_back(fibs[k-1] + fibs[k-2]);
            k++;
        }
        
        int n = k - 1; // Ilość iteracji
        
        double a0 = a;
        double b0 = b;
        
        // Punkty początkowe
        double c0 = b0 - (double)fibs[n-1] / fibs[n] * (b0 - a0);
        double d0 = a0 + (double)fibs[n-1] / fibs[n] * (b0 - a0);
        
        double fc = ff(matrix(c0), ud1, ud2)(0,0);
        double fd = ff(matrix(d0), ud1, ud2)(0,0);
        
        for (int i = 0; i < n - 1; i++) {
            if (fc < fd) {
                b0 = d0;
                d0 = c0;
                fd = fc;
                c0 = b0 - (double)fibs[n-i-2] / fibs[n-i-1] * (b0 - a0);
                fc = ff(matrix(c0), ud1, ud2)(0,0);
            } else {
                a0 = c0;
                c0 = d0;
                fc = fd;
                d0 = a0 + (double)fibs[n-i-2] / fibs[n-i-1] * (b0 - a0);
                fd = ff(matrix(d0), ud1, ud2)(0,0);
            }
        }
        
        solution Xopt;
        // Zwracamy punkt w połowie znalezionego przedziału
        Xopt.x = (a0 + b0) / 2.0;
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}
//oryginalnie solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
solution Rosen(matrix(*ff)(matrix, matrix, matrix), vector<double> x0, vector<double> s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Wielkość wektorów musi być sobie równa
		if (x0.size() != s0.size()) throw errc::invalid_argument;
		int n = x0.size();
		//Inicjalizacja tablicy kierunków
		matrix d0(n,n);
		for (int i=0; i<n; i++) {d0(i,i)=1.0;}
		//Tablica względnego wydłużenia
		vector<double> lambda(n,0.0);
		//Tablica porażek
		vector<int> p(n,0);
		while (*max_element(s0.begin(), s0.end()) < epsilon) {
			for (int j = 0; j<n; j++) {
				
			}
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
