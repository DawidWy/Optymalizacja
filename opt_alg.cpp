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

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d,
                  double alpha, int Nmax, matrix ud1, matrix ud2) {
  try {
    int i = 0;
    double *p = new double[2]{0, 0};
    solution::clear_calls();
    solution X0(x0);
    solution X1(x0 + d);
    X0.fit_fun(ff);
    X1.fit_fun(ff);
    vector<solution> x_vector;
    x_vector.push_back(X0);
    x_vector.push_back(X1);
    if (X0.y(0) == X1.y(0)) {
      p[0] = X0.x(0);
      p[1] = X1.x(0);
      return p;
    }
    if (X1.y(0) > X0.y(0)) {
      d = -d;
      X1.x(0) = X0.x(0) + d;
      X1.fit_fun(ff);
      x_vector[1] = X1;
      if (X1.y(0) >= X0.y(0)) {
        p[0] = X1.x(0);
        p[1] = X0.x(0) - d;
        return p;
      }
    }
    do {
      if (X0.f_calls > Nmax)
        break;
      i++;
      x_vector.push_back(x0 + (pow(alpha, i) * d));
      x_vector[i + 1].fit_fun(ff);
    } while (x_vector[i].y(0) >= x_vector[i + 1].y(0));

    if (d > 0) {
      p[0] = x_vector[i - 1].x(0);
      p[1] = x_vector[i + 1].x(0);
    } else {
      p[0] = x_vector[i + 1].x(0);
      p[1] = x_vector[i - 1].x(0);
    }
    return p;
  } catch (string ex_info) {
    throw("double* expansion(...):\n" + ex_info);
  }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
		solution Xopt;
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
		Xopt.f_calls++;
        double fd = ff(matrix(d0), ud1, ud2)(0,0);
		Xopt.f_calls++;
        
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
			Xopt.f_calls++;
        }
        
        // Zwracamy punkt w połowie znalezionego przedziału
        Xopt.x = (a0 + b0) / 2.0;
		Xopt.y = ff(Xopt.x,ud1,ud2);
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
		double x1 = a, x2 = b;
		double x3 = (x1 + x2) / 2.0;
		double x_min;
		double fx1 = ff(matrix(x1), ud1, ud2)(0);
		Xopt.f_calls++;
		double fx2 = ff(matrix(x2), ud1, ud2)(0);
		Xopt.f_calls++;
		double fx3 = ff(matrix(x3), ud1, ud2)(0);
		Xopt.f_calls++;

		int i = 0;

		while (true)
		{
			double numerator = ((x2 * x2 - x3 * x3) * fx1 + (x3 * x3 - x1 * x1) * fx2 + (x1 * x1 - x2 * x2) * fx3);
			double denominator = ((x2 - x3) * fx1 + (x3 - x1) * fx2 + (x1 - x2) * fx3);

			if (fabs(denominator) < 1e-12)
			{
				Xopt.flag = -1;
				break;
			}

			x_min = 0.5 * numerator / denominator;

			double f_min = ff(matrix(x_min), ud1, ud2)(0);
			Xopt.f_calls++;

			if (fabs(x_min - x3) < epsilon || i > Nmax)
			{
				Xopt.x = matrix(x_min);
				Xopt.y = f_min;
				Xopt.flag = 1;
				break;
			}

			if (x_min < x3)
			{
				if (f_min < fx3)
				{
					x2 = x3; fx2 = fx3;
					x3 = x_min; fx3 = f_min;
				}
				else
				{
					x1 = x_min; fx1 = f_min;
				}
			}
			else
			{
				if (f_min < fx3)
				{
					x1 = x3; fx1 = fx3;
					x3 = x_min; fx3 = f_min;
				}
				else
				{
					x2 = x_min; fx2 = f_min;
				}
			}

			if (fabs(x2 - x1) < gamma)
			{
				Xopt.x = matrix(x_min);
				Xopt.y = f_min;
				Xopt.flag = 1;
				break;
			}

			i++;
		}

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
        solution::clear_calls();
        solution XB(x0); 
        XB.y = XB.fit_fun(ff, ud1, ud2);
        double f_XB = m2d(XB.y);

        solution X;
        solution XB_old;
        double f_X;

        do
        {
            X = HJ_trial(ff, XB, s, ud1, ud2);
            f_X = m2d(X.y);
            if (solution::f_calls > Nmax)
                throw string("Przekroczono limit wywolan funkcji celu (Nmax).");

            if (f_X < f_XB)
            {
                do
                {
					 XB_old = XB; 
                    XB = X; 
                    f_XB = f_X;
                    matrix x_pattern_vec = 2 * XB.x - XB_old.x;
                    solution X_pattern_start(x_pattern_vec);
                    X_pattern_start.y = X_pattern_start.fit_fun(ff, ud1, ud2); 
                    if (solution::f_calls > Nmax)
                        throw string("Przekroczono limit wywolan funkcji celu (Nmax).");
                    X = HJ_trial(ff, X_pattern_start, s, ud1, ud2);
                    f_X = m2d(X.y);
                    if (solution::f_calls > Nmax)
                        throw string("Przekroczono limit wywolan funkcji celu (Nmax).");

                } while (f_X < f_XB);
            }
            else
            {
                s = alpha * s;
            }
            if (solution::f_calls > Nmax)
                throw string("Przekroczono limit wywolan funkcji celu (Nmax).");
        } while (s >= epsilon);
        return XB;
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
        solution X_current = XB; 
        double f_current = m2d(X_current.y); 
        int* dims = get_size(XB.x);
        int n_rows = dims[0];
        int m_cols = dims[1];
        int n;
        bool is_col_vector = true;

        if (m_cols == 1) {
            n = n_rows;
            is_col_vector = true;
        } else if (n_rows == 1) {
            n = m_cols;
            is_col_vector = false;
        } else {
            n = n_rows;
            is_col_vector = true;
        }
        
        for (int j = 0; j < n; j++)
        {
            matrix ej(n_rows, m_cols, 0.0); 

            if (is_col_vector) {
                ej(j, 0) = 1.0;
            } else {
                ej(0, j) = 1.0;
            }

            solution X_test(X_current.x + s * ej);
            double f_test = m2d(X_test.fit_fun(ff, ud1, ud2));

            if (f_test < f_current)
            {
                X_current = X_test; 
                f_current = f_test; 
            }
            else
            {
                X_test.x = X_current.x - s * ej; 
                f_test = m2d(X_test.fit_fun(ff, ud1, ud2));

                if (f_test < f_current)
                {
                    X_current = X_test;
                    f_current = f_test;
                }
            }
        }
        return X_current;
    }
    catch (string ex_info)
    {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), const matrix& x0, const matrix& s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    
    solution Xopt;
    // Sprawdzenie poprawności danych wejściowych
    if (get_len(x0) != get_len(s0)) 
        throw string("Rosen: x0 i s0 muszą być tych samych wymiarów");
    if (alpha <= 1) 
        throw string("Rosen: alpha musi być > 1");
    if (beta <= 0 || beta >= 1) 
        throw string("Rosen: beta musi być pomiędzy 0 i 1");
    if (epsilon <= 0) 
        throw string("Rosen: epsilon musi być > 0");
    if (Nmax <= 0) 
        throw string("Rosen: Nmax musi być > 0");

    int n = get_len(x0);
    int fcalls = 0;
    
    // Inicjalizacja kierunków - wektory jednostkowe
    vector<matrix> d(n);
    for (int i = 0; i < n; i++) {
        d[i] = matrix(n, 1, 0.0);
        d[i](i, 0) = 1.0;
    }

    matrix lambda(n, 1, 0.0);
    matrix p(n, 1, 0.0);
    matrix s = s0;
    matrix xB = x0;
    matrix x_current = x0;

    int iteration = 0;
    
    do {
        // Główna pętla po kierunkach
        for (int j = 0; j < n; j++) {
            // Sprawdzenie liczby wywołań funkcji
            if (fcalls > Nmax) {
                throw string("rosenbrock: Przekroczono max wywołania");
            }

            // Obliczenie nowego punktu: xB + s[j] * d[j]
            matrix x_new = xB + s(j, 0) * d[j];

            // Ewaluacja funkcji
            double f_new = m2d(ff(x_new, ud1, ud2));
            double f_old = m2d(ff(xB, ud1, ud2));
            fcalls += 2;

            if (f_new < f_old) {
                // Krok udany - ekspansja
                xB = x_new;
                lambda(j, 0) += s(j, 0);
                s(j, 0) *= alpha;
            } else {
                // Krok nieudany - kontrakcja
                s(j, 0) = -beta * s(j, 0);
                p(j, 0) += 1;
            }
        }

        x_current = xB;
        iteration++;

        // Sprawdzenie warunku zmiany bazy
        bool change_basis = true;
        for (int j = 0; j < n; j++) {
            if (lambda(j, 0) == 0.0 || p(j, 0) == 0) {
                change_basis = false;
                break;
            }
        }

        if (change_basis) {
            // Ortogonalizacja Grama-Schmidta względem poprzednich wektorów
            vector<matrix> new_d = d;
            
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < j; k++) {
                    double dot_product = 0.0;
                    double norm_sq = 0.0;
                    
                    for (int i = 0; i < n; i++) {
                        dot_product += new_d[j](i, 0) * new_d[k](i, 0);
                        norm_sq += new_d[k](i, 0) * new_d[k](i, 0);
                    }
                    
                    if (norm_sq > 1e-15) {
                        double scale = dot_product / norm_sq;
                        for (int i = 0; i < n; i++) {
                            new_d[j](i, 0) -= scale * new_d[k](i, 0);
                        }
                    }
                }
                
                // Normalizacja
                double norm_val = 0.0;
                for (int i = 0; i < n; i++) {
                    norm_val += new_d[j](i, 0) * new_d[j](i, 0);
                }
                norm_val = sqrt(norm_val);
                
                if (norm_val > 1e-15) {
                    for (int i = 0; i < n; i++) {
                        new_d[j](i, 0) /= norm_val;
                    }
                }
            }
            
            d = new_d;
            
            // Resetowanie parametrów
            for (int j = 0; j < n; j++) {
                lambda(j, 0) = 0.0;
                p(j, 0) = 0;
                s(j, 0) = s0(j, 0);
            }
        }

        // Sprawdzenie warunku stopu - maksymalny krok
        double max_step = 0.0;
        for (int j = 0; j < n; j++) {
            double step = fabs(s(j, 0));
            if (step > max_step) {
                max_step = step;
            }
        }
        
        if (max_step < epsilon) {
            break;
        }

    } while (true);
    
    
    // Obliczenie wartości funkcji celu w znalezionym punkcie
    double f_value = m2d(ff(x_current, ud1, ud2));
    fcalls++;
    
    // Ustawienie pól solution
    Xopt.x = x_current;
    Xopt.y = f_value;
    Xopt.flag = 0;
    
    solution::f_calls = fcalls;
    
    return Xopt;
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
