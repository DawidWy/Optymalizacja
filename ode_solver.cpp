//Ten plik nie powinien by� edytowany

#include"ode_solver.h"

std::pair<matrix,matrix> solve_ode(std::function<matrix(double, matrix, matrix, matrix)> diff, double t0, double dt, double tend, matrix Y0, matrix ud1, matrix ud2)
{
	try
	{
		int N = static_cast<int>(floor((tend - t0) / dt) + 1);
		if (N < 2)
			throw string("matrix* solve_ode(...):\nprzedzial czasu nie jest zdefiniowany poprawnie");
		std::pair<int,int> s = get_size(Y0);
		if (s.second != 1)
			throw string("matrix* solve_ode(...):\nwarunek poczatkowy musi byc wektorem pionowym");
		int n = s.first;

		std::pair<matrix, matrix> S = { matrix(N,1), matrix(n,N) };
		S.first(0) = t0;
		for (int i = 0; i < n; ++i)
			S.second(i, 0) = Y0(i);
		matrix k1(n, 1), k2(n, 1), k3(n, 1), k4(n, 1);
		for (int i = 1; i < N; ++i)
		{
			S.first(i) = S.first(i - 1) + dt;
			k1 = dt * diff(S.first(i - 1), S.second[i - 1], ud1, ud2);
			k2 = dt * diff(S.first(i - 1) + 0.5 * dt, S.second[i - 1] + 0.5 * k1, ud1, ud2);
			k3 = dt * diff(S.first(i - 1) + 0.5 * dt, S.second[i - 1] + 0.5 * k2, ud1, ud2);
			k4 = dt * diff(S.first(i - 1) + dt, S.second[i - 1] + k3, ud1, ud2);
			for (int j = 0; j < n; ++j)
				S.second(j, i) = S.second(j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
		}
		S.second = trans(S.second);
		return S;
	}
	catch (string ex_info)
	{
		throw ("matrix* solve_ode(...):\n" + ex_info);
	}
}
