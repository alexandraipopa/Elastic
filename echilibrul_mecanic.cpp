#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

#define DIM 37

ifstream fin("hexagon37.txt");
ofstream fout("hexagon37echilibru.txt");

int vecin[DIM][6];
double sol[4 * DIM];
double prev_pos[4 * DIM];
double L = 0.6;
double r[DIM];
double k = 1;

double x, y;

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();


int main(void)
{
	double mu = 10;

	gsl_odeiv2_system sys = { func, jac, 4*DIM, &mu };
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

	int i;
	int j;

	for (i = 0; i < DIM; i++)
	{
		fin >> sol[4 * i] >> sol[4 * i + 2];
		for (int j = 0; j < 6; j++)
		{
			fin >> vecin[i][j];
			vecin[i][j] -= 1; //vecinii sunt indexati de la 1
		}
		sol[4 * i + 1] = 0;
		sol[4 * i + 3] = 0;
		r[i] = 0.2;
	}

	for (i = 0; i < DIM; i++)
		fout << sol[4 * i] << ' ' << sol[4 * i + 2] << '\n';

	//molecula de index 18 va trece in HS - are coordonatele (0, 0), deci e centrala
	r[18] = 0.22;
	
	double t = 0.0, t1 = 100.0;
	bool ech=0;
	int count = 0;

	while (ech==0)
	{
		for (int j = 0; j < DIM; j += 1)
		{
			prev_pos[4 * j] = sol[4 * j];
			prev_pos[4 * j + 2] = sol[4 * j + 2];
		}
		count++;
		cout << count << '\n';
		double ti = count * t1 / 100.0;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, sol);
		if (status != GSL_SUCCESS)
		{
			printf("error, return value=%d\n", status);
			break;
		}

		ech = ajuns_la_echilibru();

		//afisare
		fout << "\n\n\n\nConfiguratia nr: " << count << "\n\n\n\n";	
		for (j = 0; j < DIM; j++)
			fout << sol[4 * j] << ' ' << sol[4 * j + 2] << '\n';

	}
	fout << "\n\n\n\nConfiguratia stabila:\n\n\n\n";
	for (i = 0; i < DIM; i++)
		fout << sol[4 * i] << ' ' << sol[4 * i + 2] << '\n';
	gsl_odeiv2_driver_free(d);
	fout.close();
	fin.close();
	return 0;
}

bool ajuns_la_echilibru()
{
	int i;
	double diffx, diffy;
	for (i = 0; i < DIM; i++)
	{
		diffx = prev_pos[i * 4] - sol[4 * i];
		diffy = prev_pos[i * 4 + 2] - sol[4 * i + 2];
		if (abs(diffx) > 0.01 || abs(diffy) > 0.01)
			return 0;
	}
	return 1;
}


int func(double t, const double sol[], double f[], void* params)
{
	(void)(t); /* avoid unused parameter warning */
	double mu = *(double*)params;

	double deltax, deltay, Fex, Fey;
	for (int i = 0; i < DIM; i++)
	{
		deltax = 0;
		deltay = 0;
		for (int j = 0; j < 6; j++)
			if (vecin[i][j] < 0)
				continue;
			else
			{
				deltax += sol[4 * i] - sol[4 * vecin[i][j]] - r[i] - r[vecin[i][j]] - L;
				//L=lungimea resortului relaxat=0.6
				//r=0.2 in LS
				//r=0.22 in HS
				deltay += sol[4 * i + 2] - sol[4 * vecin[i][j] + 2] - r[i] - r[vecin[i][j]] - L;
			}
		Fex = deltax * k;
		Fey = deltay * k;
		f[4 * i] = sol[4 * i + 1];
		f[4 * i + 1] = Fex - mu * sol[4 * i + 1];
		f[4 * i + 2] = sol[4 * i + 3];
		f[4 * i + 3] = Fey - mu * sol[4 * i + 3];
	}

	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
