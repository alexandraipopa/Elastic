#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>
#include <math.h>

using namespace std;

#define DIM 37


ifstream fin("hexagon37.dat");
ofstream fout1("nHS.txt");
ofstream fout2("configuratie_initiala.txt");
//adaug apoi fisiere pentru configuratii intermediare


int vecin[DIM][6];
double sol[4 * DIM];
double prev_pos[4 * DIM];

double r[DIM];
bool verificat[DIM];

double x, y;

double nHS;
int timp;
int verif;


void reading();
int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();
double random_number();

void sistem_in_echilibru();

void clearverificate();

double calculnHS();
double probabilityHS_LS(int dot);
double probabilityLS_HS(int dot);

double subunitar;

double t = 0.0, t1 = 100.0;
bool ech = 0;
int count = 0;

double tau = 2;
double D=1000;
double dS=20;
double T=50;
double Ea = 400;
double L = 0.6;
double mu = 10;
double k_elastic = 1450;

int total_time = 100000;

int main(void)
{
	srand(time(0));

	int dot=0;
	double P;


	reading(); 

	for (int i = 0; i < DIM; i++)
	{
		sol[4 * i] *= 1.04;
		sol[4 * i + 2] *= 1.04;
	}

	for (int i = 0; i < DIM; i++)
	{
		fout2 << sol[4 * i] << ' ' << sol[4 * i + 2] << '\n';
	}

	for (timp = 1; timp <= total_time; timp++)
	{
		nHS = calculnHS();
		fout1 << timp << ' ' << nHS << '\n';
		cout << timp << ' ' << nHS << '\n';
		if (nHS == 0)
			exit(0);
		while (verif < DIM)
		{
			dot = rand() % DIM;
			if (verificat[dot] == 1)
				continue;
			verif++;
			verificat[dot] = 1;
			subunitar = random_number();

			if (r[dot] == 0.2)
			{
				P = probabilityLS_HS(dot);
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[dot] = 0.22;
				}
			}
			else
			{
				P = probabilityHS_LS(dot);
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[dot] = 0.2;
				}
			}
		}
		sistem_in_echilibru();
		clearverificate();
	}
	sistem_in_echilibru();

	fin.close();
	return 0;
}

double probabilityHS_LS(int dot)
{
	double P=0;
	double sum_delta=0;
	int neigh;
	int part = 4 * dot;
	double radical;
	double elongation; 

	for (int j = 0; j < 6; j++)
	{
		if (vecin[dot][j] < 0)
			continue;
		neigh = vecin[dot][j] * 4;
		radical = sqrt((sol[neigh] - sol[part]) * (sol[neigh] - sol[part]) + (sol[neigh + 2] - sol[part + 2]) * (sol[neigh + 2] - sol[part + 2]));
		elongation = radical - r[dot] - r[vecin[dot][j]] - L;
		if (radical == 1.04 && r[vecin[dot][j]] == 0.22)
			elongation = 0;
		sum_delta+= elongation;
	}

	P = 1 / tau * exp((D - T * dS) / (2 * T))*exp(-((Ea + k_elastic * sum_delta) / T));
	return P;
}

double probabilityLS_HS(int dot)
{
	double P = 0;
	double sum_delta = 0;
	int neigh;
	int part = 4 * dot;
	double radical;
	double elongation;
	for (int j = 0; j < 6; j++)
	{
		if (vecin[dot][j] < 0)
			continue;
		neigh = vecin[dot][j] * 4;
		radical = sqrt((sol[neigh] - sol[part]) * (sol[neigh] - sol[part]) + (sol[neigh + 2] - sol[part + 2]) * (sol[neigh + 2] - sol[part + 2]));
		elongation = radical - r[dot] - r[vecin[dot][j]] - L;
		sum_delta += elongation;
	}

	P = 1 / tau * exp(-((D - T * dS) / (2 * T)))*exp(-((Ea - k_elastic * sum_delta) / (T)));
	return P;
}

void reading()
{
	for (int i = 0; i < DIM; i++)
	{
		fin >> sol[4 * i] >> sol[4 * i + 2];
		for (int j = 0; j < 6; j++)
		{
			fin >> vecin[i][j];
			vecin[i][j] -= 1; //vecinii sunt indexati de la 1
		}
		sol[4 * i + 1] = 0;
		sol[4 * i + 3] = 0;
		r[i] = 0.22;
	}
}

double calculnHS()
{
	double nHS = 0;
	for (int i = 0; i < DIM; i++)
		if (r[i]==0.22)
			nHS += 1;
	nHS = (double)nHS / DIM;
	return nHS;
}

double random_number()
{
	int a, b;
	double rez=0;
	a = rand();
	b = rand();
	if (a < b)
	{
		rez = (double)a / b;
	}
	else
		rez = (double)b / a;
	return rez;
}

void clearverificate()
{
	int i;
	verif = 0;
	for (i = 0; i < DIM; i++)
		verificat [i] = 0;
}

bool ajuns_la_echilibru()
{
	int i;
	double diffx, diffy;
	for (i = 0; i < DIM; i++)
	{
		diffx = prev_pos[i * 4] - sol[4 * i];
		diffy = prev_pos[i * 4 + 2] - sol[4 * i + 2];
		if (abs(diffx) > 1e-4 || abs(diffy) > 1e-4)
			return 0;
	}
	return 1;
}

void sistem_in_echilibru()
{
	gsl_odeiv2_system sys = { func, jac, 4 * DIM, &mu };
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
	int count = 0;
	while (ech == 0)
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
	}
	gsl_odeiv2_driver_free(d);
}

int func(double t, const double sol[], double f[], void* params)
{
	(void)(t); /* avoid unused parameter warning */
	double mu = *(double*)params;
	int part, neigh;
	double radical, Fe, Fex, Fey;
	for (int i = 0; i < DIM; i++)
	{
		Fex = 0;
		Fey = 0;
		part = 4 * i;
		for (int j = 0; j < 6; j++)
		{
			neigh = 4 * vecin[i][j];
			if (neigh < 0)
				continue;
			else
			{
				radical = sqrt((sol[neigh] - sol[part]) * (sol[neigh] - sol[part]) + (sol[neigh + 2] - sol[part + 2]) * (sol[neigh + 2] - sol[part + 2]));

				Fe = k_elastic * (radical - r[i] - r[vecin[i][j]] - L);
				Fex += Fe * (sol[neigh] - sol[part]) / radical;
				Fey += Fe * (sol[neigh + 2] - sol[part + 2]) / radical;
			}
		}

		f[part] = sol[part + 1];
		f[part + 1] = Fex - mu * sol[part + 1];
		f[part + 2] = sol[4 * i + 3];
		f[part + 3] = Fey - mu * sol[part + 3];
	}

	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
