#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <algorithm>
#include <random>

using namespace std;

#define DIM 2791

ifstream fin("hexagon2791.dat");
ofstream fout("nHS.txt");
ofstream fout1("1000MCS.txt");
ofstream fout2("2000MCS.txt");
ofstream fout3("3000MCS.txt");
ofstream fout4("4000MCS.txt");
ofstream fout5("5000MCS.txt");
ofstream fout6("6000MCS.txt");
ofstream fout7("7000MCS.txt");
ofstream fout8("8000MCS.txt");
ofstream fout9("9000MCS.txt");
ofstream fout10("10000MCS.txt");




//adaug apoi fisiere pentru configuratii intermediare


int vecin[DIM][6];
double sol[4 * DIM];
double prev_pos[4 * DIM];
double r[DIM];

double nHS;
int timp;

void reading();
void configuratii_intermediare();

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();
void sistem_in_echilibru();

double probabilityHS_LS(int dot);
double probabilityLS_HS(int dot);

double subunitar;

double t = 0.0, t1 = 100.0;
int count = 0;

double tau = 400;
double D=1000;
double dS=7;
double T=50;
double Ea = 400;
double L = 0.6;
double mu = 10;
double k_elastic = 1450;
int total_time = 1000000;

int index[DIM];

int main(void)
{
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	
	int part;
	double P;

	reading(); 

	for (int i = 0; i < DIM; i++)
	{
		sol[4 * i] *= 1.04;
		sol[4 * i + 2] *= 1.04;
		r[i] = 0.22;
	}

	nHS = DIM;

	for (timp = 1; timp <= total_time; timp++)
	{
		fout << timp << ' ' << (double) nHS / DIM << '\n';
		cout << timp << ' ' << (double) nHS / DIM << '\n';
		if (nHS == 0)
			exit(0);

		std::random_shuffle(std::begin(index), std::end(index));

		for (int i = 0; i < DIM; i++)
		{
			part = index[i];
			subunitar = distribution(generator);
			if (r[part] < 0.21)
			{
				P = probabilityLS_HS(part);
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[part] = 0.22;
					nHS ++;
				}
			}
			else
			{
				P = probabilityHS_LS(part);
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[part] = 0.2;
					nHS --;
				}
			}
		}
		sistem_in_echilibru();
		configuratii_intermediare();
	}

	fin.close();
	fout.close();
	return 0;
}

void configuratii_intermediare()
{
	switch (timp)
	{
	case 1000: 
		for (int i = 0; i < DIM; i++)
		{
			fout1 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout1.close();
		break;
	case 2000:
		for (int i = 0; i < DIM; i++)
		{
			fout2 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout2.close();
		break;
	case 3000:
		for (int i=0; i<DIM; i++)
		{
			fout3 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout3.close();
		break;
	case 4000:
		for (int i = 0; i < DIM; i++)
		{
			fout4 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout4.close();
		break;
	case 5000:
		for (int i = 0; i < DIM; i++)
		{
			fout5 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout5.close();
		break;
	case 6000:
		for (int i = 0; i < DIM; i++)
		{
			fout6 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout6.close();
		break;
	case 7000:
		for (int i = 0; i < DIM; i++)
		{
			fout7 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout7.close();
		break;
	case 8000:
		for (int i = 0; i < DIM; i++)
		{
			fout8 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout8.close();
		break;
	case 9000:
		for (int i = 0; i < DIM; i++)
		{
			fout9 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout9.close();
		break;
	case 10000:
		for (int i = 0; i < DIM; i++)
		{
			fout10 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
		}
		fout10.close();
		break;
	}
}

double probabilityHS_LS(int dot)
{
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
		sum_delta+= elongation;
	}

	return 1 / tau * exp((D - T * dS) / (2 * T))*exp(-((Ea + k_elastic * sum_delta) / T));
}

double probabilityLS_HS(int dot)
{
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

	return 1 / tau * exp(-((D - T * dS) / (2 * T)))*exp(-((Ea - k_elastic * sum_delta) / (T)));
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
		index[i] = i;
	}
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
	bool ech = 0;
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
