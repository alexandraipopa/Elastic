#define _CRT_SECURE_NO_WARNINGS
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

int vecin[DIM][6];
double sol[4 * DIM];
double prev_pos[4 * DIM];
double r[DIM];

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();
void sistem_in_echilibru();

double probabilityHS_LS(int dot);
double probabilityLS_HS(int dot);

double tau = 400;
double D = 1000;
double dS = 7;
double T = 50;
double Ea = 400;
double L = 0.6;
double mu = 10;
double kappa = 5;
double k_elastic = 1450;

int index[DIM];

int main(void)
{
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);

	int part;
	double P;
	double subunitar;

	double nHS;
	int timp;


	char cale_in[100] = "D:\\Simulari\\Spin\\Elastic\\Hexagoane";
	char cale_out[150] = "D:\\Simulari\\Spin\\Elastic\\Rezultate\\Relaxare\\hexagon2791\\kappa_5";

	char fisier_in[200], fisier_out1[300], fisier_out2[300];

	sprintf(fisier_in, "%s\\hexagon%d.dat", cale_in, DIM);

	ifstream fin(fisier_in);
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
	fin.close();

	for (int i = 0; i < DIM; i++)
	{
		sol[4 * i] *= 1.04;
		sol[4 * i + 2] *= 1.04;
		r[i] = 0.22;
	}

	sprintf(fisier_out1, "%s\\hexagon%d_nHS.dat", cale_out, DIM);
	ofstream fout1(fisier_out1);

	nHS = DIM;
	int total_time = 1000000;

	for (timp = 1; timp <= total_time; timp++)
	{
		fout1 << timp << ' ' << nHS / DIM << '\n';
		cout << timp << ' ' << nHS / DIM << '\n';
		if (nHS < 0.002)
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
					nHS++;
				}
			}
			else
			{
				P = probabilityHS_LS(part);
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[part] = 0.2;
					nHS--;
				}
			}
		}
		sistem_in_echilibru();
		if (double (nHS/DIM) < 0.7 && double (nHS/DIM)>0.6 )
		{
			sprintf(fisier_out2, "%s\\hexagon%d_timp%d_nHS%f.dat", cale_out, DIM, timp, nHS / DIM);
			ofstream fout2(fisier_out2);
			for (int i = 0; i < DIM; i++)
			{
				fout2 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << r[i] << '\n';
			}
			fout2.close();
		}
	}
	fout1.close();
	return 0;
}

double probabilityHS_LS(int dot)
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

	return 1 / tau * exp((D - T * dS) / (2 * T))*exp(-((Ea + kappa * k_elastic * sum_delta) / T));
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

	return 1 / tau * exp(-((D - T * dS) / (2 * T)))*exp(-((Ea - k_elastic * kappa * sum_delta) / (T)));
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
	double t = 0.0, t1 = 100.0;
	int count = 0;
	bool ech = 0;
	gsl_odeiv2_system sys = { func, jac, 4 * DIM, &mu };
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
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
		f[part + 2] = sol[part + 3];
		f[part + 3] = Fey - mu * sol[part + 3];
	}

	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
