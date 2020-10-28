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

#define DIM 10981

int vecin[DIM][6];
double sol[4 * DIM];
double prev_pos[4 * DIM];
double prev_elongation[DIM][DIM];

struct raza
{
	double x, y;
};
raza r[DIM];
double pressure[DIM];

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();
void sistem_in_echilibru();

bool spin_state[DIM];

double tau = 4000;
double D = 1000;
double dS = 7;
double T = 50;
double Ea = 400;
double L = 0.6;
double mu = 10;
double kappa = 1450;
double k = 8;
double drx = 0.025;
double dry = 0.015; 

int index[DIM];
int timp;

int main(void)
{
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);

	int part;
	double P;
	double subunitar;

	double nHS;

	char cale_in[100] = "D:\\Simulari\\Spin\\Elastic\\Hexagoane\\dreptunghi.dat";
	char cale_out[150] = "D:\\Simulari\\Spin\\Elastic\\Rezultate\\Anizotropie raze\\dreptunghi\\k=8 drx=0.025 dry=0.015";

	char fisier_out1[300], fisier_out2[300];

	ifstream fin(cale_in);
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
		r[i].x = 0.2;
		r[i].y = 0.2;
	}
	fin.close();

	for (int i = 0; i < DIM; i++)
	{
		sol[4 * i] *= 1.04;
		sol[4 * i + 2] *= 1.04;
		r[i].x +=drx;
		r[i].y += dry;
		spin_state[i] = 1;
	}
	//sistem_in_echilibru();
	/*sprintf(fisier_out2, "%s\\dreptunghi_inceput.dat", cale_out);
	ofstream fout2(fisier_out2);
	for (int i = 0; i < DIM; i++)
	{
		fout2 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << spin_state[i] << '\n';
	}
	fout2.close();*/

	nHS = DIM;

	sprintf(fisier_out1, "%s\\nHS.dat", cale_out);
	ofstream fout1(fisier_out1);

	int total_time = 1000000;

	for (timp = 1; timp <= total_time; timp++)
	{
		fout1 << timp << ' ' << nHS / DIM << '\n';
		cout << timp << ' ' << nHS / DIM << '\n';
		if (nHS / DIM < 0.005)
			break;

		std::random_shuffle(std::begin(index), std::end(index));

		for (int i = 0; i < DIM; i++)
		{

			part = index[i];
			//cout << pressure[part] << '\n';
			subunitar = distribution(generator);
			if (spin_state==0) //adica molecula e in low spin
			{
				P = 1 / tau * exp(-((D - T * dS) / (2 * T)))*exp(-((Ea - kappa * pressure[part]) / (T)));;
				cout << P << '\n';

				if (subunitar < P)
				{
					r[part].x += drx;
					r[part].y += dry;
					spin_state[part] = 1;
					nHS++;
				}
			}
			else
			{
				P = 1 / tau * exp((D - T * dS) / (2 * T))*exp(-((Ea + kappa * pressure[part]) / T));
				//cout << P << '\n';

				if (subunitar < P)
				{
					r[part].x -= drx;
					r[part].y -= dry;
					spin_state[part] = 0;
					nHS--;
				}
			}
		}
		sistem_in_echilibru();
		sprintf(fisier_out2, "%s\\timp%d.dat", cale_out, timp);
		ofstream fout3(fisier_out2);
		for (int i = 0; i < DIM; i++)
		{
			fout3 << sol[4 * i] << ' ' << sol[4 * i + 2] << ' ' << spin_state[i] << '\n';
		}
		fout3.close();
	}

	fout1.close();
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
	double Fe, Fex, Fey;
	double Fx, Fy;
	double elongation_x, elongation_y, elongation;
	double prev_dist_x, prev_dist_y;
	double dist_x, dist_y;
	for (int i = 0; i < DIM; i++)
	{
		Fx = 0;
		Fy = 0;
		part = 4 * i;
		pressure[i] = 0;
		for (int j = 0; j < 6; j++)
		{
			neigh = 4 * vecin[i][j];
			if (neigh < 0)
				continue;
			else
			{

				prev_dist_x = fabs(prev_pos[part] - prev_pos[neigh]) - r[i].x - r[vecin[i][j]].x;
				dist_x=fabs(sol[neigh]-sol[part]) - r[i].x - r[vecin[i][j]].x;
				elongation_x = dist_x - prev_dist_x;
				Fex = k * elongation_x;
				
				prev_dist_y = fabs(prev_pos[part+2] - prev_pos[neigh+2]) - r[i].y - r[vecin[i][j]].y;
				dist_y = fabs(sol[neigh+2] - sol[part+2]) - r[i].y - r[vecin[i][j]].y;
				elongation_y = dist_y - prev_dist_y;
				Fey = k * elongation_y;

				Fe = sqrt(Fex * Fex + Fey * Fey);
				elongation = sqrt(elongation_x*elongation_x + elongation_y * elongation_y);


				if (elongation-prev_elongation[i][j] < 0) //inseamna ca e comprimat
					pressure[i] += Fe * (-1);
				else
					pressure[i] += Fe;
				prev_elongation[i][j] = elongation;
				Fx += Fex;
				Fy += Fey;
			}
		}

		f[part] = sol[part + 1];
		f[part + 1] = Fx - mu * sol[part + 1];
		f[part + 2] = sol[part + 3];
		f[part + 3] = Fy - mu * sol[part + 3];
	}

	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
