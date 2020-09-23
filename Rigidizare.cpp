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
struct vec
{
	int index;
	int type;
};
struct anchor
{
	double x, y, x_init, y_init;
};
struct damping_params
{
	double translatie, rotatie;
};
damping_params damping;

anchor anch[DIM][7];
vec vecin[DIM][6];

double sol[6 * DIM]; //6 * DIM fiindca avem si theta si d theta
double prev_pos[6 * DIM];
double r[DIM];
double pressure[DIM];
double Fx[DIM], Fy[DIM], Mom[DIM];

void tip_vecini();
void update_anch_init_pos(int dot);
void anchors(int dot);

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);

void sistem_in_echilibru();
bool ajuns_la_echilibru();

const double tau = 10;
const double D = 1000;
const double dS = 7;
const double T = 50;
const double Ea = 400;
const double L = 0.6;

const double kappa = 1450;
const double k_elastic = 5;
const double rigid_translatie = 1;
const double rigid_rotatie = 1;

const double Pi = 3.1415926535897932384626433832795;
const double sinPi3 = sin(Pi / 3.0);
const double cosPi3 = cos(Pi / 3.0);
const double sinPi6 = sin(Pi / 6.0);
const double cosPi6 = cos(Pi / 6.0);

int index[DIM]; //pentru shuffle

double x_minim;
int central, molec_2;
double prev_x, prev_y, prev_y2;

int main(void)
{
	damping.translatie = 10;
	damping.rotatie = 1;
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);

	char cale_in[100] = "D:\\Simulari\\Spin\\Rigidizare\\Hexagoane";
	char cale_out[200] = "D:\\Simulari\\Spin\\Rigidizare\\Relaxare\\Hexagon 2791\\k=5 r_tr=1 r_rot=1 damp_tr=10 damp_rot=1";
	char cale_out2[200] = "D:\\Simulari\\Spin\\Rigidizare\\Relaxare\\Hexagon 2791\\k=5 r_tr=1 r_rot=1 damp_tr=10 damp_rot=1\\Config inter";
	char fisier_out[300], fisier_out2[300], fisier_in[200];

	sprintf(fisier_in, "%s\\hexagon%d.dat", cale_in, DIM);

	ifstream fin(fisier_in);
	for (int i = 0; i < DIM; i++)
	{
		fin >> sol[6 * i] >> sol[6 * i + 2];
		for (int j = 0; j < 6; j++)
		{
			fin >> vecin[i][j].index;
			vecin[i][j].index -= 1; //vecinii sunt indexati de la 1
		}
		x_minim = min(x_minim, sol[6 * i]); // pentru fixarea moleculei cu 1 grad de libertate;
		sol[6 * i + 1] = 0;
		sol[6 * i + 3] = 0;
		sol[6 * i + 4] = 0; //theta initial este 0
		sol[6 * i + 5] = 0;
		index[i] = i;
	}
	fin.close();
	for (int i = 0; i < DIM; i++)
	{
		if (abs(sol[6 * i]) < 0.01&&abs(sol[6 * i + 2]) < 0.01)
			central = i;
	}
	prev_x = sol[6 * central];
	prev_y = sol[6 * central + 2];
	for (int i = 0; i < DIM; i++)
	{
		if (abs(sol[6 * i + 2] - prev_y) < 0.01&&i != central)
			molec_2 = i;
	}
	tip_vecini();


	for (int i = 0; i < DIM; i++)
	{
		r[i] = 0.22;
		sol[6 * i] *= 1.04;
		sol[6 * i + 2] *= 1.04;
	}


	double nHS = DIM;
	sprintf(fisier_out, "%s\\Relaxare_nHS.dat", cale_out);
	ofstream fout(fisier_out);

	//*************************************RELAXARE********************************************
	double subunitar;
	double P_LS_HS, P_HS_LS;
	int dot;
	int timp = 1, timp_total = 1000000;
	for (; timp <= timp_total; timp++)
	{
		fout << timp << ' ' << (double)nHS / DIM << '\n';
		cout << timp << ' ' << (double)nHS / DIM << '\n';
		if (double(nHS / DIM) < 0.002)
			exit(0);
		std::random_shuffle(std::begin(index), std::end(index));
		for (int i = 0; i < DIM; i++)
		{
			dot = index[i];
			subunitar = distribution(generator);
			if (r[dot] < 0.21) //inseamna ca e in LS
			{
				P_LS_HS = 1 / tau * exp((D - T * dS) / (2 * T))*exp(-(Ea - kappa * pressure[dot]) / T);
				//cout << P_HS_LS << '\n';
				if (subunitar < P_LS_HS)
				{
					r[dot] = 0.22;
					nHS++;
				}
			}
			else
			{
				P_HS_LS = 1 / tau * exp((D - T * dS) / (2 * T))*exp(-(Ea + kappa * pressure[dot]) / T);
				//cout << P_HS_LS << '\n';
				if (subunitar < P_HS_LS)
				{
					r[dot] = 0.2;
					nHS--;
					//cout << "ARE LOC O TRANZITIE" << '\n';
				}
			}
		}
		sistem_in_echilibru();
		if (!(timp % 2))
		{
			sprintf(fisier_out2, "%s\\timp%d centru_fix.dat", cale_out2, timp);
			ofstream fout2(fisier_out2);
			for (int i = 0; i < DIM; i++)
			{
				fout2 << sol[6 * i] << ' ' << sol[6 * i + 2] << ' ' << pressure[i] << ' ' << r[i] << '\n';
			}
			fout2.close();
		}
	}
	fout.close();
	//ALTE DETALII
	char fisier_out3[300];
	sprintf(fisier_out3, "%s\\Unghiuri finale.dat", cale_out);
	ofstream fout3(fisier_out3);
	for (int i = 0; i < DIM; i++)
		fout3 << sol[6 * i + 4] << '\n';
	fout3.close();
}

void sistem_in_echilibru()
{
	int count = 0;
	double t = 0.0, t1 = 100.0;
	gsl_odeiv2_system sys = { func, jac, 6 * DIM, &damping };
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

	bool ech = 0;
	while (!ech)
	{
		for (int j = 0; j < DIM; j += 1)
		{
			prev_pos[6 * j] = sol[6 * j];
			prev_pos[6 * j + 2] = sol[6 * j + 2];
		}
		count++;
		double ti = count * t1 / 100.0;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, sol);
		if (status != GSL_SUCCESS)
			printf("error, return value=%d\n", status);
		sol[6 * central] = prev_x;
		sol[6 * central + 2] = prev_y;
		//********************MODIFICARE**********************************
		//sol[6 * molec_2 + 2] = prev_y;

		ech = ajuns_la_echilibru();
	}
	gsl_odeiv2_driver_free(d);
}

int func(double t, const double sol[], double f[], void* params)
{
	(void)(t);
	damping_params damping = *(damping_params*)params;
	int part;
	double radical;
	double neigh_x_centru, neigh_y_centru;
	double x_nou, y_nou;
	double theta;
	int neigh, tip;
	char ordin[] = { 0, 4, 3, 5, 2, 0, 1 };
	double Fex = 0, Fey = 0, F_rig = 0, Fe = 0;
	double own_anch_x, own_anch_y, neigh_anch_x, neigh_anch_y, diff_anch_x, diff_anch_y;
	for (int i = 0; i < DIM; i++)
	{
		update_anch_init_pos(i);
		anchors(i);
		Fx[i] = 0;
		Fy[i] = 0;
		Mom[i] = 0;
		pressure[i] = 0;
	}

	for (int dot = 0; dot < DIM; dot++)
	{
		theta = sol[6 * dot + 4];
		for (int j = 0; j < 6; j++)
		{
			Fex = 0;
			Fey = 0;
			neigh = vecin[dot][j].index;
			if (neigh < 0)
				continue;
			tip = vecin[dot][j].type;
			own_anch_x = anch[dot][tip].x;
			own_anch_y = anch[dot][tip].y;
			neigh_anch_x = anch[neigh][7 - tip].x;
			neigh_anch_y = anch[neigh][7 - tip].y;
			diff_anch_x = neigh_anch_x - own_anch_x;
			diff_anch_y = neigh_anch_y - own_anch_y;
			radical = sqrt(diff_anch_x*diff_anch_x + diff_anch_y * diff_anch_y);
			Fe = k_elastic * (radical - L);
			Fex = Fe * diff_anch_x / radical;
			Fey = Fe * diff_anch_y / radical;
			neigh_x_centru = sol[6 * neigh];
			neigh_y_centru = sol[6 * neigh + 2];
			if (tip == 1 || tip == 6)
			{
				x_nou = (neigh_anch_x - neigh_x_centru) / (neigh_anch_y - neigh_y_centru)*
					(own_anch_y - neigh_anch_y) + neigh_anch_x;
				F_rig = rigid_translatie / L * (x_nou - own_anch_x);
				Fex += F_rig;
			}
			else
			{
				y_nou = (neigh_anch_y - neigh_y_centru) / (neigh_anch_x - neigh_x_centru) *
					(own_anch_x - neigh_anch_x) + neigh_anch_y;
				F_rig = rigid_translatie / L * (y_nou - own_anch_y);
				Fey += F_rig;
			}

			Mom[dot] += Fey * r[dot] * cos(Pi / 6 + theta + ordin[tip]) - Fex * r[dot] * sin(Pi / 6 + theta + ordin[tip])
				- rigid_rotatie * r[dot] * sin(theta);
			Fx[dot] += Fex;
			Fy[dot] += Fey;
			pressure[dot] += Fe; // NU LUAM IN CALCUL FORTELE DE TORSIUNE
		}

		part = 6 * dot;
		f[part] = sol[part + 1];
		f[part + 1] = Fx[dot] - damping.translatie * sol[part + 1];
		f[part + 2] = sol[part + 3];
		f[part + 3] = Fy[dot] - damping.translatie * sol[part + 3];
		f[part + 4] = sol[part + 5];
		f[part + 5] = Mom[dot] - damping.rotatie * sol[part + 5];
	}
	return GSL_SUCCESS;
}

void tip_vecini()
{
	int neigh;
	double x_neigh, y_neigh, x, y;
	for (int i = 0; i < DIM; i++)
	{
		x = sol[6 * i];
		y = sol[6 * i + 2];
		for (int j = 0; j < 6; j++)
		{
			neigh = vecin[i][j].index;
			if (neigh < 0)
			{
				vecin[i][j].type = -100000;
				continue;
			}
			x_neigh = sol[6 * neigh];
			y_neigh = sol[6 * neigh + 2];

			if (abs(x_neigh - x) < 0.05 && y_neigh < y)
			{
				vecin[i][j].type = 1;
				continue;
			}
			if (x_neigh < x && y_neigh < y)
			{
				vecin[i][j].type = 2;
				continue;
			}
			if (x_neigh > x && y_neigh < y)
			{
				vecin[i][j].type = 3;
				continue;
			}
			if (x_neigh < x && y_neigh > y)
			{
				vecin[i][j].type = 4;
				continue;
			}
			if (x_neigh > x && y_neigh > y)
			{
				vecin[i][j].type = 5;
				continue;
			}
			if (abs(x_neigh - x) < 0.05 && y_neigh > y)
			{
				vecin[i][j].type = 6;
				continue;
			}
		}
	}
}

void update_anch_init_pos(int dot)
{
	double x_centru = sol[6 * dot];
	double y_centru = sol[6 * dot + 2];
	//OK
	anch[dot][1].x_init = x_centru;
	anch[dot][1].y_init = y_centru - r[dot];
	//OK
	anch[dot][2].x_init = x_centru - r[dot] * cosPi6;
	anch[dot][2].y_init = y_centru - r[dot] * sinPi6;
	//OK
	anch[dot][3].x_init = x_centru + r[dot] * cosPi6;
	anch[dot][3].y_init = y_centru - r[dot] * sinPi6;
	//OK
	anch[dot][4].x_init = x_centru - r[dot] * cosPi6;
	anch[dot][4].y_init = y_centru + r[dot] * sinPi6;
	//OK
	anch[dot][5].x_init = x_centru + r[dot] * cosPi6;
	anch[dot][5].y_init = y_centru + r[dot] * sinPi6;
	//OK
	anch[dot][6].x_init = x_centru;
	anch[dot][6].y_init = y_centru + r[dot];
}

void anchors(int dot)
{
	int i;
	double x_centru, y_centru, theta;
	x_centru = sol[6 * dot];
	y_centru = sol[6 * dot + 2];
	theta = sol[6 * dot + 4];
	for (i = 1; i <= 6; i++)
	{
		anch[dot][i].x = x_centru + (anch[dot][i].x_init - x_centru)*cos(theta) - (anch[dot][i].y_init - y_centru)*sin(theta);
		anch[dot][i].y = y_centru + (anch[dot][i].x_init - x_centru)*sin(theta) + (anch[dot][i].y_init - y_centru)*cos(theta);
	}
}

bool ajuns_la_echilibru()
{
	int i;
	double diffx, diffy;
	for (i = 0; i < DIM; i++)
	{
		diffx = prev_pos[i * 6] - sol[6 * i];
		diffy = prev_pos[i * 6 + 2] - sol[6 * i + 2];
		if (abs(diffx) > 1e-4 || abs(diffy) > 1e-4)
			return 0;
	}
	return 1;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
