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

#define DIM 19
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
	double translatie, rotatie, mom_inertie;
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

void calc_Fx_Fy_Mom(int dot);

int func(double t, const double sol[], double f[], void* params);
int jac(double t, const double y[], double* dfdy, double dfdt[], void* params);
bool ajuns_la_echilibru();

double tau = 400;
double D = 1000;
double dS = 7;
double T = 50;
double Ea = 400;
double L = 0.6;

double kappa = 1450;
double k_elastic = 20;

double rigid_translatie=10;
double rigid_rotatie=10;

double mom_inertie;

const double Pi = 3.1415926535897932384626433832795;

int index[DIM]; //pentru shuffle

double x_minim;
int central, molec_2;

int main(void)
{
	damping.translatie = 10;
	damping.rotatie = 10;
	
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);

	char cale_in[100] = "D:\\Simulari\\Spin\\Rigidizare\\Hexagoane";
	

	char fisier_in[200];

	sprintf(fisier_in, "%s\\hexagon%d.dat", cale_in, DIM);

	ifstream fin(fisier_in);
	for (int i = 0; i < DIM; i++)
	{
		fin >> sol[6 * i] >> sol[6 * i + 2];
		r[i] = 0.2;
		update_anch_init_pos(i);
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
		index[i] = i; //pentru shuffle
		//std::random_shuffle(std::begin(index), std::end(index));
	}
	fin.close();
	tip_vecini();

	//Moleculele care vor fi fixate:

	for (int i = 0; i < DIM; i++)
	{
		if (abs(sol[6 * i]) < 0.01&&abs(sol[6 * i + 2]) < 0.01)
			central = i;
		//cout << central << '\n';
		if (sol[6 * i]<x_minim / 2 + 0.5 && sol[6 * i] > x_minim / 2 - 0.5 && abs(sol[6 * i + 2]) < 0.01)
			molec_2 = i;
	}

	//MOLECULELE DINTR-O PARTE A HEXAGONULUI SE MODIFICA
	for (int i = 0; i < DIM; i++)
	{
		if (sol[6 * i] < -0.85)
		{
			r[i] = 0.22;
			update_anch_init_pos(i);
		}
	}
	// ECHILIBRUL MECANIC
	double t = 0.0, t1 = 100.0;

	int count = 0;

	gsl_odeiv2_system sys = { func, jac, 6 * DIM, &damping };
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

	bool ech = 0;

	double prev_x, prev_y, prev_y2;

	prev_x = sol[6 * central];
	prev_y = sol[6 * central + 2];
	prev_y2 = sol[6 * molec_2 + 2];

	while (!ech)
	{
		double ti = count * t1 / 100.0;
		count++;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, sol);
		if (status != GSL_SUCCESS)
			printf("error, return value=%d\n", status);
		gsl_odeiv2_driver_free(d);

		sol[6 * central] = prev_x;
		sol[6 * central + 2] = prev_y;
		sol[6 * molec_2 + 2] = prev_y2;

		ech = ajuns_la_echilibru();

		if (count % 2 == 1)
		{
			// ************************** AFISARE *************************************
			char fisier_out[300];
			char cale_out[150] = "D:\\Simulari\\Spin\\Rigidizare\\Rezultate\\Ech_mecanic\\Hexagon 19";
			sprintf(fisier_out, "%s\\pas %d.dat", cale_out, count);
			ofstream fout(fisier_out);
			for (int i = 0; i < DIM; i++)
			{
				fout << sol[6 * i] << ' ' << sol[6 * i + 2] << ' ' << pressure[i] << '\n';
			}
			fout.close();
		}
	}

}

void calc_Fx_Fy_Mom(int dot)
{
	char ordin[] = { 0, 4, 3, 5, 2, 0, 1 };
	double Fex = 0, Fey = 0, F_rig = 0;

	Fx[dot] = 0;
	Fy[dot] = 0;
	Mom[dot] = 0;
	double radical;
	double neigh_x_centru, neigh_y_centru;
	double x_nou, y_nou;
	double theta = sol[6 * dot + 4];
	int neigh;
	int tip;

	for (int j = 0; j < 6; j++)
	{
		neigh = vecin[dot][j].index;
		if (neigh < 0)
			continue;

		tip = vecin[dot][j].type;
		//cout << "TIPUL VECINULUI: " << tip << '\n';
		//calculez radicalul;
		radical = sqrt( (anch[dot][tip].x-anch[neigh][7-tip].x)*(anch[dot][tip].x - anch[neigh][7 - tip].x) +
			(anch[dot][tip].y - anch[neigh][7 - tip].y)*(anch[dot][tip].y - anch[neigh][7 - tip].y) );
		Fex = k_elastic * (radical - L)*(anch[neigh][7 - tip].x - anch[dot][tip].x) / radical;
		Fey = k_elastic * (radical - L)*(anch[neigh][7 - tip].y - anch[dot][tip].y) / radical;
		//cout << "Fortele elastice pe Ox: " << Fex << " si pe Oy: " << Fey << '\n';
		neigh_x_centru = sol[6 * neigh];
		neigh_y_centru = sol[6 * neigh + 2];
		if (tip == 1 || tip == 6)
		{
			x_nou = (anch[neigh][7 - tip].x - neigh_x_centru) / (anch[neigh][7 - tip].y - neigh_y_centru)*
				(anch[dot][tip].y - anch[neigh][7 - tip].y) + anch[neigh][7 - tip].x;
			F_rig = rigid_translatie / L * (x_nou - anch[dot][tip].x);
			//cout << "F_rig: " << F_rig << '\n';
			Fex += F_rig;
		}
		else
			if (tip > 1 && tip < 6)
			{
				y_nou = (anch[neigh][7 - tip].y - neigh_y_centru) / (anch[neigh][7 - tip].x - neigh_x_centru) *
					(anch[dot][tip].x - anch[neigh][7 - tip].x) + anch[neigh][7 - tip].y;
				F_rig = rigid_translatie / L * (y_nou - anch[dot][tip].y);
				//cout << "F_rig: " << F_rig << '\n';
				Fey += F_rig;
			}
			else
				cout << "PANIC!!!!!\n";

		Mom[dot] += Fey * r[dot] * cos(Pi / 6 + theta + ordin[tip]) - Fex * r[dot] * sin(Pi / 6 + theta + ordin[tip])
			- rigid_rotatie * r[dot] * sin(theta);
		Fx[dot] += Fex;
		Fy[dot] += Fey;
	}
	//**************************************************MODIFICARE****************************************************
	cout << "INDEXUL este: " << dot << '\n';
	cout << "Fx este: " << Fx[dot] << " si Fy este: " << Fy[dot] << '\n';
	cout << "Mom este: " << Mom[dot] << '\n';
	pressure[dot] = sqrt(Fx[dot] * Fx[dot] + Fy[dot] * Fy[dot]);
	cout << "PRESIUNEA: " << pressure[dot] << '\n';
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
			if (neigh == -5)
			{
				vecin[i][j].type = -100000;
				continue;
			}
			x_neigh = sol[6 * neigh];
			y_neigh = sol[6 * neigh + 2];
			
			if (x_neigh == x && y_neigh < y)
				vecin[i][j].type = 1;
			
			if (x_neigh < x && y_neigh < y)
				vecin[i][j].type = 2;
			
			if (x_neigh > x && y_neigh < y)
				vecin[i][j].type = 3;

			if (x_neigh < x && y_neigh > y)
				vecin[i][j].type = 4; 
			
			if (x_neigh > x && y_neigh > y)
				vecin[i][j].type = 5;

			if (x == x_neigh && y_neigh > y)
				vecin[i][j].type = 6;
		}
	}
}

void update_anch_init_pos(int dot)
{
	double x_centru = sol[6 * dot];
	double y_centru = sol[6 * dot + 2];

	anch[dot][1].x_init = x_centru;
	anch[dot][1].y_init = y_centru-r[dot];

	anch[dot][2].x_init = x_centru - r[dot] * sin(Pi / 3);
	anch[dot][2].y_init = y_centru - r[dot] * cos(Pi / 3);

	anch[dot][3].x_init = x_centru + r[dot] * sin(Pi / 3);
	anch[dot][3].y_init = y_centru - r[dot] * cos(Pi / 3);

	anch[dot][4].x_init = x_centru - r[dot] * sin(Pi / 3);
	anch[dot][4].y_init = y_centru + r[dot] * cos(Pi / 3);

	anch[dot][5].x_init = x_centru + r[dot] * sin(Pi / 3);
	anch[dot][5].y_init = y_centru + r[dot] * cos(Pi / 3);

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

int func(double t, const double sol[], double f[], void* params)
{
	(void)(t);
	damping_params damping = *(damping_params*)params;
	int part;

	double m=1;
	for (int j = 0; j < DIM; j += 1)
	{
		prev_pos[6 * j] = sol[6 * j];
		prev_pos[6 * j + 2] = sol[6 * j + 2];
		anchors(j);
	}
	for (int j = 0; j < DIM; j++)
		calc_Fx_Fy_Mom(j);

	for (int i = 0; i < DIM; i++)
	{
		part = 6 * i;
		f[part] = sol[part + 1];
		f[part + 1] = Fx[i] - damping.translatie * sol[part + 1];
		f[part + 2] = sol[part + 3];
		f[part + 3] = Fy[i] - damping.translatie * sol[part + 3];
		f[part + 4] = sol[part + 5];
		
		//Nu am reusit sa introduc un moment de inertie pentru rotatie
		f[part + 5] = Mom[i] - damping.rotatie * sol[part + 5];
	}
	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
