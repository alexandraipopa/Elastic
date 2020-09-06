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

double tau = 400;
double D = 1000;
double dS = 7;
double T = 50;
double Ea = 400;
double L = 0.6;

double kappa = 1000;
double k_elastic = 10;

double mom_inertie;

const double rigid_translatie = 0;
const double rigid_rotatie = 10;
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
	damping.rotatie = 0;
	std::random_device device;
	std::mt19937 generator(device());
	std::uniform_real_distribution <double> distribution(0.0, 1.0);

	char cale_in[100] = "D:\\Simulari\\Spin\\Rigidizare\\Hexagoane";
	char cale_out[150] = "D:\\Simulari\\Spin\\Rigidizare\\Rezultate\\Ech_mecanic\\Hexagon 19";
	char fisier_out[300], fisier_in[200];

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
		if (abs(sol[6*i+2]-prev_y)<0.01&&i!=central)
			molec_2 = i;
	}
	tip_vecini();

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
		//le fixez manual: 
		/*sol[6 * central] = prev_x;
		sol[6 * central + 2] = prev_y;
		sol[6 * molec_2 + 2] = prev_y;*/

		ech = ajuns_la_echilibru();

		sprintf(fisier_out, "%s\\pas %d.dat", cale_out, count);
		ofstream fout(fisier_out);
		for (int i = 0; i < DIM; i++)
		{
			fout << sol[6 * i] << ' ' << sol[6 * i + 2] << ' ' << pressure[i] << '\n';
		}
		fout.close();
	}
	gsl_odeiv2_driver_free(d);

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

			if (abs(x_neigh - x)<0.05 && y_neigh < y)
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

int func(double t, const double sol[], double f[], void* params)
{
	(void)(t);
	damping_params damping = *(damping_params*)params;
	int part;
	double radical;
	double radical2;
	double neigh_x_centru, neigh_y_centru;
	double x_nou, y_nou;
	double theta;
	int neigh;
	int tip;
	char ordin[] = { 0, 4, 3, 5, 2, 0, 1 };
	double Fex = 0, Fey = 0, F_rig = 0, Fe = 0, Fe2 = 0, Fex2=0, Fey2=0;
	double diff, diff2;
	for (int i=0; i < DIM; i++)
	{
		update_anch_init_pos(i);
		anchors(i);
		Fx[i] = 0;
		Fy[i] = 0;
		Mom[i] = 0;	
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
			radical = sqrt((anch[dot][tip].x - anch[neigh][7 - tip].x)*(anch[dot][tip].x - anch[neigh][7 - tip].x) +
				(anch[dot][tip].y - anch[neigh][7 - tip].y)*(anch[dot][tip].y - anch[neigh][7 - tip].y));
			Fe = kappa * (radical - L);
			Fex = Fe * (anch[neigh][7 - tip].x - anch[dot][tip].x) / radical;
			Fey = Fe * (anch[neigh][7 - tip].y - anch[dot][tip].y) / radical;
			
			radical2 = sqrt((sol[6 * neigh] - sol[6 * dot]) * (sol[6 * neigh] - sol[6 * dot]) + (sol[6 * neigh + 2] - sol[6 * dot + 2]) * (sol[6 * neigh + 2] - sol[6 * dot + 2]));
			Fe2 = kappa * (radical2 - r[dot] - r[neigh] - L);
			Fex2 = Fe2 * (sol[6 * neigh] - sol[6 * dot]) / radical2;
			Fey2 = Fe2 * (sol[6 * neigh + 2] - sol[6 * dot + 2]) / radical2;
			diff = (anch[neigh][7 - tip].y - anch[dot][tip].y) / radical;
			diff2 = (sol[6 * neigh + 2] - sol[6 * dot + 2]) / radical2;
			/*if (abs(diff - diff2) > 0.1)
			{
				cout << "\n\n WARNING!!!\n\n\n";
				cout << "Index: " << dot << " vecin cu indexul: " << neigh << " de tipul: " << tip << '\n';
				cout << "Coordonatele centrului: " << sol[6 * dot] << ' ' << sol[6 * dot + 2] << '\n';
				cout << "Unghiul: " << sol[6 * dot + 4] << '\n';
				cout << "Ancora - coord initiale: " << anch[dot][tip].x_init << ' ' << anch[dot][tip].y_init << '\n';
				cout << "Ancora: "<< anch[dot][tip].x << ' ' << anch[dot][tip].y << "\n\n\n";
				cout << "Coordonatele vecinului: " << sol[6 * neigh] << ' ' << sol[6 * neigh + 2] << '\n';
				cout << "Unghiul vecinului: " << sol[6 * neigh + 4] << '\n';
				cout << "Ancora vecinului - coord initiale: " << anch[neigh][7 - tip].x_init << ' ' << anch[neigh][7 - tip].y_init << '\n';
				cout << "Ancora vecinului: " << anch[neigh][7 - tip].x << ' ' << anch[neigh][7 - tip].y << '\n';
				
			}*/
			
			//cout << "Y diff: " << (anch[neigh][7 - tip].y - anch[dot][tip].y) / radical << "-----  " << (sol[6 * neigh + 2] - sol[6 * dot + 2]) / radical2 << "\n";
			//neigh_x_centru = sol[6 * neigh];
			//neigh_y_centru = sol[6 * neigh + 2];
			//if (tip == 1 || tip == 6)
			//{
			//	x_nou = (anch[neigh][7 - tip].x - neigh_x_centru) / (anch[neigh][7 - tip].y - neigh_y_centru)*
			//		(anch[dot][tip].y - anch[neigh][7 - tip].y) + anch[neigh][7 - tip].x;
			//	F_rig = rigid_translatie / L * (x_nou - anch[dot][tip].x);
			//	//cout << "F_rig: " << F_rig << '\n';
			//	Fex += F_rig;
			//}
			//else //tipul e 2, 3, 4, 5;
			//{
			//	y_nou = (anch[neigh][7 - tip].y - neigh_y_centru) / (anch[neigh][7 - tip].x - neigh_x_centru) *
			//		(anch[dot][tip].x - anch[neigh][7 - tip].x) + anch[neigh][7 - tip].y;
			//	F_rig = rigid_translatie / L * (y_nou - anch[dot][tip].y);
			//	//cout << "F_rig: " << F_rig << '\n';
			//	Fey += F_rig;
			//}

		 //Mom[dot] += Fey * r[dot] * cos(Pi / 6 + theta + ordin[tip]) - Fex * r[dot] * sin(Pi / 6 + theta + ordin[tip])
			//	- rigid_rotatie * r[dot] * sin(theta);
			Fx[dot] += Fex;
			Fy[dot] += Fey;
		}
		pressure[dot] = sqrt(Fx[dot] * Fx[dot] + Fy[dot] * Fy[dot]);
		part = 6 * dot;
		f[part] = sol[part + 1];
		f[part + 1] = Fx[dot] - damping.translatie * sol[part + 1];
		f[part + 2] = sol[part + 3];
		f[part + 3] = Fy[dot] - damping.translatie * sol[part + 3];
		f[part + 4] = 0 /*sol[part + 5]*/;
		f[part + 5] = 0 /*Mom[i] - damping.rotatie * sol[part + 5]*/;
	}
	return GSL_SUCCESS;
}

int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	return GSL_SUCCESS;
}
