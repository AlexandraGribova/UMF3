#include "LU.h"
#include <iostream>
#include <fstream>
#include <vector>
#include<locale.h>
#include<iomanip>
#include<time.h>


using namespace std;
vector<double> DI, AL, AU, V;//размерность, vec, di,ai,L,U - в таком пор€дке
vector<int>AI;
int n;
void input()
{
	string line;
	ifstream in("input.txt"); // окрываем файл дл€ чтени€
	in >> n;
	V.resize(n);
	DI.resize(n);
	AI.resize(n + 1);
	for (int i = 0; i < n; i++)
		in >> V[i];
	for (int i = 0; i < n; i++)
		in >> DI[i];
	for (int i = 0; i < n + 1; i++)
		in >> AI[i];
	AL.resize(AI[n] - 1);
	AU.resize(AI[n] - 1);
	for (int i = 0; i < AI[n] - 1; i++)
		in >> AL[i];
	for (int i = 0; i < AI[n] - 1; i++)
		in >> AU[i];
	in.close(); // закрываем файл
}

void LUdecomp()
{
	for (int i = 0; i < n; i++)
	{
		int i0 = AI[i] - 1;//первый эл-т строки
		int i1 = AI[i + 1] - 1;//последний эл-т строки
		int j = i - (i1 - i0);//столбец (номер текущий)

		long double sum_d = 0;

		for (int m = i0; m < i1; m++, j++)
		{
			long double sum_l = 0;
			long double sum_u = 0;

			int j0 = AI[j] - 1;//первый эл-т строки
			int j1 = AI[j + 1] - 1;//последний эл-т строки

			int mi = i0;
			int mj = j0;

			int kui = m - i0;//сколько
			int kuj = j1 - j0;//сколько
			int ku = kui - kuj;

			if (ku > 0)
			{
				mi += ku;
			}
			else
			{
				mj -= ku;
			}

			for (; mi < m; mi++, mj++)
			{
				sum_l += AL[mi] * AU[mj];
				sum_u += AU[mi] * AL[mj];
			}
			AL[m] = (AL[m] - sum_l) / DI[j];
			AU[m] = (AU[m] - sum_u);
			sum_d += AL[m] * AU[m];
		}


		DI[i] = DI[i] - sum_d;
	}
}

void Find()
{

	/*Ly=b*/
	for (int i = 0; i < n; i++)
	{
		int i0 = AI[i] - 1;//первый эл-т строки
		int i1 = AI[i + 1] - 2;//последний эл-т строки
		int j = i - (i1 + 1 - i0);//столбец (номер текущий)

		long double sum = 0;
		if (i0 != i1 + 1)
			for (int m = i0; m <= i1; m++, j++)
				sum += AL[m] * V[j];
		V[i] -= sum;
	}
	/*Ux=y*/


	for (int i = n - 1; i >= 0; i--)
	{
		long double sum = 0;
		for (int j = i + 1, k = 0; j < n; j++, k++)
		{
			int j1 = AI[j + 1] - 2;
			int j0 = AI[j] - 1;
			if (j0 <= j1 - k)
				sum += AU[j1 - k] * V[j];
		}
		V[i] = (V[i] - sum) / DI[i];

	}




}

void Output()
{
	cout << "Ёлементы диагонали: ";
	for (int i = 0; i < n; i++)
	{
		cout << DI[i] << " ";
	}
	cout << "\n" << "Ёлементы матрицы L: ";
	for (int i = 0; i < AI[n] - 1; i++)
	{
		cout << AL[i] << " ";
	}
	cout << "\n" << "Ёлементы матрицы U: ";
	for (int i = 0; i < AI[n] - 1; i++)
	{
		cout << AU[i] << " ";
	}

	cout.setf(ios::scientific);
	cout << "\n" << "–ешение системы:" << endl;
	for (int i = 0; i < n; i++)
		cout << setprecision(15) << V[i] << endl;

	cout << "\n" << "ѕогрешность:" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << setprecision(15) << V[i] - i - 1 << endl;
	}
}

void LU_solve(vector <double>& di, vector <double>& gl, vector <double>& gu,
	vector <double>& b, vector <double>& q, vector <int>& ig)
{
	setlocale(LC_ALL, "Russian");
	n = di.size();
	AI = ig;
	for (int i = 0; i < AI.size(); i++)
		AI[i]++;
	DI = di;
	AL = gl;
	AU = gu;
	V = b;
	clock_t tStart = clock();
	LUdecomp();
	Find();
	cout << "  ¬рем€ LU: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;
	q = V;
	//Output();
}
