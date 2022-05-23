#include <iostream>
#include <fstream>
#include <vector>
#include<locale.h>
#include<iomanip>

using namespace std;
vector<double> DI, AL, AU, V, X;//размерность, vec, di,ai,L,U - в таком порядке
vector<int>AI;
int n;
void input()
{
	string line;
	ifstream in("input.txt"); // окрываем файл для чтения
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

void Gauss() {
	ifstream GaussInput("GaussInput.txt");
	ofstream GaussOutput("GaussOutput.txt");
	GaussInput >> n;
	long double** M;
	long double* F;
	M = new long double* [n];
	for (int i = 0; i < n; i++)
		M[i] = new long double[n];
	F = new long double[n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			GaussInput >> M[i][j];
	for (int i = 0; i < n; i++)
		GaussInput >> F[i];
	for (int i = 1; i < n; i++)
		for (int j = i; j < n; j++) {
			long double m = M[j][i - 1] / M[i - 1][i - 1];
			for (int k = 0; k < n; k++)
				M[j][k] = M[j][k] - m * M[i - 1][k];
			F[j] = F[j] - m * F[i - 1];
		}
	for (int k = n - 1; k >= 0; k--) {
		double buf = 0;
		for (int j = k + 1; j < n; j++)
			buf += M[k][j] * F[j];
		F[k] = F[k] - buf;
		F[k] = F[k] / M[k][k];
	}
	GaussOutput.setf(ios::scientific);
	GaussOutput << "\n" << "Решение системы:" << endl;
	for (int i = 0; i < n; i++)
		GaussOutput << setprecision(15) << V[i] << endl;

	GaussOutput << "\n" << "Погрешность:" << endl;
	for (int i = 0; i < n; i++)
	{
		GaussOutput << setprecision(15) << V[i] - i - 1 << endl;
	}

	delete[]F;
	for (int i = 0; i < n; i++)
		delete[]M[i];
	delete[]M;
	GaussInput.close();
	GaussOutput.close();

}

void Gilbert()
{

	V.resize(n);
	DI.resize(n);
	AI.resize(n + 1);
	AL.resize((n * n - n) / 2);
	AU.resize((n * n - n) / 2);
	X.resize(n);
	int k = 0;
	for (int i = 1; i < (n + 1), k < (n * n - n) / 2; i++)
		for (int j = 1; j < i; j++, k++)
			AL[k] = 1.0 / (i + j - 1);
	k = 0;
	for (int j = 1; j < n + 1, k < (n * n - n) / 2; j++)
		for (int i = 1; i < j; i++, k++)
			AU[k] = 1.0 / (i + j - 1);

	for (int j = 1; j < n + 1; j++)
		DI[j - 1] = 1.0 / (2.0 * j - 1);

	AI[0] = AI[1] = 1;
	for (int j = 1; j < n; j++)
		AI[j + 1] = AI[j] + j;

	//умножение
	for (int i = 0; i < n; i++) {
		for (int k = 1; k <= n; k++)
			X[k - 1] = k;
		V[i] = DI[i] * X[i];
		int i0 = AI[i] - 1;
		int i1 = AI[i + 1] - 1;
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++, j++)
		{
			V[i] += AL[k] * X[j];
			V[j] += AU[k] * X[i];
		}
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
	cout << "Элементы диагонали: ";
	for (int i = 0; i < n; i++)
	{
		cout << DI[i] << " ";
	}
	cout << "\n" << "Элементы матрицы L: ";
	for (int i = 0; i < AI[n] - 1; i++)
	{
		cout << AL[i] << " ";
	}
	cout << "\n" << "Элементы матрицы U: ";
	for (int i = 0; i < AI[n] - 1; i++)
	{
		cout << AU[i] << " ";
	}

	cout.setf(ios::scientific);
	cout << "\n" << "Решение системы:" << endl;
	for (int i = 0; i < n; i++)
		cout << setprecision(15) << V[i] << endl;

	cout << "\n" << "Погрешность:" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << setprecision(15) << V[i] - i - 1 << endl;
	}
}

void LU_solve(vector <double> &di, vector <double>& gl, vector <double>& gu,
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
	LUdecomp();
	Find();
	q = V;
	//Output();
}
