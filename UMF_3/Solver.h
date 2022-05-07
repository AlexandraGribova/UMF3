#pragma once
#include <string>
#include <vector>
void r0_s();
void L_1(std::vector <double> pr, std::vector <double>& r);
void U_1(std::vector <double> pr, std::vector <double>& z);
void p0();
double scalar_mult(std::vector <double> v1, std::vector <double> v2, int size);

void X_k(double a, std::vector <double> z);
void R_k(double a, std::vector <double> p);
double Norm(std::vector<double> X);

void AVec(std::vector <double>& x, std::vector <double>& y);
void Z_k(std::vector <double> dat, double b, std::vector <double> d, std::vector <double>& z);
void vec_DI(std::vector <double>vec, std::vector <double>& res);
double _nev();
void LOS_LU();
void LUS_factorisation();
void LOS();
void LOS_solve(std::vector<double>& q, std::vector<double>& _di,
    std::vector<double>& _ggl, std::vector<double>& _ggu, std::vector<double>& b, const std::vector<int>& _ig, const std::vector<int>& _jg);