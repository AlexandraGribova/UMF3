﻿#pragma once
#include "Grid.cpp"
#include <iostream>
#include <iomanip>
class FEM
{
public:
        FEM(Grid& grid) {
        n = grid.getNx() * grid.getNy();
        q.resize(2*n);
       
        b.resize(2*n);
        L.resize(2*n - 1);
        D.resize(2*n);
        U.resize(2*n - 1);
    }

   
    void print_vector(ofstream& fout, double t2) {
        fout << "t = " << t2 << endl;
        for (int i = 0; i < n; i++)
            fout << setprecision(19) << std::scientific << q[i] << endl;
        fout << "-----------------------" << endl;
    }

private:
    vector <vector <double>> G = { {0, 0, 0, 0},
                                   {0, 0, 0, 0},
                                   {0, 0, 0, 0} ,
                                   {0, 0, 0, 0} };
    vector <vector<double>> G1 = { {2 / 6., -2 / 6., 1 / 6., -1 / 6.},
                                   {-2 / 6., 2 / 6., -1 / 6., 1 / 6.},
                                   {1 / 6., -1 / 6., 2 / 6., -2 / 6.},
                                   {-1 / 6., 1 / 6., -2 / 6., 2 / 6.} };

    vector <vector<double>> G2 = { {2 / 6., 1 / 6., -2 / 6., -1 / 6.},
                                   {1 / 6., 2 / 6., -1 / 6., -2 / 6.},
                                   {-2 / 6., -1 / 6., 2 / 6., 1 / 6.},
                                   {-1 / 6., -2 / 6., 1 / 6., 2 / 6.} };
   
    vector <vector <double>> M = { {0, 0, 0, 0},
                                   {0, 0, 0, 0},
                                   {0, 0, 0, 0} ,
                                   {0, 0, 0, 0} };

    vector <vector<double>> C = { {4 / 36., 2 / 36., 2 / 36., 1 / 36.},
                                   {2 / 36., 4 / 36., 1 / 36., 2 / 36.},
                                   {2 / 36., 1 / 36., 4 / 36., 2 / 36.},
                                   {1 / 36., 2 / 36., 2 / 36., 4 / 36.} };

    int n;  
    double eps = 1E-13;
    vector <double> q;
    vector <double> b; 
    vector <double> L;
    vector <double> U;
    vector <double> D;



    void add_local_matrix(vector<vector <double>>& a, int k) {
        
    }


    // Óìíîæåíèå ãëîáàëüíîé ìàòðèöû{l, d, r} íà âåêòîð q, ðåçóëüòàò â product
    void multMatrVect() {
       
    }


    void G_matrix_assembly(Grid& grid, double t2) {
        double lambd;
        int nx=grid.getNx(), ny= grid.getNy();
        double hx, hy;
        double x1, x2, y1, y2;
        double g1, g2;
        int wi;
        for (int j = 0; j < ny - 1; j++)
        {
            y2 = grid.getY(j + 1);
            y1 = grid.getY(j);
            hy = y2 - y1;        
            for (int i = 0; i < nx - 1; i++) 
            {
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                wi = grid.inSubArea(i, j);
                lambd = grid.lambda(wi);
                g1 = lambd*hy / hx;
                g2 = lambd * hx / hy;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        G[i][j] = g1 * G1[i][j]+ g2 * G2[i][j];
                add_local_matrix(G, i);////
            }
        }

    }


    void M_matrix_assembly(Grid& grid, double m) {        
        int nx = grid.getNx(), ny = grid.getNy();
        double hx, hy;
        double x1, x2, y1, y2;
        double g;
        int wi;  
        for (int j = 0; j < ny - 1; j++)
        {
            y2 = grid.getY(j + 1);
            y1 = grid.getY(j);
            hy = y2 - y1;
            for (int i = 0; i < nx - 1; i++) {
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                wi = grid.inSubArea(i, j);
                g = hx*hy *m;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        M[i][j] = g * C[i][j];
                add_local_matrix(M, i);////
            }
        }
    }


    void b_vector_assembly(Grid& grid) {
        int nx = grid.getNx(), ny = grid.getNy();
        double hx, hy;
        double x1, x2, y1, y2;
        double f1, f2, f3, f4;
        double g;
        int wi;
        for (int j = 0; j < ny - 1; j++)
        {
            y2 = grid.getY(j + 1);
            y1 = grid.getY(j);
            hy = y2 - y1;
            for (int i = 0; i < nx - 1; i++) {
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                f1 = grid.f(wi, x1, y1);
                f2 = grid.f(wi, x2, y1);
                f3 = grid.f(wi, x1, y2);
                f4 = grid.f(wi, x2, y2);
                wi = grid.inSubArea(i, j);
                g = hx * hy;
                for (int i = 0; i < 4; i++)
                    b[i] += g * (C[i][0] * f1 + C[i][1] * f2 + C[i][2] * f3 + C[i][3] * f4);
            }
        }
    }

   

    void boundary_condition(Grid& grid, double t2) {
        double u1_ = u1(t2);
        double u2_ = u2(t2);
        r[0] = 0;
        d[0] = 1;
        b[0] = u1_;

        l[n - 2] = 0;
        d[n - 1] = 1;
        b[n - 1] = u2_;
    }


    double vector_norm(vector <double>& v) {
        double norm = 0;
        for (int i = 0; i < n; i++)
            norm += v[i] * v[i];
        return sqrt(norm);
    }

    // ÆÓÒÊÎ ÍÅÝÔÔÅÊÒÈÂÍÎ
    void update_matrix(Grid& grid, double t1, double t2) {
        for (int i = 0; i < n - 1; i++) { // êîñòûëü
            d[i] = 0;
            l[i] = 0;
            r[i] = 0;
        }
        d[n - 1] = 0;
        M_matrix_assembly(grid, t1, t2); // êîñòûëü
        G_matrix_assembly(grid, t2);
        boundary_condition(grid, t2);  // êîñòûëü
    }
};