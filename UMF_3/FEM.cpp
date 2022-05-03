#pragma once
#include "Grid.cpp"
#include "Functi.h"
#include <iostream>
#include <iomanip>
#include<algorithm>
class FEM :Functi_1
{
public:
    FEM(Grid& grid) {
        n = grid.getNx() * grid.getNy();
        q.resize(2 * n);
        b.resize(2 * n);
        L.resize(8);
        /* L.resize(2*n - 1);
         D.resize(2*n);
         U.resize(2*n - 1);*/
        generate_portrait(grid);
        G_matrix_assembly(grid);
        M_matrix_assembly(grid, 1);
    }


    void print_vector(ofstream& fout, double t2) {
        fout << "t = " << t2 << endl;
        for (int i = 0; i < n; i++)
            fout << setprecision(19) << std::scientific << q[i] << endl;
        fout << "-----------------------" << endl;
    }

private:
    int n_jgg;            // Размерность векторов gg и jg
    // Матрица A--------------
    std::vector <double> di;
    std::vector <double> gu, gl;//????
    std::vector <int> ig;
    std::vector <int> jg;
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
    vector <int> L;
    vector <double> U;
    vector <double> D;


    // Óìíîæåíèå ãëîáàëüíîé ìàòðèöû{l, d, r} íà âåêòîð q, ðåçóëüòàò â product
    void multMatrVect() {

    }


    void G_matrix_assembly(Grid& grid) {
        double lambd;
        int nx = grid.getNx(), ny = grid.getNy();
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
                lambd = lambda(wi);
                g1 = lambd * hy / hx;
                g2 = lambd * hx / hy;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        G[i][j] = g1 * G1[i][j] + g2 * G2[i][j];
                add_local_matrix_p(G, i, j, grid);////
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
                g = hx * hy * m;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        M[i][j] = g * C[i][j];
                add_local_matrix_c(M, i, j, grid);////
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
        int l;
        for (int j = 0; j < ny - 1; j++)
        {
            y2 = grid.getY(j + 1);
            y1 = grid.getY(j);
            hy = y2 - y1;
            for (int i = 0; i < nx - 1; i++) {
                l = j * nx + i;
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                wi = grid.inSubArea(i, j);
                if (!(l % 2)) {  // если чётное
                    f1 = fc(wi, x1, y1);
                    f2 = fc(wi, x2, y1);
                    f3 = fc(wi, x1, y2);
                    f4 = fc(wi, x2, y2);
                }
                else {
                    f1 = fs(wi, x1, y1);
                    f2 = fs(wi, x2, y1);
                    f3 = fs(wi, x1, y2);
                    f4 = fs(wi, x2, y2);
                }
                g = hx * hy;
                for (int i = 0; i < 4; i++)
                    b[l] += g * (C[i][0] * f1 + C[i][1] * f2 + C[i][2] * f3 + C[i][3] * f4);
            }
        }
    }



    void boundary_condition(Grid& grid, double t2) {

    }


    double vector_norm(vector <double>& v) {
        double norm = 0;
        for (int i = 0; i < n; i++)
            norm += v[i] * v[i];
        return sqrt(norm);
    }


    void generate_portrait(Grid& grid) {

        std::vector <std::vector <int>> list(2*n);
        int g1, g2;  // Глобальные номера базисных функций
        bool not_in;
        // Цикл по конечным элементам
        for (int j = 0; j < grid.getNy() - 1; j++)
            for (int i = 0; i < grid.getNx() - 1; i++) 
            {
                L[0] = 2*grid.global_num(i, j); L[1] = 2 * grid.global_num(i, j)+1;
                L[2] = 2*grid.global_num(i + 1, j); L[3] = 2 * grid.global_num(i + 1, j)+1;
                L[4] = 2*grid.global_num(i, j + 1); L[5] = 2 * grid.global_num(i, j + 1)+1;
                L[6] = 2*grid.global_num(i + 1, j + 1); L[7] = 2 * grid.global_num(i + 1, j + 1)+1;
                // Цикл по ненулевым базисным функциям
                for (int in = 0; in < 8; in++) {
                    g1 = L[in];
                    for (int jn = in + 1; jn < 8; jn++) {
                        // g2 > g1
                        g2 = L[jn];
                        // Перед добавлением проверяем наличие элемента в списке
                        not_in = true;
                        for (int l = 0; l < list[g2].size() && not_in; l++)
                            if (g1 == list[g2][l])
                                not_in = false;

                        // Добавляем
                        if (not_in)
                        {
                            list[g2].push_back(g1);
                        }
                    }
                }
            }


        // Сортировка списков по возрастанию
        for (int i = 0; i < n; i++)
            sort(list[i].begin(), list[i].end());

        // Формирование вектора ig
        ig.resize(2 * n + 1);
        for (int i = 0; i < list.size(); i++)
            ig[i + 1] = ig[i] + list[i].size();


        n_jgg = ig[2 * n];
        jg.resize(n_jgg);
        gl.resize(n_jgg);
        gu.resize(n_jgg);
        di.resize(2 * n);

        // Формирование вектора jg
        for (int i = 1, j = 0; i < 2 * n; i++)
            for (int k = 0; k < list[i].size(); k++, j++)
                jg[j] = list[i][k];
    }

    //

    void add_local_matrix_p(std::vector<std::vector <double>>& local_matrix, int ix, int jy, Grid& grid) {
        int ibeg, iend, med;
        int k = 4;
        L[0] = grid.global_num(ix, jy);//глобальные номера узлов на каждом КЭ 
        L[1] = grid.global_num(ix + 1, jy);
        L[2] = grid.global_num(ix, jy + 1);
        L[3] = grid.global_num(ix + 1, jy + 1);

        for (int i = 0; i < k; i++)
        {

            di[2 * L[i]] += local_matrix[i][i];
            di[2 * L[i] + 1] += local_matrix[i][i];
        }

        for (int i = 0; i < k; i++) {//построене для p
            ibeg = ig[2 * L[i]];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 1] - 1;
                while (jg[ibeg] != L[j]) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < L[j])
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] += local_matrix[i][j];
                gu[ibeg] += local_matrix[j][i];
                ibeg++;
            }

            ibeg = ig[2 * L[i] + 1];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 2] - 1;
                while (jg[ibeg] != L[j] + 1) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < L[j] + 1)
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] = local_matrix[i][j];
                gu[ibeg] = local_matrix[j][i];
                ibeg++;
            }
        }

    }

    void add_local_matrix_c(std::vector<std::vector <double>>& local_matrix, int ix, int jy, Grid& grid) {
        int ibeg, iend, med;
        int k = 4;
        L[0] = grid.global_num(ix, jy);//глобальные номера узлов на каждом КЭ 
        L[1] = grid.global_num(ix + 1, jy);
        L[2] = grid.global_num(ix, jy + 1);
        L[3] = grid.global_num(ix + 1, jy + 1);

        for (int i = 0; i < k; i++) {//построене для p
            ibeg = ig[2 * L[i] + 1];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 2] - 1;
                while (jg[ibeg] != L[j] + 1) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < L[j] + 1)
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] += local_matrix[i][j];
                gu[ibeg] += local_matrix[j][i];
                ibeg++;
            }

            ibeg = ig[2 * L[i] + 2];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 3] - 1;
                while (jg[ibeg] != L[j] + 2) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < L[j] + 2)
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] = local_matrix[i][j];
                gu[ibeg] = local_matrix[j][i];
                ibeg++;
            }
        }

    }
};