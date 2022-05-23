#pragma once
#include "Grid.cpp"
#include "Solver.h"
#include "LU.h"
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
        generate_portrait(grid);
        read_bc1();
        G_matrix_assembly(grid);
        M_matrix_assembly(grid, 1);
        M_matrix_assembly(grid, 0);
        b_vector_assembly_cos(grid);
        b_vector_assembly_sin(grid);
        use_bc1(grid);
        LOS_solve(q, di, gl, gu, b, ig, jg);
        //inaccuracy_min(q);
        inaccuracy_max(q);
        regenerate_profile();
        LU_solve(di, gl_p, gu_p, b, q, ig_p);
        inaccuracy_max(q);
}


    void print_vector(ofstream& fout, double t2) {
        fout << "t = " << t2 << endl;
        for (int i = 0; i < n; i++)
            fout << setprecision(19) << std::scientific << q[i] << endl;
        fout << "-----------------------" << endl;
    }

private:
    int ns1;
    struct bc {
        int nx1, nx2;
        int ny1, ny2;
        int ni;
    };
    vector <bc> bc1;
    int n_jgg;            // Размерность векторов gg и jg
    // Матрица A--------------
    std::vector <double> di;
    std::vector <double> gu, gl;//????
    std::vector <int> ig;
    std::vector <int> jg;

    // ПРОФИЛЬНЫЙ ФОРМАТ
    std::vector <double> gu_p, gl_p;
    std::vector <int> ig_p;

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

    void inaccuracy_max(vector <double> q) {
        vector <double> q_right;
        double inaccuracy = 0, sum = 0;
        for (double k = 0; k < 1.0009; k += 0.0068965517241379309)
            for (double l = 0; l < 1.0009; l += 0.0068965517241379309)
            {
                q_right.push_back(l);
                q_right.push_back(k);
            }

        for (int i = 0; i < q.size(); i++)
        {
            inaccuracy += abs(q_right[i] - q[i]);
            sum += q_right[i];
        }

        cout << "Погрешность: " << (inaccuracy) / (sum) << endl;
    }


    void inaccuracy_min(vector <double> q) {
        vector <double> q_right;
        double inaccuracy=0, sum=0;
        for (double k = 0;k <1.009; k += 0.04)
            for (double l = 0; l<1.009; l += 0.04)
            {
                q_right.push_back(l);
                q_right.push_back(k);
            }
            
        for (int i = 0; i < q.size(); i++)
        {
            inaccuracy += abs(q_right[i] - q[i]);
            sum += q_right[i];
        }
        
        cout << "Погрешность: " << (inaccuracy)/(sum) << endl;
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
                for (int ik = 0; ik < 4; ik++)
                    for (int jk = 0; jk < 4; jk++)
                        G[ik][jk] = g1 * G1[ik][jk] + g2 * G2[ik][jk];
                add_local_matrix_p(G, i, j, grid);
            }
        }
    }


    void M_matrix_assembly(Grid& grid, bool flag) {//flag=true -> p
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
                g = hx * hy;
                if (flag) g *= -1 * w * w * chi(wi);
                else g *= w * sigma(wi);
                for (int ik = 0; ik < 4; ik++)
                    for (int jk = 0; jk < 4; jk++)
                        M[ik][jk] = g * C[ik][jk];
                if (flag) add_local_matrix_p(M, i, j, grid);
                else  add_local_matrix_c(M, i, j, grid);
            }
        }
    }

    void b_vector_assembly_sin(Grid& grid) {
        int l;
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
                wi = grid.inSubArea(i, j);
                L[0] = grid.global_num(i, j);
                L[1] = grid.global_num(i + 1, j);
                L[2] = grid.global_num(i, j + 1);
                L[3] = grid.global_num(i + 1, j + 1);
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                f1 = fs(wi, x1, y1);
                f2 = fs(wi, x2, y1);
                f3 = fs(wi, x1, y2);
                f4 = fs(wi, x2, y2);
                g = hx * hy;
                for (int ik = 0; ik < 4; ik++) {
                    int elem = 2 * L[ik];
                    b[elem] += g * (C[ik][0] * f1 + C[ik][1] * f2 + C[ik][2] * f3 + C[ik][3] * f4);
                }
            }
        }
    }


    void b_vector_assembly_cos(Grid& grid) {
        int l1, l2, l3, l4;
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
                wi = grid.inSubArea(i, j);
                L[0] = grid.global_num(i, j);
                L[1] = grid.global_num(i + 1, j);
                L[2] = grid.global_num(i, j + 1);
                L[3] = grid.global_num(i + 1, j + 1);
                x2 = grid.getX(i + 1);
                x1 = grid.getX(i);
                hx = x2 - x1;
                f1 = fc(wi, x1, y1);
                f2 = fc(wi, x2, y1);
                f3 = fc(wi, x1, y2);
                f4 = fc(wi, x2, y2);
                g = hx * hy;
                for (int ik = 0; ik < 4; ik++)
                {
                    int elem = 2 * L[ik]  + 1;
                    b[elem] += g * (C[ik][0] * f1 + C[ik][1] * f2 + C[ik][2] * f3 + C[ik][3] * f4);
                }
            }
        }
    }

    void read_bc1() {
        ifstream bc_file;
        bc_file.open("bc.txt");
        bc_file >> ns1;
        bc1.resize(ns1);
        for (int i = 0; i < ns1; i++) {
            bc_file >> bc1[i].ni >> bc1[i].nx1 >> bc1[i].nx2 >> bc1[i].ny1 >> bc1[i].ny2;
            --bc1[i].nx1;
            --bc1[i].nx2;
            --bc1[i].ny1;
            --bc1[i].ny2;
        }
        bc_file.close();
    }

    void use_bc1(Grid& grid) {
        int p;
        int i_beg, i_end, j_beg, j_end;
        double x1, y1;
        int l;
        for (int s = 0; s < ns1; s++) {
            i_beg = grid.getIXw(bc1[s].nx1);
            i_end = grid.getIXw(bc1[s].nx2);
            j_beg = grid.getIYw(bc1[s].ny1);
            j_end = grid.getIYw(bc1[s].ny2);
            p = bc1[s].ni;
            if (i_beg == i_end) {
                int i = i_beg;
                x1 = grid.getX(i);
                for (int j = j_beg; j <= j_end; j++) {
                    y1 = grid.getY(j);
                    l = grid.global_num(i, j);
                    di[2 * l] = 1;
                    for (int i = ig[2 * l]; i < ig[2 * l + 1]; i++)
                        gl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                        if (jg[i] == 2 * l)
                            gu[i] = 0;

                    di[2 * l + 1] = 1;
                    for (int i = ig[2 * l + 1]; i < ig[2 * l + 2]; i++)
                        gl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                        if (jg[i] == 2 * l + 1)
                            gu[i] = 0;
                    b[2 * l] = us(p, x1, y1);
                    b[2 * l + 1] = uc(p, x1, y1);
                }
            }
            else {
                int j = j_beg;
                y1 = grid.getY(j);
                for (int i = i_beg; i <= i_end; i++) {
                    x1 = grid.getX(i);
                    l = grid.global_num(i, j);

                    di[2 * l] = 1;
                    for (int i = ig[2 * l]; i < ig[2 * l + 1]; i++)
                        gl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                        if (jg[i] == 2 * l)
                            gu[i] = 0;

                    di[2 * l + 1] = 1;
                    for (int i = ig[2 * l + 1]; i < ig[2 * l + 2]; i++)
                        gl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                        if (jg[i] == 2 * l + 1)
                            gu[i] = 0;

                    b[2 * l] = us(p, x1, y1);
                    b[2 * l + 1] = uc(p, x1, y1);
                }
            }
        }
    }


    void generate_portrait(Grid& grid) {

        std::vector <std::vector <int>> list(2 * n);
        int g1, g2;  // Глобальные номера базисных функций
        bool not_in;
        // Цикл по конечным элементам
        for (int j = 0; j < grid.getNy() - 1; j++)
            for (int i = 0; i < grid.getNx() - 1; i++)
            {
                L[0] = 2 * grid.global_num(i, j); L[1] = 2 * grid.global_num(i, j) + 1;
                L[2] = 2 * grid.global_num(i + 1, j); L[3] = 2 * grid.global_num(i + 1, j) + 1;
                L[4] = 2 * grid.global_num(i, j + 1); L[5] = 2 * grid.global_num(i, j + 1) + 1;
                L[6] = 2 * grid.global_num(i + 1, j + 1); L[7] = 2 * grid.global_num(i + 1, j + 1) + 1;
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
        for (int i = 0; i < list.size(); i++)
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

        for (int i = 0; i < k; i++)
        {//построене для p
            int test = 2 * L[i];
            ibeg = ig[2 * L[i]];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 1] - 1;
                while (jg[ibeg] != 2 * L[j]) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < 2 * L[j])
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] += local_matrix[i][j];
                gu[ibeg] += local_matrix[j][i];
                ibeg++;
            }
        }
        for (int i = 0; i < k; i++)
        {
            int test = 2 * L[i] + 1;
            ibeg = ig[2 * L[i] + 1];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 2] - 1;
                while (jg[ibeg] != 2 * L[j] + 1) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < 2 * L[j] + 1)
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] += local_matrix[i][j];
                gu[ibeg] += local_matrix[j][i];
                ibeg++;
            }
        }

    }

    void add_local_matrix_c(std::vector<std::vector <double>>& local_matrix, int ix, int jy, Grid& grid) {
        int sign = 1;
        int ibeg, iend, med;
        int k = 4;
        L[0] = grid.global_num(ix, jy);//глобальные номера узлов на каждом КЭ 
        L[1] = grid.global_num(ix + 1, jy);
        L[2] = grid.global_num(ix, jy + 1);
        L[3] = grid.global_num(ix + 1, jy + 1);

        for (int i = 0; i < k; i++) {
            ibeg = ig[2 * L[i] + 1];
            iend = ig[2 * L[i] + 2] - 1;
            while (jg[ibeg] != 2 * L[i]) {
                med = (ibeg + iend) / 2;
                if (jg[med] < 2 * L[i])
                    ibeg = med + 1;
                else
                    iend = med;
            }
            gl[ibeg] += sign * local_matrix[i][i];
            gu[ibeg] -= sign * local_matrix[i][i];
            //sign *= -1;
        }

        for (int i = 0; i < k; i++) {//построене для c
            ibeg = ig[2 * L[i]];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[2 * L[i] + 1] - 1;
                while (jg[ibeg] != 2 * L[j] + 1) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < 2 * L[j] + 1)
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] -= sign * local_matrix[i][j];
                gu[ibeg] += sign * local_matrix[j][i];
                //sign *= -1;
                ibeg++;
            }
        }
        for (int i = 0; i < k; i++) {
            ibeg = ig[2 * L[i] + 1];
            for (int j = 0; j <= i-1; j++) {
                iend = ig[2 * L[i] + 2] - 1;
                while (jg[ibeg] != 2 * L[j]) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < 2 * L[j])
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                gl[ibeg] += sign*local_matrix[i][j];
                gu[ibeg] -= sign*local_matrix[j][i];
                //sign *= -1;
                ibeg++;
            }
        }
    }
    // Перегенерация в профильный формат
    void regenerate_profile() {
        int zeros;
        ig_p = ig;
        for (int i = 0; i < ig.size() - 1; i++) {
            for (int j = ig[i]; j < ig[i + 1] - 1; j++) {
                zeros = jg[j + 1] - jg[j];
                gl_p.push_back(gl[j]);
                gu_p.push_back(gu[j]);
                for (int z = 0; z < zeros - 1; z++) {
                    gl_p.push_back(0);
                    gu_p.push_back(0);
                }
                if (zeros > 1) {
                    for (int ss = i + 1; ss < ig_p.size(); ss++)
                        ig_p[ss] += (zeros - 1);
                }
            }
            if (i != 0) {
                gl_p.push_back(gl[ig[i + 1] - 1]);
                gu_p.push_back(gu[ig[i + 1] - 1]);
                zeros = i - jg[ig[i + 1] - 1];
                for (int z = 0; z < zeros - 1; z++) {
                    gl_p.push_back(0);
                    gu_p.push_back(0);
                }
                if (zeros > 1) {
                    for (int ss = i + 1; ss < ig_p.size(); ss++)
                        ig_p[ss] += (zeros - 1);
                }
            }
        }
    }
};