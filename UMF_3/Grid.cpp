#pragma once
#include <fstream>
#include<iostream>
#include <vector>
#include <iomanip>
using namespace std;


struct sub_area {
    int nx1; 
    int ny1;
    int nx2;  
    int ny2;
    int ni; 
};


class Area
{
protected: 
    int nw;  
    vector <double> Xw, Yw;  //координаты подобластей
    vector <sub_area> Mw;  
    void read_area() {  
        ifstream area_file;
        area_file.open("Area.txt");
        area_file >> nw;
        Xw.resize(nw + 1);
        Yw.resize(nw + 1);
        Mw.resize(nw);
        for (int i = 0; i < nw + 1; i++)
            area_file >> Xw[i]>> Yw[i];
        for (int i = 0; i < nw; i++) {
            area_file >> Mw[i].ni >> Mw[i].nx1 >> Mw[i].ny1 >> Mw[i].nx2 >> Mw[i].ny2;
            Mw[i].nx1 -= 1;
            Mw[i].nx2 -= 1;
            Mw[i].ny1 -= 1;
            Mw[i].ny2 -= 1;
        }
        area_file.close();
    }
};
class Grid : public Area
{
private:
    vector <double> X, Y;  
    vector <int> IXw, IYw;  
    vector <int> nx, ny;  
    vector <double> qx, qy;  
    int Nx, Ny; 

    void read_grid() {
        ifstream GridX, GridY;
        GridX.open("GridX.txt");
        GridY.open("GridY.txt");
        nx.resize(nw);
        ny.resize(nw);
        qx.resize(nw);
        qy.resize(nw);
        for (int i = 0; i < nw; i++) {
            GridX >> nx[i] >> qx[i];
        }
        for (int i = 0; i < nw; i++) {
            GridY >> ny[i] >> qy[i];
        }
        GridX.close();
        GridY.close();
        
    }


    void makeGrid(vector <double>& XY, double left, double right, int n, double qxy, int& i0, vector <int>& IXYw) {
        double h0;
        static int j = 0;
        if (qxy - 1 < 1E-16)
            h0 = (right - left) / n;
        else
            h0 = (right - left) * (1 - qxy) / (1 - pow(qxy, n));

        XY[i0] = left;
        IXYw[j] = i0; j++;
        for (int i = i0 + 1; i < n + i0; i++) {
            XY[i] = XY[i - 1] + h0;
            h0 *= qxy;
        }
        i0 = n + i0;
    };


public:
    Grid() {   
        read_area();
        read_grid();
        Nx = 0;
        Ny = 0;
        for (int i = 0; i < nw; i++)
        {
            Nx += nx[i];
            Ny += ny[i];
        }
        Nx++; Ny++;
        X.resize(Nx);
        Y.resize(Ny);
        IXw.resize(nw + 1);
        IYw.resize(nw + 1);
        int ix0 = 0, iy0=0;

        for (int i = 0; i < nw; i++)
        {
            makeGrid(X, Xw[i], Xw[i + 1], nx[i], qx[i], ix0, IXw);
            makeGrid(Y, Yw[i], Yw[i + 1], ny[i], qy[i], iy0, IYw);
        }
        X[Nx - 1] = Xw[nw];
        Y[Ny - 1] = Yw[nw];
        IXw[nw] = ix0;
        IYw[nw] = iy0;
    }

   
    int inSubArea(int px, int py) {
        int ixw1, ixw2, iyw1, iyw2;
        for (int i = 0; i < nw; i++) {
            ixw1 = IXw[Mw[i].nx1];
            ixw2 = IXw[Mw[i].nx2];
            iyw1 = IYw[Mw[i].ny1];
            iyw2 = IYw[Mw[i].ny2];
            bool flag1, flag2;
            flag1 = px >= ixw1 && px <= ixw2 && px + 1 >= ixw1 && px + 1 <= ixw2;
            flag2 = py >= iyw1 && py <= iyw2 && py + 1 >= iyw1 && py + 1 <= iyw2;
            if (flag1 && flag2)
                return Mw[i].ni;
        }
    }

    double chi(int wi) {//русская буква х
        return 1;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi) {
        return 1;
    }

    double fs(int wi, double x, double y) {
        return -1;
    }

    double fc(int wi, double x, double y) {
        return -1;
    }


    double getX(int i) { return X[i]; }
    double getY(int i) { return Y[i]; }
    int getNx() { return Nx; }
    int getNy() { return Ny; }
};