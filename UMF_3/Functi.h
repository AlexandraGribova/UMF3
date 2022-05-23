#pragma once
class Functi_1
{
public:

    double chi(int wi) {
        return 10E-11;
    }
    double w = 10E-4;

    double sigma(int wi) {
        return 10E+4;
    }

    double lambda(int wi) {
        return 10E+4;
    }

    double fs(int wi, double x, double y) {
        return sigma(wi)*-y * w-chi(wi)*x*w*w;
    }

    double fc(int wi, double x, double y) {
        return sigma(wi)*x * w- chi(wi) * y*w*w;
    }

    double us(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return x;
        case 2:
            return 1;
        case 3:
            return x;
        case 4:
            return 0;
        }
    }

    double uc(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return 0;
        case 2:
            return y;
        case 3:
            return 1;
        case 4:
            return y;
        }
    }
};