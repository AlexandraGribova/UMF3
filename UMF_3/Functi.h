#pragma once
class Functi_1
{
public:

    double chi(int wi) {//русская буква х
        return 0;
    }
    double w = 1;

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi) {
        return 1;
    }

    double fs(int wi, double x, double y) {
        return - y*w;
    }

    double fc(int wi, double x, double y) {
        return x*w;
    }

    double us(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return x;
        case 2:
            return 3;
        case 3:
            return x;
        case 4:
            return 1;
        }
    }

    double uc(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return 1;
        case 2:
            return y;
        case 3:
            return 3;
        case 4:
            return y;
        }
    }
};

