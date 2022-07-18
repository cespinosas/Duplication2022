///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LINEAR_H
#define LINEAR_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include "horse.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ios;

class Linear 
{
public:
    Horse hor;
    Linear();
    Linear(int semilla);
    void make(double **matriz, int hils, int cols, double intnw);
    void dot(double **matriz, int **matriz2, int *vector1, int *vector2, int tam, int n);
    void dot(double **matriz1, int **matriz2, int **matriz3, int tam1, int tam2, int n, int p);
    void dot2(double **matriz1, int **matriz2, int **matriz3, int tam1, int tam2, int n, int p, int m, int veces);
private:
    };

#endif
