///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FINDER_H
#define FINDER_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include "horse.h"
#include "linear.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ios;

class Finder 
{
public:
    Finder();
    Linear lin;
    void atractor(double **matriz1, int **trayectoria, int *vectorentrada, int tam, int& n, int& p);
    void atractor2(double **matriz1, int **trayectoria, int *vectorentrada, int tam, int& n, int& p, int& pasosmax);
//     void atractor1(double **matriz1, int **trayectoria, int *vectorentrada, int tam, int& n, int& p);
    double distfenotipica(int *ancestro, int *fenotipo, int tam); //Si el ancestro y el fenotipo tienen una hilera
    double distfenotipica(int *ancestro, int **fenotipo, int hilsf, int tam); //Si el ancestro tiene una hilera pero el fenotipo es una matriz
    double distfenotipica(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam); //Si el ancestro y el fenotipo son matrices
//     double distfenotipica2(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam);
    int numintnw(double **matriz, int hils, int cols);
    int negativos(double **matriz, int hils, int cols);
    int positivos(double **matriz, int hils, int cols);
    void media(double *Vector1, double *Vector2, int num1, int num2);
    void media(double *Vector, double **Matriz, int num1, int num2, int num3);
    void desv(double *Vector1, double *Vector2, double *Vector3, int num1, int num2);
    void desv(double *Vector1, double **Matriz1, double **Matriz2, int num1, int num2, int num3);
    double compara1(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2);
    double compara2(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2);
    void minmaxmed(double **Matriz1, double **Matriz2, int hils, int cols);
    void minmaxmedin(double **Matriz1, double **Matriz2, int hils, int cols);
    void dist_btw_esp(int per_antes, int *Per_antes, int Ns, int **Fn_antes, string arch);
    void dist_ant_vs_desp(int Nf, int Ns, int cuenta1, int **acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes);
    void dist_btw_esp2(int per_antes, int *Per_antes, int Ns, int **Fn_antes, string arch1, string arch2, int gen);
    void dist_ant_vs_desp2(int Nf, int Ns, int cuenta1, int **acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes);
//     void dist_ant_vs_desp2(int Nf, int Ns, int cuenta1, int **acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes);
    void dist_ant_vs_desp_tol(int Nf, int Ns, int cuenta1, int acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes, int gen);
private:
    Horse hor;
    double dist_aux(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam);
    double dist_aux(int *vector1, int cols, int *vector2);//distancia entre vectores, la uso para fenotipo vs deseado
    int atractores(int *vector1, int *vector2, int **matriz, int tam, int n);
    };

#endif
