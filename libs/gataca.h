///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//cd Documents/Doctorado/Resultados/libs
//g++ -Wall -c *.cc
#ifndef GATACA_H
#define GATACA_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "horse.h"
#include "linear.h"
#include "finder.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ios;
using std::string;

class GATACA
{
public:
    Horse hor;
    Linear lin;
    Finder find;
    GATACA();
    GATACA(int semilla);
    void condin(int *vector, int cols);
    void condin(int *vector, int cols, int especial);
    void escalon(int *vector, int cols);
    void walk(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int tam1, int tam2, int tolerancia, double intnw, int *vectorentrada, int num, int& n, int pfenotipo, int &pasos, int &stanby, int m);
    void cambiafenotipo(double **Matriz1, int *Vector1, int **Matriz2, int tam1, int tam2, int tam3);
    void inducirduplicacion(double** &Matriz1, double** &Matriz2, double **Matriz3, double **Matriz4, int x, int tam1, int tam2);
    void control1(double** &Matriz1, double** &Matriz2, double **Matriz3, double **Matriz4, int tam1, int tam2, double intnw);
    void control2(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2);
    void exploracion(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek);//Muta red de muestra, saca distancia fenotípica y la guarda en una matriz
    void exploraciondup(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int **Matriz5, double *Vector1, double *Vector2, double& alldists, double& notoldists_noori, int& notolrm_noori, int tam1, int tam2, int tolerancia, double intnw, int n, int pd, int p, int num1, int num2, int &cuenta1, int &cuenta2, int &cuenta3, int &Fnuevo, int &Fnuevoin, int tipo, int *valordek, int m);
    void mediadist(double **Matriz1, double *Vector1, int num1, int num2); //Toma la matriz de distancias, saca promedio por muestra y lo guarda en un vector 1por muestra
    void ceros(double **Matriz1, double *Vector1, int hils, int cols);//Toma la matriz de distancias, saca la proporción de ceros por muestra y la guarda en un vector por muestra
    void robustezfenotipica(double *Vector1, double *Vector2, int tam);//1 menos la distancia media
    double Shannon(int *valordek, int tam);
    double Shannon2(int *valordek, int tam);
    void Shannon3(int *valordek, int tam);
    double rmutarriba(double **Matriz1, int hils1, int cols1, int **Matriz2, int hils2, int p, double tolerancia, double intnw, double& PromDist);
    double rmutabajo(double **Matriz1, int hils1, int cols1, int **Matriz2, int hils2, int **Matriz3, int hils3, double tolerancia, double intnw, double& PromDist);
    void exploracionduppuntofijo(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int **Matriz5, double *Vector1, double *Vector2, double& alldists, int tam1, int tam2, int tolerancia, double intnw, int n, int pd, int p, int num1, int num2, int &cuenta1, int &cuenta2, int &cuenta3, int &Fnuevo, int &Fnuevoin, int tipo, int *valordek, int &puntosfijos);
    void exploracionpuntofijo(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek, int &puntosfijos);
    void exploracion(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek, int **Fn_notol, int per_notol, int &cuenta_mut_dup_notol);
    void interaction_lost(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *vectorentrada);
    void interaction_gain(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla);
    void interaction_change(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla);
    void interaction_lost2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *vectorentrada, int semilla);
    void interaction_gain2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla);
    void interaction_change2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla);
    void interaction_lost_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_lost_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial);
    void interaction_gain_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla);
    void interaction_change_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla);
    void interaction_lost_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_lost_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_lost_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    
    void interaction_lost_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    
    void interaction_lost_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial);
    void interaction_gain_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla);
    void interaction_change_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla);
    void interaction_lost_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_lost_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial);
    void interaction_gain_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);
    void interaction_change_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla);

    
private:
    void back(double **Matriz, int x, int y, double valor);
    void choose(double **Matriz1, int hils1, int cols1, double **Matriz2, int hils2, int cols2, int& matmut);
    void mut(double **matriz, int hils, int cols, int tolerancia, double intnw, int& x, int& y, double& valor);
};

#endif
