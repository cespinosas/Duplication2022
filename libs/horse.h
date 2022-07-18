///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HORSE_H
#define HORSE_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ifstream;
using std::string;
using std::ios;

class Horse 
{
public:
    gsl_rng * r;
    Horse();
    Horse(int semilla);
    void start_rng(int semilla);
    double randreal();
    double randgauss(int sigma);
    void close_rng();
    int find_min(int *vec, int tam); //Encuentra la distancia mínima posible
    double find_min(double *vec, int tam);
    void fillv0(int vec[], int s); //Llena con 0's el vector para la distancia
    void fillv0(double vec[], int s);
    void fillv1(int vec[], int s);//Llena con 1's el vector
    void fillv0(bool vec[], int s);
    void fillm0(double **Matriz, int tam1, int tam2);
    void fillm0(int **Matriz, int tam1, int tam2);
    int randflat(int min, int max);
    double randuniform(int min, int max);
    void spearman(double *Vector1, double *Vector2, int num, double &rho, double &pval);
    void pearson(double *Vector1, double *Vector2, int num, double &rho, double &pval);
    double coin();
    void printmat(int **matriz, int hils, int cols);
    void printmat(double **matriz, int hils, int cols);
    void printvec(int *vector, int cols);
    void printvec(double *vector, int cols);
    void espacio(double** &Matriz, int hils, int cols);
    void espacio(int **&Matriz, int hils, int cols);
    void espacio(double *&Vector, int cols);
    void espacio(int *&Vector, int cols);
    void espacio(bool *&Vector, int cols);
//     void espacio(double *&Vector, string cols);
//     void espacio(string *&Vector, int cols);
    void borrar(double** &Matriz, int hils);
    void borrar(int** &Matriz, int hils);
    void borrar(double* &Vector);
    void borrar(int* &Vector);
    void borrar(bool* &Vector);
//     void borrar(string* &Vector);
    void rellenar(double **Matriz1, double **Matriz2, double **Matriz3, int hils1, int cols1, int hils2, int cols2);
    int round(double num1, int num2);
    void rellenar(double **Matriz1, double **Matriz2, int hils, int cols);
    string inttostring(int num);
    void open_ifstream(ifstream& fe, string nomb);
    void open_ofstream(ofstream& fs, string nomb);
    string doubletostring_tex(double num);
    int char_in_string(char c, string cue, int from, int until);
    int char_in_string(char c, string cue);
    void sumarmatriz(int **Matriz, int *Vector, int hils, int cols);
    void media(double **Matriz1, double *Vector1, int hils, int cols);
    double densidadparcial(double **Matriz1, int tam1, int tam2);
    double densidadtotal(double **Matriz1, int tam1, int tam2, double **Matriz2, int tam3, int tam4);
    void mediacols(double **Matriz1, double *Vector1, int hils, int cols);
    double mediacols(double **Matriz1, int hils);
    double mediacols(double *Vector1, int hils);
    void mediacols(double *Vector1, double& media, int hils);
    void desv(double *Vector1, double *Vector2, double *Vector3, int num1, int num2);
    void desv(double *Vector1, double& media, double& desv, int num1);
    double desv(double *Vector1, int num1);
    void nonullscuenta(double *Vector1, double *Vector2, int tam, int& cuenta);
    void nonullsvecs(double *Vector1, double *Vector2, int tam, double *Vector3, double *Vector4);
    void sort(int* vec, int* nvec, int tam);
    int find_max(int* vec, int tam);
    double find_max(double* vec, int tam);
    int count_ordinals(int *Vec, int hils1);
    void ordinals_value(int *Vec, int hils1, double **Mat, int hils2);
    double get_median(int* vec, int tam);
    double get_median(double* vec, int tam);
    //void sort(int* vec, int* nvec, int tam);
    void sort(double* vec, double* nvec, int tam);
    double dist_fenotipica(int *ancestro, int *fenotipo, int tam);
    double dist_fenotipica(int *ancestro, int **fenotipo, int hilsf, int tam);
    double dist_fenotipica(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam);
    double dist_aux2(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam);
    double dist_aux2(int *vector1, int cols, int *vector2);
    void pantalla(int **M, int hils, int cols);
    void pantalla(double **M, int hils, int cols);
    void pantalla(int *V, int hils);
    void pantalla(double *V, int hils);
    double sh_index(int *valordek, int tam);
    void reduce_fenotipos(int **Fenotipos_para_reducir, int columnas, int *per_para_reducir, int cuenta_p, int *per_reducidos);
    void ofstream_valor(string root, string carpeta, string archivo, double valor);
    void ofstream_valor(string root, string archivo, string name, double valor1, double valor2, double valor3);
    void ofstream_valor(string root, string archivo, string name, double valor1, double valor2, double valor3, double valor4, double valor5, double valor6, double valor7, double valor8, double valor9);
};

#endif
