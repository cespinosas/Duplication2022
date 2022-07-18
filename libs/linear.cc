///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include "linear.h"
#include "horse.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ios;

Linear::Linear()
{
}

Linear::Linear(int semilla)
{
    hor.start_rng(semilla);
}

void Linear::make(double **matriz, int hils, int cols, double intnw){
    
    //Llenamos Matriz    myfile1.open ("matrices.txt", ios::app);

    
    int celdas = hils * cols;
    double novac = intnw * celdas;
    int vac = celdas - (int) novac;
    int sumnovac = 0;
    int sumvac = 0;
    int sumceros;
    int noceros = 1;
    
    do{ //010319
        hor.fillm0(matriz, hils, cols);
        for (int i = 0; i < hils; i++){
            for (int j = 0; j < cols; j++){
                //Si el volado dice que celda diferente de cero.
                if (hor.randreal() < intnw){
                    //Nos aseguramos que falten celdas diferentes de cero.
                    if(sumnovac != (int) novac){
                        //Si no hay suficientes, se llena con gauss.
                        matriz[i][j] = hor.randgauss(1);
                        sumnovac++;
                    }
                    else {
                        //Si ya hay suficientes, la celda será cero.
//                         matriz[i][j] = 0;
                        sumvac++;
                    }
                }
                //Si el volado dice que la celda es cero.
                else {
                    //Nos asegurames que faten celdas diferentes de cero.
                    if(sumvac != vac){
                        //Si no hay suficientes, se llena con cero.
//                         matriz[i][j] = 0;
                        sumvac++;
                    }
                    //Si ya hay suficientes ceros, será diferente de cero.
                    else{
                        matriz[i][j] = hor.randgauss(1);
                        sumnovac++;
                    }
                }
            }
        }
        //010319 Para evitar hileras de ceros
        for (int i = 0; i < hils; i++){
            sumceros = 0;
            for (int j = 0; j < cols; j++){
                sumceros = sumceros + matriz[i][j];
            }
            if(sumceros == 0){
                noceros = 0;
                break;
            }
        }
        //010319 Para evitar columnas de ceros
        if(noceros != 0){
            for (int i = 0; i < cols; i++){
                sumceros = 0;
                for (int j = 0; j < hils; j++){
                    sumceros = sumceros + matriz[j][i];
                }
                if(sumceros == 0){
                    noceros = 0;
                    break;
                }
            }
        }
    }while(noceros == 0);
}

void Linear::dot(double **matriz1, int **matriz2, int *vector1, int *vector2, int tam, int n){

    //Matriz1: factores
    //Matriz2: trayectoria
//     ofstream myfile1/*, myfile6*/;
//     myfile1.open ("matrices.txt", ios::app);
//     myfile1 << "En dot: " << endl;
//     myfile6.open ("trayectoria.txt", ios::app);
    
    //Multiplicamos una matriz continua por un vector discreto, lo que da como resultado un vector nuevo continuo.
  
    double *newvector;
    newvector = new double[tam]; 
  
    for(int i = 0; i < tam; i++){
        newvector[i] = 0;
        for(int j = 0; j < tam; j++){
            newvector[i] = newvector[i] + (vector1[j] * matriz1[i][j]);
        }
//         myfile1 << " nwv" << i << ": " << newvector[i];
    }
//     myfile1 << "\n" << endl;
    
    //Convertimos ese vector nuevo de continuo a discreto.
    
    for(int i = 0; i < tam; i++){
        if(newvector[i] < 0){ matriz2[n][i] = -1;}
        if((newvector[i] == 0) & (n == 1)){ matriz2[n][i] = -1;} //Al inicio tomamos como que el gen está apagado
        if((newvector[i] == 0) & (n > 1)){ matriz2[n][i] = matriz2[n - 1][i];}
        if(newvector[i] > 0){ matriz2[n][i] = 1;}
//          myfile1 << matriz2[n][i];
    }
// myfile1 << "\n";
    delete [] newvector;
    
    
    //Guardamos el vector resultado en la matriz de vectores
    
    for(int i = 0; i < tam; i++){
        vector2[i] = matriz2[n][i];
    }

//     myfile1.close ();
//     myfile6.close ();
}

void Linear::dot(double **matriz1, int **matriz2, int **matriz3, int tam1, int tam2, int n, int p)
{    //Para obtener fenotipo matrizNp por matrizAtractor
    //Multiplicamos una matriz continua por una matriz discreta, lo que da como resultado una matriz nueva continua.
    double **newmat;
    hor.espacio(newmat, p, tam2);

    for(int i = 0; i < p; i++){
        for(int j = 0; j < tam2; j++){
            newmat[i][j] = 0;
            for(int k = 0; k < tam1; k++){
                newmat[i][j] = newmat[i][j] + (matriz1[j][k]*matriz2[n - 1 - p + i][k]);
            }
        }
    }
    //Convertimos cada elemento de la matriz nueva de continuo a discreto.
    
    for(int i = 0; i < p; i++){
        for(int j = 0; j < tam2; j++){
        if(newmat[i][j] < 0){ matriz3[i][j] = -1;}
        if((newmat[i][j] == 0) & (i == 0)){ matriz3[i][j] = -1;} //Al inicio tomamos como que el gen está apagado
        if((newmat[i][j] == 0) & (i != 0)){ matriz3[i][j] = matriz3[i - 1][j];}
        if(newmat[i][j] > 0){ matriz3[i][j] = 1;}
        }
    }
  for (int i = 0; i < p; i++)
      delete [] newmat[i];
  delete [] newmat;
}

void Linear::dot2(double **matriz1, int **matriz2, int **matriz3, int tam1, int tam2, int n, int p, int m, int veces)
{    
//     Matriz1: Marcadores
//     Matriz2: Trayectoria
//     Matriz3: Fenotipo
    //Multiplicamos una matriz continua por una matriz discreta, lo que da como resultado una matriz nueva continua.
    if((m == 37) & (veces == 7)){cout << "dot 1\n";}
    double **newmat;
//     if(m == 32 & veces == 119){cout << "dot 1.1\n";}
    hor.espacio(newmat, p, tam2);
//     if(m == 32 & veces == 119){cout << "dot 1.2\n";}

    for(int i = 0; i < p; i++){
//         if(m == 32 & veces == 119){cout << "i: " << i << endl;}
        for(int j = 0; j < tam2; j++){
//             if(m == 32 & veces == 119){cout << "j: " << j << endl;}
            newmat[i][j] = 0;
            for(int k = 0; k < tam1; k++){
                if((m == 37) & (veces == 7)){
                cout << "k: " << k << endl;
                cout << "nm: " << newmat[i][j] << endl;
                cout << "m1: " << matriz1[j][k] << endl;
                cout << "hils de m2: " << n - 1 - p + i << endl;
                cout << "m2: " << matriz2[n - 1 - p + i][k] << endl;
            }
                newmat[i][j] = newmat[i][j] + (matriz1[j][k]*matriz2[n - 1 - p + i][k]);
//                 if(m == 32 & veces == 119){cout << "nm i: " << i << " j: " << j << " k: " << k << endl;}
            }
        }
    }
if((m == 37) & (veces == 7)){cout << "dot 2\n";}
    //Convertimos cada elemento de la matriz nueva de continuo a discreto.
    
    for(int i = 0; i < p; i++){
        for(int j = 0; j < tam2; j++){
        if(newmat[i][j] < 0){ matriz3[i][j] = -1;}
        if((newmat[i][j] == 0) & (i == 0)){ matriz3[i][j] = -1;} //Al inicio tomamos como que el gen está apagado
        if((newmat[i][j] == 0) & (i != 0)){ matriz3[i][j] = matriz3[i - 1][j];}
        if(newmat[i][j] > 0){ matriz3[i][j] = 1;}
        }
    }
// if(m == 32 & veces == 119){cout << "dot 3\n";}    
  for (int i = 0; i < p; i++)
      delete [] newmat[i];
  delete [] newmat;
}
