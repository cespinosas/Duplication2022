///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///cd Documents/Doctorado/Resultados/libs
//g++ -Wall -c *.cc
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gataca.h"
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


GATACA::GATACA()
{
}

GATACA::GATACA(int semilla) //No sé si está de más.
{
  hor.start_rng(semilla);
}

void GATACA::condin(int *vector, int cols, int especial)
{
  for(int i = 0; i < cols; i++){
      if(i == especial){vector[i] = 1;}
      else{vector[i] = -1;}
  }
}
void GATACA::condin(int *vector, int cols)
{
  for(int i = 0; i < cols; i++){
      vector[i] = -1;
  }
}

void GATACA::escalon(int *vector, int cols)
{
  for(int i = 0; i < cols; i++)
    if(hor.randgauss(1) < 0)
      vector[i] = -1;
    else{vector[i] = 1;}
}

void GATACA::mut(double **matriz, int hils, int cols, int tolerancia, double intnw, int& x, int& y, double& valor)
{
  double coin;
  bool res;
  int conex;
  
  //      Horse hor(semilla);
  //     Linear lin(semilla);
  //     Finder find;
  int minconex = (intnw * (hils * cols)) - tolerancia;
  int maxconex = (intnw * (hils * cols)) + tolerancia;
    
  do{
    conex = find.numintnw(matriz, hils, cols);
    x = hor.randflat(0, hils);
    y = hor.randflat(0, cols);
    valor = matriz[x][y];
    res = true;
    coin = hor.coin();
    
    if(coin == valor){
      res = false;
    } //Si el nuevo valor es igual al anterior, no se vale y se intenta de nuevo.
    if(coin == 0)
    {
      if(conex - 1 >= minconex){
        if(matriz[x][y] != 0){
          matriz[x][y] = coin; //Si coin = 0, matriz[x][y] no, y se vale perderlo; entonces matriz[x][y] = 0 (coin) [res = true].
        }
        else{
          res = false;
        } //Si coin = 0 y matriz[x][y] = 0, aunque se vale perderlo, no nos sirve y se sale.
      }
      else{
        res = false;} //Si no se vale perderlo, se sale.
    }
    else{
      if(matriz[x][y] == 0){
        if(conex + 1 <= maxconex){
          matriz[x][y] = hor.randgauss(1); //Si no había, coin dice que haya y se puede, matriz[x][y] ya /*no será cero.
        }
        else{
          res = false;} //Si ya había y se conserva, pero no se vale. Se sale.
      }
      else{
        if(matriz[x][y] > 0){
          matriz[x][y] = coin; //Si se cambia el valor de matriz[x][y], quien además es positiva. Se conserva el signo +.
        }
        else{matriz[x][y] = coin * -1;
        } //Si matriz[x][y] era negativa, se conserva el signo -.
      }
    }
  }
  while(res == false);
  
}

void GATACA::choose(double **Matriz1, int hils1, int cols1, double **Matriz2, int hils2, int cols2, int& matmut){
  /*
   Horse hor(semilla);
   GATACA gat(semilla);*/
  
  //     ofstream myfile1;
  
  int cont1 = hils1 * cols1; //Cuenta número de elementos en la matriz1
  int cont2 = hils2 * cols2; //Cuenta número de elementos en la matriz2
  double elementos = cont1 + cont2; //Elementos totales
  double peso1 = cont1 / elementos; //Peso de la matriz1
  double peso2 = cont2 / elementos; //Peso de la matriz2
  double ruleta = hor.randuniform(0,1);
  
  if(peso1 == peso2){ //Si los pesos son iguales
    if(ruleta < 0.5)
      matmut = 1;
    else{matmut = 2;}
  }
  if(peso1 < peso2){
    if(ruleta < peso1)
      matmut = 1;
    else{matmut = 2;}
  }
  else{if(ruleta < peso2)
    matmut = 2;
  else{matmut = 1;}
  }
}

void GATACA::back(double **Matriz, int x, int y, double valor){
  Matriz[x][y] = valor;
}

void GATACA::walk(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int tam1, int tam2, int tolerancia, double intnw, int *vectorentrada, int num, int& n, int pfenotipo, int &pasos, int &stanby, int m)
{
  //     Matriz1: Factores
  //     Matriz2: Marcadores
  //     Matriz3: Trayectoria inicial
  //     Matriz4: Fenotipo
  //     num: Máximo de steps
  //     tam1 = hils1 = cols1 = cols2 = Nf
  //     tam2 = hils2 = Np
    
    double dist, valor;
    int per, npend, x, y;
    double **red; //Una matriz que gobierne a todas
    hor.espacio(red, tam1 + tam2, tam1);
    int **trayectoria;
    hor.espacio(trayectoria, 200, tam1);
    int **keeptrayectoria;
    hor.espacio(keeptrayectoria, 200, tam1); //Para guardar trayectoria en caso de regresarse (como en red)
    int **Fenotipo2;
    
    for(int i = 0; i < n; i++){
            for(int j = 0; j < tam1; j++){
                trayectoria[i][j] = Matriz3[i][j]; //Por si los primeros pasos se salen de la red neutral fenotípica
            }
        }
    
    do{
        //llenamos matriz madre
        for(int i = 0; i < tam1; i++){
            for(int j = 0; j < tam1; j++){
                red[i][j] = Matriz1[i][j];
            }
        }
        
        for(int i = tam1; i < tam1 + tam2; i++){
            for(int j = 0; j < tam1; j++){
                red[i][j] = Matriz2[i - tam1][j];
            }
        }
        
        for(int i = 0; i < n; i++){
            for(int j = 0; j < tam1; j++){
                keeptrayectoria[i][j] = trayectoria[i][j]; //Por si los primeros pasos se salen de la red neutral fenotípica
            }
        }
        
        //Mutamos un elemento de la red
        mut(red, tam1 + tam2, tam1, tolerancia, intnw, x, y, valor);
        
        //Separamos a la red para obtener fenotipo
        for(int i = 0; i < tam1; i++){
            for(int j = 0; j < tam1; j++){
                Matriz1[i][j] = red[i][j];
            }
        }
        
        for(int i = tam1; i < tam1 + tam2; i++){
            for(int j = 0; j < tam1; j++){
                Matriz2[i - tam1][j] = red[i][j];
            }
        }
        
        //Obtenemos atractor
        npend = 1;
        for(int i = 0; i < tam1; i++){
            vectorentrada[i] = Matriz3[0][i]; //Que trabaje desde la condición incial, sino n = 1 no sirve
        }
            find.atractor(Matriz1, trayectoria, vectorentrada, tam1, npend, per);//Actualiza el valor de n y p
            hor.espacio(Fenotipo2, per, tam2);
        //Obtenemos Fenotipo2, puede salir de más de una hilera
        lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, npend, per);
        //Comparamos fenotipos
        dist = find.distfenotipica(Matriz4, pfenotipo, Fenotipo2, per, tam2);
        
        //Aceptación o rechazo de la mutación
        if(dist != 0 || pfenotipo != per){ //010319
            pasos++;
            stanby++;
            back(red, x, y, valor); //Regresamos a la última red que dio el fenotipo
            for(int i = 0; i < tam1; i++){
                for(int j = 0; j < tam1; j++){
                    Matriz1[i][j] = red[i][j];
                }
            }
            for(int i = tam1; i < tam1 + tam2; i++){
                for(int j = 0; j < tam1; j++){
                    Matriz2[i - tam1][j] = red[i][j];
                }
            }
            for(int i = 0; i < n; i++){
                for(int j = 0; j < tam1; j++){
                    trayectoria[i][j] = keeptrayectoria[i][j]; //Si el paso se sale de la red neutral fenotípica
                }
            }
        }
        else{
            pasos++;
            n = npend; //Es importante que n salga correcta porque luego entra a exploración
        }

        for (int i = 0; i < per; i++)
            delete [] Fenotipo2[i];
        delete [] Fenotipo2;

    }while(pasos < num);
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < tam1; j++){
            Matriz3[i][j] = trayectoria[i][j]; //Actualizamos trayectoria de la red antes de salir
        }
    }
    
    for (int i = 0; i < (tam1 + tam2); i++)
        delete [] red[i];
    delete [] red;
    
    for (int i = 0; i < 200; i++)
        delete [] trayectoria[i];
    delete [] trayectoria;
    
    for (int i = 0; i < 200; i++)
        delete [] keeptrayectoria[i];
    delete [] keeptrayectoria;

  ofstream guardapasos;
//   guardapasos.open ("Condin/outs/Pasos/muestrapasos.txt", ios::app);
  guardapasos.open ("outs/Pasos/muestrapasos.txt", ios::app);
  guardapasos << n - 1 << " ";
  guardapasos.close ();
}

void GATACA::cambiafenotipo(double **Matriz1, int *Vector1, int **Matriz2, int tam1, int tam2, int tam3)
{
  //     Matriz1: Matrizp
  //     Vector1: Wish (Fenotipo deseado)
  //     Matriz2: Fenotipo
  //     Tam1: Nf
  //     Tam2: Np
  //     Tam3: p
    for(int i = 0; i < tam3; i++){
        for(int j = 0; j < tam2; j++){
//             cout << i << " " << j << " " << Vector1[j] << " " << Matriz2[i][j];
//             if(Vector1[j] == Matriz2[i][j]){cout << " Igual\n";}
            if(Vector1[j] != Matriz2[i][j]){ /*cout << " Diferente\n";*/
                for(int k = 0; k < tam1; k++)
                    Matriz1[j][k] = Matriz1[j][k] * -1;
            }
        }
    }
    
}

void GATACA::inducirduplicacion(double** &Matriz1, double** &Matriz2, double **Matriz3, double **Matriz4, int x, int tam1, int tam2)
{
  //Matriz1: Factores duplicadas
  //Matriz2: Marcadores duplicadas
  //Matriz3: Factores
  //Matriz4: Marcadores
  //Tam1: Nf
  //Tam2: Np
  
  //Rellenamos matrices con las no duplicadas
  for(int i = 0; i < tam1; i++)
    for(int j = 0; j < tam1; j++)
      Matriz1[i][j] = Matriz3[i][j];
  
  for(int i = 0; i < tam2; i++)
    for(int j = 0; j < tam1; j++)
      Matriz2[i][j] = Matriz4[i][j];
  
  //Extendemos las matrices 1 columna
  for(int i = 0; i < tam1; i++)
    Matriz1[i][tam1] = Matriz1[i][x];
  for(int i = 0; i < tam2; i++)
    Matriz2[i][tam1] = Matriz2[i][x];
  //Extendamos la matriz de factores 1 hilera
  for(int i = 0; i < tam1; i++)
    Matriz1[tam1][i] = Matriz3[x][i];
  //Rellenamos el espacio faltante en factores
  Matriz1[tam1][tam1] = Matriz1[x][x];
}

void GATACA::control1(double** &Matriz1, double** &Matriz2, double **Matriz3, double **Matriz4, int tam1, int tam2, double intnw)
{
  //Matriz1: Factores extendida
  //Matriz2: Marcadores extendida
  //Matriz3: Factores
  //Matriz4: Marcadores
  //Tam1: Nf
  //Tam2: Np
  //intnw
  //Rellenamos matrices con las existentes y el resto con ceros
  for(int i = 0; i < tam1; i++)
    for(int j = 0; j < tam1; j++)
      Matriz1[i][j] = Matriz3[i][j];
  
  for(int i = 0; i < tam2; i++)
    for(int j = 0; j < tam1; j++)
      Matriz2[i][j] = Matriz4[i][j];
  
  //     for(int i = 0; i < tam1 + 1; i++){
  //         Matriz1[i][tam1] = 0;
  //         Matriz1[tam1][i] = 0;
  //     }
  //     for(int i = 0; i < tam2, i++)
  //         Matriz2[i][tam1] = 0;
  
  //Extendemos matriz de factores
  int celdas = (tam1 + 1) * (tam1 + 1);
  int novac = hor.round(intnw, celdas);
  //     cout << "novac: " << novac << endl;
  int vac = celdas - novac;
  //     cout << "vac: " << vac << endl;
  int sumnovac = find.numintnw(Matriz1, tam1, tam1);
  //     cout << "sumnovac: " << sumnovac << endl;
  int sumvac = (tam1 * tam1) - sumnovac;
  //     cout << "sumvac: " << sumvac << endl;
  double real;
  
  //Hileras de M1
  for (int i = 0; i < tam1 + 1; i++){
    real = hor.randreal();
    // //         cout << "*****i: " << i << endl;
    //         cout << "real: " << real << endl;
    //Si el volado dice que celda diferente de cero.
    if (real < intnw){
      //Nos aseguramos que falten celdas diferentes de cero.
      if(sumnovac < novac){
        //Si no hay suficientes, se llena con gauss.
        Matriz1[i][tam1] = hor.randgauss(1);
        sumnovac++;
      }
      else {
        //Si ya hay suficientes, la celda será cero.
        Matriz1[i][tam1] = 0;
        sumvac++;
      }
    }
    //Si el volado dice que la celda es cero.
    else {
      //Nos asegurames que faten celdas diferentes de cero.
      if(sumvac < vac){
        //Si no hay suficientes, se llena con cero.
        Matriz1[i][tam1] = 0;
        sumvac++;
      }
      //Si ya hay suficientes ceros, será diferente de cero.
      else{
        Matriz1[i][tam1] = hor.randgauss(1);
        sumnovac++;
      }
    }
    //             cout << "m1" << i << tam1 << ": " << Matriz1[i][tam1] << endl;
    //             cout << "sumnovac: " << sumnovac << endl;
    //             cout << "sumvac: " << sumvac << endl;
  }
  
  //Columnas de M1
  for (int i = 0; i < tam1; i++){ //Hasta tam1 porque m[tam1 + 1][tam1 + 1] ya se creó arriba.
    real = hor.randreal();
    //         cout << "*****i: " << i << endl;
    //         cout << "real: " << real << endl;
    //Si el volado dice que celda diferente de cero.
    if (real < intnw){
      //Nos aseguramos que falten celdas diferentes de cero.
      if(sumnovac < novac){
        //Si no hay suficientes, se llena con gauss.
        Matriz1[tam1][i] = hor.randgauss(1);
        sumnovac++;
      }
      else {
        //Si ya hay suficientes, la celda será cero.
        Matriz1[tam1][i] = 0;
        sumvac++;
      }
    }
    //Si el volado dice que la celda es cero.
    else {
      //Nos asegurames que faten celdas diferentes de cero.
      if(sumvac < vac){
        //Si no hay suficientes, se llena con cero.
        Matriz1[tam1][i] = 0;
        sumvac++;
      }
      //Si ya hay suficientes ceros, será diferente de cero.
      else{
        Matriz1[tam1][i] = hor.randgauss(1);
        sumnovac++;
      }
    }
    //             cout << "m1" << i << tam1 << ": " << Matriz1[i][tam1] << endl;
    //             cout << "sumnovac: " << sumnovac << endl;
    //             cout << "sumvac: " << sumvac << endl;
  }
  
  
  //Extendemos matriz de marcadores
  int celdas2 = tam2 * (tam1 + 1);
  int novac2 = hor.round(intnw, celdas2);
  //     cout << "novac2: " << novac2 << endl;
  int vac2 = celdas2 - novac2;
  //     cout << "vac2: " << vac2 << endl;
  int sumnovac2 = find.numintnw(Matriz2, tam2, tam1);
  //     cout << "sumnovac2: " << sumnovac2 << endl;
  int sumvac2 = (tam2 * tam1) - sumnovac2;
  //     cout << "sumvac2: " << sumvac2 << endl;
  double real2;
  
  //Hileras de M2
  for (int i = 0; i < tam2; i++){
    real2 = hor.randreal();
    //         cout << "*****i: " << i << endl;
    //         cout << "real2: " << real2 << endl;
    //Si el volado dice que celda diferente de cero.
    if (real2 < intnw){
      //Nos aseguramos que falten celdas diferentes de cero.
      if(sumnovac2 < novac2){
        //Si no hay suficientes, se llena con gauss.
        Matriz2[i][tam1] = hor.randgauss(1);
        sumnovac++;
      }
      else {
        //Si ya hay suficientes, la celda será cero.
        Matriz2[i][tam1] = 0;
        sumvac++;
      }
    }
    //Si el volado dice que la celda es cero.
    else {
      //Nos asegurames que faten celdas diferentes de cero.
      if(sumvac2 < vac2){
        //Si no hay suficientes, se llena con cero.
        Matriz2[i][tam1] = 0;
        sumvac++;
      }
      //Si ya hay suficientes ceros, será diferente de cero.
      else{
        Matriz2[i][tam1] = hor.randgauss(1);
        sumnovac++;
      }
    }
    //             cout << "m2" << i << tam1 << ": " << Matriz2[i][tam1] << endl;
    //             cout << "sumnovac2: " << sumnovac2 << endl;
    //             cout << "sumvac2: " << sumvac2 << endl;
  }
}

void GATACA::control2(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2)//Misma cantidad de interacciones que el gen duplicado
{
  //     Matriz1: CMatrizf2
  //     Matriz2: CMatrizp2
  //     Matriz3: Matrizfd
  //     Matriz4: Matrizpd
  //     num1 y 2: Nf, Np
  //Rellenamos matrices con las existentes
  for(int i = 0; i < num1; i++)
    for(int j = 0; j < num1; j++)
      Matriz1[i][j] = Matriz3[i][j];
  
  for(int i = 0; i < num2; i++)
    for(int j = 0; j < num1; j++)
      Matriz2[i][j] = Matriz4[i][j];
  
  for(int i = 0; i < num1 + 1; i++){
    Matriz1[i][num1] = 0;
    Matriz1[num1][i] = 0;
  }
  
  for(int i = 0; i < num2; i++)
    Matriz2[i][num1] = 0;
  
  int conex1 = find.numintnw(Matriz3, num1 + 1, num1 + 1) + find.numintnw(Matriz4, num2, num1 + 1);
  int conex2 = find.numintnw(Matriz1, num1 + 1, num1 + 1) + find.numintnw(Matriz2, num2, num1 + 1);
  int sum = 0;
  double real1, real2;
  int comodin1, comodin2;
  
  //Hileras de factores
  do{
    real1 = hor.randgauss(1);
    real2 = hor.randgauss(1);
    //Si el volado dice que celda diferente de cero.
    if (real1 < 0){
      comodin1 = hor.randflat(0,num1 + 1);
      if(real2 < 0){
        Matriz1[comodin1][num1] = hor.randgauss(1);
        sum++;
      }
      else {
        Matriz1[num1][comodin1] = hor.randgauss(1);
        sum++;
      }
    }
    else {
      comodin2 = hor.randflat(0,num2);
      Matriz2[comodin2][num1] = hor.randgauss(1);
      sum++;
    }
    conex1 = find.numintnw(Matriz3, num1 + 1, num1 + 1) + find.numintnw(Matriz4, num2, num1 + 1);
    conex2 = find.numintnw(Matriz1, num1 + 1, num1 + 1) + find.numintnw(Matriz2, num2, num1 + 1);
  }while(conex2 < conex1);
}

void GATACA::exploracion(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek)
{
  //     Matriz1: Factores
  //     Matriz2: Marcadores
  //     Matriz3: Trayectoria
  //     Matriz4: Fenotipo original
  //     Matriz5: Distancias
  //     num: Máximo de steps
  //     tam1 = hils1 = cols1 = cols2 = Nf
  //     tam2 = hils2 = Np
  //     num1: Número de muestra para saber en qué hilera va la distancia
  //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
    
    int **Fenotipo2;
    int *vectorentrada;
    hor.espacio(vectorentrada, tam1);
    ///////////////////////////////////////////////////
    int **Fenotipos;
    hor.espacio(Fenotipos, 10000, tam2);
    int *Periodos;
    hor.espacio(Periodos, num2);
    ofstream vecpasos, phenotypes, nuevo, diferente/*, acces*/, index;
    //   int cuentafenotipos= 0;//Cuenta de número de fenotipos diferentes
    for(int i = 0; i < 10000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos[i][j] = 0;
        }
    }
    for(int i = 0; i < num2; i++){Periodos[i] = 0;}
    int per1, per2, n2, hileras = 0;
    Fenotiponuevo = 0;
    double distf, dist;
    int **Comparafenotipo;//Para tomar fenotipo por fenotipo de la matriz de fenotipos para comparar
    //   double *distf; //Distancia entre fenotipo pasado y acutal
    //   hor.espacio(distf, num2);
    ///////////////////////////////////////////////////
    int **trayectoria; //Trayectoria para cuando se muten Nf's
    hor.espacio(trayectoria, 1000, tam1);
    double valor = 0;
    int x = 0, y = 0, matmut = 0; //Nos dicen cuál matriz y en qué elemento fue mutada y su valor anterior para poder
    int per;
    int veces = 0;
//     double valor2;
    //Robusto mut no cambia
//     int cuenta_gen = 0;
//     int cuenta_cero_nocero = 0;
//     int cuenta_nocero_cero = 0;
//     int cuenta_cambiosigno = 0;
//     int cuenta_cambiovalor = 0;
    //No robusto mut cambia
//     int cuenta_gen2 = 0;
//     int cuenta_cero_nocero2 = 0;
//     int cuenta_nocero_cero2 = 0;
//     int cuenta_cambiosigno2 = 0;
//     int cuenta_cambiovalor2 = 0;
    //Redundancia
    ifstream fi;
    int cuenta_redundancia_dup = 0;
    int cuenta_noredundancia_dup = 0;
    double dist_redundancia_dup = 0;
    double dist_noredundancia_dup = 0;
    double distsum_redundancia_dup = 0;
    double distsum_noredundancia_dup = 0;
    int cuenta_weack_dup = 0;
    int cuenta_moderate_dup = 0;
    int cuenta_strong_dup = 0;
    int cuenta_lethal_dup = 0;
    int cuenta_weack2_dup = 0;
    int cuenta_moderate2_dup = 0;
    int cuenta_strong2_dup = 0;
    int cuenta_lethal2_dup = 0;
    
    int cuenta_redundancia_c1 = 0;
    int cuenta_noredundancia_c1 = 0;
    double dist_redundancia_c1 = 0;
    double dist_noredundancia_c1 = 0;
    double distsum_redundancia_c1 = 0;
    double distsum_noredundancia_c1 = 0;
    int cuenta_weack_c1 = 0;
    int cuenta_moderate_c1 = 0;
    int cuenta_strong_c1 = 0;
    int cuenta_lethal_c1 = 0;
    int cuenta_weack2_c1 = 0;
    int cuenta_moderate2_c1 = 0;
    int cuenta_strong2_c1 = 0;
    int cuenta_lethal2_c1 = 0;
    
    int cuenta_redundancia_c2 = 0;
    int cuenta_noredundancia_c2 = 0;
    double dist_redundancia_c2 = 0;
    double dist_noredundancia_c2 = 0;
    double distsum_redundancia_c2 = 0;
    double distsum_noredundancia_c2 = 0;
    int cuenta_weack_c2 = 0;
    int cuenta_moderate_c2 = 0;
    int cuenta_strong_c2 = 0;
    int cuenta_lethal_c2 = 0;
    int cuenta_weack2_c2 = 0;
    int cuenta_moderate2_c2 = 0;
    int cuenta_strong2_c2 = 0;
    int cuenta_lethal2_c2 = 0;
    
    int *Vec_eventos_dup;
    hor.espacio(Vec_eventos_dup, 1000);
    hor.open_ifstream(fi, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/acces/eventos_dup.txt");
    for(int i = 0; i < 1000; i++){fi >> Vec_eventos_dup[i];}
    fi.close();
    
    int *Vec_eventos_c1;
    hor.espacio(Vec_eventos_c1, 1000);
    hor.open_ifstream(fi, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/acces/eventos_c1.txt");
    for(int i = 0; i < 1000; i++){fi >> Vec_eventos_c1[i];}
    fi.close();
    
    int *Vec_eventos_c2;
    hor.espacio(Vec_eventos_c2, 1000);
    hor.open_ifstream(fi, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/acces/eventos_c2.txt");
    for(int i = 0; i < 1000; i++){fi >> Vec_eventos_c2[i];}
    fi.close();
    
    
  //////////////////////Rmutparciales/////////////////////////////////

  ofstream fo;
  do{
      choose(Matriz1, tam1, tam1, Matriz2, tam2, tam1, matmut);
      //Si choose fue para la matriz de factores de transcripción
      if(matmut == 1){
//           cuentaup++; //Contamos que la mutación está arriba
          n2 = 1; //Tiene su propio tamaño de trayectoria
          mut(Matriz1, tam1, tam1, tolerancia, intnw, x, y, valor);
//           valor2 = Matriz1[x][y];
          //Buscamos atractor, puede salir de más de una hilera
          for(int i = 0; i < tam1; i++){
              vectorentrada[i] = Matriz3[0][i];
              trayectoria[0][i] = Matriz3[0][i];
          }//Que trabaje desde la condición incial, sino n = 1 no sirve
          find.atractor(Matriz1, trayectoria, vectorentrada, tam1, n2, per); //Tiene su propia trayectoria
          hor.espacio(Fenotipo2, per, tam2);
          //Obtenemos Fenotipo2, puede salir de más de una hilera
          lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, n2, per); //Saca su fenotipo vecino
          
          //###################Robustez redundante
      //Veces que la mutación cayó en gen con duplicado
//      if(y == Vec_eventos_dup[num1]){
//          cuenta_redundancia_dup++;
//          distsum_redundancia_dup += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_redundancia_dup = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_redundancia_dup == 0){cuenta_lethal_dup++;}
//          else{if(dist_redundancia_dup < 0.66){cuenta_strong_dup++;}
//              else{if(dist_redundancia_dup < 0.83){cuenta_moderate_dup++;}
//              else{cuenta_weack_dup++;}
//              }
//          }
//      }
//      else{
//          cuenta_noredundancia_dup++;
//          distsum_noredundancia_dup += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_noredundancia_dup = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_noredundancia_dup == 0){cuenta_lethal2_dup++;}
//          else{if(dist_noredundancia_dup < 0.66){cuenta_strong2_dup++;}
//              else{if(dist_noredundancia_dup < 0.83){cuenta_moderate2_dup++;}
//              else{cuenta_weack2_dup++;}
//              }
//          }
//      }
//      //Veces que la mutación cayó en gen con CI
//      if(y == Vec_eventos_c1[num1]){
//          cuenta_redundancia_c1++;
//          distsum_redundancia_c1 += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_redundancia_c1 = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_redundancia_c1 == 0){cuenta_lethal_c1++;}
//          else{if(dist_redundancia_c1 < 0.66){cuenta_strong_c1++;}
//              else{if(dist_redundancia_c1 < 0.83){cuenta_moderate_c1++;}
//              else{cuenta_weack_c1++;}
//              }
//          }
//      }
//      else{
//          cuenta_noredundancia_c1++;
//          distsum_noredundancia_c1 += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_noredundancia_c1 = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_noredundancia_c1 == 0){cuenta_lethal2_c1++;}
//          else{if(dist_noredundancia_c1 < 0.66){cuenta_strong2_c1++;}
//              else{if(dist_noredundancia_c1 < 0.83){cuenta_moderate2_c1++;}
//              else{cuenta_weack2_c1++;}
//              }
//          }
//      }
//      //Veces que la mutación cayó en gen con CII
//      if(y == Vec_eventos_c2[num1]){
//          cuenta_redundancia_c2++;
//          distsum_redundancia_c2 += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_redundancia_c2 = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_redundancia_c2 == 0){cuenta_lethal_c2++;}
//          else{if(dist_redundancia_c2 < 0.66){cuenta_strong_c2++;}
//              else{if(dist_redundancia_c2 < 0.83){cuenta_moderate_c2++;}
//              else{cuenta_weack_c2++;}
//              }
//          }
//      }
//      else{
//          cuenta_noredundancia_c2++;
//          distsum_noredundancia_c2 += 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          dist_noredundancia_c2 = 1 - find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//          if(dist_noredundancia_c2 == 0){cuenta_lethal2_c2++;}
//          else{if(dist_noredundancia_c2 < 0.66){cuenta_strong2_c2++;}
//              else{if(dist_noredundancia_c2 < 0.83){cuenta_moderate2_c2++;}
//              else{cuenta_weack2_c2++;}
//              }
//          }
//      }
      //######################################################
          
          //Regresamos a la matriz muestra
          back(Matriz1, x, y, valor);
      }
      //Si choose fue para la matriz de marcadores
      if(matmut == 2){
//           cuentadown++; //Contamos que la mutación está abajo
          mut(Matriz2, tam2, tam1, tolerancia, intnw, x, y, valor);
          //Obtenemos Fenotipo2
          per = p; //Ya tiene la trayectoria y atractor definidos
          hor.espacio(Fenotipo2, per, tam2);
          lin.dot(Matriz2, Matriz3, Fenotipo2, tam1, tam2, n, per); //Saca su fenotipo con la trayectoria del espécimen
//           distparcial = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//           if(distparcial != 0){down++;}
          //Regresamos a la matriz muestra
          back(Matriz2, x, y, valor);
      }
      dist = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//       cout << "dist" << dist << endl;
      if(dist != 0){//Nos aseguramos de no contar el fenotipo original
          cuenta++;
          //TIPO DE ROBUSTEZ 140720
//          if(matmut == 1){
//              if((y == num1)||(y==tam1)){cuenta_gen++;}
//              if(valor == 0){}
//          }
          //
          //###############Guardamos fenotipos###############270420
          hor.open_ofstream(fo, "outs/acces/fn/fn_" + hor.inttostring(num1) + ".txt");
          for(int j = 0; j < per; j++){
              for(int k = 0; k < tam2; k++){
                  fo << Fenotipo2[j][k] << " ";
              }fo << endl;
          }fo << endl;
          fo.close();
          hor.open_ofstream(fo, "outs/acces/fn/per_" + hor.inttostring(num1) + ".txt");
          fo << per << endl;
          fo.close();
          //#################################################
//           cout << "Fdif = " << cuenta << endl;
          if(cuenta == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
              Fenotiponuevo++;
//               cout << "Fnuevo1 = " << Fenotiponuevo << endl;
              valordek[0] = 1;
//               cout << "1 valordek[0]= " << valordek[0] << endl;
              Periodos[0] = per; //Almacenamos el primer periodo
              for(int i = 0; i < per; i++){
                  for(int j = 0; j < tam2; j++){
                      Fenotipos[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo diferente
                  }
              }
          }
          else{//Si ya tenemos un fenotipo para comparar
              Periodos[cuenta - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
              hileras += Periodos[cuenta - 2];
              per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
              for(int i = hileras; i < (hileras + per); i++){
                  for(int j = 0; j < tam2; j++){
                      Fenotipos[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                  }
              }
              for(int k = 0; k < (cuenta - 1); k++){//Comparar todos los fenotipos uno a uno
                  per2 = Periodos[k];//sacamos periodo de cada fenotipo a comparar
                  hor.espacio(Comparafenotipo, per2, tam2);
                  for(int i = per1; i < (per1 + per2); i++){
                      for(int j = 0; j < tam2; j++){
                          Comparafenotipo[i - per1][j] = Fenotipos[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                      }
                  }
                  distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                  hor.borrar(Comparafenotipo, per2);
                  if(distf == 0){
                      valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
//                       cout << "valordek[" << k << "]= " << valordek[k] << endl;
                      break;
                  }//Break porque ya sabemos que no es nuevo y ya no necesita seguir comparando
                  //                 else{hor.borrar(Comparafenotipo, per2);}
                  if(k == (cuenta - 2)){
                      Fenotiponuevo++;
//                       cout << "Fnuevo = " << Fenotiponuevo << endl;
                      valordek[k + 1] = 1; //Suma 1 al fenotipo k que se repite //051019
//                       cout << "valordeknuevo[" << k + 1 << "]= " << valordek[k + 1] << endl;
                  }
                  per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
              }
          }
      }
      
      //Guardamos la distancia en la matriz de distancias
    Matriz5[veces] = 1 - dist;
    
    hor.borrar(Fenotipo2, per);
    
    veces++;
    
//     vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
//     vecpasos << n -1 << " ";
//     vecpasos.close ();
  }while(veces != num2);
  
  //Guardamos info de robustez por redundancia DUP
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_redundancia_dup.txt");
//  fo << cuenta_redundancia_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_redundancia_dup.txt");
//  fo << distsum_redundancia_dup / cuenta_redundancia_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_noredundancia_dup.txt");
//  fo << cuenta_noredundancia_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_noredundancia_dup.txt");
//  fo << distsum_noredundancia_dup / cuenta_noredundancia_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_red_dup.txt");
//  fo << cuenta_lethal_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_red_dup.txt");
//  fo << cuenta_strong_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_red_dup.txt");
//  fo << cuenta_moderate_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_red_dup.txt");
//  fo << cuenta_weack_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_nored_dup.txt");
//  fo << cuenta_lethal2_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_nored_dup.txt");
//  fo << cuenta_strong2_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_nored_dup.txt");
//  fo << cuenta_moderate2_dup << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_nored_dup.txt");
//  fo << cuenta_weack2_dup << endl;
//  fo.close();
//
//  //Guardamos info de robustez por redundancia CI
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_redundancia_c1.txt");
//  fo << cuenta_redundancia_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_redundancia_c1.txt");
//  fo << distsum_redundancia_c1 / cuenta_redundancia_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_noredundancia_c1.txt");
//  fo << cuenta_noredundancia_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_noredundancia_c1.txt");
//  fo << distsum_noredundancia_c1 / cuenta_noredundancia_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_red_c1.txt");
//  fo << cuenta_lethal_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_red_c1.txt");
//  fo << cuenta_strong_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_red_c1.txt");
//  fo << cuenta_moderate_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_red_c1.txt");
//  fo << cuenta_weack_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_nored_c1.txt");
//  fo << cuenta_lethal2_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_nored_c1.txt");
//  fo << cuenta_strong2_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_nored_c1.txt");
//  fo << cuenta_moderate2_c1 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_nored_c1.txt");
//  fo << cuenta_weack2_c1 << endl;
//  fo.close();
//
//  //Guardamos info de robustez por redundancia CII
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_redundancia_c2.txt");
//  fo << cuenta_redundancia_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_redundancia_c2.txt");
//  fo << distsum_redundancia_c2 / cuenta_redundancia_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Cuenta_noredundancia_c2.txt");
//  fo << cuenta_noredundancia_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/Dist_noredundancia_c2.txt");
//  fo << distsum_noredundancia_c2 / cuenta_noredundancia_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_red_c2.txt");
//  fo << cuenta_lethal_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_red_c2.txt");
//  fo << cuenta_strong_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_red_c2.txt");
//  fo << cuenta_moderate_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_red_c2.txt");
//  fo << cuenta_weack_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/lethal_nored_c2.txt");
//  fo << cuenta_lethal2_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/strong_nored_c2.txt");
//  fo << cuenta_strong2_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/moderate_nored_c2.txt");
//  fo << cuenta_moderate2_c2 << endl;
//  fo.close();
//
//  hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/EffectAntes/weack_nored_c2.txt");
//  fo << cuenta_weack2_c2 << endl;
//  fo.close();
//
      //#############################################
  
//   vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
//   vecpasos << "\n";
//   vecpasos.close ();

  hor.borrar(Fenotipos, 10000);
  hor.borrar(trayectoria, 1000);
  hor.borrar(Periodos);
  hor.borrar(vectorentrada);
}

void GATACA::exploraciondup(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int **Matriz5, double *Vector1, double *Vector2, double& alldists, double& notoldists_noori, int& notolrm_noori, int tam1, int tam2, int tolerancia, double intnw, int n, int pd, int p, int num1, int num2, int &cuenta1, int &cuenta2, int &cuenta3, int &Fnuevo, int &Fnuevoin, int tipo, int *valordek, int m)
{
  //     Matriz1: Factores
  //     Matriz2: Marcadores
  //     Matriz3: Trayectoria
  //     Matriz4: Fenotipo del duplicado
  //     Matriz5: Fenotipo original
  //     Vector1: Distancia con el fenotipo del duplicado
  //     Vector2: Distancia con el fenotipo original
  //     alldists: Distancia del fmutar con f ori
  //     num: Máximo de steps
  //     tam1 = hils1 = cols1 = cols2 = Nf + 1
  //     tam2 = hils2 = Np
  //     num1: Número de gen para saber en qué hilera va la distancia
  //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
  //     cuenta1: Redes duplicadas que conservan f original, y que al mutar No conservaron fenotipo.
  //     cuenta2: Redes duplicadas que no conservan f original, y que al mutar No conserva fenotipo.
  //     cuenta3: Los que regresan.
  
//     ofstream abc;

    int **Fenotipo2;
    int **Fenotipos1;
    hor.espacio(Fenotipos1, 100000, tam2);
    int **trayectoria; //Trayectoria para cuando se muten Nf's
    hor.espacio(trayectoria, 1000, tam1);
    int *Periodos1;
    hor.espacio(Periodos1, num2);
    int **Fenotipos2;
    hor.espacio(Fenotipos2, 100000, tam2);
    int *Periodos2;
    hor.espacio(Periodos2, num2);
    ofstream dvecpasos, c1vecpasos, c2vecpasos, phenotypes;
//   string tipoa = tipo;
//   int cuentafenotipos= 0;//Cuenta de número de fenotipos diferentes

    for(int i = 0; i < 100000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos1[i][j] = 0;
        }
    }
    for(int i = 0; i < 100000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos2[i][j] = 0;
        }
    }
    for(int i = 0; i < num2; i++){Periodos1[i] = 0;}
    for(int i = 0; i < num2; i++){Periodos2[i] = 0;}
    int per1, per2, hileras = 0, n2;
    Fnuevo = 0;
    Fnuevoin = 0;
    double distf;

    int **Comparafenotipo;//Para tomar fenotipo por fenotipo de la matriz de fenotipos para comparar
    double valor = 0;
    int x = 0, y = 0, matmut = 0; //Nos dicen cuál matriz y en qué elemento fue mutada y su valor anterior para poder
    int per;
    int veces = 0;
    double dist0, dist1, dist2, dist3;
    int *vectorentrada;
    hor.espacio(vectorentrada, tam1);
    ofstream fo;
    int cuenta_redundancia = 0;
    int cuenta_noredundancia = 0;
    double dist_redundancia = 0;
    double dist_noredundancia = 0;
    double distsum_redundancia = 0;
    double distsum_noredundancia = 0;
    int cuenta_weack = 0;
    int cuenta_moderate = 0;
    int cuenta_strong = 0;
    int cuenta_lethal = 0;
    int cuenta_weack2 = 0;
    int cuenta_moderate2 = 0;
    int cuenta_strong2 = 0;
    int cuenta_lethal2 = 0;
    double dist_no_ori;
  do{
  choose(Matriz1, tam1, tam1, Matriz2, tam2, tam1, matmut);
  //Si choose fue para la matriz de factores de transcripción
  if(matmut == 1){
      n2 = 1;
      mut(Matriz1, tam1, tam1, tolerancia, intnw, x, y, valor);
      //Buscamos atractor, puede salir de más de una hilera
      for(int i = 0; i < tam1; i++){
          vectorentrada[i] = Matriz3[0][i]; //Que trabaje desde la condición incial, sino n = 1 no sirve
          trayectoria[0][i] = Matriz3[0][i];
      }
      find.atractor(Matriz1, trayectoria, vectorentrada, tam1, n2, per);
      hor.espacio(Fenotipo2, per, tam2);
      //Obtenemos Fenotipo2, puede salir de más de una hilera
      lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, n2, per);
      
      //###################Robustez redundante
      //Veces que la mutación cayó en gen con duplicado
//      if(y == num1){
//          cuenta_redundancia++;
//          distsum_redundancia += 1 - find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
//          dist_redundancia = 1 - find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
//          if(dist_redundancia == 0){cuenta_lethal++;}
//          else{if(dist_redundancia < 0.66){cuenta_strong++;}
//              else{if(dist_redundancia < 0.83){cuenta_moderate++;}
//              else{cuenta_weack++;}
//              }
//          }
//      }
//      else{
//          cuenta_noredundancia++;
//          distsum_noredundancia += 1 - find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
//          dist_noredundancia = 1 - find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
//          if(dist_noredundancia == 0){cuenta_lethal2++;}
//          else{if(dist_noredundancia < 0.66){cuenta_strong2++;}
//              else{if(dist_noredundancia < 0.83){cuenta_moderate2++;}
//              else{cuenta_weack2++;}
//              }
//          }
//      }
      //######################################################
      
      //Regresamos a la matriz muestra
      back(Matriz1, x, y, valor);
      
  }
    //Si choose fue para la matriz de marcadores
    if(matmut == 2){
        mut(Matriz2, tam2, tam1, tolerancia, intnw, x, y, valor);
        //Obtenemos Fenotipo2
        per = pd;
        hor.espacio(Fenotipo2, per, tam2);
        lin.dot(Matriz2, Matriz3, Fenotipo2, tam1, tam2, n, per);
        
        //Regresamos a la matriz muestra
        back(Matriz2, x, y, valor);
    }
    //Tomamos distancia con el fenotipo de la muestra
    
    //Únicamente si el fenotipo del duplicado es diferente al fenotipo original
    
    //////////////////Robustez a mutaciones después de duplicar, de redes que conservaron el fenotipo original//////////////////
    //////////////////Accesibilidad después de duplicar: fenotipos nuevos de redes que conservaron el fenotipo original después de duplicar//////////////////
    dist3 = find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
    alldists += dist3;
    
    dist0 = find.distfenotipica(Matriz4, pd, Matriz5, p, tam2); //Comparamos fenotipo de duplicar vs original
    if(dist0 == 0){//Si al duplicar mantuvieron el fenotipo original
        dist1 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Comparamos el fenotipo del duplicado vs el de después de mutar
        if(dist1 != 0){//El fenotipo al mutar no se parece al original (fmdup = Fori)
            //###############Guardamos fenotipos###############270420
            hor.open_ofstream(fo, "outs/acces/fn/fn_" + hor.inttostring(m) + "_"+ hor.inttostring(tipo) +"_" + hor.inttostring(num1) + ".txt");
            for(int j = 0; j < per; j++){
                for(int k = 0; k < tam2; k++){
                    fo << Fenotipo2[j][k] << " ";
                }fo << endl;
            }fo << endl;
            fo.close();
            hor.open_ofstream(fo, "outs/acces/fn/per_" + hor.inttostring(m) + "_"+ hor.inttostring(tipo) +"_" + hor.inttostring(num1) + ".txt");
            fo << per << endl;
            fo.close();
            cuenta1++; //Para la robustez a mutaciones después de duplicar, contamos fenotipos que no se mantuvieron [Rob = (1 - cuenta)/vecinos]
            if(cuenta1 == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                Fnuevo++;
                valordek[0] = 1;
                Periodos1[0] = per; //Almacenamos el primer periodo
                for(int i = 0; i < per; i++){
                    for(int j = 0; j < tam2; j++){
                        Fenotipos1[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo
                    }
                }
            }
            else{//Si es diferente del original
                Periodos1[cuenta1 - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                hileras += Periodos1[cuenta1 - 2];
                per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                for(int i = hileras; i < (hileras + per); i++){
                    for(int j = 0; j < tam2; j++){
                        Fenotipos1[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                    }
                }
                for(int k = 0; k < (cuenta1 - 1); k++){
                    //Comparar todos los fenotipos uno a uno
                    per2 = Periodos1[k];//sacamos periodo de cada fenotipo a comparar
                    hor.espacio(Comparafenotipo, per2, tam2);
                    for(int i = per1; i < (per1 + per2); i++){
                        for(int j = 0; j < tam2; j++){
                            Comparafenotipo[i - per1][j] = Fenotipos1[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                        }
                    }
                    distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                    hor.borrar(Comparafenotipo, per2);
                    if(distf == 0){
                        valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
                        break;
                    }//Break porque ya sabemos que no es nuevo y ya no necesita seguir comparando
//                     else{hor.borrar(Comparafenotipo, per2);}
                    if(k == (cuenta1 - 2)){
                        Fnuevo++;
                        valordek[k + 1] = 1; //Suma 1 al fenotipo k que se repite //051019
                    }
                    per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                }
            }
        }
    //Guardamos la distancia del f de mutar cuando no se parece ni al fori ni al fdup (fori = fdup)
    Vector2[veces] = 1 - dist1;
    Vector1[veces] = 2;//Para identificar las celdas que no se deben tomar en cuenta porque el fdup es igual al fori
    }
    else{//Si el fdup no se parece al fenotipo ori (M4!=M5)
        Vector2[veces] = 2;//Para identificar las celdas que no se deben tomar en cuenta porque el fdup no es igual al fori
        //////////////////Robustez a mutaciones in después de duplicar, de redes que no conservaron el f original//////////////////
        ////////////////Accesibilidad después de duplicar: fenotipos nuevos de redes que NO conservaron el fori después de duplicar NI el fdup////////////////
        //Tomamos distancia con el fenotipo de la muestra
        dist3 = find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);//Buscamos si al mutar regresa al fenotipo original.
        if(dist3 == 0){
            dist2 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Vemos la distancia entre el fenotipo original y el fenotipo del no tolerante a la duplicación.
            notoldists_noori += dist2;
            notolrm_noori++;
//            dist2 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Fdup vs mutar. No se mantuvo el fenotipo del duplicado.
            cuenta3++;//if(dist2 != 0) 260120 PORQUE SON LOS QUE REGRESAN
        }//Si regresa al fenotipo original
        else{//Si el fenotipo al mutar es diferente al original y el de duplicar es diferente al original también
            dist2 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Fdup vs mutar. No se mantuvo el fenotipo del duplicado.
            notoldists_noori += dist2; //18 abril 22, para tener similaridad de no tolerantes a su nuevo fenotipo
            if(dist2 != 0){
                notolrm_noori++; //26Jun22, 1 - robustez
                //Cuando el fenotipo que resulta de mutar no se parece al original
                //###############Guardamos fenotipos###############270420
                hor.open_ofstream(fo, "outs/acces/fn/fn_" + hor.inttostring(m) + "_"+ hor.inttostring(tipo) +"_" + hor.inttostring(num1) + ".txt");
                for(int j = 0; j < per; j++){
                    for(int k = 0; k < tam2; k++){
                        fo << Fenotipo2[j][k] << " ";
                    }fo << endl;
                }fo << endl;
                fo.close();
                hor.open_ofstream(fo, "outs/acces/fn/per_" + hor.inttostring(m) + "_"+ hor.inttostring(tipo) +"_" + hor.inttostring(num1) + ".txt");
                fo << per << endl;
                fo.close();
                //#################################################
                cuenta2++;//Para la robustez a mutaciones in después de duplicar, contamos fenotipos que no se mantuvieron [Robin = (1 - cuenta)/vecinos]
//                 cout << "Encontró fdif= " << cuenta2 << endl;
                if(cuenta2 == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                    Fnuevoin++;
//                     cout << "1 Encontró Fnuevoin= " << Fnuevoin << endl;
                    valordek[0] = 1;
//                     cout << "1 valordek[0]= " << valordek[0] << endl;
                    Periodos2[0] = per; //Almacenamos el primer periodo
                    for(int i = 0; i < per; i++){
                        for(int j = 0; j < tam2; j++){
                            Fenotipos2[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo
                        }
                    }
                }
                else{//Si ya tenemos un fenotipo para comparar
//                     if(m == 32 & veces == 124){cout << "n\n";}
                    Periodos2[cuenta2 - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                    hileras += Periodos2[cuenta2 - 2];
                    per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                    for(int i = hileras; i < (hileras + per); i++){
                        for(int j = 0; j < tam2; j++){
                            Fenotipos2[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                        }
                    }
                    for(int k = 0; k < (cuenta2 - 1); k++){//Comparar todos los fenotipos uno a uno
                        per2 = Periodos2[k];//sacamos periodo de cada fenotipo a comparar
                        hor.espacio(Comparafenotipo, per2, tam2);
                        for(int i = per1; i < (per1 + per2); i++){
                            for(int j = 0; j < tam2; j++){
                                Comparafenotipo[i - per1][j] = Fenotipos2[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                            }
                        }
                        distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                        hor.borrar(Comparafenotipo, per2);
                        if(distf == 0){//Si se encontró otro fenotipo igual en la lista de fenotipos
                            valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
//                             cout << "valordek["<<k<<"]= " << valordek[k] << endl;
                            break;
                        }//Break porque ya sabemos que es nuevo y ya no necesita seguir comparando
//                         else{hor.borrar(Comparafenotipo, per2);}
                        if(k == (cuenta2 - 2)){
                            Fnuevoin++;
//                             cout << "Encontró Fnuevoin= " << Fnuevoin << endl; //051019
                            valordek[k + 1] = 1; //Suma 1 al fenotipo k que se repite //051019
//                             cout << "valordek["<<k + 1<<"]= " << valordek[k + 1] << endl; //051019
                        }
                        per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                    }
                }
            }
        }
        //Guardamos la distancia en la matriz de distancias
        Vector1[veces] = 1 - dist3;//Rm con distancia al original
    }
    
    hor.borrar(Fenotipo2, per);
    veces++;
//    if(tipo == 1){
////         dvecpasos.open ("Condin/outs/Pasos/dvecpasos.txt", ios::app);
//        dvecpasos.open ("outs/Pasos/dvecpasos.txt", ios::app);
//        dvecpasos << n -1 << " ";
//        dvecpasos.close ();
//    }
//    if(tipo == 2){
////         c1vecpasos.open ("Condin/outs/Pasos/c1vecpasos.txt", ios::app);
//        c1vecpasos.open ("outs/Pasos/c1vecpasos.txt", ios::app);
//        c1vecpasos << n -1 << " ";
//        c1vecpasos.close ();
//    }
//    if(tipo == 3){
////         c2vecpasos.open ("Condin/outs/Pasos/c2vecpasos.txt", ios::app);
//        c2vecpasos.open ("outs/Pasos/c2vecpasos.txt", ios::app);
//        c2vecpasos << n -1 << " ";
//        c2vecpasos.close ();
//    }
//     if(m == 32 & veces == 124){cout << "r\n";}
  }while(veces != num2);
    

//  if(tipo == 1){
//         dvecpasos.open ("Condin/outs/Pasos/dvecpasos.txt", ios::app);
//        dvecpasos.open ("outs/Pasos/dvecpasos.txt", ios::app);
//        dvecpasos << "\n";
//        dvecpasos.close ();
      //#############################################
      //Guardamos info de robustez por redundancia
//     hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_redundancia_dup.txt");
//     if(num1 == 11){
//         fo << cuenta_redundancia << endl;
//     }
//     else{fo << cuenta_redundancia << " ";}
//     fo.close();
//
//     hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_redundancia_dup.txt");
//     if(num1 == 11){
//         fo << distsum_redundancia / cuenta_redundancia << endl;
//     }
//     else{fo << distsum_redundancia / cuenta_redundancia << " ";}
//     fo.close();
//
//     hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_noredundancia_dup.txt");
//     if(num1 == 11){
//         fo << cuenta_noredundancia << endl;
//     }
//     else{fo << cuenta_noredundancia << " ";}
//     fo.close();
//
//     hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_noredundancia_dup.txt");
//     if(num1 == 11){
//         fo << distsum_noredundancia / cuenta_noredundancia << endl;
//     }
//     else{fo << distsum_noredundancia / cuenta_noredundancia << " ";}
//     fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_red_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_lethal << endl;
//      }
//      else{fo << cuenta_lethal << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_red_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_strong << endl;
//      }
//      else{fo << cuenta_strong << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_red_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_moderate << endl;
//      }
//      else{fo << cuenta_moderate << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_red_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_weack << endl;
//      }
//      else{fo << cuenta_weack << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_nored_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_lethal2 << endl;
//      }
//      else{fo << cuenta_lethal2 << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_nored_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_strong2 << endl;
//      }
//      else{fo << cuenta_strong2 << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_nored_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_moderate2 << endl;
//      }
//      else{fo << cuenta_moderate2 << " ";}
//      fo.close();
//
//      hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_nored_dup.txt");
//      if(num1 == 11){
//          fo << cuenta_weack2 << endl;
//      }
//      else{fo << cuenta_weack2 << " ";}
//      fo.close();
//
//      //#############################################
//    }
//    if(tipo == 2){
////         c1vecpasos.open ("Condin/outs/Pasos/c1vecpasos.txt", ios::app);
////        c1vecpasos.open ("outs/Pasos/c1vecpasos.txt", ios::app);
////        c1vecpasos << "\n";
////        c1vecpasos.close ();
//        //#############################################
//        //Guardamos info de robustez por redundancia
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_redundancia_c1.txt");
//       if(num1 == 11){
//           fo << cuenta_redundancia << endl;
//       }
//       else{fo << cuenta_redundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_redundancia_c1.txt");
//       if(num1 == 11){
//           fo << distsum_redundancia / cuenta_redundancia << endl;
//       }
//       else{fo << distsum_redundancia / cuenta_redundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_noredundancia_c1.txt");
//       if(num1 == 11){
//           fo << cuenta_noredundancia << endl;
//       }
//       else{fo << cuenta_noredundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_noredundancia_c1.txt");
//       if(num1 == 11){
//           fo << distsum_noredundancia / cuenta_noredundancia << endl;
//       }
//       else{fo << distsum_noredundancia / cuenta_noredundancia << " ";}
//       fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_red_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_lethal << endl;
//        }
//        else{fo << cuenta_lethal << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_red_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_strong << endl;
//        }
//        else{fo << cuenta_strong << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_red_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_moderate << endl;
//        }
//        else{fo << cuenta_moderate << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_red_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_weack << endl;
//        }
//        else{fo << cuenta_weack << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_nored_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_lethal2 << endl;
//        }
//        else{fo << cuenta_lethal2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_nored_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_strong2 << endl;
//        }
//        else{fo << cuenta_strong2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_nored_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_moderate2 << endl;
//        }
//        else{fo << cuenta_moderate2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_nored_c1.txt");
//        if(num1 == 11){
//            fo << cuenta_weack2 << endl;
//        }
//        else{fo << cuenta_weack2 << " ";}
//        fo.close();
//
//        //#############################################
//    }
//    if(tipo == 3){
////         c2vecpasos.open ("Condin/outs/Pasos/cnum12vecpasos.txt", ios::app);
////        c2vecpasos.open ("outs/Pasos/cnum12vecpasos.txt", ios::app);
////        c2vecpasos << "\n";
////        c2vecpasos.close ();
//        //#############################################
//        //Guardamos info de robustez por redundancia
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_redundancia_c2.txt");
//       if(num1 == 11){
//           fo << cuenta_redundancia << endl;
//       }
//       else{fo << cuenta_redundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_redundancia_c2.txt");
//       if(num1 == 11){
//           fo << distsum_redundancia / cuenta_redundancia << endl;
//       }
//       else{fo << distsum_redundancia / cuenta_redundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Cuenta_noredundancia_c2.txt");
//       if(num1 == 11){
//           fo << cuenta_noredundancia << endl;
//       }
//       else{fo << cuenta_noredundancia << " ";}
//       fo.close();
//
//       hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Dist_noredundancia_c2.txt");
//       if(num1 == 11){
//           fo << distsum_noredundancia / cuenta_noredundancia << endl;
//       }
//       else{fo << distsum_noredundancia / cuenta_noredundancia << " ";}
//       fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_red_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_lethal << endl;
//        }
//        else{fo << cuenta_lethal << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_red_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_strong << endl;
//        }
//        else{fo << cuenta_strong << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "//home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_red_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_moderate << endl;
//        }
//        else{fo << cuenta_moderate << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_red_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_weack << endl;
//        }
//        else{fo << cuenta_weack << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/lethal_nored_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_lethal2 << endl;
//        }
//        else{fo << cuenta_lethal2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/strong_nored_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_strong2 << endl;
//        }
//        else{fo << cuenta_strong2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/moderate_nored_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_moderate2 << endl;
//        }
//        else{fo << cuenta_moderate2 << " ";}
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/Rm/Redundancia/Effect/weack_nored_c2.txt");
//        if(num1 == 11){
//            fo << cuenta_weack2 << endl;
//        }
//        else{fo << cuenta_weack2 << " ";}
//        fo.close();
        
        //#############################################
//    }
    
    
    
    
  hor.borrar(Fenotipos1, 100000);
  hor.borrar(trayectoria, 1000);
  hor.borrar(Fenotipos2, 100000);
  hor.borrar(Periodos1);
  hor.borrar(Periodos2);
}

void GATACA::mediadist(double **Matriz1, double *Vector1, int num1, int num2)
{
  //     Matriz1: Distancias, una hilera por muestra
  //     Vector1: Promedios de distancias, columna por muestra
  //     num1: Cantidad de muestras para saber número de hileras
  //     num2: Número de exploraciones, número de columnas
  double sum;
  for(int i = 0; i < num1; i++){
    sum = 0;
    for(int j = 0; j < num2; j++){
      sum = sum + Matriz1[i][j];
    }
    Vector1[i] = sum / num2;
  }
}

void GATACA::ceros(double **Matriz1, double *Vector1, int hils, int cols)
{
  //     Matriz1: Distancias, una hilera por muestra
  //     Vector1: Proporción de ceros por muestra
  //     hils: Número de muestras
  //     cols: Número de exploraciones
  int ceros;
  
  for(int i = 0; i < hils; i++){
    ceros = 0;
    for(int j = 0; j < cols; j++){
      if(Matriz1[i][j] == 0)
        ceros++;
    }
    Vector1[i] = (double)ceros / cols;
  }
}

void GATACA::robustezfenotipica(double *Vector1, double *Vector2, int tam)
{
  for(int i = 0; i < tam; i++)
    Vector2[i] = 1 - Vector1[i];
}

double GATACA::Shannon(int *valordek, int tam)
{
    double nlogn, fi, sum, H;
    //Cantidad de fenotipos
    int tamk = 0;
    int *vecesxfenotipo;
    int n = 0; //Número total de mutaciones
    int k = 0; //Número de fenotipos diferentes
    for(int i = 0; i < tam; i++){
                if(valordek[i] != 0){
                    tamk++;
                }
            }
     hor.espacio(vecesxfenotipo, tamk);
     //Quitamos ceros
     for(int i = 0; i < tam; i++){
         if(valordek[i] != 0){
             vecesxfenotipo[k] = valordek[i];
             n = n + vecesxfenotipo[k];
             k++;
//              cout << "valordek: " << valordek[i] << endl;
        }
     }
     //Ahora sí sacamos el índice
     nlogn = n*log(n); //n = suma de los valores de los datos, log base e
     sum = 0;
     for(int i = 0; i < tamk; i++){
         fi = vecesxfenotipo[i];
         sum = sum + (fi*log(fi));
     }
     H = (nlogn - sum)/n;
     
//      cout << "n: " << n << " nlogn: " << nlogn << " sum: " << sum << " H: " << H << endl;
     
     if(H==0){cout << "Sh1 n: " << n <<  " k: " << k << " nlogn: " << nlogn << " sum: " << sum << " H: " << H << endl;}
     
     return H;
}

double GATACA::Shannon2(int *valordek, int tam)
{
//     ofstream datadup2;
    double nlogn, fi, sum, H;
    //Cantidad de fenotipos
    int tamk = 1; //Por el ori
    int *vecesxfenotipo;
    int n = tam;//Por el ori
    int N = 0; //No ori
    int k = 1;
    for(int i = 0; i < tam; i++){
                if(valordek[i] != 0){
                    tamk++;
                }
            }
     hor.espacio(vecesxfenotipo, tamk);
     //Quitamos ceros
     for(int i = 0; i < tam; i++){
         if(valordek[i] != 0){
             vecesxfenotipo[k] = valordek[i];
             N = N + vecesxfenotipo[k];
             k++;
        }
     }
     vecesxfenotipo[0] = 432 - N; //Fori
     
//      if(redes == 0){
//          datadup2.open("outs/Shannon/datadup2.txt", ios::app);
//          for(int i = 0; i < tamk; i++){
//              datadup2 << vecesxfenotipo[i] << endl;
//          }
//          datadup2.close();
//      }
     
     //Ahora sí sacamos el índice
     nlogn = (n)*log(n);
     sum = 0;
     for(int i = 0; i < tamk; i++){
         fi = vecesxfenotipo[i];
         sum = sum + (fi*log(fi));
     }
     H = (nlogn - sum)/(n);
     
     if(H==0){cout << "Sh2 n: " << n << " nlogn: " << nlogn << " sum: " << sum << " H: " << H << endl;}     
     
     return H;
}

void GATACA::Shannon3(int *valordek, int tam)
{
//     ofstream datadup2;
    double nlogn, fi, sum, H;
    //Cantidad de fenotipos
    int tamk = 1; //Por el ori
    int *vecesxfenotipo;
    int n = tam;//Por el ori
    int N = 0; //No ori
    int k = 1;
    for(int i = 0; i < tam; i++){
                if(valordek[i] != 0){
                    tamk++;
                }
            }
     hor.espacio(vecesxfenotipo, tamk);
     //Quitamos ceros
     for(int i = 0; i < tam; i++){
         if(valordek[i] != 0){
             vecesxfenotipo[k] = valordek[i];
             N = N + vecesxfenotipo[k];
             k++;
        }
     }
     vecesxfenotipo[0] = 432 - N; //Fori
//      
     //Ahora sí sacamos el índice
     nlogn = (n)*log(n);
     sum = 0;
     for(int i = 0; i < tamk; i++){
         fi = vecesxfenotipo[i];
         sum = sum + (fi*log(fi));
     }
     H = (nlogn - sum)/(n);
     
//      cout << "N: " << N << " nlogn: " << nlogn << " sum: " << sum << " H: " << H << endl;
     
}

double GATACA::rmutarriba(double **Matriz1, int hils1, int cols1, int **Matriz2, int hils2, int p, double tolerancia, double intnw, double& PromDist)
{
    //Matriz1: Matriz a mutar
    //Matriz2: Trayectoria de la Matriz1 [n][nf]
    //Vector1: Condición inicial [Nf]
    //Nf: hils1, cols1, cols2
    //n original: hils2
    int veces = 0;
    int x, y, n, per;
    int cuenta = 0;
    double valor, dist, sum;
    int celdas = hils1 * cols1;
    int mutaciones = 2 * celdas;
    int *vectorentrada, **trayectoria;
    hor.espacio(vectorentrada, hils1);
    hor.espacio(trayectoria, 1000, hils1);
    hor.fillm0(trayectoria, 1000, hils1);
    int **atractor1, **atractor2;
    hor.espacio(atractor1, p, hils1);
    for(int i = 0; i < p; i++){
        for(int j = 0; j < hils1; j++){
            atractor1[i][j] = Matriz2[hils2 - 1 - p + i][j]; //Guardamos el atractor de la trayectoria original
//             cout << atractor1[i][j] << endl;
        }
    }
//     hor.printmat(Matriz2, hils2, cols1);
//     hor.printmat(atractor1, p, hils1);
    sum = 0;
    do{
        veces++;
        x = hor.randflat(0, hils1);
        y = hor.randflat(0, cols1);
        
        n = 1; //Tiene su propio tamaño de trayectoria
        mut(Matriz1, hils1, cols1, tolerancia, intnw, x, y, valor);
        //Buscamos atractor, puede salir de más de una hilera
        for(int i = 0; i < cols1; i++){
            vectorentrada[i] = Matriz2[0][i];
            trayectoria[0][i] = Matriz2[0][i];
        }//Que trabaje desde la condición incial, sino n = 1 no sirve
        find.atractor(Matriz1, trayectoria, vectorentrada, hils1, n, per); //Tiene su propia trayectoria
        hor.espacio(atractor2, per, hils1);
        for(int i = 0; i < per; i++){
            for(int j = 0; j < hils1; j++){
                atractor2[i][j] = trayectoria[n - 1 - per + i][j]; //Guardamos el atractor de la trayectoria nueva
            }
        }
//         hor.printmat(trayectoria, n, cols1);
//         hor.printmat(atractor2, per, cols1);
        dist = find.distfenotipica(atractor1, p, atractor2, per, hils1);
        sum += dist;
        if(dist != 0){cuenta++;}
//         cout << "veces: " << veces << " dist: " << dist << " cuenta: " << cuenta << endl;
        //Regresamos a la matriz muestra
        back(Matriz1, x, y, valor);
    }
    while(veces < mutaciones);
    
    PromDist = 1 - (sum / mutaciones);
    
    hor.borrar(vectorentrada);
    hor.borrar(trayectoria, 1000);
    hor.borrar(atractor1, p);
    hor.borrar(atractor2, per);
    
    return cuenta;
}



double GATACA::rmutabajo(double **Matriz1, int hils1, int cols1, int **Matriz2, int hils2, int **Matriz3, int hils3, double tolerancia, double intnw, double& PromDist)
{
    //Matriz1: Matriz a mutar
    //Matriz2: Fenotipo [p][np]
    //Matriz3: Trayectoria de la Matriz1 [n][nf]
    //Nf: cols1, cols3
    //Np: hils1, cols2
    //p: hils2
    //n: hils3
    
    int veces = 0;
    int x, y;
    int cuenta = 0;
    double valor, dist, sum;
    int celdas = hils1 * cols1;
    int mutaciones = 2 * celdas;
    
    int **Fenotipo2;
    hor.espacio(Fenotipo2, hils2, hils1);
    sum = 0;
    do{
        veces++;
        x = hor.randflat(0, hils1);
        y = hor.randflat(0, cols1);
        
        mut(Matriz1, hils1, cols1, tolerancia, intnw, x, y, valor);
        //Obtenemos Fenotipo2
        lin.dot(Matriz1, Matriz3, Fenotipo2, cols1, hils1, hils3, hils2);
        //lin.dot(Matriz2, Matriz3, Fenotipo2, Nf, Np, n, per); //Saca su fenotipo con la trayectoria del espécimen
        dist = find.distfenotipica(Matriz2, hils2, Fenotipo2, hils2, hils1);
        sum += dist;
        if(dist != 0){cuenta++;}
        //Regresamos a la matriz muestra
        back(Matriz1, x, y, valor);
    }
    while(veces < mutaciones);
    
    PromDist = 1 - (sum / mutaciones);
    
    hor.borrar(Fenotipo2, hils2);
               
    return cuenta;
}

void GATACA::exploracionduppuntofijo(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, int **Matriz5, double *Vector1, double *Vector2, double& alldists, int tam1, int tam2, int tolerancia, double intnw, int n, int pd, int p, int num1, int num2, int &cuenta1, int &cuenta2, int &cuenta3, int &Fnuevo, int &Fnuevoin, int tipo, int *valordek, int &puntosfijos)
{
    //     Matriz1: Factores
    //     Matriz2: Marcadores
    //     Matriz3: Trayectoria
    //     Matriz4: Fenotipo del duplicado
    //     Matriz5: Fenotipo original
    //     Vector1: Distancia con el fenotipo del duplicado
    //     Vector2: Distancia con el fenotipo original
    //     alldists: Distancia del fmutar con f ori
    //     num: Máximo de steps
    //     tam1 = hils1 = cols1 = cols2 = Nf + 1
    //     tam2 = hils2 = Np
    //     num1: Número de gen para saber en qué hilera va la distancia
    //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
    //     cuenta1: Redes duplicadas que conservan f original, y que al mutar No conservaron fenotipo.
    //     cuenta2: Redes duplicadas que no conservan f original, y que al mutar No conserva fenotipo.
    //     cuenta3: Los que regresan.
    
    //     ofstream abc;
    
    int **Fenotipo2;
    int **Fenotipos1;
    hor.espacio(Fenotipos1, 100000, tam2);
    int **trayectoria; //Trayectoria para cuando se muten Nf's
    hor.espacio(trayectoria, 1000, tam1);
    int *Periodos1;
    hor.espacio(Periodos1, num2);
    int **Fenotipos2;
    hor.espacio(Fenotipos2, 100000, tam2);
    int *Periodos2;
    hor.espacio(Periodos2, num2);
    ofstream dvecpasos, c1vecpasos, c2vecpasos, phenotypes;
    
    double distpuntofijo;
    bool res;
    int *pormientras;
    hor.espacio(pormientras, tam2);
    
    for(int i = 0; i < 100000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos1[i][j] = 0;
        }
    }
    for(int i = 0; i < 100000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos2[i][j] = 0;
        }
    }
    for(int i = 0; i < num2; i++){Periodos1[i] = 0;}
    for(int i = 0; i < num2; i++){Periodos2[i] = 0;}
    int per1, per2, hileras = 0, n2;
    Fnuevo = 0;
    Fnuevoin = 0;
    double distf;
    int **Comparafenotipo;//Para tomar fenotipo por fenotipo de la matriz de fenotipos para comparar
    double valor = 0;
    int x = 0, y = 0, matmut = 0; //Nos dicen cuál matriz y en qué elemento fue mutada y su valor anterior para poder
    int per;
    int veces = 0;
    double dist0, dist1, dist2, dist3;
    int *vectorentrada;
    hor.espacio(vectorentrada, tam1);
    ofstream fo;
    do{
        res = false;
        choose(Matriz1, tam1, tam1, Matriz2, tam2, tam1, matmut);
        //Si choose fue para la matriz de factores de transcripción
        if(matmut == 1){
            n2 = 1;
            mut(Matriz1, tam1, tam1, tolerancia, intnw, x, y, valor);
            //Buscamos atractor, puede salir de más de una hilera
            for(int i = 0; i < tam1; i++){
                vectorentrada[i] = Matriz3[0][i]; //Que trabaje desde la condición incial, sino n = 1 no sirve
                trayectoria[0][i] = Matriz3[0][i];
            }
            find.atractor(Matriz1, trayectoria, vectorentrada, tam1, n2, per);
            hor.espacio(Fenotipo2, per, tam2);
            //Obtenemos Fenotipo2, puede salir de más de una hilera
            lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, n2, per);
            
            //Regresamos a la matriz muestra
            back(Matriz1, x, y, valor);
            
        }
        //Si choose fue para la matriz de marcadores
        if(matmut == 2){
            mut(Matriz2, tam2, tam1, tolerancia, intnw, x, y, valor);
            //Obtenemos Fenotipo2
            per = pd;
            hor.espacio(Fenotipo2, per, tam2);
            lin.dot(Matriz2, Matriz3, Fenotipo2, tam1, tam2, n, per);
            
            //Regresamos a la matriz muestra
            back(Matriz2, x, y, valor);
        }
        //Tomamos distancia con el fenotipo de la muestra
        
        //Únicamente si el fenotipo del duplicado es diferente al fenotipo original
        
        //////////////////Robustez a mutaciones después de duplicar, de redes que conservaron el fenotipo original//////////////////
        
        //////////////////////////////////PARA PUNTOS FIJOS//////////////////////////////
        if(per > 1){
            for(int a = 0; a < tam2; a++){
                pormientras[a] = Fenotipo2[0][a];
            }
            distpuntofijo = find.distfenotipica(pormientras, Fenotipo2, per, tam2);//Por si a = aaa
            if(distpuntofijo == 0){res = true;}
        }
        ////////////////////////////////////////////////////////////////////////////////
        if((per == 1)||(res == true)){
            puntosfijos++;
            dist3 = find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);
            alldists += dist3;
            
            dist0 = find.distfenotipica(Matriz4, pd, Matriz5, p, tam2); //Comparamos fenotipo de duplicar vs original
            if(dist0 == 0){//Si al duplicar mantuvieron el fenotipo original
                dist1 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Comparamos el fenotipo del duplicado vs el de después de mutar
                if(dist1 != 0){//El fenotipo al mutar no se parece al original (fmdup = Fori)
                    cuenta1++; //Para la robustez a mutaciones después de duplicar, contamos fenotipos que no se mantuvieron [Rob = (1 - cuenta)/vecinos]
                    if(cuenta1 == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                        Fnuevo++;
                        Periodos1[0] = per; //Almacenamos el primer periodo
                        for(int i = 0; i < per; i++){
                            for(int j = 0; j < tam2; j++){
                                Fenotipos1[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo
                            }
                        }
                    }
                    else{//Si es diferente del original
                        Periodos1[cuenta1 - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                        hileras += Periodos1[cuenta1 - 2];
                        per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                        for(int i = hileras; i < (hileras + per); i++){
                            for(int j = 0; j < tam2; j++){
                                Fenotipos1[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                            }
                        }
                        for(int k = 0; k < (cuenta1 - 1); k++){
                            //Comparar todos los fenotipos uno a uno
                            per2 = Periodos1[k];//sacamos periodo de cada fenotipo a comparar
                            hor.espacio(Comparafenotipo, per2, tam2);
                            for(int i = per1; i < (per1 + per2); i++){
                                for(int j = 0; j < tam2; j++){
                                    Comparafenotipo[i - per1][j] = Fenotipos1[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                                }
                            }
                            distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                            hor.borrar(Comparafenotipo, per2);
                            if(distf == 0){
                                valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
                                break;
                            }//Break porque ya sabemos que no es nuevo y ya no necesita seguir comparando
                            //                     else{hor.borrar(Comparafenotipo, per2);}
                            if(k == (cuenta1 - 2)){
                                Fnuevo++;
                            }
                            per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                        }
                    }
                }
                //Guardamos la distancia del f de mutar cuando no se parece ni al fori ni al fdup (fori = fdup)
                Vector2[veces] = 1 - dist1;
                Vector1[veces] = 2;//Para identificar las celdas que no se deben tomar en cuenta porque el fdup es igual al fori
            }
            else{//Si el fdup no se parvecesfenotipoece al fori (M4!=M5)
                Vector2[veces] = 2;//Para identificar las celdas que no se deben tomar en cuenta porque el fdup no es igual al fori
                //////////////////Robustez a mutaciones in después de duplicar, de redes que no conservaron el f original//////////////////
                ////////////////Accesibilidad después de duplicar: fenotipos nuevos de redes que NO conservaron el fori después de duplicar NI el fdup////////////////
                //Tomamos distancia con el fenotipo de la muestra
                dist3 = find.distfenotipica(Matriz5, p, Fenotipo2, per, tam2);//Buscamos si al mutar regresa al fenotipo original.
                if(dist3 == 0){
                    dist2 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Fdup vs mutar. No se mantuvo el fenotipo del duplicado.
                    if(dist2 != 0){cuenta3++;}
                }//Si regresa al fenotipo original
                else{//Si el fenotipo al mutar es diferente al original y al de duplicar
                    dist2 = find.distfenotipica(Matriz4, pd, Fenotipo2, per, tam2);//Fdup vs mutar. No se mantuvo el fenotipo del duplicado.
                    if(dist2 != 0){
                        //Cuando el fenotipo que resulta de mutar no se parece al original
                        cuenta2++;//Para la robustez a mutaciones in después de duplicar, contamos fenotipos que no se mantuvieron [Robin = (1 - cuenta)/vecinos]
                        if(cuenta2 == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                            Fnuevoin++;
                            Periodos2[0] = per; //Almacenamos el primer periodo
                            for(int i = 0; i < per; i++){
                                for(int j = 0; j < tam2; j++){
                                    Fenotipos2[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo
                                }
                            }
                        }
                        else{//Si ya tenemos un fenotipo para comparar
                            //                     if(m == 32 & veces == 124){cout << "n\n";}
                            Periodos2[cuenta2 - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                            hileras += Periodos2[cuenta2 - 2];
                            per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                            for(int i = hileras; i < (hileras + per); i++){
                                for(int j = 0; j < tam2; j++){
                                    Fenotipos2[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                                }
                            }
                            for(int k = 0; k < (cuenta2 - 1); k++){//Comparar todos los fenotipos uno a uno
                                per2 = Periodos2[k];//sacamos periodo de cada fenotipo a comparar
                                hor.espacio(Comparafenotipo, per2, tam2);
                                for(int i = per1; i < (per1 + per2); i++){
                                    for(int j = 0; j < tam2; j++){
                                        Comparafenotipo[i - per1][j] = Fenotipos2[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                                    }
                                }
                                distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                                hor.borrar(Comparafenotipo, per2);
                                if(distf == 0){//Si se encontró otro fenotipo igual en la lista de fenotipos
                                    valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
                                    break;
                                }//Break porque ya sabemos que es nuevo y ya no necesita seguir comparando
                                //                         else{hor.borrar(Comparafenotipo, per2);}
                                if(k == (cuenta2 - 2)){
                                    Fnuevoin++;
                                }
                                per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                            }
                        }
                    }
                }
                //Guardamos la distancia en la matriz de distancias
                Vector1[veces] = 1 - dist2;
            }
        }
        else{Vector1[veces] = 2; Vector2[veces] = 2;}
        
        hor.borrar(Fenotipo2, per);
        veces++;
        if(tipo == 1){
            //         dvecpasos.open ("Condin/outs/Pasos/dvecpasos.txt", ios::app);
            dvecpasos.open ("outs/Pasos/dvecpasos.txt", ios::app);
            dvecpasos << n -1 << " ";
            dvecpasos.close ();
        }
        if(tipo == 2){
            //         c1vecpasos.open ("Condin/outs/Pasos/c1vecpasos.txt", ios::app);
            c1vecpasos.open ("outs/Pasos/c1vecpasos.txt", ios::app);
            c1vecpasos << n -1 << " ";
            c1vecpasos.close ();
        }
        if(tipo == 3){
            //         c2vecpasos.open ("Condin/outs/Pasos/c2vecpasos.txt", ios::app);
            c2vecpasos.open ("outs/Pasos/c2vecpasos.txt", ios::app);
            c2vecpasos << n -1 << " ";
            c2vecpasos.close ();
        }
    }while(veces != num2);
    
    if(tipo == 1){
        //         dvecpasos.open ("Condin/outs/Pasos/dvecpasos.txt", ios::app);
        dvecpasos.open ("outs/Pasos/dvecpasos.txt", ios::app);
        dvecpasos << "\n";
        dvecpasos.close ();
    }
    if(tipo == 2){
        //         c1vecpasos.open ("Condin/outs/Pasos/c1vecpasos.txt", ios::app);
        c1vecpasos.open ("outs/Pasos/c1vecpasos.txt", ios::app);
        c1vecpasos << "\n";
        c1vecpasos.close ();
    }
    if(tipo == 3){
        //         c2vecpasos.open ("Condin/outs/Pasos/cnum12vecpasos.txt", ios::app);
        c2vecpasos.open ("outs/Pasos/cnum12vecpasos.txt", ios::app);
        c2vecpasos << "\n";
        c2vecpasos.close ();
    }
    
    hor.borrar(Fenotipos1, 100000);
    hor.borrar(trayectoria, 1000);
    hor.borrar(Fenotipos2, 100000);
    hor.borrar(Periodos1);
    hor.borrar(Periodos2);
    hor.borrar(pormientras);
}

void GATACA::exploracionpuntofijo(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek, int &puntosfijos)
{
    //     Matriz1: Factores
    //     Matriz2: Marcadores
    //     Matriz3: Trayectoria
    //     Matriz4: Fenotipo original
    //     Matriz5: Distancias
    //     num: Máximo de steps
    //     tam1 = hils1 = cols1 = cols2 = Nf
    //     tam2 = hils2 = Np
    //     num1: Número de muestra para saber en qué hilera va la distancia
    //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
    
    int **Fenotipo2;
    int *vectorentrada;
    hor.espacio(vectorentrada, tam1);
    ///////////////////////////////////////////////////
    int **Fenotipos;
    hor.espacio(Fenotipos, 10000, tam2);
    int *Periodos;
    hor.espacio(Periodos, num2);
    ofstream vecpasos, phenotypes, nuevo, diferente/*, acces*/, index;
    //   int cuentafenotipos= 0;//Cuenta de número de fenotipos diferentes
    for(int i = 0; i < 10000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos[i][j] = 0;
        }
    }
    for(int i = 0; i < num2; i++){Periodos[i] = 0;}
    int per1, per2, n2, hileras = 0;
    Fenotiponuevo = 0;
    double distf, dist;
    int **Comparafenotipo;//Para tomar fenotipo por fenotipo de la matriz de fenotipos para comparar
    //   double *distf; //Distancia entre fenotipo pasado y acutal
    //   hor.espacio(distf, num2);
    ///////////////////////////////////////////////////
    int **trayectoria; //Trayectoria para cuando se muten Nf's
    hor.espacio(trayectoria, 1000, tam1);
    double valor = 0;
    int x = 0, y = 0, matmut = 0; //Nos dicen cuál matriz y en qué elemento fue mutada y su valor anterior para poder
    int per;
    int veces = 0;
    
    double distpuntofijo;
    bool res;
    int *pormientras;
    hor.espacio(pormientras, tam2);
    
    ofstream fo;
    do{
        res = false;
        choose(Matriz1, tam1, tam1, Matriz2, tam2, tam1, matmut);
        //Si choose fue para la matriz de factores de transcripción
        if(matmut == 1){
            //           cuentaup++; //Contamos que la mutación está arriba
            n2 = 1; //Tiene su propio tamaño de trayectoria
            mut(Matriz1, tam1, tam1, tolerancia, intnw, x, y, valor);
            //Buscamos atractor, puede salir de más de una hilera
            for(int i = 0; i < tam1; i++){
                vectorentrada[i] = Matriz3[0][i];
                trayectoria[0][i] = Matriz3[0][i];
            }//Que trabaje desde la condición incial, sino n = 1 no sirve
            find.atractor(Matriz1, trayectoria, vectorentrada, tam1, n2, per); //Tiene su propia trayectoria
            hor.espacio(Fenotipo2, per, tam2);
            //Obtenemos Fenotipo2, puede salir de más de una hilera
            lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, n2, per); //Saca su fenotipo vecino
            
            //Regresamos a la matriz muestra
            back(Matriz1, x, y, valor);
        }
        //Si choose fue para la matriz de marcadores
        if(matmut == 2){
            //           cuentadown++; //Contamos que la mutación está abajo
            mut(Matriz2, tam2, tam1, tolerancia, intnw, x, y, valor);
            //Obtenemos Fenotipo2
            per = p; //Ya tiene la trayectoria y atractor definidos
            hor.espacio(Fenotipo2, per, tam2);
            lin.dot(Matriz2, Matriz3, Fenotipo2, tam1, tam2, n, per); //Saca su fenotipo con la trayectoria del espécimen
            //           distparcial = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
            //           if(distparcial != 0){down++;}
            //Regresamos a la matriz muestra
            back(Matriz2, x, y, valor);
        }
        
        //////////////////////////////////PARA PUNTOS FIJOS//////////////////////////////
        if(per > 1){
            for(int a = 0; a < tam2; a++){
                pormientras[a] = Fenotipo2[0][a];
            }
            distpuntofijo = find.distfenotipica(pormientras, Fenotipo2, per, tam2);//Por si a = aaa
            if(distpuntofijo == 0){res = true;}
        }
        ////////////////////////////////////////////////////////////////////////////////
        if((per == 1)||(res == true)){
            puntosfijos++;
            dist = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
            //       cout << "dist" << dist << endl;
            if(dist != 0){//Nos aseguramos de no contar el fenotipo original
                cuenta++;  
                if(cuenta == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                    Fenotiponuevo++;
                    Periodos[0] = per; //Almacenamos el primer periodo
                    for(int i = 0; i < per; i++){
                        for(int j = 0; j < tam2; j++){
                            Fenotipos[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo
                        }
                    }
                }
                else{//Si ya tenemos un fenotipo para comparar
                    Periodos[cuenta - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                    hileras += Periodos[cuenta - 2];
                    per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                    for(int i = hileras; i < (hileras + per); i++){
                        for(int j = 0; j < tam2; j++){
                            Fenotipos[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                        }
                    }
                    for(int k = 0; k < cuenta; k++){//Comparar todos los fenotipos uno a uno
                        per2 = Periodos[k];//sacamos periodo de cada fenotipo a comparar
                        hor.espacio(Comparafenotipo, per2, tam2);
                        for(int i = per1; i < (per1 + per2); i++){
                            for(int j = 0; j < tam2; j++){
                                Comparafenotipo[i - per1][j] = Fenotipos[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                            }
                        }
                        distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                        hor.borrar(Comparafenotipo, per2);
                        if(distf == 0){
                            valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
                            break;
                        }//Break porque ya sabemos que no es nuevo y ya no necesita seguir comparando
                        //                 else{hor.borrar(Comparafenotipo, per2);}
                        if(k == (cuenta - 2)){
                            Fenotiponuevo++;
                        }
                        per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                    }
                }
            }
            
            //Guardamos la distancia en la matriz de distancias
            Matriz5[veces] = 1 - dist;
        }
        else{Matriz5[veces] = 2;}
        hor.borrar(Fenotipo2, per);
        
        veces++;
        
        vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
        vecpasos << n -1 << " ";
        vecpasos.close ();
    }while(veces != num2);
    
    vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
    vecpasos << "\n";
    vecpasos.close ();
    
    hor.borrar(Fenotipos, 10000);
    hor.borrar(trayectoria, 1000);
    hor.borrar(Periodos);
    hor.borrar(vectorentrada);
    hor.borrar(pormientras);
}
void GATACA::exploracion(double **Matriz1, double **Matriz2, int **Matriz3, int **Matriz4, double *Matriz5, int tam1, int tam2, int tolerancia, double intnw, int n, int p, int num1, int num2, int &cuenta, int &Fenotiponuevo, int *valordek, int **Fn_notol, int per_notol, int &cuenta_mut_dup_notol)
{
  //     Matriz1: Factores
  //     Matriz2: Marcadores
  //     Matriz3: Trayectoria
  //     Matriz4: Fenotipo original
  //     Matriz5: Distancias
  //     num: Máximo de steps
  //     tam1 = hils1 = cols1 = cols2 = Nf
  //     tam2 = hils2 = Np
  //     num1: Número de muestra para saber en qué hilera va la distancia
  //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
    
    int **Fenotipo2;
    int *vectorentrada;
    hor.espacio(vectorentrada, tam1);
    ///////////////////////////////////////////////////
    int **Fenotipos;
    hor.espacio(Fenotipos, 10000, tam2);
    int *Periodos;
    hor.espacio(Periodos, num2);
    ofstream vecpasos, phenotypes, nuevo, diferente/*, acces*/, index;
    //   int cuentafenotipos= 0;//Cuenta de número de fenotipos diferentes
    for(int i = 0; i < 10000; i++){
        for(int j = 0; j < tam2; j++){
            Fenotipos[i][j] = 0;
      }
  }
  for(int i = 0; i < num2; i++){Periodos[i] = 0;}
  int per1, per2, n2, hileras = 0;
  Fenotiponuevo = 0;
  double distf, dist, dist_notol;
  int **Comparafenotipo;//Para tomar fenotipo por fenotipo de la matriz de fenotipos para comparar
//   double *distf; //Distancia entre fenotipo pasado y acutal
//   hor.espacio(distf, num2);
  ///////////////////////////////////////////////////
  int **trayectoria; //Trayectoria para cuando se muten Nf's
  hor.espacio(trayectoria, 1000, tam1);
  double valor = 0;
  int x = 0, y = 0, matmut = 0; //Nos dicen cuál matriz y en qué elemento fue mutada y su valor anterior para poder
  int per;
  int veces = 0;
  
  //////////////////////Rmutparciales/////////////////////////////////

  ofstream fo;
  do{
      choose(Matriz1, tam1, tam1, Matriz2, tam2, tam1, matmut);
      //Si choose fue para la matriz de factores de transcripción
      if(matmut == 1){
//           cuentaup++; //Contamos que la mutación está arriba
          n2 = 1; //Tiene su propio tamaño de trayectoria
          mut(Matriz1, tam1, tam1, tolerancia, intnw, x, y, valor);
          //Buscamos atractor, puede salir de más de una hilera
          for(int i = 0; i < tam1; i++){
              vectorentrada[i] = Matriz3[0][i];
              trayectoria[0][i] = Matriz3[0][i];
          }//Que trabaje desde la condición incial, sino n = 1 no sirve
          find.atractor(Matriz1, trayectoria, vectorentrada, tam1, n2, per); //Tiene su propia trayectoria
          hor.espacio(Fenotipo2, per, tam2);
          //Obtenemos Fenotipo2, puede salir de más de una hilera
          lin.dot(Matriz2, trayectoria, Fenotipo2, tam1, tam2, n2, per); //Saca su fenotipo vecino
          
          //Regresamos a la matriz muestra
          back(Matriz1, x, y, valor);
      }
      //Si choose fue para la matriz de marcadores
      if(matmut == 2){
//           cuentadown++; //Contamos que la mutación está abajo
          mut(Matriz2, tam2, tam1, tolerancia, intnw, x, y, valor);
          //Obtenemos Fenotipo2
          per = p; //Ya tiene la trayectoria y atractor definidos
          hor.espacio(Fenotipo2, per, tam2);
          lin.dot(Matriz2, Matriz3, Fenotipo2, tam1, tam2, n, per); //Saca su fenotipo con la trayectoria del espécimen
//           distparcial = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
//           if(distparcial != 0){down++;}
          //Regresamos a la matriz muestra
          back(Matriz2, x, y, valor);
      }
      dist_notol = find.distfenotipica(Fn_notol, per_notol, Fenotipo2, per, tam2);
      if(dist_notol == 0){//Quitamos fenotipos nuevos de no tolerantes
          cuenta_mut_dup_notol++;//Número de mutaciones que dan fenotipo no tolerante
          Matriz5[veces] = 2;
      }
      else{
          dist = find.distfenotipica(Matriz4, p, Fenotipo2, per, tam2);
          //       cout << "dist" << dist << endl;
          if(dist != 0){//Nos aseguramos de no contar el fenotipo original
              cuenta++;
              //###############Guardamos fenotipos###############270420
//              hor.open_ofstream(fo, "outs/acces/fn/fn_" + hor.inttostring(num1) + ".txt");
//              for(int j = 0; j < per; j++){
//                  for(int k = 0; k < tam2; k++){
//                      fo << Fenotipo2[j][k] << " ";
//                  }fo << endl;
//              }fo << endl;
//              fo.close();
//              hor.open_ofstream(fo, "outs/acces/fn/per_" + hor.inttostring(num1) + ".txt");
//              fo << per << endl;
//              fo.close();
              //#################################################
              //           cout << "Fdif = " << cuenta << endl;
              if(cuenta == 1){//Si es el primer fenotipo, sólo lo contamos, no lo comparamos.
                  Fenotiponuevo++;
                  //               cout << "Fnuevo1 = " << Fenotiponuevo << endl;
                  valordek[0] = 1;
                  //               cout << "1 valordek[0]= " << valordek[0] << endl;
                  Periodos[0] = per; //Almacenamos el primer periodo
                  for(int i = 0; i < per; i++){
                      for(int j = 0; j < tam2; j++){
                          Fenotipos[i][j] = Fenotipo2[i][j]; //Almacenamos primer fenotipo diferente
                      }
                  }
              }
              else{//Si ya tenemos un fenotipo para comparar
                  Periodos[cuenta - 1] = per; //Almacenamos los periodos de los fenotipos siguientes
                  hileras += Periodos[cuenta - 2];
                  per1 = 0; //Para que inicie en la hilera 0, per1 lleva la suma de periodos para hubicar el fenotipo a comparar
                  for(int i = hileras; i < (hileras + per); i++){
                      for(int j = 0; j < tam2; j++){
                          Fenotipos[i][j] = Fenotipo2[i - hileras][j]; //Almacenamos fenotipo
                      }
                  }
                  for(int k = 0; k < (cuenta - 1); k++){//Comparar todos los fenotipos uno a uno
                      per2 = Periodos[k];//sacamos periodo de cada fenotipo a comparar
                      hor.espacio(Comparafenotipo, per2, tam2);
                      for(int i = per1; i < (per1 + per2); i++){
                          for(int j = 0; j < tam2; j++){
                              Comparafenotipo[i - per1][j] = Fenotipos[i][j]; //Sacamos al fenotipo a comparar del saco de fenotipos
                          }
                      }
                      distf = find.distfenotipica(Comparafenotipo, per2, Fenotipo2, per, tam2); //Almacenamos primer fenotipo
                      hor.borrar(Comparafenotipo, per2);
                      if(distf == 0){
                          valordek[k] = valordek[k] + 1; //Suma 1 al fenotipo k que se repite
                          //                       cout << "valordek[" << k << "]= " << valordek[k] << endl;
                          break;
                      }//Break porque ya sabemos que no es nuevo y ya no necesita seguir comparando
                      //                 else{hor.borrar(Comparafenotipo, per2);}
                      if(k == (cuenta - 2)){
                          Fenotiponuevo++;
                          //                       cout << "Fnuevo = " << Fenotiponuevo << endl;
                          valordek[k + 1] = 1; //Suma 1 al fenotipo k que se repite //051019
                          //                       cout << "valordeknuevo[" << k + 1 << "]= " << valordek[k + 1] << endl;
                      }
                      per1 += per2; //Actualizamos la hilera en la que vamos conforme avanzamos en los fenotipos
                  }
              }
          }
          
          //Guardamos la distancia en la matriz de distancias
          Matriz5[veces] = 1 - dist;
      }
    hor.borrar(Fenotipo2, per);
    
    veces++;
    
    vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
    vecpasos << n -1 << " ";
    vecpasos.close ();
  }while(veces != num2);
  
  vecpasos.open ("outs/Pasos/vecpasos.txt", ios::app);
  vecpasos << "\n";
  vecpasos.close ();

  hor.borrar(Fenotipos, 10000);
  hor.borrar(trayectoria, 1000);
  hor.borrar(Periodos);
  hor.borrar(vectorentrada);
}

// void GATACA::muta_una_celda(int veces, double Matriz, int num_gen, int evento){
//     double coin;
//     int trayectoria;
//     hor.espacio(trayectoria, 1000, Nf);
//     for(int num_mut = 0; num_mut < veces; num_mut++){
//         cout << "num_mut: " << num_mut << endl;
//         coin = hor.randgauss(std);
//         if(valor == 0){ Matrizf[num_gen][evento] = coin;}
//         else{
//             if(coin < 0){ Matrizf[num_gen][evento]] = 0;}
//             else{
//                 if(valor < 0){ Matrizf[num_gen][evento] = abs(hor.randgauss(std)) * (-1);}
//                 else{ Matriz[num_gen][evento] = abs(hor.randgauss(std));}
//             }
//         }
//         //Revisamos si cambia de fenotipo
//         cout << "FTs mut" << endl;
//         hor.pantalla(Matrizf, Nf, Nf);
//         p = 0;
//         n = 1;
//         find.atractor(Matrizf, trayectoria, vectorentrada, Nf, n, p);
//         hor.espacio(Fenotipo, p, Np);
//         lin.dot(Matrizp, trayectoria, Fenotipo, Nf, Np, n, p);
//         dist = find.distfenotipica(Wish, Fenotipo, p, Np);
//         hor.borrar(Fenotipo, p);
//         sum_sim_antes += 1 - dist;
//         if(dist != 0){
//             cuenta_dif_antes++;
//         }
//         //Regresamos el valor original
//         Matriz[num_gen][evento] = valor;
//         cout << "dist: " << dist << " sum_sim_antes: " << sum_sim_antes << " cuenta_dif_antes: " << cuenta_dif_antes << endl;
//     }
//     hor.borrar(trayectoria, 1000);
// }


void GATACA::interaction_lost(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial){
//     cout << "1.1" << endl;
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    //Matrizf
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][evento] != 0){
//             cout << "1.2" << endl;
            contador1++;
            Matrizf[num_gen][evento] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum1 += 1 - dist;
            if(dist != 0){
                count1++;
            }
        }
        //Regresamos el valor original
        Matrizf[num_gen][evento] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[evento][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][num_gen] != 0){
//             cout << "1.3" << endl;
            contador2++;
            Matrizf[evento][num_gen] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum2 += 1 - dist;
            if(dist != 0){
                count2++;
            }
        }
        //Regresamos el valor original
        Matrizf[evento][num_gen] = valor;
    }
//     cout << "1.4" << endl;
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
//         cout << "1.5" << endl;
        //Interacciones de entrada
        valor = Matrizp[num_gen][evento];
//         cout << "valor: " << valor << endl;
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][evento] != 0){
//             cout << "1.5.1" << endl;
            contador3++;
            Matrizp[num_gen][evento] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum3 += 1 - dist;
            if(dist != 0){
                count3++;
            }
        }
//         cout << "1.5.2" << endl;
        //Regresamos el valor original
        Matrizp[num_gen][evento] = valor;
    }
//     cout << "1.7" << endl;
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}


void GATACA::interaction_gain(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla){
//      cout << "1.1" << endl;
     Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    //Matrizf
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][evento] == 0){
//                          cout << "1.2" << endl;
            contador1++;
            for(int i = 0; i < 5; i++){
                Matrizf[num_gen][evento] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
            for(int i = 0; i < 5; i++){
                Matrizf[num_gen][evento] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[num_gen][evento] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[evento][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][num_gen] == 0){
//             cout << "1.3" << endl;
            contador2++;
            for(int i = 0; i < 5; i++){
//                 cout << "1.3.1" << endl;
//                 cout << "hor.randgauss(1): " << endl;
                Matrizf[evento][num_gen] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
//                 cout << "Matrizf[evento][num_gen]: " << Matrizf[evento][num_gen] << endl;
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
//                 cout << "1.3.2" << endl;
            }
//             cout << "1.3.3" << endl;
            for(int i = 0; i < 5; i++){
                Matrizf[evento][num_gen] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[evento][num_gen] = valor;
    }
//     if(contador1 != 0){sum1 = sum1 / (10 * contador1);}
//     if(contador2 != 0){sum2 = sum2 / (10 * contador2);}
//      cout << "1.4" << endl;
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
//          cout << "1.5" << endl;
        //Interacciones de entrada
        valor = Matrizp[num_gen][evento];
//         cout << "valor: " << valor << endl;
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][evento] == 0){
            //             cout << "1.5.1" << endl;
            contador3++;
            for(int i = 0; i < 5; i++){
                Matrizp[num_gen][evento] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
            for(int i = 0; i < 5; i++){
                Matrizp[num_gen][evento] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
        }
        //         cout << "1.5.2" << endl;
        //Regresamos el valor original
        Matrizp[num_gen][evento] = valor;
    }
//     if(contador3 != 0){sum3 = sum3 / (10 * contador3);}
//     cout << "1.7" << endl;
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}


void GATACA::interaction_change(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla){
//      cout << "1.1" << endl;
     Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    //Matrizf
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][evento] != 0){
            contador1++;
            for(int i = 0; i < 10; i++){
                Matrizf[num_gen][evento] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[num_gen][evento] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[evento][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][num_gen] != 0){
//             cout << "1.3" << endl;
            contador2++;
            for(int i = 0; i < 10; i++){
                Matrizf[evento][num_gen] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[evento][num_gen] = valor;
    }
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
        //Interacciones de entrada
        valor = Matrizp[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][evento] != 0){
            contador3++;
            for(int i = 0; i < 10; i++){
                Matrizp[num_gen][evento] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
        }
        //Regresamos el valor original
        Matrizp[num_gen][evento] = valor;
    }
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}

void GATACA::interaction_lost2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla){
//Eliminar interacciones en un singleton de la red
    double valor, dist;
    int n, p, singleton;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    //Matrizf
    Horse hor(semilla);
    
    do{singleton = hor.randflat(0, hilsf);}
    while(singleton == evento);
//     cout << "evento: " << evento << " singleton: " << singleton << endl;
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][singleton];
//         cout << "valor: " << valor << endl;
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][singleton] != 0){
//             cout << "Entro a eliminar" << endl;
            contador1++;
            Matrizf[num_gen][singleton] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum1 += 1 - dist;
            if(dist != 0){
                count1++;
            }
//             cout << "dist: " << dist << endl;
        }
        //Regresamos el valor original
        Matrizf[num_gen][singleton] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[singleton][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[singleton][num_gen] != 0){
//             cout << "1.3" << endl;
            contador2++;
            Matrizf[singleton][num_gen] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum2 += 1 - dist;
            if(dist != 0){
                count2++;
            }
        }
        //Regresamos el valor original
        Matrizf[singleton][num_gen] = valor;
    }
//     cout << "1.4" << endl;
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
//         cout << "1.5" << endl;
        //Interacciones de entrada
        valor = Matrizp[num_gen][singleton];
//         cout << "valor: " << valor << endl;
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][singleton] != 0){
//             cout << "1.5.1" << endl;
            contador3++;
            Matrizp[num_gen][singleton] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum3 += 1 - dist;
            if(dist != 0){
                count3++;
            }
        }
//         cout << "1.5.2" << endl;
        //Regresamos el valor original
        Matrizp[num_gen][singleton] = valor;
    }
//     cout << "1.7" << endl;
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}

void GATACA::interaction_gain2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla){
//      cout << "1.1" << endl;
     Horse hor(semilla);
    double valor, dist;
    int n, p, singleton;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    
    do{singleton = hor.randflat(0, hilsf);}
    while(singleton == evento);
    
    //Matrizf
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][singleton];
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][singleton] == 0){
            contador1++;
            for(int i = 0; i < 5; i++){
                Matrizf[num_gen][singleton] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
            for(int i = 0; i < 5; i++){
                Matrizf[num_gen][singleton] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[num_gen][singleton] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[singleton][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[singleton][num_gen] == 0){
            contador2++;
            for(int i = 0; i < 5; i++){
                Matrizf[singleton][num_gen] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
            }
            for(int i = 0; i < 5; i++){
                Matrizf[singleton][num_gen] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[singleton][num_gen] = valor;
    }
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
        //Interacciones de entrada
        valor = Matrizp[num_gen][singleton];
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][singleton] == 0){
            contador3++;
            for(int i = 0; i < 5; i++){
                Matrizp[num_gen][singleton] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
            for(int i = 0; i < 5; i++){
                Matrizp[num_gen][singleton] = abs(hor.randgauss(1)) * (-1);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
        }
        //Regresamos el valor original
        Matrizp[num_gen][singleton] = valor;
    }
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}

void GATACA::interaction_change2(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2, double &sum3, double &count3, double &contador3, int *Wish, int *condinicial, int semilla){
//      cout << "1.1" << endl;
    Horse hor(semilla);
    double valor, dist;
    int n, p, singleton;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    
    do{singleton = hor.randflat(0, hilsf);}
    while(singleton == evento);
    evento = singleton;
    //Matrizf
    for(int num_gen = 0; num_gen < hilsf; num_gen++){
        //Interacciones de entrada
        valor = Matrizf[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizf[num_gen][evento] != 0){
            contador1++;
            for(int i = 0; i < 10; i++){
                Matrizf[num_gen][evento] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum1 += 1 - dist;
                if(dist != 0){
                    count1++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[num_gen][evento] = valor;
        //////////////////////////////////////////////////////////////////
        //Interacciones de salida
        valor = Matrizf[evento][num_gen];
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][num_gen] != 0){
//             cout << "1.3" << endl;
            contador2++;
            for(int i = 0; i < 10; i++){
                Matrizf[evento][num_gen] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum2 += 1 - dist;
                if(dist != 0){
                    count2++;
                }
            }
        }
        //Regresamos el valor original
        Matrizf[evento][num_gen] = valor;
    }
    //Matrizp
    for(int num_gen = 0; num_gen < hilsp; num_gen++){
        //Interacciones de entrada
        valor = Matrizp[num_gen][evento];
        //Robustez ante perder una conección de entrada
        if(Matrizp[num_gen][evento] != 0){
            contador3++;
            for(int i = 0; i < 10; i++){
                Matrizp[num_gen][evento] = hor.randgauss(valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum3 += 1 - dist;
                if(dist != 0){
                    count3++;
                }
            }
        }
        //Regresamos el valor original
        Matrizp[num_gen][evento] = valor;
    }
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}

// void GATACA::interaction_lost_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial){
//     double valor, dist;
//     int n, p;
//     int **Fenotipo, **trayectoria, *vectorentrada;    
//     hor.espacio(trayectoria, 1000, hilsf);
//     hor.espacio(vectorentrada, hilsf);
//     double *vectorsum, *contador, *cuentacambio;
//     int *vechils;
//     hor.espacio(vectorsum, 2);//entrada 1, salida 2
//     hor.espacio(contador, 2);
//     hor.espacio(cuentacambio, 2);
//     hor.espacio(vechils,2);
//     vechils[0] = evento;
//     vechils[1] = 12;
//     hor.fillv0(vectorsum, 2);
//     hor.fillv0(contador, 2);
//     hor.fillv0(cuentacambio, 2);
//        
//     for(int num = 0; num < 2; num++){
//         //Entradas
//         valor = Matrizf[evento][vechils[num]];
//         //Robustez ante perder una conección de entrada
//         if(Matrizf[evento][vechils[num]] != 0){
//             contador[0]++;
//             Matrizf[evento][vechils[num]] = 0;
//             //Revisamos si cambia de fenotipo
//             for(int i = 0; i < hilsf; i++){
//                 vectorentrada[i] = condinicial[i];
//                 trayectoria[0][i] = condinicial[i];
//             }
//             p = 0;
//             n = 1;
//             find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//             hor.espacio(Fenotipo, p, hilsp);
//             lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//             dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//             hor.borrar(Fenotipo, p);
//             vectorsum[0] += 1 - dist;
//             if(dist != 0){
//                 cuentacambio[0]++;
//             }
//         }
//         //Regresamos el valor original
//         Matrizf[evento][vechils[num]] = valor;
//         //Salidas
//         valor = Matrizf[vechils[num]][evento];
//         //Robustez ante perder una conección de entrada
//         if(Matrizf[vechils[num]][evento] != 0){
//             contador[1]++;
//             Matrizf[vechils[num]][evento] = 0;
//             //Revisamos si cambia de fenotipo
//             for(int i = 0; i < hilsf; i++){
//                 vectorentrada[i] = condinicial[i];
//                 trayectoria[0][i] = condinicial[i];
//             }
//             p = 0;
//             n = 1;
//             find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//             hor.espacio(Fenotipo, p, hilsp);
//             lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//             dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//             hor.borrar(Fenotipo, p);
//             vectorsum[1] += 1 - dist;
//             if(dist != 0){
//                 cuentacambio[1]++;
//             }
//         }
//         //Regresamos el valor original
//         Matrizf[vechils[num]][evento] = valor;
//     }  
//     sum1 = vectorsum[0];
//     count1 = contador[0];
//     contador1 = cuentacambio[0];
//     sum2 = vectorsum[1];
//     count2 = contador[1];
//     contador2 = cuentacambio[1];
//       
//     hor.borrar(trayectoria, 1000);
//     hor.borrar(vectorentrada);
//     hor.borrar(vectorsum);//entrada, salida
//     hor.borrar(contador);
//     hor.borrar(cuentacambio);
//     hor.borrar(vechils);
// }

void GATACA::interaction_lost_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int *vechils;
    hor.espacio(vechils,2);
    vechils[0] = evento;
    vechils[1] = 12;
       
    for(int num = 0; num < 2; num++){
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][vechils[num]] != 0){
//             cout << "Dup>Dup [" << evento << "][" << vechils[num] << "]: " << Matrizf[evento][vechils[num]] << endl;
            valor = Matrizf[evento][vechils[num]];
            contador++;
            Matrizf[evento][vechils[num]] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
            //Regresamos el valor original
            Matrizf[evento][vechils[num]] = valor;
        }
        //Robustez ante perder una conección de entrada
        if(Matrizf[vechils[1]][vechils[num]] != 0){
//             cout << "Dup>Dup [" << vechils[1] << "][" << vechils[num] << "]: " << Matrizf[vechils[1]][vechils[num]] << endl;
            valor = Matrizf[vechils[1]][vechils[num]];
            contador++;
            Matrizf[vechils[1]][vechils[num]] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
            //Regresamos el valor original
            Matrizf[vechils[1]][vechils[num]] = valor;
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vechils);
}

// void GATACA::interaction_gain_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
//     Horse hor(semilla);
//     double valor, dist;
//     int n, p;
//     int **Fenotipo, **trayectoria, *vectorentrada;    
//     hor.espacio(trayectoria, 1000, hilsf);
//     hor.espacio(vectorentrada, hilsf);
//     double *vectorsum, *contador, *cuentacambio;
//     int *vechils;
//     hor.espacio(vectorsum, 2);//entrada, salida
//     hor.espacio(contador, 2);
//     hor.espacio(cuentacambio, 2);
//     hor.espacio(vechils, 2);
//     vechils[0] = evento;
//     vechils[1] = 12;
//     hor.fillv0(vectorsum, 2);
//     hor.fillv0(contador, 2);
//     hor.fillv0(cuentacambio, 2);
//     
//     for(int veces = 0; veces < 5; veces++){
//         for(int num = 0; num < 2; num++){
//             //Entradas
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[evento][vechils[num]] == 0){
//                 valor = Matrizf[evento][vechils[num]];
//                 contador[0]++;
//                 Matrizf[evento][vechils[num]] = abs(hor.randgauss(0));
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[0] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[0]++;
//                 }
//                 //Regresamos el valor original
//                 Matrizf[evento][vechils[num]] = valor;
//             }
//             //Salidas
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[vechils[num]][evento] == 0){
//                 valor = Matrizf[vechils[num]][evento];
//                 contador[1]++;
//                 Matrizf[vechils[num]][evento] = abs(hor.randgauss(0));
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[1] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[1]++;
//                 }
//                 //Regresamos el valor original
//                 Matrizf[vechils[num]][evento] = valor;
//             }
//         }
//     }
//     
//     for(int veces = 0; veces < 5; veces++){
//         for(int num = 0; num < 2; num++){
//             //Entradas
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[evento][vechils[num]] == 0){
//                 valor = Matrizf[evento][vechils[num]];
//                 contador[0]++;
//                 Matrizf[evento][vechils[num]] = abs(hor.randgauss(0)) * -1;
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[0] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[0]++;
//                 }
//                 //Regresamos el valor original
//                 Matrizf[evento][vechils[num]] = valor;
//             }
//             //Salidas
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[vechils[num]][evento] == 0){
//                 valor = Matrizf[vechils[num]][evento];
//                 contador[1]++;
//                 Matrizf[vechils[num]][evento] = abs(hor.randgauss(0)) * -1;
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[1] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[1]++;
//                 }
//                 //Regresamos el valor original
//                 Matrizf[vechils[num]][evento] = valor;
//             }
//         }
//     }
//     
//     sum1 = vectorsum[0];
//     count1 = contador[0];
//     contador1 = cuentacambio[0];
//     sum2 = vectorsum[1];
//     count2 = contador[1];
//     contador2 = cuentacambio[1];
//     
//     hor.borrar(trayectoria, 1000);
//     hor.borrar(vectorentrada);
//     hor.borrar(vectorsum);//entrada, salida
//     hor.borrar(contador);
//     hor.borrar(cuentacambio);
//     hor.borrar(vechils);
//     
//     hor.close_rng();
// }

void GATACA::interaction_gain_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int *vechils;
    hor.espacio(vechils, 2);
    vechils[0] = evento;
    vechils[1] = 12;
    
    for(int veces = 0; veces < 5; veces++){
        for(int num = 0; num < 2; num++){
            //Robustez ante perder una conección de entrada
            if(Matrizf[evento][vechils[num]] == 0){
//                 cout << "Dup>Dup [" << evento << "][" << vechils[num] << "]: " << Matrizf[evento][vechils[num]] << endl;
                valor = Matrizf[evento][vechils[num]];
                contador++;
                Matrizf[evento][vechils[num]] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[evento][vechils[num]] = valor;
            }
            //Robustez ante perder una conección de entrada
            if(Matrizf[vechils[1]][vechils[num]] == 0){
//                 cout << "Dup>Dup [" << vechils[1] << "][" << vechils[num] << "]: " << Matrizf[vechils[1]][vechils[num]] << endl;
                valor = Matrizf[vechils[1]][vechils[num]];
                contador++;
                Matrizf[vechils[1]][vechils[num]] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[vechils[1]][vechils[num]] = valor;
            }
        }
    }
    
    for(int veces = 0; veces < 5; veces++){
        for(int num = 0; num < 2; num++){
            //Robustez ante perder una conección de entrada
            if(Matrizf[evento][vechils[num]] == 0){
//                 cout << "Dup>Dup [" << evento << "][" << vechils[num] << "]: " << Matrizf[evento][vechils[num]] << endl;
                valor = Matrizf[evento][vechils[num]];
                contador++;
                Matrizf[evento][vechils[num]] = abs(hor.randgauss(1)) * -1;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[evento][vechils[num]] = valor;
            }
            //Robustez ante perder una conección de entrada
            if(Matrizf[vechils[1]][vechils[num]] == 0){
//                 cout << "Dup>Dup [" << vechils[1] << "][" << vechils[num] << "]: " << Matrizf[vechils[1]][vechils[num]] << endl;
                valor = Matrizf[vechils[1]][vechils[num]];
                contador++;
                Matrizf[vechils[1]][vechils[num]] = abs(hor.randgauss(1)) * -1;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[vechils[1]][vechils[num]] = valor;
            }
        }
    }
    
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vechils);
    
    hor.close_rng();
}

// void GATACA::interaction_change_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
//     Horse hor(semilla);
//     double valor, dist;
//     int n, p;
//     int **Fenotipo, **trayectoria, *vectorentrada;    
//     hor.espacio(trayectoria, 1000, hilsf);
//     hor.espacio(vectorentrada, hilsf);
//     double *vectorsum, *contador, *cuentacambio;
//     int *vechils;
//     hor.espacio(vectorsum, 2);//entrada, salida
//     hor.espacio(contador, 2);
//     hor.espacio(cuentacambio, 2);
//     hor.espacio(vechils, 2);
//     vechils[0] = evento;
//     vechils[1] = 12;
//     hor.fillv0(vectorsum, 2);
//     hor.fillv0(contador, 2);
//     hor.fillv0(cuentacambio, 2);
//     
//     for(int veces = 0; veces < 10; veces++){
//         for(int num = 0; num < 2; num++){
//             //Entradas
//             valor = Matrizf[evento][vechils[num]];
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[evento][vechils[num]] != 0){
//                 contador[0]++;
//                 Matrizf[evento][vechils[num]] = abs(hor.randgauss(valor));
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[0] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[0]++;
//                 }
//             }
//             //Regresamos el valor original
//             Matrizf[evento][vechils[num]] = valor;
//             //Salidas
//             valor = Matrizf[vechils[num]][evento];
//             //Robustez ante perder una conección de entrada
//             if(Matrizf[vechils[num]][evento] != 0){
//                 contador[1]++;
//                 Matrizf[vechils[num]][evento] = abs(hor.randgauss(valor));
//                 //Revisamos si cambia de fenotipo
//                 for(int i = 0; i < hilsf; i++){
//                     vectorentrada[i] = condinicial[i];
//                     trayectoria[0][i] = condinicial[i];
//                 }
//                 p = 0;
//                 n = 1;
//                 find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                 hor.espacio(Fenotipo, p, hilsp);
//                 lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                 dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                 hor.borrar(Fenotipo, p);
//                 vectorsum[1] += 1 - dist;
//                 if(dist != 0){
//                     cuentacambio[1]++;
//                 }
//             }
//             //Regresamos el valor original
//             Matrizf[vechils[num]][evento] = valor;
//         }
//     }
//     
//     sum1 = vectorsum[0];
//     count1 = contador[0];
//     contador1 = cuentacambio[0];
//     sum2 = vectorsum[1];
//     count2 = contador[1];
//     contador2 = cuentacambio[1];
//     
//     hor.borrar(trayectoria, 1000);
//     hor.borrar(vectorentrada);
//     hor.borrar(vectorsum);//entrada, salida
//     hor.borrar(contador);
//     hor.borrar(cuentacambio);
//     hor.borrar(vechils);
//     
//     hor.close_rng();
// }

void GATACA::interaction_change_dup_dup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int *vechils;
    hor.espacio(vechils, 2);
    vechils[0] = evento;
    vechils[1] = 12;
    
    for(int veces = 0; veces < 10; veces++){
        for(int num = 0; num < 2; num++){
            //Entradas
            if(Matrizf[evento][vechils[num]] != 0){
//                 cout << "Dup>Dup [" << evento << "][" << vechils[num] << "]: " << Matrizf[evento][vechils[num]] << endl;
                valor = Matrizf[evento][vechils[num]];
                contador++;
                Matrizf[evento][vechils[num]] = abs(hor.randgauss(1) + valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[evento][vechils[num]] = valor;
            }
            //Robustez ante perder una conección de entrada
            if(Matrizf[vechils[1]][vechils[num]] != 0){
//                 cout << "Dup>Dup [" << vechils[1] << "][" << vechils[num] << "]: " << Matrizf[vechils[1]][vechils[num]] << endl;
                valor = Matrizf[vechils[1]][vechils[num]];
                contador++;
                Matrizf[vechils[1]][vechils[num]] = abs(hor.randgauss(1) + valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                //Regresamos el valor original
                Matrizf[vechils[1]][vechils[num]] = valor;
            }
        }
    }
        
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vechils);
    
    hor.close_rng();
}

void GATACA::interaction_lost_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
    
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            //Ndup a dup
            //Robustez ante perder una conección de entrada
            if(Matrizf[evento][num] != 0){
//                 cout << "Ndup>Dup [" << evento << "][" << num << "]: " << Matrizf[evento][num] << endl;
                valor = Matrizf[evento][num];
                contador[0]++;
                Matrizf[evento][num] = 0;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[0] += 1 - dist;
                if(dist != 0){
                    cuentacambio[0]++;
                }
                Matrizf[evento][num] = valor;
            }
            //Dup a Ndup
            //Robustez ante perder una conección de entrada
            if(Matrizf[num][evento] != 0){
//                 cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                valor = Matrizf[num][evento];
                contador[1]++;
                Matrizf[num][evento] = 0;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[1] += 1 - dist;
                if(dist != 0){
                    cuentacambio[1]++;
                }
                Matrizf[num][evento] = valor;
            }
        }  
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
}

void GATACA::interaction_gain_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
    
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            for(int veces = 0; veces < 5; veces++){
                //Ndup a dup
                //Robustez ante perder una conección de entrada
                if(Matrizf[evento][num] == 0){
                    valor = Matrizf[evento][num];
                    contador[0]++;
                    Matrizf[evento][num] = abs(hor.randgauss(1));
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[0] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[0]++;
                    }
                    Matrizf[evento][num] = valor;
                }
                //Dup a Ndup
                //Robustez ante perder una conección de entrada
                if(Matrizf[num][evento] == 0){
//                     cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                    valor = Matrizf[num][evento];
                    contador[1]++;
                    Matrizf[num][evento] = abs(hor.randgauss(1));
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[1] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[1]++;
                    }
                    Matrizf[num][evento] = valor;
                }
            } 
//             cout << "8\n";
            for(int veces = 0; veces < 5; veces++){
                //Ndupt a dup
                //Robustez ante perder una conección de entrada
//                 cout << "9\n";
                if(Matrizf[evento][num] == 0){
//                     cout << "Ndup>Dup [" << evento << "][" << num << "]: " << Matrizf[evento][num] << endl;
                    valor = Matrizf[evento][num];
                    contador[0]++;
                    Matrizf[evento][num] = abs(hor.randgauss(1)) * -1;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[0] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[0]++;
                    }
//                     cout << "ND_D Matrizf[" << evento << "][" << num << "]= " << Matrizf[evento][num] << " dist= " << dist << endl;
                    Matrizf[evento][num] = valor;
                }
                //Dup a Ndup
                //Robustez ante perder una conección de entrada
                if(Matrizf[num][evento] == 0){
//                     cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                    valor = Matrizf[num][evento];
                    contador[1]++;
                    Matrizf[num][evento] = abs(hor.randgauss(1)) * -1;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[1] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[1]++;
                    }
                    Matrizf[num][evento] = valor;
                }
            }
        }
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
    
    hor.close_rng();
}

void GATACA::interaction_change_dup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
       
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            for(int veces = 0; veces < 10; veces++){
                //Ndup a dup
                //Robustez ante perder una conección de entrada
                if(Matrizf[evento][num] != 0){
//                     cout << "Ndup>Dup [" << evento << "][" << num << "]: " << Matrizf[evento][num] << endl;
                    valor = Matrizf[evento][num];
                    contador[0]++;
                    Matrizf[evento][num] = abs(hor.randgauss(1) + valor);
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[0] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[0]++;
                    }
                    Matrizf[evento][num] = valor;
                }
                //Dup a Ndup
                //Robustez ante perder una conección de entrada
                if(Matrizf[num][evento] != 0){
//                     cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                    valor = Matrizf[num][evento];
                    contador[1]++;
                    Matrizf[num][evento] = abs(hor.randgauss(1) + valor);
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    vectorsum[1] += 1 - dist;
                    if(dist != 0){
                        cuentacambio[1]++;
                    }
                    Matrizf[num][evento] = valor;
                }
            } 
        }
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
    
    hor.close_rng();
}

void GATACA::interaction_lost_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
           
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            for(int num2 = 0; num2 < (hilsf - 1); num2++){
                if(num2 != evento){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizf[num][num2] != 0){
//                         cout << "Ndup>Ndup [" << num << "][" << num2 << "]: " << Matrizf[num][num2] << endl;
                        valor = Matrizf[num][num2];
                        contador++;
                        Matrizf[num][num2] = 0;
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizf[num][num2] = valor;
                    }
                    //Salidas
                    //Robustez ante perder una conección de entrada
//                     if(Matrizf[num2][num] != 0){
//                         cout << "Ndup>Ndup [" << num2 << "][" << num << "]: " << Matrizf[num2][num] << endl;
//                         valor = Matrizf[num2][num];
//                         contador++;
//                         Matrizf[num2][num] = 0;
//                         //Revisamos si cambia de fenotipo
//                         for(int i = 0; i < hilsf; i++){
//                             vectorentrada[i] = condinicial[i];
//                             trayectoria[0][i] = condinicial[i];
//                         }
//                         p = 0;
//                         n = 1;
//                         find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                         hor.espacio(Fenotipo, p, hilsp);
//                         lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                         dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                         hor.borrar(Fenotipo, p);
//                         sum += 1 - dist;
//                         if(dist != 0){
//                             cambio++;
//                         }
//                         Matrizf[num2][num] = valor;
//                     }
                }  
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

void GATACA::interaction_gain_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            for(int num2 = 0; num2 < (hilsf - 1); num2++){
                if(num2 != evento){
                    for(int veces = 0; veces < 5; veces++){
                        //Entradas
                        //Robustez ante perder una conección de entrada
                        if(Matrizf[num][num2] == 0){
//                             cout << "Ndup>Ndup [" << num << "][" << num2 << "]: " << Matrizf[num][num2] << endl;
                            valor = Matrizf[num][num2];
                            contador++;
                            Matrizf[num][num2] = abs(hor.randgauss(1));
                            //Revisamos si cambia de fenotipo
                            for(int i = 0; i < hilsf; i++){
                                vectorentrada[i] = condinicial[i];
                                trayectoria[0][i] = condinicial[i];
                            }
                            p = 0;
                            n = 1;
                            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                            hor.espacio(Fenotipo, p, hilsp);
                            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                            hor.borrar(Fenotipo, p);
                            sum += 1 - dist;
                            if(dist != 0){
                                cambio++;
                            }
                            Matrizf[num][num2] = valor;
                        }
                        //Salidas
                        //Robustez ante perder una conección de entrada
//                         if(Matrizf[num2][num] == 0){
//                             cout << "Ndup>Ndup [" << num2 << "][" << num << "]: " << Matrizf[num2][num] << endl;
//                             valor = Matrizf[num2][num];
//                             contador++;
//                             Matrizf[num2][num] = abs(hor.randgauss(0));
//                             //Revisamos si cambia de fenotipo
//                             for(int i = 0; i < hilsf; i++){
//                                 vectorentrada[i] = condinicial[i];
//                                 trayectoria[0][i] = condinicial[i];
//                             }
//                             p = 0;
//                             n = 1;
//                             find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                             hor.espacio(Fenotipo, p, hilsp);
//                             lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                             dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                             hor.borrar(Fenotipo, p);
//                             sum += 1 - dist;
//                             if(dist != 0){
//                                 cambio++;
//                             }
//                             Matrizf[num2][num] = valor;
//                         }
                    } 
                    for(int veces = 0; veces < 5; veces++){
                        //Entradas
                        //Robustez ante perder una conección de entrada
                        if(Matrizf[num][num2] == 0){
//                             cout << "Ndup>Ndup [" << num << "][" << num2 << "]: " << Matrizf[num][num2] << endl;
                            valor = Matrizf[num][num2];
                            contador++;
                            Matrizf[num][num2] = abs(hor.randgauss(1)) * -1;
                            //Revisamos si cambia de fenotipo
                            for(int i = 0; i < hilsf; i++){
                                vectorentrada[i] = condinicial[i];
                                trayectoria[0][i] = condinicial[i];
                            }
                            p = 0;
                            n = 1;
                            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                            hor.espacio(Fenotipo, p, hilsp);
                            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                            hor.borrar(Fenotipo, p);
                            sum += 1 - dist;
                            if(dist != 0){
                                cambio++;
                            }
                            Matrizf[num][num2] = valor;
                        }
                        //Salidas
                        //Robustez ante perder una conección de entrada
//                         if(Matrizf[num2][num] == 0){
//                             cout << "Ndup>Ndup [" << num2 << "][" << num << "]: " << Matrizf[num2][num] << endl;
//                             valor = Matrizf[num2][num];
//                             contador++;
//                             Matrizf[num2][num] = abs(hor.randgauss(0)) * -1;
//                             //Revisamos si cambia de fenotipo
//                             for(int i = 0; i < hilsf; i++){
//                                 vectorentrada[i] = condinicial[i];
//                                 trayectoria[0][i] = condinicial[i];
//                             }
//                             p = 0;
//                             n = 1;
//                             find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                             hor.espacio(Fenotipo, p, hilsp);
//                             lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                             dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                             hor.borrar(Fenotipo, p);
//                             sum += 1 - dist;
//                             if(dist != 0){
//                                 cambio++;
//                             }
//                             Matrizf[num2][num] = valor;
//                         }
                    }
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_change_Ndup_Ndup(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < (hilsf - 1); num++){
        if(num != evento){
            for(int num2 = 0; num2 < (hilsf - 1); num2++){
                if(num2 != evento){
                    for(int veces = 0; veces < 10; veces++){
                        //Entradas
                        //Robustez ante perder una conección de entrada
                        if(Matrizf[num][num2] != 0){
//                             cout << "Ndup>Ndup [" << num << "][" << num2 << "]: " << Matrizf[num][num2] << endl;
                            valor = Matrizf[num][num2];
                            contador++;
                            Matrizf[num][num2] = abs(hor.randgauss(1) + valor);
                            //Revisamos si cambia de fenotipo
                            for(int i = 0; i < hilsf; i++){
                                vectorentrada[i] = condinicial[i];
                                trayectoria[0][i] = condinicial[i];
                            }
                            p = 0;
                            n = 1;
                            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                            hor.espacio(Fenotipo, p, hilsp);
                            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                            hor.borrar(Fenotipo, p);
                            sum += 1 - dist;
                            if(dist != 0){
                                cambio++;
                            }
                            Matrizf[num][num2] = valor;
                        }
                        //Salidas
                        //Robustez ante perder una conección de entrada
//                         if(Matrizf[num2][num] == 0){
//                             cout << "Ndup>Ndup [" << num2 << "][" << num << "]: " << Matrizf[num2][num] << endl;
//                             valor = Matrizf[num2][num];
//                             contador++;
//                             Matrizf[num2][num] = abs(hor.randgauss(valor));
//                             //Revisamos si cambia de fenotipo
//                             for(int i = 0; i < hilsf; i++){
//                                 vectorentrada[i] = condinicial[i];
//                                 trayectoria[0][i] = condinicial[i];
//                             }
//                             p = 0;
//                             n = 1;
//                             find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
//                             hor.espacio(Fenotipo, p, hilsp);
//                             lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
//                             dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
//                             hor.borrar(Fenotipo, p);
//                             sum += 1 - dist;
//                             if(dist != 0){
//                                 cambio++;
//                             }
//                             Matrizf[num2][num] = valor;
//                         }
                    } 
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_lost_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
           
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizp[num][num2] != 0){
//                     cout << "Ndup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                    valor = Matrizp[num][num2];
                    contador++;
                    Matrizp[num][num2] = 0;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizp[num][num2] = valor;
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

void GATACA::interaction_gain_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                for(int veces = 0; veces < 5; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizp[num][num2] == 0){
//                         cout << "Ndup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                        valor = Matrizp[num][num2];
                        contador++;
                        Matrizp[num][num2] = abs(hor.randgauss(1));
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizp[num][num2] = valor;
                    }
                } 
                for(int veces = 0; veces < 5; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizp[num][num2] == 0){
//                         cout << "Ndup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                        valor = Matrizp[num][num2];
                        contador++;
                        Matrizp[num][num2] = abs(hor.randgauss(1)) * -1;
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizp[num][num2] = valor;
                    }
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_change_Ndup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                for(int veces = 0; veces < 10; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizp[num][num2] != 0){
//                         cout << "Ndup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                        valor = Matrizp[num][num2];
                        contador++;
                        Matrizp[num][num2] = abs(hor.randgauss(1) + valor);
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizp[num][num2] = valor;
                    }
                } 
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_lost_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int num2 = evento;
    
    for(int num = 0; num < hilsp; num++){
        //Entradas
        //Robustez ante perder una conección de entrada
        if(Matrizp[num][num2] != 0){
//             cout << "Dup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
            valor = Matrizp[num][num2];
            contador++;
            Matrizp[num][num2] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
            Matrizp[num][num2] = valor;
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

void GATACA::interaction_gain_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;    
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int num2 = evento;
    
    for(int num = 0; num < hilsp; num++){
        for(int veces = 0; veces < 5; veces++){
            //Entradas
            //Robustez ante perder una conección de entrada
            if(Matrizp[num][num2] == 0){
//                 cout << "Dup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                valor = Matrizp[num][num2];
                contador++;
                Matrizp[num][num2] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                Matrizp[num][num2] = valor;
            }
        } 
        for(int veces = 0; veces < 5; veces++){
            //Entradas
            //Robustez ante perder una conección de entrada
            if(Matrizp[num][num2] == 0){
//                 cout << "Dup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                valor = Matrizp[num][num2];
                contador++;
                Matrizp[num][num2] = abs(hor.randgauss(1)) * -1;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                Matrizp[num][num2] = valor;
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_change_Dup_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    int num2 = evento;
    
    for(int num = 0; num < hilsp; num++){
        for(int veces = 0; veces < 10; veces++){
            //Entradas
            //Robustez ante perder una conección de entrada
            if(Matrizp[num][num2] != 0){
                //                 cout << "Dup>ST [" << num << "][" << num2 << "]: " << Matrizp[num][num2] << endl;
                valor = Matrizp[num][num2];
                contador++;
                Matrizp[num][num2] = abs(hor.randgauss(1) + valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                sum += 1 - dist;
                if(dist != 0){
                    cambio++;
                }
                Matrizp[num][num2] = valor;
            }
        }
    }
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

//Además mide auto-interacciones NoControl 16Mar22
void GATACA::interaction_lost_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    evento = 12; //Es el gen control
    
    if(Matrizf[evento][evento] != 0){//Si existe autorregulación
            valor = Matrizf[evento][evento];//Guarda interacción
            contador++;
            Matrizf[evento][evento] = 0;//Elimina interacción
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
            Matrizf[evento][evento] = valor;//Regresa el valor original de la interacciión
        }
   
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

//Además mide auto-interacciones NoControl 16Mar22
void GATACA::interaction_gain_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    evento = 12; //Es el gen control
    
    if(Matrizf[evento][evento] == 0){//Si no existe autorregulación
        for(int veces = 0; veces < 5; veces++){
            contador++;
            Matrizf[evento][evento] = abs(hor.randgauss(1));//Crea interacción
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
        }
        for(int veces = 0; veces < 5; veces++){
            contador++;
            Matrizf[evento][evento] = abs(hor.randgauss(1)) * -1;//Crea interacción
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
        }
        Matrizf[evento][evento] = 0;//Regresa el valor original de la interacciión
    }
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}

//Además mide auto-interacciones NoControl 16Mar22
void GATACA::interaction_change_Control_Control(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    evento = 12; //Es el gen control
    
    if(Matrizf[evento][evento] != 0){//Si existe autorregulación
        valor = Matrizf[evento][evento];
        for(int veces = 0; veces < 20; veces++){
            contador++;
            Matrizf[evento][evento] = abs(hor.randgauss(1) + valor);//Cambia interacción
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            sum += 1 - dist;
            if(dist != 0){
                cambio++;
            }
        }
        Matrizf[evento][evento] = valor;//Regresa el valor original de la interacciión
    }
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.close_rng();
}


//Para los controles 121121
//Elimina una a una las interacciones entre el gen control y los demás FTs
void GATACA::interaction_lost_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
    
    for(int num = 0; num < (hilsf - 1); num++){
        //NoC a C
        //Robustez ante perder una conección de entrada
        if(Matrizf[evento][num] != 0){//Si existe interacción
            valor = Matrizf[evento][num];//Guarda interacción
            contador[0]++;
            Matrizf[evento][num] = 0;//Elimina interacción
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            vectorsum[0] += 1 - dist;
            if(dist != 0){
                cuentacambio[0]++;
            }
            Matrizf[evento][num] = valor;//Regresa el valor original de la interacciión
        }
        //C a NC
        //Robustez ante perder una conección de entrada
        if(Matrizf[num][evento] != 0){
            //                 cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
            valor = Matrizf[num][evento];
            contador[1]++;
            Matrizf[num][evento] = 0;
            //Revisamos si cambia de fenotipo
            for(int i = 0; i < hilsf; i++){
                vectorentrada[i] = condinicial[i];
                trayectoria[0][i] = condinicial[i];
            }
            p = 0;
            n = 1;
            find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
            hor.espacio(Fenotipo, p, hilsp);
            lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
            dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
            hor.borrar(Fenotipo, p);
            vectorsum[1] += 1 - dist;
            if(dist != 0){
                cuentacambio[1]++;
            }
            Matrizf[num][evento] = valor;
        }
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
}

void GATACA::interaction_gain_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
    
    for(int num = 0; num < (hilsf - 1); num++){
        for(int veces = 0; veces < 5; veces++){
            //NC a C
            //Robustez ante ganar una conección de entrada
            if(Matrizf[evento][num] == 0){//Si no hay interacción
                valor = Matrizf[evento][num];
                contador[0]++;
                Matrizf[evento][num] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[0] += 1 - dist;
                if(dist != 0){
                    cuentacambio[0]++;
                }
                Matrizf[evento][num] = valor;
            }
            //Dup a Ndup
            //Robustez ante perder una conección de entrada
            if(Matrizf[num][evento] == 0){
                //                     cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                valor = Matrizf[num][evento];
                contador[1]++;
                Matrizf[num][evento] = abs(hor.randgauss(1));
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[1] += 1 - dist;
                if(dist != 0){
                    cuentacambio[1]++;
                }
                Matrizf[num][evento] = valor;
            }
        }
        //             cout << "8\n";
        for(int veces = 0; veces < 5; veces++){
            //Ndupt a dup
            //Robustez ante ganar una conección de entrada
            if(Matrizf[evento][num] == 0){
                valor = Matrizf[evento][num];
                contador[0]++;
                Matrizf[evento][num] = abs(hor.randgauss(1)) * -1;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[0] += 1 - dist;
                if(dist != 0){
                    cuentacambio[0]++;
                }
                Matrizf[evento][num] = valor;
            }
            //Dup a Ndup
            //Robustez ante perder una conección de entrada
            if(Matrizf[num][evento] == 0){
                valor = Matrizf[num][evento];
                contador[1]++;
                Matrizf[num][evento] = abs(hor.randgauss(1)) * -1;
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[1] += 1 - dist;
                if(dist != 0){
                    cuentacambio[1]++;
                }
                Matrizf[num][evento] = valor;
            }
        }
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
    
    hor.close_rng();
}

void GATACA::interaction_change_Control_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum1, double &count1, double &contador1, double &sum2, double &count2, double &contador2,int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    double *vectorsum, *contador, *cuentacambio;
    hor.espacio(vectorsum, 2);//entrada 1, salida 2
    hor.espacio(contador, 2);
    hor.espacio(cuentacambio, 2);
    hor.fillv0(vectorsum, 2);
    hor.fillv0(contador, 2);
    hor.fillv0(cuentacambio, 2);
    
    for(int num = 0; num < (hilsf - 1); num++){
        for(int veces = 0; veces < 10; veces++){
            //NC a C
            //Robustez ante perder una conección de entrada
            if(Matrizf[evento][num] != 0){
                valor = Matrizf[evento][num];
                contador[0]++;
                Matrizf[evento][num] = abs(hor.randgauss(1) + valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[0] += 1 - dist;
                if(dist != 0){
                    cuentacambio[0]++;
                }
                Matrizf[evento][num] = valor;
            }
            //Dup a Ndup
            //Robustez ante perder una conección de entrada
            if(Matrizf[num][evento] != 0){
                //                     cout << "Dup>Ndup [" << num << "][" << evento << "]: " << Matrizf[num][evento] << endl;
                valor = Matrizf[num][evento];
                contador[1]++;
                Matrizf[num][evento] = abs(hor.randgauss(1) + valor);
                //Revisamos si cambia de fenotipo
                for(int i = 0; i < hilsf; i++){
                    vectorentrada[i] = condinicial[i];
                    trayectoria[0][i] = condinicial[i];
                }
                p = 0;
                n = 1;
                find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                hor.espacio(Fenotipo, p, hilsp);
                lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                hor.borrar(Fenotipo, p);
                vectorsum[1] += 1 - dist;
                if(dist != 0){
                    cuentacambio[1]++;
                }
                Matrizf[num][evento] = valor;
            }
        }
    }
    sum1 = vectorsum[0];
    count1 = contador[0];
    contador1 = cuentacambio[0];
    sum2 = vectorsum[1];
    count2 = contador[1];
    contador2 = cuentacambio[1];
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    hor.borrar(vectorsum);//entrada, salida
    hor.borrar(contador);
    hor.borrar(cuentacambio);
    
    hor.close_rng();
}

void GATACA::interaction_lost_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    
    for(int num = 0; num < (hilsf - 1); num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizf[num][num2] != 0){
                    valor = Matrizf[num][num2];
                    contador++;
                    Matrizf[num][num2] = 0;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizf[num][num2] = valor;
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

void GATACA::interaction_gain_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
    
    for(int num = 0; num < (hilsf - 1); num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                for(int veces = 0; veces < 5; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizf[num][num2] == 0){
                        valor = Matrizf[num][num2];
                        contador++;
                        Matrizf[num][num2] = abs(hor.randgauss(1));
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizf[num][num2] = valor;
                    }
                }
                for(int veces = 0; veces < 5; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizf[num][num2] == 0){
                        valor = Matrizf[num][num2];
                        contador++;
                        Matrizf[num][num2] = abs(hor.randgauss(1)) * -1;
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizf[num][num2] = valor;
                    }
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_change_NControl_NControl(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < (hilsf - 1); num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            if(num2 != evento){
                for(int veces = 0; veces < 10; veces++){
                    //Entradas
                    //Robustez ante perder una conección de entrada
                    if(Matrizf[num][num2] != 0){
                        valor = Matrizf[num][num2];
                        contador++;
                        Matrizf[num][num2] = abs(hor.randgauss(1) + valor);
                        //Revisamos si cambia de fenotipo
                        for(int i = 0; i < hilsf; i++){
                            vectorentrada[i] = condinicial[i];
                            trayectoria[0][i] = condinicial[i];
                        }
                        p = 0;
                        n = 1;
                        find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                        hor.espacio(Fenotipo, p, hilsp);
                        lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                        dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                        hor.borrar(Fenotipo, p);
                        sum += 1 - dist;
                        if(dist != 0){
                            cambio++;
                        }
                        Matrizf[num][num2] = valor;
                    }
                }
            }
        }
    }
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_lost_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial){
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
           
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizp[num][num2] != 0){
                    valor = Matrizp[num][num2];
                    contador++;
                    Matrizp[num][num2] = 0;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizp[num][num2] = valor;
                }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
}

void GATACA::interaction_gain_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            for(int veces = 0; veces < 5; veces++){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizp[num][num2] == 0){
                    valor = Matrizp[num][num2];
                    contador++;
                    Matrizp[num][num2] = abs(hor.randgauss(1));
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizp[num][num2] = valor;
                }
            }
            for(int veces = 0; veces < 5; veces++){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizp[num][num2] == 0){
                    valor = Matrizp[num][num2];
                    contador++;
                    Matrizp[num][num2] = abs(hor.randgauss(1)) * -1;
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizp[num][num2] = valor;
                }
            }
        }
    }
      
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}

void GATACA::interaction_change_NControl_ST(double **Matrizf, int hilsf, double **Matrizp, int hilsp, int evento, double &sum, double &contador, double &cambio, int *Wish, int *condinicial, int semilla){
    Horse hor(semilla);
    double valor, dist;
    int n, p;
    int **Fenotipo, **trayectoria, *vectorentrada;
    hor.espacio(trayectoria, 1000, hilsf);
    hor.espacio(vectorentrada, hilsf);
       
    for(int num = 0; num < hilsp; num++){
        for(int num2 = 0; num2 < (hilsf - 1); num2++){
            for(int veces = 0; veces < 10; veces++){
                //Entradas
                //Robustez ante perder una conección de entrada
                if(Matrizp[num][num2] != 0){
                    valor = Matrizp[num][num2];
                    contador++;
                    Matrizp[num][num2] = abs(hor.randgauss(1) + valor);
                    //Revisamos si cambia de fenotipo
                    for(int i = 0; i < hilsf; i++){
                        vectorentrada[i] = condinicial[i];
                        trayectoria[0][i] = condinicial[i];
                    }
                    p = 0;
                    n = 1;
                    find.atractor(Matrizf, trayectoria, vectorentrada, hilsf, n, p);
                    hor.espacio(Fenotipo, p, hilsp);
                    lin.dot(Matrizp, trayectoria, Fenotipo, hilsf, hilsp, n, p);
                    dist = find.distfenotipica(Wish, Fenotipo, p, hilsp);
                    hor.borrar(Fenotipo, p);
                    sum += 1 - dist;
                    if(dist != 0){
                        cambio++;
                    }
                    Matrizp[num][num2] = valor;
                }
            }
        }
    }
    
    hor.borrar(trayectoria, 1000);
    hor.borrar(vectorentrada);
    
    hor.close_rng();
}
