///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////FUNCIONES UTILIZADAS EN MAESTRÍA, PCI, 2016-2018/////////////////////////
//////////////////////////////////YURIDIA SELENE POSADAS GARCIA////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//g++ -Wall -c *.cc
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include "finder.h"
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

Finder::Finder()
{
}

void Finder::atractor(double **matriz1, int **trayectoria, int *vectorentrada, int tam, int& n, int& p)
{
    int *vectorsalida;
    hor.espacio(vectorsalida, tam);
    ofstream myfile2;

    do{

        lin.dot(matriz1, trayectoria, vectorentrada, vectorsalida, tam, n);//Genera siguiente paso en tray y vectorsal
        n = n + 1;
        p = atractores(vectorentrada, vectorsalida, trayectoria, tam, n);//Compara tray con vecsal, si no hay uno igual p = -1 y genera vectorentrada
    }
    while(p == -1);//El ciclo termina cuando se encuentra un atractor de periodo p
    delete [] vectorsalida;
}

void Finder::atractor2(double **matriz1, int **trayectoria, int *vectorentrada, int tam, int& n, int& p, int& pasosmax)
{
    int *vectorsalida;
    hor.espacio(vectorsalida, tam);
    ofstream myfile2;

    do{

        lin.dot(matriz1, trayectoria, vectorentrada, vectorsalida, tam, n);//Genera siguiente paso en tray y vectorsal
        n = n + 1;
        p = atractores(vectorentrada, vectorsalida, trayectoria, tam, n);//Compara tray con vecsal, si no hay uno igual p = -1 y genera vectorentrada
        if(n > 200){
            pasosmax = 1;
            break;
        }
    }
    while(p == -1);//El ciclo termina cuando se encuentra un atractor de periodo p
    delete [] vectorsalida;
}

int Finder::atractores(int *vectorentrada, int *vectorsalida, int **trayectoria, int tam, int n) //Compara vector salida (resultado del estado por la matriz) con los vectores de la trayectoria para ver si encuentra uno igual e indica el periodo del atractor.
{
    int k = 0;
    int res;
    int periodo = 0;
    
    for(int j = n - 2; j >= 0; j--){ //Prueba el vector n-1 de la trayectoria, al inicio n = 2, n = 0 es cond inic
        periodo = periodo + 1;
        for(int i = 0; i < tam; i++){
            if(trayectoria[j][i] != vectorsalida[i]){ //Prueba cada elemento del vector
                break;//Si encuentra elementos diferentes, se sale de este ciclo for
            }
                k = i + 1; //Si todos los elementos son iguales, K = tam
        }
        if(k == tam){
            res = periodo;
            break; //Si al salir del otro ciclo se encontró un vector igual, sale de este ciclo for
        }
        if((j == 0) & (k < tam)){ //Si no encontró vectores similares en la trayectoria, res falso y sale del ciclo for
            res = -1;
        }
    }
            
    for (int i = 0; i < tam; i++){
        vectorentrada[i] = vectorsalida[i]; //Cambiar de nombre el vectorsalida para que regrese al do-while como vectorentrante.
    }
    
    return res; //si el tamaño = 1, es punto fijo; de lo contrario es ciclo.
}

double Finder::distfenotipica(int *ancestro, int *fenotipo, int tam) //Si el ancestro y el fenotipo tienen una hilera
{
  int dif = 0,i;
  for (i=0; i<tam; i++)
    if (ancestro[i] != fenotipo[i])
      dif++; //Suma 1 cada vez que no coincide
//   ofstream myfile2;
//   myfile2.open ("infosystem.txt", ios::app);
//   myfile2 << "Distancia 1: " << double(dif)/double(tam);
//   myfile2 << "\n" << endl;
//   myfile2.close ();
  
  return double(dif)/double(tam);
}

double Finder::distfenotipica(int *ancestro, int **fenotipo, int hilsf, int tam) //Si el ancestro tiene una hilera pero el fenotipo es una matriz
{
  int i,j;
  double dif = 0;
  for (i=0; i<hilsf; i++)
    for (j=0; j < tam; j++)
      if (ancestro[j] != fenotipo[i][j])
        dif = dif + (1.0/double(hilsf)); //Compara el ancestro con cada hilera del fenotipo y si difieren suma 1 y saca el promedio

  return dif/double(tam);
}

double Finder::distfenotipica(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam) //Si el ancestro y el fenotipo son matrices
{ /*cout << "dist 1\n";*/
  int i,j,k;
  double dis;
  if (hilsa != hilsf) {
    int mulng, mulna;
    if (hilsa < hilsf) {
      if ((hilsf%hilsa)==0) {
        mulna = 1;
        mulng = hilsf/hilsa;
      }
      else {
        mulng = hilsf;
        mulna = hilsa;
      }
    }
    else {
      if ((hilsa%hilsf)==0) {
        mulng = 1;
        mulna = hilsa/hilsf;
      }
      else {
        mulng = hilsf;
        mulna = hilsa;
      }
    }
//     cout << "dist 2\n";
//     cout << " hilsa*mulng: " << hilsa*mulng << " tam: " << tam << endl;
    int **nugoal, **nuatr;
    nugoal = new int*[hilsa*mulng];
    for (i= 0; i<(hilsa*mulng); i++){
//         cout << "nugoal i: " << i << endl;
      nugoal[i] = new int[tam];
    }
//     cout << "dist 2.1\n";
//     cout << " hilsf*mulna: " << hilsf*mulna << " tam: " << tam << endl;
    nuatr = new int*[hilsf*mulna];
    for (i=0; i <(hilsf*mulna); i++){
//         cout << "nuatr i: " << i << endl;
      nuatr[i] = new int[tam];
    }
//     /*cout*/ << "dist 2.2\n";
    for (i=0; i < mulng; i++)
      for (j=0; j < hilsa; j++)
        for (k=0; k < tam; k++)
          nugoal[(hilsa*i)+j][k] = ancestro[j][k];
//     cout << "dist 2.3\n";
    for (i=0; i < mulna; i++)
      for (j=0; j < hilsf; j++)
        for (k=0; k < tam;k++){
          nuatr[(hilsf*i)+j][k] = fenotipo[j][k];}
//     cout << "dist 2.4\n";
    dis = dist_aux(nugoal, (hilsa*mulng), nuatr, (hilsf*mulna), tam);
    for (i= 0; i<(hilsa*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(hilsf*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
  }
  else
    dis = dist_aux(ancestro, hilsa, fenotipo, hilsf, tam);
//   cout << "dist 3" << endl;

  return dis;
}

double Finder::dist_aux(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam)
{
  if (hilsa != hilsf) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists,d;
  dists = new double[hilsa];
  hor.fillv0(dists, hilsa);
  for (i=0; i<hilsa; i++) {
    dists[i] =0;
    for (j=0; j<hilsa; j++) {
      for (k=0;k<tam;k++){
        if (ancestro[j][k] != fenotipo[(j+i)%hilsa][k]){
          dists[i] = dists[i] + (1.0/double(hilsa));
        }
      }
    }
  }
  
  d = hor.find_min(dists, hilsa);

  d = d/double(tam);
  delete [] dists;
  return d;
}

double Finder::dist_aux(int *vector1, int cols, int *vector2)
{
  int i,j;
  double *dists,d;
  dists = new double[cols];
  hor.fillv0(dists, cols);
  for (i=0; i < cols; i++) {
    dists[i] =0;
    for (j=0; j < cols; j++) {
        cout << "W" << i << ": " << vector1[i] << endl;
        cout << "F" << i << ": " << vector2[(j+i)%cols] << endl;
        if (vector1[i] != vector2[(j+i)%cols]){
          dists[i] = dists[i] + (1.0/double(cols));
          cout << vector1[i] << " != " << vector2[(j+i)%cols] << endl;
          cout << "dists" << i << ": " << dists[i] << endl;
        }
      }
    }
  
  d = hor.find_min(dists, cols);
  cout << "d1: " << d << endl;
  d = d/double(cols);
  cout << "d2: " << d << endl;
  delete [] dists;
  return d;
}

int Finder::numintnw(double **matriz, int hils, int cols)
{
    int num = 0;
    for(int i = 0; i < hils; i++){
        for(int j = 0; j < cols; j++){
            if(matriz[i][j] != 0){
                num++;
            }
        }
    }
    return num;
}

int Finder::negativos(double **matriz, int hils, int cols)
{
    int num = 0;
    for(int i = 0; i < hils; i++){
        for(int j = 0; j < cols; j++){
            if(matriz[i][j] < 0){
                num++;
            }
        }
    }
    return num;
}

int Finder::positivos(double **matriz, int hils, int cols)
{
    int num = 0;
    for(int i = 0; i < hils; i++){
        for(int j = 0; j < cols; j++){
            if(matriz[i][j] > 0){
                num++;
            }
        }
    }
    return num;
}

void Finder::media(double *Vector1, double *Vector2, int num1, int num2)
{
    //Vector: Distancia
    //Matriz: Media
    //num1: límite de j
    //num2: muestra en la que va
//     cout << "num1: " << num1 << endl;
    double sum = 0;
    int cuenta = 0;
          for(int j = 0; j < num1; j++){
              if(Vector1[j] != 2){sum += Vector1[j];}
              if(Vector1[j] == 2){cuenta++;}
          }        
          Vector2[num2] = sum / (num1 - cuenta);
}

void Finder::media(double *Vector1, double **Matriz1, int num1, int num2, int num3)
{
    //Vector: Distancia
    //Matriz: Media
    //num1: límite de j
    //num2 y num3: muestra en la que va
    double sum = 0;
    int cuenta = 0;
    for(int j = 0; j < num1; j++){
        if(Vector1[j] != 2){sum += Vector1[j];}
        if(Vector1[j] == 2){cuenta++;}
    }
    Matriz1[num2][num3] = sum / (num1 - cuenta);
}

void Finder::desv(double *Vector1, double *Vector2, double *Vector3, int num1, int num2)
{
    //Vector: Distancia
    //Vector2: Media
    //Vector3: Desvst
    //num1: Cantidad de valores, límite de i
    //num2: muestra en la que va
    double Var; //Varianza de las distancias por muestra
    
    int cuenta = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] == 2) cuenta++;
    }
    
    double sum = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] != 2){
            sum += (Vector1[i] - Vector2[num2]) * (Vector1[i] - Vector2[num2]);
        }
    }
          
          
          Var = sum / (cuenta - 1);
          Vector3[num2] = sqrt(Var);
}

void Finder::desv(double *Vector1, double **Matriz1, double **Matriz2, int num1, int num2, int num3)
{
    //Vector: Distancia
    //Matriz: Media
    //Matriz2: Dsv
    //num1: Cantidad de valores, límite de i
    //num2 y 3: muestra en la que va
    double Var; //Varianza de las distancias por duplicado por muestra
    
    int cuenta = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] == 2) cuenta++;
    }
    
    double sum = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] != 2){
            sum += (Vector1[i] - Matriz1[num2][num3]) * (Vector1[i] - Matriz1[num2][num3]);
        }
    }
    
    Var = sum / (cuenta - 1);
    Matriz2[num2][num3] = sqrt(Var);
}

double Finder::compara1(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2)
{
//Número de mutaciones necesarias para llegar de una red a otra -Ojo: no necesariamente el número de mutaciones que sucedieron.
//Se compara a la red(0) con las otras 1,000 y se toma el valor más grande.
//Matriz1: factores1
//Matriz2: marcadores1
//Matriz3: factores2 (después de caminata)
//Matriz4: marcadores2 (después de caminata)
//num1: Nf
//num2: Np
    int sum = 0;
    double prop, valor;
    int novac1 = numintnw(Matriz1, num1, num1);
    int novac2 = numintnw(Matriz2, num2, num1);
    int novac3 = numintnw(Matriz3, num1, num1);
    int novac4 = numintnw(Matriz4, num2, num1);
    int novac = novac1 + novac2 + novac3 + novac4;
    
    
    for(int i = 0; i < num1; i++){
        for(int j = 0; j < num1; j++){
            if(Matriz1[i][j] != Matriz3[i][j]){
                valor = Matriz1[i][j] * Matriz3[i][j]; //Si son de diferente signo es igual a 2 mutaciones posibles
                if(valor < 0){sum += 2;}
                if(valor == 0){
                    if(Matriz1[i][j] != 0 || Matriz3[i][j] != 0){sum++;} //Si pasó de algo a 0 o de 0 a algo es 1 mutación
                }
            }
        }
    }
        
        for(int i = 0; i < num2; i++){
            for(int j = 0; j < num1; j++){
                if(Matriz2[i][j] != Matriz4[i][j]){
                    valor = Matriz2[i][j] * Matriz4[i][j];
                    if(valor < 0){sum += 2;}
                    if(valor == 0){
                        if(Matriz2[i][j] != 0 || Matriz4[i][j] != 0){sum++;}
                    }
                }
            }
        }

            prop = (double)sum / novac; //Total de mutaciones posibles entre el número de regulaciones (celdas diferentes de cero) totales de ambas redes.
            
return prop;
}

double Finder::compara2(double **Matriz1, double **Matriz2, double **Matriz3, double **Matriz4, int num1, int num2)
{
    //Diferencia en pesos en la red entre el número de regulaciones (celdas diferentes de cero) de ambas redes.
//Matriz1: factores1
//Matriz2: marcadores1
//Matriz3: factores2
//Matriz4: marcadores2
//num1: Nf
//num2: Np

    double valor = 0;
    double prop;
    int novac1 = numintnw(Matriz1, num1, num1);
    int novac2 = numintnw(Matriz2, num2, num1);
    int novac3 = numintnw(Matriz3, num1, num1);
    int novac4 = numintnw(Matriz4, num2, num1);
    int novac = novac1 + novac2 + novac3 + novac4;

for(int i = 0; i < num1; i++)
for(int j = 0; j < num1; j++){
valor += fabs(Matriz1[i][j] - Matriz3[i][j]);
}

for(int i = 0; i < num2; i++)
for(int j = 0; j < num1; j++){
valor += fabs(Matriz2[i][j] - Matriz4[i][j]);
}

prop = valor / novac;

return prop;
}

void Finder::minmaxmed(double **Matriz1, double **Matriz2, int hils, int cols)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Matriz1: entrada
    //Matriz2: salida
    double sum;
    
    for(int i = 0; i < hils; i++){
        sum = 0;
        Matriz2[i][0] = Matriz1[i][0];
        Matriz2[i][1] = Matriz1[i][0];
        for(int j = 0; j < cols; j++){
              sum += Matriz1[i][j];//Para el promedio
            if(Matriz1[i][j] < Matriz2[i][0]){Matriz2[i][0] = Matriz1[i][j];}
            if(Matriz1[i][j] > Matriz2[i][1]){Matriz2[i][1] = Matriz1[i][j];}
        }
        Matriz2[i][2] = sum / cols;
    }
}

void Finder::minmaxmedin(double **Matriz1, double **Matriz2, int hils, int cols)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Matriz1: entrada
    //Matriz2: salida
    double sum;
    int cuenta;
    
    for(int i = 0; i < hils; i++){
        sum = 0;
        cuenta = 0;
        for(int j = 0; j < cols; j++){
            if(Matriz1[i][j] != 2){
                sum += Matriz1[i][j];
                cuenta++;
            }//Para el promedio
        }
        Matriz2[i][2] = sum / cuenta;//sacamos promedio
        if(cuenta == 0){Matriz2[i][2] = -1;}
    }
    
    int k;
    for(int i = 0; i < hils; i++){
        //Encontramos el primer número que no sea 2
        for(k = 0; k < cols; k++){
            if(Matriz1[i][k] != 2){
                Matriz2[i][0] = Matriz1[i][k];//Mínimo inicial
                Matriz2[i][1] = Matriz1[i][k];//Máximo inicial
                break;
            }
        }
        if((k + 1) < cols)
            for(int j = k + 1; j < cols; j++)
                if(Matriz1[i][j] != 2){
                    if(Matriz1[i][j] < Matriz2[i][0]){Matriz2[i][0] = Matriz1[i][j];}
                    if(Matriz1[i][j] > Matriz2[i][1]){Matriz2[i][1] = Matriz1[i][j];}
                }
    }
}

void Finder::dist_btw_esp(int per_antes, int *Per_antes, int Ns, int **Fn_antes, string arch)
{
    ofstream fo;
    int sum_esp1, sum_esp2, **nuevoantes, **nuevodespues;
    double dist_fn_esp/*, sumdist_fn_esp*/;
    
    sum_esp1 = 0;
    dist_fn_esp = 0;
    for(int compara = 0; compara < (per_antes - 1); compara++){
        hor.espacio(nuevoantes, Per_antes[compara], Ns);
        for(int j = 0; j < Per_antes[compara]; j++){
            for(int k = 0; k < Ns; k++){
                nuevoantes[j][k] = Fn_antes[j + sum_esp1][k];
            }
        }
        
        sum_esp2 = 0;
        for(int n = 0; n < (compara + 1); n++){sum_esp2 += Per_antes[n];}
        for(int compara2 = 0; compara2 < (per_antes - (compara + 1)); compara2++){ //cada vez compara con uno menos para evitar repetición
            hor.espacio(nuevodespues, Per_antes[compara + compara2 + 1], Ns);
            for(int l = 0; l < Per_antes[compara + compara2 + 1]; l++){ //Fn a comparar
//                 hor.espacio(nuevodespues, Per_antes[compara + compara2 + 1], Ns);
                for(int m = 0; m < Ns; m++){
                    nuevodespues[l][m] = Fn_antes[l + sum_esp2][m];
                }
            }
            dist_fn_esp += distfenotipica(nuevoantes, Per_antes[compara], nuevodespues, Per_antes[compara + compara2 + 1], Ns);

            hor.borrar(nuevodespues, Per_antes[compara + compara2 + 1]);
            sum_esp2 += Per_antes[compara + compara2 + 1];
        }
        sum_esp1 += Per_antes[compara];
        hor.borrar(nuevoantes, Per_antes[compara]);
    }
    hor.open_ofstream(fo, ""+ arch +".txt");
    fo << dist_fn_esp / ((per_antes * (per_antes - 1))/2) << endl;
    fo.close();
}

void Finder::dist_btw_esp2(int per_antes, int *Per_antes, int Ns, int **Fn_antes, string arch1, string arch2, int gen)
{
    ofstream fo;
    int sum_esp1, sum_esp2, **nuevoantes, **nuevodespues;
    double dist_fn_esp/*, sumdist_fn_esp*/;
    
    sum_esp1 = 0;
    dist_fn_esp = 0;
    for(int compara = 0; compara < (per_antes - 1); compara++){
        hor.espacio(nuevoantes, Per_antes[compara], Ns);
        for(int j = 0; j < Per_antes[compara]; j++){
            for(int k = 0; k < Ns; k++){
                nuevoantes[j][k] = Fn_antes[j + sum_esp1][k];
            }
        }
        
        sum_esp2 = 0;
        for(int n = 0; n < (compara + 1); n++){sum_esp2 += Per_antes[n];}
        for(int compara2 = 0; compara2 < (per_antes - (compara + 1)); compara2++){ //cada vez compara con uno menos para evitar repetición
            hor.espacio(nuevodespues, Per_antes[compara + compara2 + 1], Ns);
            for(int l = 0; l < Per_antes[compara + compara2 + 1]; l++){ //Fn a comparar
//                 hor.espacio(nuevodespues, Per_antes[compara + compara2 + 1], Ns);
                for(int m = 0; m < Ns; m++){
                    nuevodespues[l][m] = Fn_antes[l + sum_esp2][m];
                }
            }
            dist_fn_esp += distfenotipica(nuevoantes, Per_antes[compara], nuevodespues, Per_antes[compara + compara2 + 1], Ns);

            hor.borrar(nuevodespues, Per_antes[compara + compara2 + 1]);
            sum_esp2 += Per_antes[compara + compara2 + 1];
        }
        sum_esp1 += Per_antes[compara];
        hor.borrar(nuevoantes, Per_antes[compara]);
    }
    hor.open_ofstream(fo, ""+ arch1 +"_"+ arch2 +"_"+ hor.inttostring(gen) +".txt");
    fo << dist_fn_esp / ((per_antes * (per_antes - 1))/2) << endl;
    fo.close();
}

void Finder::dist_ant_vs_desp(int Nf, int Ns, int cuenta1, int **acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes)
{
    //arch1: /home/rdaneel/Documentos/Doctorado/Resultados/outs/acces/
    //arch2: 1, 2 ó 3 dependiendo si es dup, c1 ó c2
    //arch3: dup, c1 ó c2
    //i: número de espécimen
//     cout << "Entra a dist_ant_vs_desp muestra: " << i << endl;
    int cuenta2, per_despues, sum1, sum2, fn_igual, *Per_despues, **nuevoantes, **nuevodespues, **Fn_despues;
    double dist;
    ifstream fi;
    ofstream fo;
    string basura;
    
    for(int gen = 0; gen < Nf; gen++){
        if(acces1[i][gen] != 0){
            cuenta2 = 0;
            hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
            while(fi >> basura)
                cuenta2++;
            fi.close();
            hor.espacio(Fn_despues, cuenta2/Ns, Ns);
            hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
            for(int j = 0; j < cuenta2/Ns; j++){
                for(int k = 0; k < Ns; k++){
                    fi >> Fn_despues[j][k];
                }
            }
            fi.close();
            hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
            
            per_despues = 0;
            while(fi >> basura)
                per_despues++;
            fi.close();
            //                 per_despues = per_despues - 1;
            hor.espacio(Per_despues, per_despues);
            hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
            for(int j = 0; j < per_despues; j++){
                fi >> Per_despues[j];
            }
            fi.close();
//             cout << "Entra a dist_btw_esp2 muestra: " << i << " gen: " << gen << endl;
            dist_btw_esp2(per_despues, Per_despues, Ns, Fn_despues, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/acces/dist/dist_fn", arch3, gen);
//             dist_btw_esp2(per_despues, Per_despues, Ns, Fn_despues, "/home/yuri/Documents/Doctorado/Resultados/outs/acces/dist/dist_fn", arch3, gen);
//             cout << "Sale de dist_btw_esp2 muestra: " << i << " gen: " << gen << endl;
//             cout << "out dist_btw_esp2: "<< endl;
            fn_igual = 0;
            sum1 = 0;
            for(int m = 0; m < per_antes; m++){ //Para todos los fns del espécimen i en curso
                hor.espacio(nuevoantes, Per_antes[m], Ns);
                for(int j = 0; j < Per_antes[m]; j++){ //Llenar vector de fn antes
                    for(int k = 0; k < Ns; k++){
                        nuevoantes[j][k] = Fn_antes[j + sum1][k];
                    }
                }
                sum2 = 0;
                
                for(int n = 0; n < per_despues; n++){ //Llenar vector de fn después
                    hor.espacio(nuevodespues, Per_despues[n], Ns);
                    for(int o = 0; o < Per_despues[n]; o++){
                        for(int q = 0; q < Ns; q++){
                            nuevodespues[o][q] = Fn_despues[o + sum2][q];
                        }
                    }
                    sum2 += Per_despues[n];
                    dist = distfenotipica(nuevoantes, Per_antes[m], nuevodespues, Per_despues[n], Ns);
//                     cout << "out distfenotipica: "<< endl;
                    if(dist == 0){
                        fn_igual++;
                        break;
                    }
                    hor.borrar(nuevodespues, Per_despues[n]);
                }
                sum1 += Per_antes[m];
                hor.borrar(nuevoantes, Per_antes[m]);
            }
            hor.borrar(Fn_despues, cuenta2/Ns);
            hor.borrar(Per_despues);
            hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
            fo << fn_igual << " ";
            fo.close();
            if(gen == (Nf - 1)){
                hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
                fo << endl;
                fo.close();
            }
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
            fo << "NA ";
            fo.close();
            if(gen == (Nf - 1)){
                hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
                fo << endl;
                fo.close();
            }
        }
    }
//     cout << "Sale de dist_ant_vs_desp muestra: " << i << endl;
}


void Finder::dist_ant_vs_desp_tol(int Nf, int Ns, int cuenta1, int acces1, int i, string arch1, string arch2, string arch3, int **Fn_antes, int per_antes, int *Per_antes, int gen)
{
    //arch1: /home/rdaneel/Documentos/Doctorado/Resultados/outs/acces/
    //arch2: 1, 2 ó 3 dependiendo si es dup, c1 ó c2
    //arch3: dup, c1 ó c2
    //i: número de espécimen
    //     cout << "Entra a dist_ant_vs_desp muestra: " << i << endl;
    int cuenta2, per_despues, sum1, sum2, fn_igual, *Per_despues, **nuevoantes, **nuevodespues, **Fn_despues;
    double dist;
    ifstream fi;
    ofstream fo;
    string basura;
    
    if(acces1 != 0){
        cuenta2 = 0;
        hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
        while(fi >> basura)
            cuenta2++;
        fi.close();
        hor.espacio(Fn_despues, cuenta2/Ns, Ns);
        hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
        for(int j = 0; j < cuenta2/Ns; j++){
            for(int k = 0; k < Ns; k++){
                fi >> Fn_despues[j][k];
            }
        }
        fi.close();
        hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
        
        per_despues = 0;
        while(fi >> basura)
            per_despues++;
        fi.close();
        //                 per_despues = per_despues - 1;
        hor.espacio(Per_despues, per_despues);
        hor.open_ifstream(fi, ""+ arch1 +"fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(gen) + ".txt");//020122 para rehacer dist_fenotipos_recurrentes261221.cc
        for(int j = 0; j < per_despues; j++){
            fi >> Per_despues[j];
        }
        fi.close();
        //             cout << "Entra a dist_btw_esp2 muestra: " << i << " gen: " << gen << endl;
        dist_btw_esp2(per_despues, Per_despues, Ns, Fn_despues, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/acces/dist/dist_fn", arch3, gen);
        
        fn_igual = 0;
        sum1 = 0;
        for(int m = 0; m < per_antes; m++){ //Para todos los fns del espécimen i en curso
            hor.espacio(nuevoantes, Per_antes[m], Ns);
            for(int j = 0; j < Per_antes[m]; j++){ //Llenar vector de fn antes
                for(int k = 0; k < Ns; k++){
                    nuevoantes[j][k] = Fn_antes[j + sum1][k];
                }
            }
            sum2 = 0;
            
            for(int n = 0; n < per_despues; n++){ //Llenar vector de fn después
                hor.espacio(nuevodespues, Per_despues[n], Ns);
                for(int o = 0; o < Per_despues[n]; o++){
                    for(int q = 0; q < Ns; q++){
                        nuevodespues[o][q] = Fn_despues[o + sum2][q];
                    }
                }
                sum2 += Per_despues[n];
                dist = distfenotipica(nuevoantes, Per_antes[m], nuevodespues, Per_despues[n], Ns);
                //                     cout << "out distfenotipica: "<< endl;
                if(dist == 0){
                    fn_igual++;
                    break;
                }
                hor.borrar(nuevodespues, Per_despues[n]);
            }
            sum1 += Per_antes[m];
            hor.borrar(nuevoantes, Per_antes[m]);
        }
        hor.borrar(Fn_despues, cuenta2/Ns);
        hor.borrar(Per_despues);
        hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
        fo << i << " " << gen << " " << fn_igual << " " << endl;
        fo.close();
//        if(gen == (Nf - 1)){
//            hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
//            fo << endl;
//            fo.close();
//        }
    }
    else{
        hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
        fo << "NA ";
        fo.close();
        if(gen == (Nf - 1)){
            hor.open_ofstream(fo, ""+ arch1 +"dist/dist_fn_esp_"+ arch3 +".txt");
            fo << endl;
            fo.close();
        }
    }
    //     cout << "Sale de dist_ant_vs_desp muestra: " << i << endl;
}
