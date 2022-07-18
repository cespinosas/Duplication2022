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
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
#include <iostream>
#include "horse.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::string;

Horse::Horse()
{
}

Horse::Horse(int semilla)
{
    start_rng(semilla);
}

void Horse::start_rng(int semilla)
{
    r = gsl_rng_alloc(gsl_rng_default); // Media en 0
    gsl_rng_set(r, semilla);
}

void Horse::espacio(double** &Matriz, int hils, int cols)
{
    Matriz = new double*[hils];
    for (int i = 0; i < hils; i++)
        Matriz[i] = new double[cols];
    return;
}

void Horse::borrar(double** &Matriz, int hils)
{
    for (int i = 0; i < hils; i++)
            delete [] Matriz[i];
        delete [] Matriz;
    return;
}

void Horse::borrar(int** &Matriz, int hils)
{
    for (int i = 0; i < hils; i++)
            delete [] Matriz[i];
        delete [] Matriz;
    return;
}


void Horse::borrar(double* &Vector)
{
    delete [] Vector;
    return;
}

void Horse::borrar(int* &Vector)
{
    delete [] Vector;
    return;
}

void Horse::borrar(bool* &Vector)
{
    delete [] Vector;
    return;
}


void Horse::espacio(int **&Matriz, int hils, int cols)
{
    Matriz = new int*[hils];
    for (int i = 0; i < hils; i++)
        Matriz[i] = new int[cols];
    return;
}

void Horse::espacio(int *&Vector, int cols)
{
    Vector = new int[cols];
    return;
}

void Horse::espacio(bool *&Vector, int cols)
{
    Vector = new bool[cols];
    return;
}

void Horse::espacio(double *&Vector, int cols)
{
    Vector = new double[cols];//double
    return;
}

// void Horse::espacio(double *&Vector, string cols)
// {
//     int columnas = stoi(cols);
//     Vector = new double[columnas];//double
//     return;
// } g++ -std=c++11 -Wall -c *.cc

// void Horse::espacio(string *&Vector, int cols)
// {
//     Vector = new int[cols];
//     return;
// }

double Horse::randreal()
{
  //returns random double in the [0, 1) interval
  double res = gsl_rng_uniform (r);
  return res;
}

double Horse::randgauss(int sigma)
{
  double res = gsl_ran_gaussian(r, sigma);
//   cout << "Horse res= " << res << endl;
  return res;
}

int Horse::randflat(int min, int max)
{
  int res = (int) gsl_ran_flat(r, min, max);
  return res;
}

double Horse::randuniform(int min, int max)
{
  double res = gsl_ran_flat(r, min, max);
  return res;
}

void Horse::spearman(double *Vector1, double *Vector2, int num, double &rho, double &pval)
{
    double *workspace = new double[2 * num];
    
    rho = gsl_stats_spearman(Vector1, 1, Vector2,1, num, workspace);

    double t = rho * (sqrt((num - 2) / (1 - (rho * rho))));
    
    pval = gsl_cdf_tdist_Q(t, num - 2);
    if (pval < 0.5)
        pval = 2*pval;
    else
        pval = 2*(1-pval);
    
    delete[]workspace;
    return;
}

void Horse::pearson(double *Vector1, double *Vector2, int num, double &rho, double &pval)
{
    
    rho = gsl_stats_correlation(Vector1, 1, Vector2, 1, num);

    double t = rho * (sqrt((num - 2) / (1 - (rho * rho))));
    
    pval = gsl_cdf_tdist_Q(t, num - 2);
    if (pval < 0.5)
        pval = 2*pval;
    else
        pval = 2*(1-pval);
    
    return;
}

void Horse::close_rng()
{
  gsl_rng_free (r);
}

int Horse::find_min(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

double Horse::find_min(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

void Horse::fillv0(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Horse::fillv1(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Horse::fillv0(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Horse::fillm0(double **Matriz, int tam1, int tam2)
{
  for(int i = 0; i < tam1; i++){
     for(int j = 0; j < tam2; j++){
         Matriz[i][j] = 0;
     }
  }
}

void Horse::fillm0(int **Matriz, int tam1, int tam2)
{
  for(int i = 0; i < tam1; i++){
     for(int j = 0; j < tam2; j++){
         Matriz[i][j] = 0;
     }
  }
}

void Horse::fillv0(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = false;
}

double Horse::coin()
{
    double coin;
    double random = randgauss(1);
    if(random > 0){coin = random;}
    else{coin = 0.0;}
    
    return coin;
}

void Horse::printmat(double **matriz, int hils, int cols){
    ofstream myfile1;
    myfile1.open ("matrices.txt", ios::app);
    for (int i = 0; i < hils; i++){
        for (int j = 0; j < cols; j++){
            //myfile1 << setw(8) << fixed << setprecision(4) << "m" << i << j << ": " << matriz[i][j];
myfile1 << setprecision(4) << matriz[i][j] << " ";
        }
        myfile1 << endl;
    }
    myfile1 << "\n" << endl;
    myfile1.close ();
}

void Horse::printmat(int **matriz, int hils, int cols){
    ofstream myfile1;
    myfile1.open ("matrices.txt", ios::app);
    for (int i = 0; i < hils; i++){
        for (int j = 0; j < cols; j++){
            //myfile1 << setw(8) << fixed << setprecision(4) << "m" << i << j << ": " << matriz[i][j];
myfile1 << setprecision(4) << matriz[i][j] << " ";
        }
        myfile1 << endl;
    }
    myfile1 << "\n" << endl;
    myfile1.close ();
}

void Horse::printvec(int *vector, int cols)
{
    ofstream myfile1;
    myfile1.open ("matrices.txt", ios::app);
    for(int i = 0; i < cols; i++){
        //myfile1 << setw(8) << fixed << setprecision(4) << "v" << i << ": " << vector[i];
myfile1 << setprecision(4) << " v" << i << ": " << vector[i];
    }
    myfile1 << "\n" << endl;
    myfile1.close ();
}

void Horse::printvec(double *vector, int cols)
{
    ofstream myfile1;
    myfile1.open ("matrices.txt", ios::app);
    for(int i = 0; i < cols; i++){
        //myfile1 << setw(8) << fixed << setprecision(4) << "v" << i << ": " << vector[i];
myfile1 << setprecision(4) << " v" << i << ": " << vector[i];
    }
    myfile1 << "\n" << endl;
    myfile1.close ();
}

void Horse::rellenar(double **Matriz1, double **Matriz2, double **Matriz3, int hils1, int cols1, int hils2, int cols2)
{
    //M1: factores
    //M2: marcadores
    //M3: Original (factores + marcadores) o muestra
    for(int i = 0; i < hils1; i++)
        for(int j = 0; j < cols1; j++){
            Matriz3[i][j] = Matriz1[i][j];
        }
    for(int i = hils1; i < (hils1 + hils2); i++)
        for(int j = 0; j < cols2; j++)
            Matriz3[i][j] = Matriz2[i - hils1][j];
}

int Horse::round(double num1, int num2)
{
    int res;
    double val = num1 * num2;
    double val2 = val * 2;
    int val3 = int(val) + int(val) + 1;
    if(val2 < val3)
        res = (int)val;
    else{res = (int)val + 1;}
    return res;
}

void Horse::rellenar(double **Matriz1, double **Matriz2, int hils, int cols)
{
    //M1: factores
    //M2: marcadores
    //M3: Original (factores + marcadores) o muestra
    for(int i = 0; i < hils; i++)
        for(int j = 0; j < cols; j++){
            Matriz2[i][j] = Matriz1[i][j];
        }
}

string Horse::inttostring(int num)
{
	char buff[50];
	int j;
	j = sprintf(buff, "%d", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to int was not possible using Horse::inttostring.\n";
		exit(1);
	}
	string res(buff);
	return res;
}

void Horse::open_ifstream(ifstream& fe, string nomb)
{
	fe.open(nomb.c_str());
	if (!fe.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Horse::open_ifstream.\n";
		exit(1);
	}
}

void Horse::open_ofstream(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str(), ios::app);
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Horse::open_ofstream.\n";
		exit(1);
	}
}

string Horse::doubletostring_tex(double num)
{
	char buff[50];
	int j,k,l;
	string esp,sig,resf;
	j = sprintf(buff, "%g", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to double was not possible using Horse::doubletostring.\n";
		exit(1);
	}
	string res(buff);
	j = char_in_string('e', res);
	if (j==(-1))
		j = char_in_string('E', res);
	if (j>=0) {
		l = res.size();
		k= char_in_string('-', res, 1, l);
		esp = res.substr(j+1, l-(j+1));
		if (esp[0]=='-')
			sig = "-";
		else {
			if (esp[0]!='+') {
				cout << "[Error]: No + sign in exponent. The error was found while performing Horse::doubletostring.\n";
				exit(1);
			}
			sig = "";
		}
		esp = esp.substr(1, (l-(j+2)));
		resf = "$"+res.substr(0, j) +"\\times10^{"+sig+esp+"}$";
	}
	else {
		k=char_in_string('.', res);
		if (k==(-1))
			resf = res;
		else {
			l = res.size();
			while ((l>1) &&(res[l-1]==0)) {
				res = res.substr(0, (l-1));
				l = res.size();
			}
			resf = res;
		}
	}
	return resf;
}

int Horse::char_in_string(char c, string cue, int from, int until)
{
	int j,i = -1;
	for (j=from; j<until; j++) {
		if (cue[j]==c) {
			i = j;
			break;
		}
	}
	return i;
}

int Horse::char_in_string(char c, string cue)
{
	int u, res;
	u = cue.length();
	res = char_in_string(c, cue, 0, u);
	return res;
}

void Horse::sumarmatriz(int **Matriz, int *Vector, int hils, int cols)
{
    for(int k = 0; k < hils; k++){
        Vector[k] = 0;
        for(int j = 0; j < cols; j++){
            Vector[k] += Matriz[k][j];
        }
    }
}

void Horse::media(double **Matriz1, double *Vector1, int hils, int cols)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Matriz1: entrada
    //Vector1: salida
    double sum;
    
    for(int i = 0; i < hils; i++){
        sum = 0;
        for(int j = 0; j < cols; j++){
            sum += Matriz1[i][j];
        }
        Vector1[i]= sum / cols;
    }
}

void Horse::mediacols(double **Matriz1, double *Vector1, int hils, int cols)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Matriz1: entrada
    //Vector1: salida
    double sum;
    int cuenta = 0;
    for(int i = 0; i < cols; i++){
        sum = 0;
        for(int j = 0; j < hils; j++){
            if(Matriz1[i][j] != -1){
                sum += Matriz1[i][j];
                cuenta++;
            }
        }
        Vector1[i]= sum / cuenta;
    }
}
double Horse::mediacols(double *Vector1, int hils)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Vector1: entrada
    //Vector2: salida
    double sum;
    int cuenta = 0;
    sum = 0;
    for(int j = 0; j < hils; j++){
        if(Vector1[j] != -1){
            sum += Vector1[j];
            cuenta++;
        }
    }
    
    return (sum / cuenta);
}

void Horse::mediacols(double *Vector1, double& media, int hils)
{
    //Saca el mínimo, máximo y promedio de una matriz (en ese orden)
    //Vector1: entrada
    //Vector2: salida
    double sum;
    int cuenta;
    cuenta = 0;
    sum = 0;
    for(int j = 0; j < hils; j++){
        if(Vector1[j] != -1){
        sum += Vector1[j];
        cuenta++;
        }
    }
    media = sum / cuenta;
}

void Horse::desv(double *Vector1, double *Vector2, double *Vector3, int num1, int num2)
{
    
    /////////////////FUNCIÓN MAL!!!!!!!!!!!!!
    //Vector1: Entrada
    //Vector2: Media
    //Vector3: Salida
    //num1: Cantidad de valores, límite de i
    //num2: muestra en la que va
    double Var; //Varianza de las distancias por muestra
    
    double sum = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] != -1){
            sum += (Vector1[i] - Vector2[num2]) * (Vector1[i] - Vector2[num2]);
        }
    }
          
          
          Var = sum / (num1 - 1);
          Vector3[num2] = sqrt(Var);
}

double Horse::desv(double *Vector1, int num1)
{
    //Vector1: Entrada
    //Vector2: Media
    //Vector3: Salida
    //num1: Cantidad de valores, límite de i
    //num2: muestra en la que va
    
    double sum, media, Var;
    int cuenta = 0;
    
    sum = 0;
    for(int j = 0; j < num1; j++){
        if(Vector1[j] != -1){
            sum += Vector1[j];
            cuenta++;
        }
    }
    
    media = sum / cuenta;
    
    sum = 0;
    for(int i = 0; i < num1; i++){
        if(Vector1[i] != -1){
            sum += (Vector1[i] - media) * (Vector1[i] - media);
        }
    }
    
    
    Var = sum / (cuenta - 1);
    return sqrt(Var);
}

double Horse::densidadparcial(double **Matriz1, int tam1, int tam2)
{
    int celdas, cuenta;
    double res;
    cuenta = 0;
    celdas = tam1 * tam2;
    
    for(int i = 0; i < tam1; i++){
        for(int j = 0; j < tam2; j++){
            if(Matriz1[i][j] != 0){cuenta++;}
        }
    }
            
    res = (double)cuenta/celdas;
    
    return res;
}

double Horse::densidadtotal(double **Matriz1, int tam1, int tam2, double **Matriz2, int tam3, int tam4)
{
    int celdas, cuenta;
    double res;
    cuenta = 0;
    celdas = (tam1 * tam2) + (tam2 * tam3);
    
    for(int i = 0; i < tam1; i++){
        for(int j = 0; j < tam2; j++){
            if(Matriz1[i][j] != 0){cuenta++;}
        }
    }
    
    for(int i = 0; i < tam3; i++){
        for(int j = 0; j < tam4; j++){
            if(Matriz2[i][j] != 0){cuenta++;}
        }
    }
            
    res = (double)cuenta/celdas;
    
    return res;
}

void Horse::nonullscuenta(double *Vector1, double *Vector2, int tam, int& cuenta)
{
    cuenta = 0;
    
    for(int i = 0; i < tam; i++){
        if(Vector1[i] != -1){
            if(Vector2[i] != -1){
                cuenta++;
            }
        }
    }
}

void Horse::nonullsvecs(double *Vector1, double *Vector2, int tam, double *Vector3, double *Vector4)
{
    for(int i = 0; i < tam; i++){
        if(Vector1[i] != -1){
            if(Vector2[i] != -1){
                Vector3[i] = Vector1[i];
                Vector4[i] = Vector2[i];
            }
        }
    }
}

void Horse::sort(int* vec, int* nvec, int tam)
{
	int min, max;
	int i, j, quedan, k;
	int* wov;
	wov = new int[tam];
// 	check_creation(wov, "working vector in Horse::sort");
	int* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
    
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
        
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new int[quedan-1];
// 		check_creation(wov2, "working vector2 in Horse::sort");
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new int[quedan];
// 		check_creation(wov, "working vector in Horse::sort");
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
// 		cout << "[Error]: sort does not converge in Horse::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

int Horse::find_max(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

double Horse::find_max(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

int Horse::count_ordinals(int *Vec, int hils1)
{
    int cuenta_ordinals = 1;
    for(int hils = 0; hils < (hils1 - 1); hils++){
        if(Vec[hils] != Vec[hils + 1]){
            cuenta_ordinals++;
        }
    }
    return cuenta_ordinals;
}

void Horse::ordinals_value(int *Vec, int hils1, double **Mat, int hils2)
{
    int cuenta_ordinals2 = 0;
    int cuenta_ordinals3 = 0;
    for(int hils = 0; hils < (hils1 - 1); hils++){
        cuenta_ordinals3++;
        //cout << "Vec_orden_antes[hils]: " << Vec_orden_antes[hils] << endl;
        //cout << "Vec_orden_antes[hils + 1]: " << Vec_orden_antes[hils + 1] << endl;
        if(Vec[hils] != Vec[hils + 1]){
            //cout << "hils: " << hils << " hils2: " << hils2 << endl;
            cuenta_ordinals2++;
            //cout << "cuenta_ordinals2: " << cuenta_ordinals2 << endl;
            //cout << "cuenta_ordinals3: " << cuenta_ordinals3 << endl;
            Mat[cuenta_ordinals2 - 1][0] = (double)Vec[hils];
            if(cuenta_ordinals3 > 1){
                Mat[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2 + 1/(double)cuenta_ordinals3;
            }
            else{
                Mat[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2;
            }
            //cout << Mat_ordinals_antes[cuenta_ordinals2][0] << " " << Mat_ordinals_antes[cuenta_ordinals2][1] << endl;
            //hils2 = hils;
            if(hils < (hils1 - 2)){cuenta_ordinals3 = 0;}
        }
    }
    
    Mat[hils2 - 1][0] = (double)Vec[hils1 - 1];
    if(cuenta_ordinals2 != hils2){
        //El ùltimo nùmero de la lista es diferente
        Mat[hils2 - 1][1] = (double)(hils2);
    }
    else{
        //El ùltimo nùmero de la lista es igual al o a los anteriores
        Mat[hils2 - 1][1] = (double)hils2 + 1/(double)(cuenta_ordinals3 + 1);
    }
}

double Horse::get_median(int* vec, int tam)
{
    double med;
    int *vt;
    vt = new int[tam];
    //check_creation(vt, "temp vec in Basics::get_median");
    sort(vec, vt, tam);
    if ((tam%2) == 0)
        med = (vt[tam/2] + vt[(tam/2)-1])/2.0;
    else
        med = vt[(tam-1)/2];
    delete [] vt;
    return med;
}

double Horse::get_median(double* vec, int tam)
{
    double med;
    double *vt;
    vt = new double[tam];
    //check_creation(vt, "temp vec in Basics::get_median");
    sort(vec, vt, tam);
    if ((tam%2) == 0)
        med = (vt[tam/2] + vt[(tam/2)-1])/2.0;
    else
        med = vt[(tam-1)/2];
    delete [] vt;
    return med;
}

void Horse::sort(double* vec, double* nvec, int tam)
{
    double min, max;
    int i, j, quedan, k;
    double* wov;
    wov = new double[tam];
    //check_creation(wov, "working vector in Basics::sort");
    double* wov2;
    for (i=0; i < tam; i++)
        wov[i] = vec[i];
    max = find_max(wov, tam);
    quedan = tam;
    i = 0;
    while (quedan > 0) {
        min = find_min(wov, quedan);
        for (j=0; j< quedan; j++) {
            if (wov[j] == min) {
                nvec[i] = wov[j];
                i++;
                break;
            }
        }
        wov2 = new double[quedan-1];
        //check_creation(wov2, "working vector2 in Basics::sort");
        for (k = 0; k < j; k++)
            wov2[k] = wov[k];
        for (k=j; k < (quedan-1); k++)
            wov2[k] = wov[k+1];
        delete [] wov;
        quedan--;
        wov = new double[quedan];
        //check_creation(wov, "working vector in Basics::sort");
        for (k=0; k < quedan; k++)
            wov[k] = wov2[k];
        delete [] wov2;
    }
    if (min != max) {
        cout << "[Error]: sort does not converge in Basics::sort.\n";
        exit(1);
    }
    delete [] wov;
    return;
}

double Horse::dist_fenotipica(int *ancestro, int *fenotipo, int tam) //Si el ancestro y el fenotipo tienen una hilera
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

double Horse::dist_fenotipica(int *ancestro, int **fenotipo, int hilsf, int tam) //Si el ancestro tiene una hilera pero el fenotipo es una matriz
{
  int i,j;
  double dif = 0;
  for (i=0; i<hilsf; i++)
    for (j=0; j < tam; j++)
      if (ancestro[j] != fenotipo[i][j])
        dif = dif + (1.0/double(hilsf)); //Compara el ancestro con cada hilera del fenotipo y si difieren suma 1 y saca el promedio

  return dif/double(tam);
}

double Horse::dist_fenotipica(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam) //Si el ancestro y el fenotipo son matrices
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
    dis = dist_aux2(nugoal, (hilsa*mulng), nuatr, (hilsf*mulna), tam);
    for (i= 0; i<(hilsa*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(hilsf*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
  }
  else
    dis = dist_aux2(ancestro, hilsa, fenotipo, hilsf, tam);
//   cout << "dist 3" << endl;

  return dis;
}

double Horse::dist_aux2(int **ancestro, int hilsa, int **fenotipo, int hilsf, int tam)
{
  if (hilsa != hilsf) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists,d;
  dists = new double[hilsa];
  fillv0(dists, hilsa);
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
  
  d = find_min(dists, hilsa);

  d = d/double(tam);
  delete [] dists;
  return d;
}

double Horse::dist_aux2(int *vector1, int cols, int *vector2)
{
  int i,j;
  double *dists,d;
  dists = new double[cols];
  fillv0(dists, cols);
  for (i=0; i < cols; i++) {
    dists[i] =0;
    for (j=0; j < cols; j++) {
//        cout << "W" << i << ": " << vector1[i] << endl;
//        cout << "F" << i << ": " << vector2[(j+i)%cols] << endl;
        if (vector1[i] != vector2[(j+i)%cols]){
          dists[i] = dists[i] + (1.0/double(cols));
//          cout << vector1[i] << " != " << vector2[(j+i)%cols] << endl;
//          cout << "dists" << i << ": " << dists[i] << endl;
        }
      }
    }
  
  d = find_min(dists, cols);
//  cout << "d1: " << d << endl;
  d = d/double(cols);
//  cout << "d2: " << d << endl;
  delete [] dists;
  return d;
}

void Horse::pantalla(int **M, int hils, int cols)
{
    for(int i = 0; i < hils; i++){
        for(int j = 0; j < cols; j++){
            cout << M[i][j] << " " ;
        }
        cout << endl;
    }
    cout << endl;
}

void Horse::pantalla(double **M, int hils, int cols)
{
    for(int i = 0; i < hils; i++){
        for(int j = 0; j < cols; j++){
            cout << M[i][j] << " " ;
        }
        cout << endl;
    }
    cout << endl;
}

void Horse::pantalla(int *V, int hils)
{
    for(int i = 0; i < hils; i++){
        cout << V[i] << " ";
    }
    cout << endl;
}

void Horse::pantalla(double *V, int hils)
{
    for(int i = 0; i < hils; i++){
        cout << V[i] << " ";
    }
    cout << endl;
}

double Horse::sh_index(int *valordek, int tam)
{
    double nlogn, fi, sum, H;
    //Cantidad de fenotipos
    int tamk = 0;
    int *vecesxfenotipo;
    int n = 0; //Número total de fenotipos
    int k = 0; //Número de fenotipos diferentes
    for(int i = 0; i < tam; i++){
                if(valordek[i] != 0){
                    tamk++;
                }
            }
     espacio(vecesxfenotipo, tamk);
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
    borrar(vecesxfenotipo);
     if(H==0){cout << "Sh1 n: " << n <<  " k: " << k << " nlogn: " << nlogn << " sum: " << sum << " H: " << H << endl;}
     
     return H;
}

//void Horse::reduce_fenotipos(int **Fenotipos_notol, int ST, int *perd, int cuenta_pnotol, int **Fenotipos_notol_reducidos, int *perd_nuevo)
void Horse::reduce_fenotipos(int **Fenotipos_para_reducir, int columnas, int *per_para_reducir, int cuenta_p, int *per_reducidos)
{
    //Reduce fenotipos que podrían ser puntos fijos en ciclos o ciclos más pequeños.
    int sum1, sum2;
    int **para_reducir, *hil_compara, **hils_compara;
    double dist;
    ifstream fi;
    
    sum1 = 0;
    sum2 = 0;
    dist = 0;
    for(int cuenta_fenotipo = 0; cuenta_fenotipo < cuenta_p; cuenta_fenotipo++){
        per_reducidos[cuenta_fenotipo] = per_para_reducir[cuenta_fenotipo];
        if(per_para_reducir[cuenta_fenotipo] > 1){
            espacio(para_reducir, per_para_reducir[cuenta_fenotipo], columnas);
            for(int i = 0; i < per_para_reducir[cuenta_fenotipo]; i++){ //Llenar vector de fn antes
                for(int j = 0; j < columnas; j++){
                    para_reducir[i][j] = Fenotipos_para_reducir[i + sum1][j];
                }
            }
            //Punto Fijo
            espacio(hil_compara, columnas);
            for(int i = 0; i < columnas; i++){
                hil_compara[i] = para_reducir[0][i];
            }
            dist = dist_fenotipica(hil_compara, para_reducir, per_para_reducir[cuenta_fenotipo], columnas);
            borrar(hil_compara);
            if(dist == 0){
                per_reducidos[cuenta_fenotipo] = 1;
            }
            else{//Ciclo
                if(per_para_reducir[cuenta_fenotipo] > 3){//Si es mayor de 3 hileras
                    if(per_para_reducir[cuenta_fenotipo]%2 == 0){//Si es par
                        for(int i = 2; i < per_para_reducir[cuenta_fenotipo]; i++){
                            if(per_para_reducir[cuenta_fenotipo]%i == 0){//Que sea divisor para no perder tiempo
                                //Lleno una fracción para comparar
                                espacio(hils_compara, i, columnas);
                                for(int hils = 0; hils < i; hils++){
                                    for(int cols = 0; cols < columnas; cols++){
                                        hils_compara[hils][cols] = para_reducir[hils][cols];
                                    }
                                }

                                dist = dist_fenotipica(hils_compara, i, para_reducir, per_para_reducir[cuenta_fenotipo], columnas);
                                borrar(hils_compara, i);
                                if(dist == 0){
                                    per_reducidos[cuenta_fenotipo] = i;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        sum1 += per_para_reducir[cuenta_fenotipo];
    }
    
    //ESto debe ir afuera
//    suma_perd_nuevo = 0;
//    for(int i = 0; i < cuenta_p; i++){
//        suma_per_nuevo += per_reducidos[i];
//    }
//    //Llenamos con fenotipos reducidos
//    sum1 = 0;
//    for(int i = 0; i < cuenta_p; i++){
//        for(int j = 0; j < per_reducidos[i]; j++){
//            for(int k = 0; k < cols; k++){
//                Fenotipos_notol_reducidos[j][k] = Fenotipos_para_reducir[j + sum1][k];
//            }
//        }
//        sum1 += per_para_reducir[i];
//    }
}


//double GraphI::eval_mod(double **mamod, int *vop) {
//  int i, j;
//  double Q = 0;
//  for (i = 0; i < size; i++)
//    for (j = 0; j < size; j++)
//      if (vop[i] == vop[j])
//        Q += mamod[i][j];
//  Q /= (2.0*number_of_edges());
//  return Q;
//}


void Horse::ofstream_valor(string root, string archivo, string name, double valor){
    ofstream fo;
    open_ofstream(fo, ""+ root +""+ archivo +""+ name +".txt");
    fo << valor << endl;
    fo.close();
    
    
}

void Horse::ofstream_valor(string root, string archivo, string name, double valor1, double valor2, double valor3){
    ofstream fo;
    open_ofstream(fo, ""+ root +""+ archivo +""+ name +".txt");
    fo << valor1 << " " << valor2 << " " << valor3 << endl;
    fo.close();
    
    
}

void Horse::ofstream_valor(string root, string archivo, string name, double valor1, double valor2, double valor3, double valor4, double valor5, double valor6, double valor7, double valor8, double valor9){
    ofstream fo;
    open_ofstream(fo, ""+ root +""+ archivo +""+ name +".txt");
    fo << valor1 << " " << valor2 << " " << valor3 << " " << valor4 << " " << valor5 << " " << valor6 << " " << valor7 << " " << valor8 << " " << valor9 << endl;
    fo.close();
    
    
}








