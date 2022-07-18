// Similarity per kind of mutation and gene interaction
// g++ -Wall simialrity_permut_pergene.cc libs/horse.o libs/linear.o libs/finder.o libs/gataca.o -lgsl -lgslcblas -lm /root dup especimenes_dup 1740

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "libs/horse.h"
#include "libs/linear.h"
#include "libs/finder.h"
#include "libs/gataca.h"

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

int main(int argc, char **argv)
{
    system("date");
    time_t rawtime;
    time ( &rawtime );
    
    int muestra = atoi(argv[1]); //Número de redes de genotipos a estudiar
    int Nf = atoi(argv[2]); //Número de genes de transcripción
    int Np = atoi(argv[3]); //Número de genes marcadores del fenotipo
    string arch1 = argv[4];//Ruta
    string arch2 = argv[5];//esp, dup, c1, c2
    string arch3 = argv[6];//archivo de especímenes: especimenes, especimenes_dup, especimenes_c1, especimenes_c2
    int semilla = atoi(argv[7]);

    Horse hor(semilla);
    Finder find;
    Linear lin;
    GATACA gat;
    
    ifstream fi;
    ofstream fo;
    
    //Variables
    int n, p;
    double dist;
    double cuenta_in, cambio_in, cuenta_out, cambio_out;
    double suma_in, suma_out;
    
    //Vectores
    int *condinicial, *vectorentrada, *eventos, *Wish;
    hor.espacio(condinicial, Nf);
    hor.espacio(vectorentrada, Nf);
    hor.espacio(eventos, muestra);
    
    //Matrices
    double **bichos, **Matrizf, **Matrizp, **mat_di;
    int **trayectoria, **Fenotipo;
    hor.espacio(bichos, muestra*(Nf + Np), Nf);
    hor.espacio(Matrizf, Nf, Nf);
    hor.espacio(Matrizp, Np, Nf);
    hor.espacio(trayectoria, 1000, Nf);
    hor.espacio(mat_di, muestra, Nf);
    
    //Condición inicial (Medio ambiente)
    for(int i = 0; i < Nf; i++){
        condinicial[i] = -1;
    }
    
    for(int i = 0; i < Nf; i++){
        vectorentrada[i] = condinicial[i];
        trayectoria[0][i] = condinicial[i];
    }
    
    //Fenotipo original
    hor.espacio(Wish, Np);
    
    Wish[0] = 1;
    Wish[1] = -1;
    Wish[2] = -1;
    Wish[3] = -1;
    Wish[4] = -1;
    Wish[5] = 1;
    
//     cout << "1" << endl; 
    
    //Llamamos a las redes muestra a la matriz Bichos
    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch3 +".txt");
    for(int i = 0; i < (muestra*(Nf + Np)); i++){
        for(int j = 0; j < Nf; j++){
            fi >> bichos[i][j];
        }
    }
    fi.close();
    
    
    //Llamamos la lista de genes a duplicar y mutar.
    hor.open_ifstream(fi, ""+ arch1 +"/outs070120/acces/eventos_dup.txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos[i];
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/outs070120/"+ arch2 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < (Nf -1); j++){
           fi >> mat_di[i][j];
        }
    }
    fi.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    
    
    hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
    fo << "Dup_Dup Ndup_Dup Dup_Ndup Ndup_Ndup Ndup_ST Dup_ST\n";
    fo.close();
    
    for(int m = 0; m < muestra; m++){
        
        cout << "Muestra: " << m << endl;
        
        for(int i = 0; i < Nf; i++){
            for(int j = 0; j < Nf; j++){
                Matrizf[i][j] = bichos[(m)*(Nf + Np) + i][j];
            }
        }
        
        for(int i = 0; i < Np; i++){
            for(int j = 0; j < Nf; j++){
                Matrizp[i][j] = bichos[((m + 1)*Nf + m*Np) + i][j];
            }
        }
        
        //Fenotipo
        for(int i = 0; i < Nf; i++){
            vectorentrada[i] = condinicial[i];
            trayectoria[0][i] = condinicial[i];
        }
        
        
        p = 0;
        n = 1;
        find.atractor(Matrizf, trayectoria, vectorentrada, Nf, n, p);
        hor.espacio(Fenotipo, p, Np);
        lin.dot(Matrizp, trayectoria, Fenotipo, Nf, Np, n, p);
        
        dist = find.distfenotipica(Wish, Fenotipo, p, Np);
        hor.borrar(Fenotipo, p);
        if(dist != 0){
             cout << "No coincide fenotipo\n";
        }
        
        ////////////////
        //Dup con Dup//
        ///////////////
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_lost_dup_dup\n";
        
        gat.interaction_lost_dup_dup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial);
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_lost.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_gain_dup_dup\n";
        
        gat.interaction_gain_dup_dup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
                
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_gain.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_change_dup_dup\n";
        
        gat.interaction_change_dup_dup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
                
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_change.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        /////////////////
        //Dup con NDup//
        ////////////////
        suma_in = 0;
        cuenta_in = 0;
        cambio_in = 0;
        suma_out = 0;
        cuenta_out = 0;
        cambio_out = 0;
        
//         cout << "interaction_lost_dup_Ndup\n";
        
        gat.interaction_lost_dup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, suma_out, cuenta_out, cambio_out, Wish, condinicial);
                
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_lost.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " " << suma_out << " " << cuenta_out << " " << cambio_out << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA" << endl;
        }
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_gral.txt");
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_out!= 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_out != 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        suma_out = 0;//suma de similitud (1 - dist)
        cuenta_out = 0;//cuántas autoregulaciones hay
        cambio_out = 0;//cambio en el fenotipo
        //Gain
        
//         cout << "interaction_gain_dup_Ndup" << endl;
        
        gat.interaction_gain_dup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, suma_out, cuenta_out, cambio_out, Wish, condinicial, hor.randflat(1000,1000000));
                
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_gain.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " " << suma_out << " " << cuenta_out << " " << cambio_out << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA" << endl;
        }
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_gral.txt");
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_out!= 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_out != 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        suma_out = 0;//suma de similitud (1 - dist)
        cuenta_out = 0;//cuántas autoregulaciones hay
        cambio_out = 0;//cambio en el fenotipo
        //Change
        
//         cout << "interaction_change_dup_Ndup" << endl;
        
        gat.interaction_change_dup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, suma_out, cuenta_out, cambio_out, Wish, condinicial, hor.randflat(1000,1000000));
        
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_change.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " " << suma_out << " " << cuenta_out << " " << cambio_out << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA" << endl;
        }
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_gral.txt");
        if(cuenta_out != 0){
            fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Dup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_Ndup_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " " << 1 - (cambio_out / cuenta_out) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_out!= 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_out != 0){
                fo << suma_out / cuenta_out << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_out!= 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_out != 0){
                fo << 1 - (cambio_out / cuenta_out) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_lost_Ndup_Ndup" << endl;
        
        gat.interaction_lost_Ndup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial);
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_lost.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipp
        
//         cout << "interaction_gain_Ndup_Ndup" << endl;
        
        gat.interaction_gain_Ndup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
                
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_gain.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_change_Ndup_Ndup" << endl;
        
        gat.interaction_change_Ndup_Ndup(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
                
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_change.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_Ndup_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_lost_Ndup_ST" << endl;
        
        gat.interaction_lost_Ndup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial);
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_lost.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_gain_Ndup_ST" << endl;
        
        gat.interaction_gain_Ndup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_gain.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_change_Ndup_ST" << endl;
        
        gat.interaction_change_Ndup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_change.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Ndup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA ";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_lost_Dup_ST" << endl;
        
        gat.interaction_lost_Dup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial);
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_lost.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/lostP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_gain_Dup_ST" << endl;
        
        gat.interaction_gain_Dup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_gain.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
        }
        else{
            fo << "NA NA ";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << " ";
            }
            else{
                fo << "NA NA ";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/gainP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        
        suma_in = 0;//suma de similitud (1 - dist)
        cuenta_in = 0;//cuántas autoregulaciones hay
        cambio_in = 0;//cambio en el fenotipo
        
//         cout << "interaction_change_Dup_ST" << endl;
        
        gat.interaction_change_Dup_ST(Matrizf, Nf, Matrizp, Np, eventos[m], suma_in, cuenta_in, cambio_in, Wish, condinicial, hor.randflat(1000,1000000));
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_change.txt");
        fo << m << " " <<  mat_di[m][eventos[m]] << " " << suma_in << " " << cuenta_in << " " << cambio_in << " ";
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << endl;
        }
        else{
            fo << "NA NA" << endl;
        }
        fo.close();
        
        hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_gral.txt");
        if(cuenta_in != 0){
            fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
        }
        else{
            fo << "NA NA\n";
        }
        fo.close();
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/Dup_ST_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << " " << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_tol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeS_ntol.txt");
            if(cuenta_in != 0){
                fo << suma_in / cuenta_in << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        
        if(mat_di[m][eventos[m]] == 1){
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_tol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/changeP_ntol.txt");
            if(cuenta_in != 0){
                fo << 1 - (cambio_in / cuenta_in) << "\n";
            }
            else{
                fo << "NA\n";
            }
            fo.close();
        }
    }
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    hor.borrar(condinicial);
    hor.borrar(vectorentrada);
    hor.borrar(eventos);
    hor.borrar(bichos, muestra*(Nf + Np));
    hor.borrar(Matrizf, Nf);
    hor.borrar(Matrizp, Np);
    hor.borrar(trayectoria, 1000);
    hor.borrar(mat_di, muestra);
    hor.borrar(Wish);
    
    hor.close_rng();
    
    return 0;
    
}
