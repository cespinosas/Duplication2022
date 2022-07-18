//Program to count the number of connections of the gene to be duplicated.
// Distance between the duplicate and the SGs
// g++ -Wall num_conex.cc libs/horse.o -lgsl -lgslcblas -lm
// ./a.out 1000 13 6 /root dup especimenes_dup

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
    string arch2 = argv[5]; //dup, c1, c2
    string arch3 = argv[6];//archivo de especímenes: especimenes, especimenes_dup, especimenes_c1, especimenes_c2
    
    ifstream fi;
    ofstream fo;
    
    Horse hor;
    
    //Variables
    int entrada, salida;
    
    //Matrices
    double **bichos, **Matrizf, **Matrizp;
    int **mat_di, *eventos;
    
    hor.espacio(bichos, muestra*(Nf + Np), Nf);
    hor.espacio(Matrizf, Nf, Nf);
    hor.espacio(Matrizp, Np, Nf);
    hor.espacio(mat_di, muestra, Nf);
    hor.espacio(eventos, muestra);
    
    //Llamamos a las redes muestra a la matriz Bichos
    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch3 +".txt");
    for(int i = 0; i < (muestra*(Nf + Np)); i++){
        for(int j = 0; j < Nf; j++){
            fi >> bichos[i][j];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/outs/acces/eventos_"+ arch2 +".txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos[i];
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/outs/"+ arch2 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < (Nf -1); j++){
            fi >> mat_di[i][j];
        }
    }
    fi.close();
    
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
        
        //Interacciones de entrada
        entrada = 0;
        for(int j = 0; j < Nf; j++){
            if(Matrizf[Nf - 1][j] != 0){
                entrada++;
            }
        }
        
        //Interacciones de salida
        salida = 0;
        for(int i = 0; i < Nf; i++){
            if(Matrizf[i][Nf - 1] != 0){salida++;}
        }
        for(int i = 0; i < Np; i++){
            if(Matrizp[i][Nf - 1] != 0){salida++;}
        }
        
        hor.open_ofstream(fo, ""+ arch1 +"/"+ arch3 +"_added_num_conex_.txt");
        fo << entrada << " " << salida << " " << entrada + salida << endl;
        fo.close();
        
        if(mat_di[m][eventos[m]] != 1){
            hor.open_ofstream(fo, ""+ arch1 +"/"+ arch3 +"_added_num_conex_notol.txt");
            fo << entrada << " " << salida << " " << entrada + salida << endl;
            fo.close();
        }
        else{
            hor.open_ofstream(fo, ""+ arch1 +"/"+ arch3 +"_added_num_conex_tol.txt");
            fo << entrada << " " << salida << " " << entrada + salida << endl;
            fo.close();
        }
    }

                                  
    hor.borrar(bichos, muestra*(Nf + Np));
    hor.borrar(Matrizf, Nf);
    hor.borrar(Matrizp, Np);
    hor.borrar(mat_di, muestra);
    hor.borrar(eventos);
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    return 0;
    
}




