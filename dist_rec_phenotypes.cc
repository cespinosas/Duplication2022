//Programa para medir las distancias entre fenotipos nuevos. Va antes de Spearmyphenotype251221

// g++ -Wall dist_fenotipos_recurrentes261221.cp libs/horse.o libs/finder.o libs/linear.o -lgsl -lgslcblas -lm
// ./a.out 1000 6 12 /Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120

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
#include "libs/finder.h"

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

int main(int argc, char **argv){
    
    system("date");
    
    int muestra = atoi(argv[1]);
    int Ns = atoi(argv[2]);
    int Nf = atoi(argv[3]);
    string arch1 = argv[4];//Ruta

    int cuenta1, per_antes, cuenta_dup, cuenta_c1, cuenta_c2;
    string basura;
    int *eventos_dup, *eventos_c1, *eventos_c2;
    int **mat_di_dup, **mat_di_c1, **mat_di_c2;
    
    Horse hor;
    Finder find;
    
    int **Fn_antes, *Per_antes, **acces1, **acces2, **acces3;
    
    hor.espacio(acces1, muestra, Nf); //dup
    hor.espacio(acces2, muestra, Nf); //c1
    hor.espacio(acces3, muestra, Nf); //c2
    
    hor.espacio(mat_di_dup, muestra, Nf);
    hor.espacio(mat_di_c1, muestra, Nf);
    hor.espacio(mat_di_c2, muestra, Nf);
    
    ofstream fo;
    ifstream fi;
    
    hor.espacio(eventos_dup, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_dup.txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos_dup[i];
    }
    fi.close();
    
    hor.espacio(eventos_c1, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_c1.txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos_c1[i];
    }
    fi.close();
    
    hor.espacio(eventos_c2, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_c2.txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos_c2[i];
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/acces/accesdup.txt");
    for(int j = 0; j < muestra; j++){
        for(int k = 0; k < Nf; k++){
            fi >> acces1[j][k];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/duppasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> mat_di_dup[i][j];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/c1pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> mat_di_c1[i][j];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/c2pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> mat_di_c2[i][j];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/acces/accesc1.txt");
    for(int j = 0; j < muestra; j++){
        for(int k = 0; k < Nf; k++){
            fi >> acces2[j][k];
        }
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/acces/accesc2.txt");
    for(int j = 0; j < muestra; j++){
        for(int k = 0; k < Nf; k++){
            fi >> acces3[j][k];
        }
    }
    fi.close();
    
    cuenta_dup = 0;
    cuenta_c1 = 0;
    cuenta_c2 = 0;
    for(int i = 0; i < muestra; i++){
        
        cout << "muestra: " << i << endl;
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + ".txt");
        cuenta1 = 0;
        while(fi >> basura)
            cuenta1++;
        fi.close();
        
        hor.espacio(Fn_antes, (cuenta1/Ns), Ns);
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + ".txt");
        for(int j = 0; j < cuenta1/Ns; j++){
            for(int k = 0; k < Ns; k++){
                fi >> Fn_antes[j][k];
            }
        }
        fi.close();
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/per_" + hor.inttostring(i) + ".txt");
        per_antes = 0;
        while(fi >> basura)
            per_antes++;
        fi.close();
        
        hor.espacio(Per_antes, per_antes);
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/per_" + hor.inttostring(i) + ".txt");
        for(int j = 0; j < per_antes; j++){
            fi >> Per_antes[j];
        }
        fi.close();
        
        ///////////////////////COMPARAR CONJUNTO DE FN DE ESPÉCIMEN
        ////////////Mide distancia entre fnuevos de los especímienes
        find.dist_btw_esp(per_antes, Per_antes, Ns, Fn_antes, ""+ arch1 +"/acces/dist/dist_fn_esp");
        
        
        ////////////Mide distancias entre fenotipos nuevos después de agregar un gen y antes contra después
        
        if(mat_di_dup[i][eventos_dup[i]] != 0){
            cuenta_dup++;
            find.dist_ant_vs_desp_tol(Nf, Ns, cuenta1, acces1[i][eventos_dup[i]], i, ""+ arch1 +"/acces/", "1", "dup", Fn_antes, per_antes, Per_antes, eventos_dup[i]);
        }
        
        if(mat_di_c1[i][eventos_c1[i]] != 0){
            cuenta_c1++;
            find.dist_ant_vs_desp_tol(Nf, Ns, cuenta1, acces2[i][eventos_c1[i]], i, ""+ arch1 +"/acces/", "2", "c1", Fn_antes, per_antes, Per_antes, eventos_c1[i]);
        }
        
        if(mat_di_c2[i][eventos_c2[i]] != 0){
            cuenta_c2++;
            find.dist_ant_vs_desp_tol(Nf, Ns, cuenta1, acces3[i][eventos_c2[i]], i, ""+ arch1 +"/acces/", "3", "c2", Fn_antes, per_antes, Per_antes, eventos_c2[i]);
        }
        
        hor.borrar(Fn_antes, (cuenta1/Ns));
        hor.borrar(Per_antes);
        
    }

    hor.borrar(acces1, muestra); //dup
    hor.borrar(acces2, muestra); //c1
    hor.borrar(acces3, muestra); //c2
    
    hor.borrar(eventos_dup);
    hor.borrar(eventos_c1);
    hor.borrar(eventos_c2);
    
    hor.borrar(mat_di_dup, muestra);
    hor.borrar(mat_di_c1, muestra);
    hor.borrar(mat_di_c2, muestra);
    
    //#########################################################################################################//
        
    cout << "------------Ejecución terminada!!-----------\n";
    
    //     cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    return 0;
}
