//Program to remove recurrently accessible phenotypes.
// g++ -Wall norec_phenotypes.cc libs/horse.o libs/linear.o libs/finder.o -lgsl -lgslcblas -lm

// ./a.out 1000 12 6 /outs dup
// ./a.out 1000 12 6 /outs c1
// ./a.out 1000 12 6 /outs c2

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <time.h>
#include "libs/horse.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
    
    Horse hor;
    
    //Sólo para los tolerantes
    int muestra = atoi(argv[1]); //Número de muestra
    int Nf = atoi(argv[2]); //Número de genes de transcripción
    int Ns = atoi(argv[3]); //Número de genes estructurales
    string arch1 = argv[4];//Ruta
    string arch2 = argv[5];//dup, c1, c2
    
    int mut = 432;
    
    ifstream fi;
    ofstream fo;
    
    int *per, **fn, **valordek;
    int **mat_di, *eventos;
    hor.espacio(eventos, muestra);
    hor.espacio(valordek, mut, muestra);
    hor.espacio(mat_di, muestra, Nf);
    int cuenta_per, cuenta_fn, sum;
    int contador, tipo, sum_per;
    string basura;
    
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_"+ arch2 +".txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos[i];
    }
    fi.close();

    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch2 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> mat_di[i][j];
        }
    }
    fi.close();
    
    if(arch2 == "dup"){tipo = 1;}
    if(arch2 == "c1"){tipo = 2;}
    if(arch2 == "c2"){tipo = 3;}
    
    contador = 0;
    for(int m = 0; m < muestra; m++){
        cout << "m: " << m << endl;
        
        if(mat_di[m][eventos[m]] == 1){
            cout << "m: " << m << " evento: " << eventos[m] << endl;
            contador++;
            
            cuenta_per = 0;
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/per_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
            while(fi >> basura)
                cuenta_per++;
            fi.close();
            
            hor.espacio(per, cuenta_per);
            
            sum_per = 0;
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/per_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
            for(int i = 0; i < cuenta_per; i++){
                fi >> per[i];
                sum_per += per[i];
            }
            fi.close();
            
            cuenta_fn = 0;
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
            while(fi >> basura)
                cuenta_fn++;
            fi.close();
            
            cuenta_fn = cuenta_fn / Ns;
            if(sum_per != cuenta_fn){
                cout << "Cuidado! No concuerda " << endl;
                cout << "cuenta_fn = " << cuenta_fn << " sum_per = " << sum_per << endl;
            }
            
            hor.espacio(fn, cuenta_fn, Ns);
            
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
            for(int j = 0; j < cuenta_fn; j++){
                for(int k = 0; k < Ns; k++){
                    fi >> fn[j][k];
                }
            }
            fi.close();
                        
            hor.open_ifstream(fi, ""+ arch1 +"/acces/valordek/valordek_"+ arch2 +"_" + hor.inttostring(m) + "de1000.txt");
            for(int k = 0; k < mut; k++){
                for(int j = 0; j < Nf; j++){
                    fi >> valordek[k][j];
                }
            }
            fi.close();
                        
            sum = 0;
            for(int i = 0; i < cuenta_per; i++){
                if(valordek[i][eventos[m]] != 0){
                    hor.open_ofstream(fo, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
                                        
                    for(int o = 0; o < per[i]; o++){
                        for(int q = 0; q < Ns; q++){
                            fo << fn[o + sum][q] << " ";
                        }
                        fo << endl;
                    }
                    fo << endl;
                    sum += per[i];
                    
                    fo.close();
                                        
                    hor.open_ofstream(fo, ""+ arch1 +"/acces/fn_norepetido/per_" + hor.inttostring(m) + "_" + hor.inttostring(tipo) + "_" + hor.inttostring(eventos[m]) + ".txt");
                    fo << per[i] << endl;
                                  fo.close();
                }
                else{sum += per[i];}
            }
                                                                    
            hor.borrar(per);
            hor.borrar(fn, cuenta_fn);
            
        }
    }
    
    hor.borrar(eventos);
    hor.borrar(mat_di, muestra);
    hor.borrar(valordek, mut);
    
    cout << "contador: " << contador << endl;
    cout << "------------Ejecución terminada!!-----------\n";
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    return 0;
}
