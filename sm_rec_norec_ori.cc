//Distance between accessible phenotypes and the original phenotype.
// g++ -Wall sm_rec_norec_ori.cc libs/horse.o libs/linear.o libs/finder.o -lgsl -lgslcblas -lm

// ./a.out /outs 1 dup
// ./a.out /outs 2 c1
// ./a.out /outs 3 c2

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

int main(int argc, char **argv)
{
    system("date");
    time_t rawtime;
    time ( &rawtime );

    string arch1 = argv[1];//Ruta
    string arch2 = argv[2];//1, 2 ó 3 (dup, c1 ó c2)
    string arch3 = argv[3];//dup, c1 ó c2
    
    int muestra = 1000;
    int Nf = 12;
    int Ns = 6;
    int cuenta1, cuenta2, per_antes, per_despues, evento, cuentavector, sum1, sum2;

    string basura;

    Horse hor;
    Finder find;

    ifstream fi;
    ofstream fo;

    double dist, dist_r, dist_todo_antes, dist_todo_desp;
    int hils_fn_antes;
    int **Fn_antes, *Per_antes, **Fn_despues, *Per_despues, **Mat_fn_desp, **Mat_conserv, **nuevoantes, **nuevodespues, *Vec_fn_antes, **Mat_dif_igual, *Vec_eventos;
    bool *boolean_desp;

    hor.espacio(Mat_conserv, muestra, Nf);
    hor.espacio(Vec_fn_antes, muestra);
    hor.espacio(Mat_fn_desp, muestra, Nf);
    
    int *Wish;
    hor.espacio(Wish, Ns);
    Wish[0] = 1;
    Wish[1] = -1;
    Wish[2] = -1;
    Wish[3] = -1;
    Wish[4] = -1;
    Wish[5] = 1;
    
    hor.espacio(Mat_dif_igual, muestra, Nf);
    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch3 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> Mat_dif_igual[i][j];
        }
    }
    fi.close();
    
    hor.espacio(Vec_eventos, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_"+ arch3 +".txt");
    for(int i = 0; i < muestra; i++){fi >> Vec_eventos[i]; /*cout << Vec_eventos[i] << endl;*/}
    fi.close();

    hor.open_ifstream(fi, ""+ arch1 +"/acces/acces"+ arch3 +".txt");
    for(int j = 0; j < muestra; j++){
        for(int k = 0; k < Nf; k++){
            fi >> Mat_fn_desp[j][k];
        }
    }
    fi.close();

    hor.open_ifstream(fi, ""+ arch1 +"/acces/accesesp.txt");
    for(int j = 0; j < muestra; j++){
        fi >> Vec_fn_antes[j];
    }
    fi.close();

//     for(int i = 0; i < 3; i++){
    for(int i = 0; i < muestra; i++){
        cout << "muestra: " << i << endl;
        dist_r = 0;
        dist_todo_antes = 0;
        dist_todo_desp = 0;
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(i) + ".txt");

        cuenta1 = 0;
        while(fi >> basura)
            cuenta1++;
        fi.close();
        hils_fn_antes = cuenta1/Ns;
        hor.espacio(Fn_antes, hils_fn_antes, Ns);
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(i) + ".txt");
        for(int j = 0; j < cuenta1/Ns; j++){
            for(int k = 0; k < Ns; k++){
                fi >> Fn_antes[j][k];
            }
        }
        fi.close();
        
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/Per_" + hor.inttostring(i) + ".txt");
        per_antes = 0;
        while(fi >> basura)
            per_antes++;
        fi.close();

        hor.espacio(Per_antes, per_antes);
        hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/Per_" + hor.inttostring(i) + ".txt");
        for(int j = 0; j < per_antes; j++){
            fi >> Per_antes[j];
        }
        fi.close();

        evento = Vec_eventos[i];

        if(Mat_fn_desp[i][evento] != 0){//Comprueva que sí haya fenotipos nuevos generados después
            hor.open_ifstream(fi, ""+ arch1 +"/acces/dist/dist_fn_esp_"+ arch3 +".txt");
            for(int i = 0; i < muestra; i++){
                for(int j = 0; j < Nf; j++){
                    fi >> Mat_conserv[i][j];
                }
            }
            fi.close();

            cuenta2 = 0;
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(i) + "_ "+ arch2 +"_" + hor.inttostring(evento) + ".txt");
            while(fi >> basura)
                cuenta2++;
            fi.close();
            hor.espacio(Fn_despues, cuenta2/Ns, Ns);

            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/fn_" + hor.inttostring(i) + "_ "+ arch2 +"_" + hor.inttostring(evento) + ".txt");
            for(int j = 0; j < cuenta2/Ns; j++){
                for(int k = 0; k < Ns; k++){
                    fi >> Fn_despues[j][k];
                }
            }
            fi.close();

            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/Per_" + hor.inttostring(i) + "_ "+ arch2 +"_" + hor.inttostring(evento) + ".txt");
            per_despues = 0;
            while(fi >> basura)
                per_despues++;
            fi.close();

            hor.espacio(Per_despues, per_despues);

            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn/Per_" + hor.inttostring(i) + "_ "+ arch2 +"_" + hor.inttostring(evento) + ".txt");
            for(int j = 0; j < per_despues; j++){
                fi >> Per_despues[j];
            }
            fi.close();

            hor.espacio(boolean_desp, per_despues);
            hor.fillv0(boolean_desp, per_despues);

            //////////////////////////Conservados//////////////////////////////////////////////
            sum1 = 0;
            cuentavector = 0; //Para ver cuántos fns antes y después son parecidos
            for(int m = 0; m < per_antes; m++){ //Para todos los fns del espécimen i en curso
                hor.espacio(nuevoantes, Per_antes[m], Ns);
                for(int j = 0; j < Per_antes[m]; j++){ //Llenar vector de fn antes
                    for(int k = 0; k < Ns; k++){
                        nuevoantes[j][k] = Fn_antes[j + sum1][k];
                    }
                }
                sum2 = 0;
                for(int n = 0; n < per_despues; n++){ //Llenar vector de fn después
                    if(boolean_desp[n] != true){//Para no volver a comparar a los que ya se les encontró igual
                        hor.espacio(nuevodespues, Per_despues[n], Ns);
                        for(int o = 0; o < Per_despues[n]; o++){
                            for(int q = 0; q < Ns; q++){
                                nuevodespues[o][q] = Fn_despues[o + sum2][q];
                            }
                        }              

                        dist = find.distfenotipica(nuevoantes, Per_antes[m], nuevodespues, Per_despues[n], Ns);
                        
                        if(dist == 0){//Si encontramos un fenotipo igual antes y después
                            boolean_desp[n] = true;
                            dist_r += 1 - find.distfenotipica(Wish, nuevoantes, Per_antes[m], Ns);
                            hor.borrar(nuevodespues, Per_despues[n]);//Para que borre el espacio antes de romper el ciclo
                            break;
                        }
                        hor.borrar(nuevodespues, Per_despues[n]);
                    }
                    sum2 += Per_despues[n];
                }//Fin n
                sum1 += Per_antes[m];
                hor.borrar(nuevoantes, Per_antes[m]);
            }//Fin m
            ///////////////////////////////////////////////////////////////////////////
            //////////////////////////DISTANCIA DE TODOS//////////////////////////////////////////////
            sum1 = 0;
            cuentavector = 0; //Para ver cuántos fns antes y después son parecidos
            for(int m = 0; m < per_antes; m++){ //Para todos los fns del espécimen i en curso
                hor.espacio(nuevoantes, Per_antes[m], Ns);
                for(int j = 0; j < Per_antes[m]; j++){ //Llenar vector de fn antes
                    for(int k = 0; k < Ns; k++){
                        nuevoantes[j][k] = Fn_antes[j + sum1][k];
                    }
                }
                sum1 += Per_antes[m];
                dist_todo_antes += 1 - find.distfenotipica(Wish, nuevoantes, Per_antes[m], Ns);
                hor.borrar(nuevoantes, Per_antes[m]);
            }//Fin m
            sum2 = 0;
            for(int n = 0; n < per_despues; n++){ //Llenar vector de fn después
                hor.espacio(nuevodespues, Per_despues[n], Ns);
                for(int o = 0; o < Per_despues[n]; o++){
                    for(int q = 0; q < Ns; q++){
                        nuevodespues[o][q] = Fn_despues[o + sum2][q];
                    }
                }
                
                dist_todo_desp += 1 - find.distfenotipica(Wish, nuevodespues, Per_despues[n], Ns);
                
                hor.borrar(nuevodespues, Per_despues[n]);
                sum2 += Per_despues[n];
            }//Fin n
            
            ///////////////////////////////////////////////////////////////////////////
            
            hor.borrar(Fn_despues, cuenta2/Ns);
            hor.borrar(Per_despues);
            hor.borrar(boolean_desp); 
        }//Fin Mat_fn_desp[i][evento] != 0
        
        hor.borrar(Per_antes);
        hor.borrar(Fn_antes, hils_fn_antes);
//    cout << evento << endl;
    hor.open_ofstream(fo, ""+ arch1 +"/dist/dist_o/dist_o_"+ arch3 +".txt");
        fo << i << " " << evento << " " << Mat_dif_igual[i][evento] << " " << dist_r << " " << dist_todo_antes << " " << dist_todo_desp << " " << Vec_fn_antes[i] << " " << Mat_fn_desp[i][evento] << " " << Mat_conserv[i][evento] << endl;
        fo.close();
    }//Fin muestra
    
    hor.borrar(Vec_fn_antes);
    hor.borrar(Mat_fn_desp, muestra);
    hor.borrar(Mat_conserv, muestra);
    hor.borrar(Wish);
    hor.borrar(Mat_dif_igual, muestra);
    hor.borrar(Vec_eventos);
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");

//    hor.close_rng();
    return 0;
}
