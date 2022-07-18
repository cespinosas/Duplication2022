//This program arrange number of new accessible phenotypes and the number of mutations leading to them. This allows perform spearman's correlation on them.
// g++ -Wall Spearmyphenotype.cc libs/horse.o libs/linear.o libs/finder.o -lgsl -lgslcblas -lm
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

int main(int argc, char **argv)
{
    system("date");
    time_t rawtime;
    time ( &rawtime );

    string arch1 = argv[1];//Ruta
    string arch2 = argv[2];//1, 2 ó 3 (dup, c1 ó c2)
    string arch3 = argv[3];//dup, c1 ó c2

    int muestra = 1000;
    int mut = 432;
    int Nf = 12;
    int Ns = 6;
    int cuenta1, cuenta2, per_antes, per_despues, evento, cuentavector, cuentavector2, cuentavector3, sum1, sum2;

    string basura;

    Horse hor;
    Finder find;

    ifstream fi;
    ofstream fo;

    double dist;
    int hils_spear, hils_fn_antes;
    int **Fn_antes, *Per_antes, **Fn_despues, *Per_despues, **Mat_spear, **Mat_numuts_antes, **Mat_numuts_desp, **Mat_fn_desp, **Mat_conserv, **nuevoantes, **nuevodespues, *Vec_numuts_antes, *Vec_numuts_desp, *Vec_fn_antes;
    bool *boolean_antes, *boolean_desp, *listademuestras;

    hor.espacio(Mat_conserv, muestra, Nf);
    hor.espacio(Vec_fn_antes, muestra);
    hor.espacio(Mat_fn_desp, muestra, Nf);
    hor.espacio(Mat_numuts_antes, mut, muestra);
    hor.espacio(Mat_numuts_desp, mut, Nf);
    hor.espacio(listademuestras, muestra);
    hor.fillv0(listademuestras, muestra);
    
    int *eventos, **mat_di;
    hor.espacio(eventos, muestra);
    hor.espacio(mat_di, muestra, Nf);
    
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_"+ arch3 +".txt");
    for(int i = 0; i < muestra; i++){
        fi >> eventos[i];
    }
    fi.close();
    
    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch3 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> mat_di[i][j];
        }
    }
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

    hor.open_ifstream(fi, ""+ arch1 +"/acces/valordek/valordek_esp.txt");
    for(int k = 0; k < mut; k++){
        for(int j = 0; j < muestra; j++){
            fi >> Mat_numuts_antes[k][j];
        }
    }
    fi.close();

    for(int i = 0; i < muestra; i++){
        cout << "muestra: " << i << endl;
        if(mat_di[i][eventos[i]] != 0){
            cout << "muestra: " << i << " evento: " << eventos[i] << endl;
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + ".txt");
            cuenta1 = 0;
            while(fi >> basura)
                cuenta1++;
            fi.close();
                        
            hils_fn_antes = cuenta1/Ns;
            hor.espacio(Fn_antes, hils_fn_antes, Ns);
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + ".txt");
            for(int j = 0; j < cuenta1/Ns; j++){
                for(int k = 0; k < Ns; k++){
                    fi >> Fn_antes[j][k];
                }
            }
            fi.close();
                        
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/Per_" + hor.inttostring(i) + ".txt");
            per_antes = 0;
            while(fi >> basura)
                per_antes++;
            fi.close();
            
            hor.espacio(boolean_antes, per_antes);
            hor.fillv0(boolean_antes, per_antes);
                        
            hor.espacio(Per_antes, per_antes);
            hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/Per_" + hor.inttostring(i) + ".txt");
            for(int j = 0; j < per_antes; j++){
                fi >> Per_antes[j];
            }
            fi.close();
            
            hor.open_ifstream(fi, ""+ arch1 +"/acces/valordek/valordek_"+ arch3 +"_" + hor.inttostring(i) + "de1000.txt");
            for(int k = 0; k < mut; k++){
                for(int j = 0; j < Nf; j++){
                    fi >> Mat_numuts_desp[k][j];
                }
            }
            fi.close();
                        
            hor.espacio(Vec_numuts_antes, per_antes);
            cuentavector = 0;
            for(int k = 0; k < mut; k++){
                if(Mat_numuts_antes[k][i] != 0){
                    Vec_numuts_antes[cuentavector] = Mat_numuts_antes[k][i];
                    cuentavector++;
                }
            }
                        
            evento = eventos[i];
            
            if(Mat_fn_desp[i][evento] != 0){//Comprueva que sí haya fenotipos nuevos generados después
                hor.open_ifstream(fi, ""+ arch1 +"/acces/dist/dist_fn_esp_"+ arch3 +".txt");
                for(int i = 0; i < muestra; i++){
                    for(int j = 0; j < Nf; j++){
                        fi >> Mat_conserv[i][j];
                    }
                }
                fi.close();
                
                hils_spear = Vec_fn_antes[i] + Mat_fn_desp[i][evento] - Mat_conserv[i][evento];
                
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/speardata/"+ arch3 +"/fns_evento_"+ arch3 +".txt");
                fo << Vec_fn_antes[i] << " " << Mat_fn_desp[i][evento] << " " << Mat_conserv[i][evento] << " " << hils_spear << endl;
                fo.close();
                
                hor.espacio(Mat_spear, hils_spear, 2);
                hor.fillm0(Mat_spear, hils_spear, 2);
                
                hor.espacio(Vec_numuts_desp, Mat_fn_desp[i][evento]);
                for(int hils = 0; hils < Mat_fn_desp[i][evento]; hils++){
                    if(Mat_numuts_desp[hils][evento] != 0){Vec_numuts_desp[hils] = Mat_numuts_desp[hils][evento];}
                }

                cuenta2 = 0;
                hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(evento) + ".txt");
                while(fi >> basura)
                    cuenta2++;
                fi.close();
                hor.espacio(Fn_despues, cuenta2/Ns, Ns);
                
                hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/fn_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(evento) + ".txt");
                for(int j = 0; j < cuenta2/Ns; j++){
                    for(int k = 0; k < Ns; k++){
                        fi >> Fn_despues[j][k];
                    }
                }
                fi.close();

                hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(evento) + ".txt");
                per_despues = 0;
                while(fi >> basura)
                    per_despues++;
                fi.close();

                hor.espacio(Per_despues, per_despues);
                
                hor.open_ifstream(fi, ""+ arch1 +"/acces/fn_norepetido/Per_" + hor.inttostring(i) + "_"+ arch2 +"_" + hor.inttostring(evento) + ".txt");
                for(int j = 0; j < per_despues; j++){
                    fi >> Per_despues[j];
                }
                fi.close();
                
                hor.espacio(boolean_desp, per_despues);
                hor.fillv0(boolean_desp, per_despues);
                                
                cuentavector = 0;
                for(int k = 0; k < mut; k++){
                    if(Mat_numuts_desp[k][evento] != 0){
                        Vec_numuts_desp[cuentavector] = Mat_numuts_desp[k][evento];
                        cuentavector++;
                    }
                }
                //////////////////////////Conservados//////////////////////////////////////////////
                sum1 = 0;
                cuentavector = 0; //Para ver cuántos fns antes y después son parecidos
                for(int m = 0; m < per_antes; m++){ //Para todos los fns del espécimen i en curso
                    
                    cout << "4" << endl;
                    
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
                            //                        cout << "1.5" << endl;
                            dist = find.distfenotipica(nuevoantes, Per_antes[m], nuevodespues, Per_despues[n], Ns);
                            if(dist == 0){//Si encontramos un fenotipo igual antes y después
                                boolean_antes[m] = true;
                                boolean_desp[n] = true;
                                Mat_spear[cuentavector][0] = Vec_numuts_antes[m];
                                Mat_spear[cuentavector][1] = Vec_numuts_desp[n];
                                listademuestras[i] = true;
                                
                                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/speardata/Cons/"+ arch3 +"/mat_spear_"+ arch3 +"_"+ hor.inttostring(i) +"_cons.txt");
                                fo << Mat_spear[cuentavector][0] << " " << Mat_spear[cuentavector][1] << endl;
                                fo.close();
                                cuentavector++;
                                hor.borrar(nuevodespues, Per_despues[n]);//Para que borre el espacio antes de romper el ciclo
                                break;
                            }
                            hor.borrar(nuevodespues, Per_despues[n]);
                        }
                        sum2 += Per_despues[n];
                    }
                    sum1 += Per_antes[m];
                    hor.borrar(nuevoantes, Per_antes[m]);
                }
                ///////////////////////////////////////////////////////////////////////////7
                hor.borrar(Fn_despues, cuenta2/Ns);
                hor.borrar(Per_despues);
                                
                ////////////////////////////Llenar el resto de la matriz/////////////////////////////////////////////
                cuentavector2 = 0;
                for(int hils = 0; hils < per_antes; hils++){
                    if(boolean_antes[hils] != true){//Llenamos con #mut por fnuevo antes y con ceros después
                        Mat_spear[cuentavector + cuentavector2][0] = Vec_numuts_antes[hils];
                        Mat_spear[cuentavector + cuentavector2][1] = 0;
                        cuentavector2++;
                    }
                }
                                
                cuentavector3 = 0;
                for(int hils = 0; hils < per_despues; hils++){
                    if(boolean_desp[hils] == false){//Llenamos con 0 por fnuevo antes y con #mut después
                        Mat_spear[cuentavector + cuentavector2 + cuentavector3][0] = 0;
                        Mat_spear[cuentavector + cuentavector2 + cuentavector3][1] = Vec_numuts_desp[hils];
                        cuentavector3++;
                    }
                }
                
                ///////////////////////////////////////////////////////////////////////////////////
                
                //Sacamos matriz para después comparar con spearman
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/speardata/"+ arch3 +"/mat_spear_"+ arch3 +"_"+ hor.inttostring(i) +".txt");
                for(int hils = 0; hils < hils_spear; hils++){
                    fo << Mat_spear[hils][0] << " " << Mat_spear[hils][1] << endl;
                }
                fo.close();
                hor.borrar(Mat_spear, hils_spear);
                hor.borrar(Vec_numuts_desp);
                hor.borrar(boolean_desp);
                
            }
            else{cout << "No hay fenotipos nuevo qué comparar... muestra " << i << " evento " << evento << endl;}
            hor.borrar(Per_antes);
            hor.borrar(Fn_antes, hils_fn_antes);
            hor.borrar(boolean_antes);
            hor.borrar(Vec_numuts_antes);
        }
    }
        
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/speardata/Cons/"+ arch3 +"/cons_"+ arch3 +".txt");
    for(int i = 0; i < muestra; i++){
        if(listademuestras[i] == true){fo << i << endl;}
    }
    fo.close();
    
    hor.borrar(Vec_fn_antes);
    hor.borrar(Mat_fn_desp, muestra);
    hor.borrar(Mat_conserv, muestra);
    hor.borrar(Mat_numuts_antes, mut);
    hor.borrar(Mat_numuts_desp, mut);
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    return 0;
}
