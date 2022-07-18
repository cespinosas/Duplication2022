//Program to obtain the number of mutations leading to recurrently accessible phenotypes
// g++ -Wall mut_lead_rec.cc libs/horse.o -lgsl -lgslcblas -lm
// ./a.out /outs dup

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
    string arch3 = argv[2];//dup, c1 ó c2
    int muestra = 1000;
    int Nf = 12;
    
    Horse hor;
    
    ifstream fi, fi2, fi3;
    ofstream fo, fo2;
    int *Vec_antes_todo, *Vec_antes_cons, *Vec_antes_nocons, *Vec_desp_todo, *Vec_desp_cons, *Vec_desp_nocons, *Vec_eventos, **Mat_todo, **Mat_dif_igual, *Vec_orden_antes, *Vec_orden_desp, *Vec_ordinals_antes, *Vec_ordinals_desp;
    int cuenta, cuenta_antes_todo, cuenta_antes_cons, cuenta_antes_nocons, cuenta_desp_todo, cuenta_desp_cons, cuenta_desp_nocons, sum_antes_cons, sum_antes_nocons, sum_desp_cons, sum_desp_nocons, cuenta_ordinals, cuenta_ordinals2, cuenta_ordinals3, evento, dif_igual, cuenta_ceros;
    string basura;
    double **Mat_ordinals_antes, **Mat_ordinals_desp, **Mat_fn_sorted, *V1;
    double median, sum, min, max;
    double shannon_antes_todo, shannon_desp_todo, shannon_antes_cons, shannon_desp_cons, shannon_antes_nocons, shannon_desp_nocons;
    int *Vec_acces_esp, **Mat_acces_desp;
    int fntotal, frec, fnorec, muttotal, mutrec, mutnorec, suma_antes_todo, suma_desp_todo;
    int max_antes_cons, min_antes_cons, max_antes_nocons, min_antes_nocons, max_desp_cons, min_desp_cons, max_desp_nocons, min_desp_nocons, median_antes_cons, median_antes_nocons, median_desp_cons, median_desp_nocons;
    
    hor.espacio(Vec_eventos, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/eventos_"+ arch3 +".txt");
    for(int i = 0; i < muestra; i++){fi >> Vec_eventos[i]; /*cout << Vec_eventos[i] << endl;*/}
    fi.close();
    
    hor.espacio(Mat_dif_igual, muestra, Nf);
    hor.open_ifstream(fi, ""+ arch1 +"/"+ arch3 +"pasosdi.txt");
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            fi >> Mat_dif_igual[i][j];
        }
    }
    fi.close();
    
    //hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/Rwilcox_antes_"+ arch3 +".sh");
    //fo << "setwd(\""+ arch1 +"\")\n";
    //fo.close();
    
    //ANTES - MUT
    //Sum gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Sum tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Sum notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //DESP - MUT
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Sum tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Sum notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Min notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Max notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_gral.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_tol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    //Median notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_notol.txt");
    fo << "MFR MFNR" << endl;
    fo.close();
    
    //ANTES
    //Sum gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Sum tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Sum notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //DESPUÉS
    //Sum gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Sum tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Sum notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Min notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Max notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median gral
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_gral.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median tol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_tol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    //Median notol
    hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_notol.txt");
    fo << "FR FNR" << endl;
    fo.close();
    
//    cout << "1" << endl;
    hor.espacio(Vec_acces_esp, muestra);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/accesesp.txt");
    for(int j = 0; j < muestra; j++){
        fi >> Vec_acces_esp[j];
    }
    fi.close();
    hor.espacio(Mat_acces_desp, muestra, Nf);
    hor.open_ifstream(fi, ""+ arch1 +"/acces/acces"+ arch3 +".txt");
    for(int hils = 0; hils < muestra; hils++){
        for(int cols = 0; cols < Nf; cols++){
            fi >> Mat_acces_desp[hils][cols];
        }
    }
    fi.close();
    
    for(int i = 0; i < muestra; i++){
    //for(int i = 0; i < 1; i++){
        cout << "muestra: " << i << endl;
        evento = Vec_eventos[i];
        dif_igual = Mat_dif_igual[i][evento];
        
        cuenta = 0;
        hor.open_ifstream(fi, ""+ arch1 +"/acces/dist/speardata/"+ arch3 +"/mat_spear_"+ arch3 +"_"+ hor.inttostring(i) +".txt");
        while(fi >> basura)
            cuenta++;
        fi.close();
        
        hor.espacio(Mat_todo, cuenta/2, 2);
        
        //cout << "1" << endl;
        hor.open_ifstream(fi, ""+ arch1 +"/acces/dist/speardata/"+ arch3 +"/mat_spear_"+ arch3 +"_"+ hor.inttostring(i) +".txt");
        for(int j = 0; j < (cuenta/2); j++){
            fi >> Mat_todo[j][0] >> Mat_todo[j][1];
            //cout << Mat_todo[j][0] << " " << Mat_todo[j][1] << endl;
            //Antes num mut - Desp num mut por fenotipo nuevo
        }
        fi.close();
        //cout << endl;
        ///////////////////////////////////////////////////////
        cuenta_antes_todo = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][0] != 0){cuenta_antes_todo++;}
        }
        
        hor.espacio(Vec_antes_todo, cuenta_antes_todo);
        
        //cout << "Vec_antes_todo" << endl;
        suma_antes_todo = 0;
        cuenta_antes_todo = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][0] != 0){
                Vec_antes_todo[cuenta_antes_todo] = Mat_todo[j][0];
//                cout << Vec_antes_todo[cuenta_antes_todo] << endl;
                suma_antes_todo += Vec_antes_todo[cuenta_antes_todo];
                cuenta_antes_todo++;
            }
        }
        if(cuenta_antes_todo > 0){shannon_antes_todo = hor.sh_index(Vec_antes_todo, cuenta_antes_todo);}
//        cout << endl;
//        cout << suma_antes_todo << endl;
        //////////////////////////////////////////////////////
        cuenta_desp_todo = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][1] != 0){cuenta_desp_todo++;}
        }
        
        hor.espacio(Vec_desp_todo, cuenta_desp_todo);
        
        //cout << "Vec_desp_todo" << endl;
        suma_desp_todo = 0;
        cuenta_desp_todo = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][1] != 0){
                Vec_desp_todo[cuenta_desp_todo] = Mat_todo[j][1];
                //cout << Vec_desp_todo[cuenta_desp_todo] << endl;
                suma_desp_todo += Vec_desp_todo[cuenta_desp_todo];
                cuenta_desp_todo++;
            }
        }
        if(cuenta_desp_todo > 0){shannon_desp_todo = hor.sh_index(Vec_desp_todo, cuenta_desp_todo);}
//        cout << endl;
        ///////////////////////////////////////////////////////
        cuenta_antes_cons = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][0] != 0){
                if(Mat_todo[j][1] != 0){cuenta_antes_cons++;}
            }
        }
        
        if(cuenta_antes_cons > 0){
            hor.espacio(Vec_antes_cons, cuenta_antes_cons);
            //         cout << "Vec_antes_cons" << endl;
            cuenta_antes_cons = 0;
            sum_antes_cons = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_todo[j][0] != 0){
                    if(Mat_todo[j][1] != 0){
                        Vec_antes_cons[cuenta_antes_cons] = Mat_todo[j][0];
                        sum_antes_cons += Mat_todo[j][0];
                        //                      cout << Vec_antes_cons[cuenta_antes_cons] << endl;
                        cuenta_antes_cons++;
                    }
                }
            }/*cout << endl;*/
            shannon_antes_cons = hor.sh_index(Vec_antes_cons, cuenta_antes_cons);
            min_antes_cons = hor.find_min(Vec_antes_cons, cuenta_antes_cons);
            max_antes_cons = hor.find_max(Vec_antes_cons, cuenta_antes_cons);
            median_antes_cons = hor.get_median(Vec_antes_cons, cuenta_antes_cons);
            hor.borrar(Vec_antes_cons);
        }
        
        //////////////////////////////////////////////////////
        cuenta_desp_cons = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][1] != 0){
                if(Mat_todo[j][0] != 0){cuenta_desp_cons++;}
            }
        }
        
        if(cuenta_desp_cons > 0){
            hor.espacio(Vec_desp_cons, cuenta_desp_cons);
            //         cout << "Vec_desp_cons" << endl;
            cuenta_desp_cons = 0;
            sum_desp_cons = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_todo[j][1] != 0){
                    if(Mat_todo[j][0] != 0){
                        Vec_desp_cons[cuenta_desp_cons] = Mat_todo[j][1];
                        sum_desp_cons += Mat_todo[j][1];
                        //                      cout << Vec_desp_cons[cuenta_desp_cons] << endl;
                        cuenta_desp_cons++;
                    }
                }
            }/*cout << endl;*/
            shannon_desp_cons = hor.sh_index(Vec_desp_cons, cuenta_desp_cons);
            min_desp_cons = hor.find_min(Vec_desp_cons, cuenta_desp_cons);
            max_desp_cons = hor.find_max(Vec_desp_cons, cuenta_desp_cons);
            median_desp_cons = hor.get_median(Vec_desp_cons, cuenta_desp_cons);
            hor.borrar(Vec_desp_cons);
        }
        //////////////////////////////////////////////////////
        cuenta_antes_nocons = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][0] != 0){
                if(Mat_todo[j][1] == 0){cuenta_antes_nocons++;}
            }
        }
        //cout << "2" << endl;
        if(cuenta_antes_nocons > 0){
            hor.espacio(Vec_antes_nocons, cuenta_antes_nocons);
            //          cout << "Vec_antes_nocons" << endl;
            cuenta_antes_nocons = 0;
            sum_antes_nocons = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_todo[j][0] != 0){
                    if(Mat_todo[j][1] == 0){
                        Vec_antes_nocons[cuenta_antes_nocons] = Mat_todo[j][0];
                        sum_antes_nocons += Mat_todo[j][0];
                        //                      cout << Vec_antes_nocons[cuenta_antes_nocons] << endl;
                        cuenta_antes_nocons++;
                    }
                }
            }/*cout << endl;*/
            shannon_antes_nocons = hor.sh_index(Vec_antes_nocons, cuenta_antes_nocons);
            min_antes_nocons = hor.find_min(Vec_antes_nocons, cuenta_antes_nocons);
            max_antes_nocons = hor.find_max(Vec_antes_nocons, cuenta_antes_nocons);
            median_antes_nocons = hor.get_median(Vec_antes_nocons, cuenta_antes_nocons);
            hor.borrar(Vec_antes_nocons);
        }
        //////////////////////////////////////////////////////
        cuenta_desp_nocons = 0;
        for(int j = 0; j < (cuenta/2); j++){
            if(Mat_todo[j][1] != 0){
                if(Mat_todo[j][0] == 0){cuenta_desp_nocons++;}
            }
        }
        
        if(cuenta_desp_nocons > 0){
            hor.espacio(Vec_desp_nocons, cuenta_desp_nocons);
            //          cout << "Vec_desp_nocons" << endl;
            cuenta_desp_nocons = 0;
            sum_desp_nocons = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_todo[j][1] != 0){
                    if(Mat_todo[j][0] == 0){
                        Vec_desp_nocons[cuenta_desp_nocons] = Mat_todo[j][1];
                        sum_desp_nocons += Mat_todo[j][1];
                        //                      cout << Vec_desp_nocons[cuenta_desp_nocons] << endl;
                        cuenta_desp_nocons++;
                    }
                }
            }/*cout << endl;*/
            shannon_desp_nocons = hor.sh_index(Vec_desp_nocons, cuenta_desp_nocons);
            min_desp_nocons = hor.find_min(Vec_desp_nocons, cuenta_desp_nocons);
            max_desp_nocons = hor.find_max(Vec_desp_nocons, cuenta_desp_nocons);
            median_desp_nocons = hor.get_median(Vec_desp_nocons, cuenta_desp_nocons);
            hor.borrar(Vec_desp_nocons);
        }
        
        //cout << "\n-----------Antes TODO\n";
        //for(int hils = 0; hils < cuenta_antes_todo; hils++){
        //cout << Vec_antes_todo[hils] << endl;
        //}
//        cout << "2" << endl;
        //-----------------------SORTING-------------------------
        //ANTES
        hor.espacio(Vec_orden_antes, cuenta_antes_todo);
        hor.espacio(Vec_ordinals_antes, cuenta_antes_todo);
        hor.sort(Vec_antes_todo, Vec_orden_antes, cuenta_antes_todo);
        //Volteo el orden de Vec_orden_desp para que los ordinales vayan de forma ascendente 090520
        //cout << "Ordinales" << endl;
        for(int hils = 0; hils < cuenta_antes_todo; hils++){
            //cout << Vec_orden_antes[hils] << endl;
        }
        
        cuenta_ordinals = hor.count_ordinals(Vec_orden_antes, cuenta_antes_todo);
        //cout << "cuenta_ordinalS: " << cuenta_ordinals << endl;
        //cout << "cuenta_antes_todo: " << cuenta_antes_todo << endl;
        //cout << endl;
        //cout << "cuenta_ordinals: " << cuenta_ordinals << endl;
        hor.espacio(Mat_ordinals_antes, cuenta_ordinals, 2);
        hor.fillm0(Mat_ordinals_antes, cuenta_ordinals, 2);
        
        cuenta_ordinals2 = 0;
        cuenta_ordinals3 = 0;
        for(int hils = 0; hils < (cuenta_antes_todo - 1); hils++){
            cuenta_ordinals3++;
            if(Vec_orden_antes[hils] != Vec_orden_antes[hils + 1]){
                cuenta_ordinals2++;
                Mat_ordinals_antes[cuenta_ordinals2 - 1][0] = (double)Vec_orden_antes[hils];
                if(cuenta_ordinals3 > 1){
                    Mat_ordinals_antes[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2 + 1/(double)cuenta_ordinals3;
                }
                else{
                    Mat_ordinals_antes[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2;
                }
                if(hils < (cuenta_antes_todo - 2)){cuenta_ordinals3 = 0;}
            }
        }
        
        Mat_ordinals_antes[cuenta_ordinals - 1][0] = (double)Vec_orden_antes[cuenta_antes_todo - 1];
        if(Vec_orden_antes[cuenta_antes_todo - 2] != Vec_orden_antes[cuenta_antes_todo - 1]){
            //El ùltimo nùmero de la lista es diferente
            Mat_ordinals_antes[cuenta_ordinals - 1][1] = (double)(cuenta_ordinals);
        }
        else{
            //El ùltimo nùmero de la lista es igual al o a los anteriores
            Mat_ordinals_antes[cuenta_ordinals - 1][1] = (double)cuenta_ordinals + 1/(double)(cuenta_ordinals3 + 1);
        }
        
        hor.pantalla(Mat_ordinals_antes, cuenta_ordinals, 2);
        
        hor.borrar(Vec_antes_todo);
        hor.borrar(Vec_orden_antes);
        hor.borrar(Vec_ordinals_antes);
        
        //Unificamos Ordinales
        //cout << "cuenta/2: " << cuenta/2 << endl;
        hor.espacio(Mat_fn_sorted, 4, cuenta/2);
        hor.fillm0(Mat_fn_sorted, 4, cuenta/2);
        //cout << "3" << endl;
        //cout << "Mat sorted antes" << endl;
        for(int hils = 0; hils < cuenta/2; hils++){
            Mat_fn_sorted[0][hils] = Mat_todo[hils][0];
            for(int hils2 = 0; hils2 < cuenta_ordinals; hils2++){
                if(Mat_fn_sorted[0][hils] == Mat_ordinals_antes[hils2][0]){
                    Mat_fn_sorted[1][hils] = Mat_ordinals_antes[hils2][1];
                    break;
                }
            }
        }
        
        hor.borrar(Mat_ordinals_antes, cuenta_ordinals);
//        cout << "3" << endl;
        //DESPUÈS
        hor.espacio(Vec_orden_desp, cuenta_desp_todo);
        hor.espacio(Vec_ordinals_desp, cuenta_desp_todo);
        hor.sort(Vec_desp_todo, Vec_orden_desp, cuenta_desp_todo);
        
        cuenta_ordinals = hor.count_ordinals(Vec_orden_desp, cuenta_desp_todo);
        
        hor.espacio(Mat_ordinals_desp, cuenta_ordinals, 2);
        hor.fillm0(Mat_ordinals_desp, cuenta_ordinals, 2);
        //Mat_ordinals_desp[][0] -> ordinales de Mat_todo[][1], mat_spears segunda fila de lista de mut por fn (desp)
        //Mat_ordinals_desp[][1] ->
        cuenta_ordinals2 = 0;
        cuenta_ordinals3 = 0;
        for(int hils = 0; hils < (cuenta_desp_todo - 1); hils++){
            cuenta_ordinals3++;
            if(Vec_orden_desp[hils] != Vec_orden_desp[hils + 1]){
                cuenta_ordinals2++;
                Mat_ordinals_desp[cuenta_ordinals2 - 1][0] = (double)Vec_orden_desp[hils];
                if(cuenta_ordinals3 > 1){
                    Mat_ordinals_desp[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2 + 1/(double)cuenta_ordinals3;
                }
                else{
                    Mat_ordinals_desp[cuenta_ordinals2 - 1][1] = (double)cuenta_ordinals2;
                }
                if(hils < (cuenta_desp_todo - 2)){cuenta_ordinals3 = 0;}
            }
        }
        
        Mat_ordinals_desp[cuenta_ordinals - 1][0] = (double)Vec_orden_desp[cuenta_desp_todo - 1];
        if(Vec_orden_desp[cuenta_desp_todo - 2] != Vec_orden_desp[cuenta_desp_todo - 1]){
            //Los dos ùltimos nùmeros son diferentes
            Mat_ordinals_desp[cuenta_ordinals - 1][1] = (double)(cuenta_ordinals);
        }
        else{
            //El ùltimo nùmero de la lista es igual al o a los anteriores
            Mat_ordinals_desp[cuenta_ordinals - 1][1] = (double)cuenta_ordinals + 1/(double)(cuenta_ordinals3 + 1);
        }
        
        hor.borrar(Vec_desp_todo);
        hor.borrar(Vec_orden_desp);
        hor.borrar(Vec_ordinals_desp);
        
        //Unificamos Ordinales
        for(int hils = 0; hils < cuenta/2; hils++){
            Mat_fn_sorted[2][hils] = Mat_todo[hils][1];
            for(int hils2 = 0; hils2 < cuenta_ordinals; hils2++){
                if(Mat_fn_sorted[2][hils] == Mat_ordinals_desp[hils2][0]){
                    Mat_fn_sorted[3][hils] = Mat_ordinals_desp[hils2][1];
                    break;
                }
            }
        }
        hor.pantalla(Mat_fn_sorted, 4, cuenta/2);
        hor.borrar(Mat_ordinals_desp, cuenta_ordinals);
//        cout << "4" << endl;
        //SUM(Wilcoxon), MAX, MIN y MEDIANA de los valores cons y no cons
        
        //Max, min y mediana cons antes
        cuenta_ceros = 0;
        if(cuenta_antes_cons > 0){
            hor.espacio(V1, cuenta_antes_cons);
            hor.fillv0(V1, cuenta_antes_cons);
            cuenta_antes_cons = 0;
            sum = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_fn_sorted[0][j] != 0){
                    if(Mat_fn_sorted[2][j] != 0){
                        V1[cuenta_antes_cons] = Mat_fn_sorted[1][j];
                        sum += Mat_fn_sorted[1][j];
                        cuenta_antes_cons++;
                    }
                }
            }
            
            min = hor.find_min(V1, cuenta_antes_cons);
            max = hor.find_max(V1, cuenta_antes_cons);
            median = hor.get_median(V1, cuenta_antes_cons);
            
            hor.borrar(V1);
            
            fntotal = Vec_acces_esp[i];
            frec = cuenta_antes_cons;
            muttotal = suma_antes_todo;
            mutrec = sum_antes_cons;
//            cout << "5" << endl;
            //Todo
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_all.txt");
            fo << i << " " << evento << " " << dif_igual << " " << fntotal <<  " " << frec << " " << muttotal << " " << mutrec << " " <<  sum << " " << min << " " << max << " " << median << " " << min_antes_cons << " " << max_antes_cons << " " << median_antes_cons << " ";
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << i << " " << evento << " " << dif_igual << " " << shannon_antes_todo << " " << shannon_antes_cons << " ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_gral.txt");
            fo << mutrec << " ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_tol.txt");
                fo << mutrec << " ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_notol.txt");
                fo << mutrec << " ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_gral.txt");
            fo << min_antes_cons << " ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_tol.txt");
                fo << min_antes_cons << " ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_notol.txt");
                fo << min_antes_cons << " ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_gral.txt");
            fo << max_antes_cons << " ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_tol.txt");
                fo << max_antes_cons << " ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
                fo << max_antes_cons << " ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_gral.txt");
            fo << median_antes_cons << " ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_tol.txt");
                fo << median_antes_cons << " ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_notol.txt");
                fo << median_antes_cons << " ";
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_gral.txt");
            fo << sum << " ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_tol.txt");
                fo << sum << " ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_notol.txt");
                fo << sum << " ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_gral.txt");
            fo << min << " ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_tol.txt");
                fo << min << " ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_notol.txt");
                fo << min << " ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_gral.txt");
            fo << max << " ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_tol.txt");
                fo << max << " ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
                fo << max << " ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_gral.txt");
            fo << median << " ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_tol.txt");
                fo << median << " ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_notol.txt");
                fo << median << " ";
                fo.close();
            }
        }
        else{
            cuenta_ceros++;
            fntotal = Vec_acces_esp[i];
            muttotal = suma_antes_todo;
            //Todo
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_all.txt");
            fo << i << " " << evento << " " << dif_igual << " " << fntotal <<  " 0 " << muttotal << " 0 0 0 0 0 0 0 0 ";
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << i << " " << evento << " " << dif_igual << " 0 0 ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_gral.txt");
            fo << "0 ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_gral.txt");
            fo << "0 ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_gral.txt");
            fo << "0 ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_gral.txt");
            fo << "0 ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_gral.txt");
            fo << "0 ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_gral.txt");
            fo << "0 ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_gral.txt");
            fo << "0 ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_gral.txt");
            fo << "0 ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_notol.txt");
                fo << "0 ";
                fo.close();
            }
        }
//        cout << "Cuenta_ceros antes cons: " << cuenta_ceros << endl;
        //Max, min y mediana nocons antes
        cuenta_ceros = 0;
        if(cuenta_antes_nocons > 0){
            hor.espacio(V1, cuenta_antes_nocons);
            hor.fillv0(V1, cuenta_antes_nocons);
            cuenta_antes_nocons = 0;
            sum = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_fn_sorted[0][j] != 0){
                    if(Mat_fn_sorted[2][j] == 0){
                        V1[cuenta_antes_nocons] = Mat_fn_sorted[1][j];
                        sum += Mat_fn_sorted[1][j];
                        cuenta_antes_nocons++;
                    }
                }
            }
            
            min = hor.find_min(V1, cuenta_antes_nocons);
            max = hor.find_max(V1, cuenta_antes_nocons);
            median = hor.get_median(V1, cuenta_antes_nocons);
            
            hor.borrar(V1);
            
            fnorec = cuenta_antes_nocons;
            mutnorec = sum_antes_nocons;
            
            //TODO
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_all.txt");
            fo << fnorec << " " << mutnorec << " " << sum << " " << min << " " << max << " " << median << " " << min_antes_nocons << " " << max_antes_nocons << " " << median_antes_nocons << endl;
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << shannon_antes_nocons << " ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_gral.txt");
            fo << mutnorec << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_tol.txt");
                fo << mutnorec << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_notol.txt");
                fo << mutnorec << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_gral.txt");
            fo << min_antes_nocons << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_tol.txt");
                fo << min_antes_nocons << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_notol.txt");
                fo << min_antes_nocons << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_gral.txt");
            fo << max_antes_nocons << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_tol.txt");
                fo << max_antes_nocons << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
                fo << max_antes_nocons << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_gral.txt");
            fo << median_antes_nocons << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_tol.txt");
                fo << median_antes_nocons << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_notol.txt");
                fo << median_antes_nocons << endl;
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_gral.txt");
            fo << sum << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_tol.txt");
                fo << sum << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_notol.txt");
                fo << sum << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_gral.txt");
            fo << min << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_tol.txt");
                fo << min << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_notol.txt");
                fo << min << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_gral.txt");
            fo << max << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_tol.txt");
                fo << max << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
                fo << max << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_gral.txt");
            fo << median << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_tol.txt");
                fo << median << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_notol.txt");
                fo << median << endl;
                fo.close();
            }
        }
        else{ cuenta_ceros++;
            //Todo
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_all.txt");
            fo << "0 0 0 0 0 0 0 0 0" << endl;
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << "0 ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_sum_notol.txt");
                fo << "0"<< endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_min_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_tol.txt");
                fo << "0"<< endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_median_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_sum_notol.txt");
                fo << "0"<< endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_min_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_tol.txt");
                fo << "0"<< endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_median_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
        }
//        cout << "Cuenta_ceros antes nocons: " << cuenta_ceros << endl;
        //Max, min y mediana cons despuès
        cuenta_ceros = 0;
        if(cuenta_desp_cons > 0){
            hor.espacio(V1, cuenta_desp_cons);
            hor.fillv0(V1, cuenta_desp_cons);
            cuenta_desp_cons = 0;
            sum = 0;
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_fn_sorted[2][j] != 0){
                    if(Mat_fn_sorted[0][j] != 0){
                        V1[cuenta_desp_cons] = Mat_fn_sorted[3][j];
                        sum += Mat_fn_sorted[3][j];
                        cuenta_desp_cons++;
                    }
                }
            }
            
            min = hor.find_min(V1, cuenta_desp_cons);
            max = hor.find_max(V1, cuenta_desp_cons);
            median = hor.get_median(V1, cuenta_desp_cons);
            
            hor.borrar(V1);
            
            fntotal = Mat_acces_desp[i][evento];
            frec = cuenta_desp_cons;
            muttotal = suma_desp_todo;
            mutrec = sum_desp_cons;
            //TODO
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_all.txt");
            fo << i << " " << evento << " " << dif_igual << " " << fntotal <<  " " << frec << " " << muttotal << " " << mutrec << " " << sum << " " << min << " " << max << " " << median << " " << min_desp_cons << " " << max_desp_cons << " " << median_desp_cons << " ";
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << shannon_desp_todo << " " << shannon_desp_cons << " ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_gral.txt");
            fo << mutrec << " ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_tol.txt");
                fo << mutrec << " ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_notol.txt");
                fo << mutrec << " ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_gral.txt");
            fo << min_desp_cons << " ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_tol.txt");
                fo << min_desp_cons << " ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_notol.txt");
                fo << min_desp_cons << " ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_gral.txt");
            fo << max_desp_cons << " ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_tol.txt");
                fo << max_desp_cons << " ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_notol.txt");
                fo << max_desp_cons << " ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_gral.txt");
            fo << median_desp_cons << " ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_tol.txt");
                fo << median_desp_cons << " ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_notol.txt");
                fo << median_desp_cons << " ";
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_gral.txt");
            fo << sum << " ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_tol.txt");
                fo << sum << " ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_notol.txt");
                fo << sum << " ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_gral.txt");
            fo << min << " ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_tol.txt");
                fo << min << " ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_notol.txt");
                fo << min << " ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_gral.txt");
            fo << max << " ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_tol.txt");
                fo << max << " ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_notol.txt");
                fo << max << " ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_gral.txt");
            fo << median << " ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_tol.txt");
                fo << median << " ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_notol.txt");
                fo << median << " ";
                fo.close();
            }
        }
        else{
            cuenta_ceros++;
            fntotal = Mat_acces_desp[i][evento];
            muttotal = suma_desp_todo;
            //Todo
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_all.txt");
            fo << i << " " << evento << " " << dif_igual << " " << fntotal <<  " 0 " << muttotal << " 0 0 0 0 0 0 0 0 ";
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << "0 0 ";
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_gral.txt");
            fo << "0 ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_gral.txt");
            fo << "0 ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_gral.txt");
            fo << "0 ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_gral.txt");
            fo << "0 ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_gral.txt");
            fo << "0 ";
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_gral.txt");
            fo << "0 ";
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_gral.txt");
            fo << "0 ";
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Max notol         REVISAR ESTE ERROR
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_notol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_gral.txt");
            fo << "0 ";
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_tol.txt");
                fo << "0 ";
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_notol.txt");
                fo << "0 ";
                fo.close();
            }
        }
//        cout << "Cuenta_ceros desp cons: " << cuenta_ceros << endl;
        
        //Max, min y mediana nocons despuès
        cuenta_ceros = 0;
        if(cuenta_desp_nocons > 0){
            hor.espacio(V1, cuenta_desp_nocons);
            hor.fillv0(V1, cuenta_desp_nocons);
            cuenta_desp_nocons = 0;
            sum = 0;
            
            for(int j = 0; j < (cuenta/2); j++){
                if(Mat_fn_sorted[2][j] != 0){
                    if(Mat_fn_sorted[0][j] == 0){
                        V1[cuenta_desp_nocons] = Mat_fn_sorted[3][j];
                        sum += Mat_fn_sorted[3][j];
                        cuenta_desp_nocons++;
                    }
                }
            }
            
            min = hor.find_min(V1, cuenta_desp_nocons);
            max = hor.find_max(V1, cuenta_desp_nocons);
            median = hor.get_median(V1, cuenta_desp_nocons);
            
            hor.borrar(V1);
            
            fnorec = cuenta_desp_nocons;
            mutnorec = sum_desp_nocons;
            
            //TODO
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_all.txt");
            fo << fnorec << " " << mutnorec << " " << sum << " " << min << " " << max << " " << median << " " << min_desp_nocons << " " << max_desp_nocons << " " << median_desp_nocons << endl;
            fo.close();
            //Shannon
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << shannon_desp_nocons << endl;
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_gral.txt");
            fo << mutnorec << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_tol.txt");
                fo << mutnorec << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_notol.txt");
                fo << mutnorec << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_gral.txt");
            fo << min_desp_nocons << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_tol.txt");
                fo << min_desp_nocons << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_notol.txt");
                fo << min_desp_nocons << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_gral.txt");
            fo << max_desp_nocons << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_tol.txt");
                fo << max_desp_nocons << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_notol.txt");
                fo << max_desp_nocons << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_gral.txt");
            fo << median_desp_nocons << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_tol.txt");
                fo << median_desp_nocons << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_notol.txt");
                fo << median_desp_nocons << endl;
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_gral.txt");
            fo << sum << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_tol.txt");
                fo << sum << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_notol.txt");
                fo << sum << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_gral.txt");
            fo << min << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_tol.txt");
                fo << min << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_notol.txt");
                fo << min << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_gral.txt");
            fo << max << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_tol.txt");
                fo << max << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_notol.txt");
                fo << max << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_gral.txt");
            fo << median << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_tol.txt");
                fo << median << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_notol.txt");
                fo << median << endl;
                fo.close();
            }
        }
        else{
            cuenta_ceros++;
            //Todo
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_all.txt");
            fo << "0 0 0 0 0 0 0 0 0" << endl;
            fo.close();
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/shannon_"+ arch3 +".txt");
            fo << "0" << endl;
            fo.close();
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_sum_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_min_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_max_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_antes_"+ arch3 +"_max_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_mut_desp_"+ arch3 +"_median_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //###########################################################################################
            //Sum gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Sum tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Sum notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_sum_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Min tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Min notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_min_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Max tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_max_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Max notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_antes_"+ arch3 +"_max_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median gral
            hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_gral.txt");
            fo << "0" << endl;
            fo.close();
            //Median tol
            if(dif_igual == 1){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_tol.txt");
                fo << "0" << endl;
                fo.close();
            }
            //Median notol
            if(dif_igual == 0){
                hor.open_ofstream(fo, ""+ arch1 +"/acces/dist/Stats_"+ arch3 +"/stats_desp_"+ arch3 +"_median_notol.txt");
                fo << "0" << endl;
                fo.close();
            }
        }
        //        cout << "Cuenta_ceros desp nocons: " << cuenta_ceros << endl;
        
        hor.borrar(Mat_fn_sorted, 4);
        hor.borrar(Mat_todo, cuenta/2);
        
    }
    
    hor.borrar(Mat_dif_igual, muestra);
    hor.borrar(Vec_eventos);
    hor.borrar(Vec_acces_esp);
    hor.borrar(Mat_acces_desp, muestra);
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    return 0;
}
