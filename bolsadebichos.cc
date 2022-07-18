//Program to obtein specimenes

// g++ -Wall bolsadebichos.cc libs/horse.o libs/linear.o libs/finder.o libs/gataca.o -lgsl -lgslcblas -lm
// ./a.out 1 1000 12 6

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
    //Datos a introucir:
    double intnw = 0.25; // Interconectividad
    int Nf = atoi(argv[3]); //Número de genes de transcripción
    int Np = atoi(argv[4]); //Número de genes marcadores del fenotipo
    int semilla = atoi(argv[1]); //Para números aleatorios
    int tolerancia = 2; //Cuánto puede variar el número de Conexiones
    int muestra = atoi(argv[2]); //Número de redes de genotipos a estudiar
    int caminata = 20 * ((Nf * Nf) + (Nf * Np)); // Número de redes a explorar antes de llegar al genotipo que será parte de la muestra
    int explorar = 2 * ((Nf * Nf) + (Nf * Np)); //Mutaciones para medir de la muestra
    //Lo anterior significa que se tomarán 3 muestras separadas por 10 pasos y se mutarán 3 veces, se les hará una duplicación y se mutarán 3 veces (para cada gen).
    
    int p, perd, Cp, Cp2; //Periodo del fenotipo
    int n, nd, Cn, Cn2; //Contador para lin.dot
    int pasos, standby;
    
    Horse hor(semilla);
    Finder find;
    Linear lin(semilla);
    GATACA gat(semilla);
    
    hor.start_rng(semilla);
    //Vectores y matrices para crear red de determiado fenotipo fenotipo
    
    double **Matrizf; //Matriz de red de genes factores de transcripción
    hor.espacio(Matrizf, Nf, Nf);
    double **Matrizp; //Matriz de red de genes reguladores del fenotipo
    hor.espacio(Matrizp, Np, Nf);
    int **trayectoria; //Trayectoria de condición inicial a atractor
    hor.espacio(trayectoria, 1000, Nf);
    int *condinicial;
    hor.espacio(condinicial, Nf);
    int *vectorinicio;
    hor.espacio(vectorinicio, Nf);
    int *vectorentrada;
    hor.espacio(vectorentrada, Nf);
    int *atractor;
    hor.espacio(atractor, Nf);
    int **Fenotipo; //Fenotipo no porque debe ser flexible en el periodo
    int *Wish; //Matriz del fenotipo deseado
    double *comp1; //Compara muestras con la primera y sólo se almacena el mayor valor. Esto para ver que no estemos dando vuelta en un solo lugar.
    hor.espacio(comp1, 2);
    double *comp2;//Compara magnitudes de las regulaciones de las muestras.
    hor.espacio(comp2, 2);
    double **Matriz1f; //Matriz de red de genes factores de transcripción
    hor.espacio(Matriz1f, Nf, Nf);
    double **Matriz1p; //Matriz de red de genes reguladores del fenotipo
    hor.espacio(Matriz1p, Np, Nf);
    
    //Vectores y matrices de la red con una duplicación
    
    double **Matrizfd; //Matriz de red de genes factores de transcripción duplicada
    hor.espacio(Matrizfd, Nf + 1, Nf + 1);
    double **Matrizpd; //Matriz de red de genes reguladores del fenotipo duplicada
    hor.espacio(Matrizpd, Np, Nf + 1);
    int **trayectoriad;
    hor.espacio(trayectoriad, 1000, Nf + 1);
    int *vectoriniciod;
    hor.espacio(vectoriniciod, Nf + 1);
    int *vectorentradad;
    hor.espacio(vectorentradad, Nf + 1);
    int **Fenotipod; //Fenotipo no porque debe ser flexible en el periodo
    
    //////////Análisis//////////
    
    //Para guardar datos de la caminata
    double *pasosfuera; //Vector que guarda el promedio de los pasos dados fuera del espacio durante la caminata por muestra (standby después de la caminata / pasos totales de la caminata).
    hor.espacio(pasosfuera, muestra);
    
    //Para guardar datos de los vecinos de la muestra
    double *fdiferente; //Guardamos el promedio de vecinos con fenotipo diferente por muestra.
    hor.espacio(fdiferente, muestra);
    int fdif; //Para llevar la cuenta de vecinos diferentes por muestra.
    double *Distvecinos; //Guardamos las distancias entre el fenotipo de la muestra y de las mutantes
    hor.espacio(Distvecinos, explorar);
    double *Mediavecinos; //Promedio de distancias por muestra
    hor.espacio(Mediavecinos, muestra);
    double *Desvst; //Desviación estándar de las distancias por muestra
    hor.espacio(Desvst, muestra);
    
    //Datos de redes con gen duplicado
    //Cambio en el fenotipo por el hecho de duplicar
    double *Distdup; //Distancias vs el f original de los duplicados
    hor.espacio(Distdup, Nf); //Hay Nf duplicados por muestra
    double *Mediadup; //Promedio de distancias del fenotipo resultado de duplicar vs el f original
    hor.espacio(Mediadup, muestra);
    double *Desvstdup; //Desviación estándar de las distancias de las resdes duplicadas por muestra
    hor.espacio(Desvstdup, muestra);
    double *fdifdup; //Cuántas redes duplicados dieron fenotipo diferente
    hor.espacio(fdifdup, muestra);
    int fddup; //Para la suma de f diferentes
    
    //Vecinos de las redes duplicadas vs el f original
    double *Distvecdup; //Distancias de los vecinos por red duplicada
    hor.espacio(Distvecdup, Nf * explorar);
    double **Mediavecdup; //Medias de vecinos por red duplicada
    hor.espacio(Mediavecdup, muestra, Nf);
    double **Desvstvecdup; //Desviación estándar por duplicado
    hor.espacio(Desvstvecdup, muestra, Nf);
    double **fdifvecdup; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red duplicada
    hor.espacio(fdifvecdup, muestra, Nf);
    int fdvecdup;
    
    //Vecinos de las redes duplicadas vs el f resultante de duplicar
    double *Distvecdupin; //Distancias de los vecinos (mutadas) de las redes duplicadas con la duplicada
    //     hor.espacio(Distvecdupin, Nf * explorar);
    hor.espacio(Distvecdupin, explorar);
    double **Mediavecdupin; //Medias de vecinos por red duplicada
    hor.espacio(Mediavecdupin, muestra, Nf);
    double **Desvstvecdupin; //Desviación estándar por duplicado
    hor.espacio(Desvstvecdupin, muestra, Nf);
    double **fdifvecdupin; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red duplicada
    hor.espacio(fdifvecdupin, muestra, Nf);
    int fdvecdupin;//Vecinos con fenotipo diferente al del duplicado
//     int fddori; //Redes con duplicado que tienen fenotipo diferente al original
    
    //Control 1
    double **CMatrizf; //Matriz de red de genes factores de transcripción extendida
    hor.espacio(CMatrizf, Nf + 1, Nf + 1);
    double **CMatrizp; //Matriz de red de genes reguladores del fenotipo extendida
    hor.espacio(CMatrizp, Np, Nf + 1);
    int **Ctrayectoria;
    hor.espacio(Ctrayectoria, 1000, Nf + 1);
    int *Cvectorinicio;
    hor.espacio(Cvectorinicio, Nf + 1);
    int *Cvectorentrada;
    hor.espacio(Cvectorentrada, Nf + 1);
    int **CFenotipo; //Fenotipo no porque debe ser flexible en el periodo
    //Cambio en el fenotipo por el hecho de expandir
    double *CDist; //Distancias vs el f original de los extendidos
    hor.espacio(CDist, Nf); //Hay Nf extendidos por muestra
    double *CMedia; //Promedio de distancias del fenotipo resultado de extenderse vs el f original
    hor.espacio(CMedia, muestra);
    double *CDesvst; //Desviación estándar de las distancias de las resdes extendidas por muestra
    hor.espacio(CDesvst, muestra);
    double *Cfdif; //Cuántas redes extendidas dieron fenotipo diferente
    hor.espacio(Cfdif, muestra);
    int Cfd; //Para la suma de f diferentes
    //Vecinos de las redes extendidas vs el f original
    double *CDistvec; //Distancias de los vecinos por red extendida
    hor.espacio(CDistvec, Nf * explorar);
    double **CMediavec; //Medias de vecinos por red extendida
    hor.espacio(CMediavec, muestra, Nf);
    double **CDesvstvec; //Desviación estándar por red extendida
    hor.espacio(CDesvstvec, muestra, Nf);
    double **Cfdifvec; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red extendida
    hor.espacio(Cfdifvec, muestra, Nf);
    int Cfdvec;
    //Vecinos de las redes extendidas vs el f resultante de extender
    double *CDistvecin; //Distancias de los vecinos (mutadas) de las redes extendidas con la extendida
    //     hor.espacio(CDistvecin, Nf * explorar);
    hor.espacio(CDistvecin, explorar);
    double **CMediavecin; //Medias de vecinos por red extendida
    hor.espacio(CMediavecin, muestra, Nf);
    double **CDesvstvecin; //Desviación estándar por red extendida
    hor.espacio(CDesvstvecin, muestra, Nf);
    double **Cfdifvecin; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red extendida
    hor.espacio(Cfdifvecin, muestra, Nf);
    int Cfdvecin;
//     int Cfddori; //Redes con duplicado que tienen fenotipo diferente al original
    
    //Control 2 para pares
    double **CMatrizf2; //Matriz de red de genes factores de transcripción extendida
    hor.espacio(CMatrizf2, Nf + 1, Nf + 1);
    double **CMatrizp2; //Matriz de red de genes reguladores del fenotipo extendida
    hor.espacio(CMatrizp2, Np, Nf + 1);
    int **Ctrayectoria2;
    hor.espacio(Ctrayectoria2, 1000, Nf + 1);
    int *Cvectorinicio2;
    hor.espacio(Cvectorinicio2, Nf + 1);
    int *Cvectorentrada2;
    hor.espacio(Cvectorentrada2, Nf + 1);
    int **CFenotipo2; //Fenotipo no porque debe ser flexible en el periodo
    //Cambio en el fenotipo por el hecho de expandir
    double *CDist2; //Distancias vs el f original de los extendidos
    hor.espacio(CDist2, Nf); //Hay Nf extendidos por muestra
    double *CMedia2; //Promedio de distancias del fenotipo resultado de extenderse vs el f original
    hor.espacio(CMedia2, muestra);
    double *CDesvst2; //Desviación estándar de las distancias de las resdes extendidas por muestra
    hor.espacio(CDesvst2, muestra);
    double *Cfdif2; //Cuántas redes extendidas dieron fenotipo diferente
    hor.espacio(Cfdif2, muestra);
    int Cfd2; //Para la suma de f diferentes
    //Vecinos de las redes extendidas vs el f original
    double *CDistvec2; //Distancias de los vecinos por red extendida
    hor.espacio(CDistvec2, Nf * explorar);
    double **CMediavec2; //Medias de vecinos por red extendida
    hor.espacio(CMediavec2, muestra, Nf);
    double **CDesvstvec2; //Desviación estándar por red extendida
    hor.espacio(CDesvstvec2, muestra, Nf);
    double **Cfdifvec2; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red extendida
    hor.espacio(Cfdifvec2, muestra, Nf);
    int Cfdvec2;
    //Vecinos de las redes extendidas vs el f resultante de extender
    double *CDistvecin2; //Distancias de los vecinos (mutadas) de las redes extendidas con la extendida
//     hor.espacio(CDistvecin2, Nf * explorar);
    hor.espacio(CDistvecin2, explorar);
    double **CMediavecin2; //Medias de vecinos por red extendida
    hor.espacio(CMediavecin2, muestra, Nf);
    double **CDesvstvecin2; //Desviación estándar por red extendida
    hor.espacio(CDesvstvecin2, muestra, Nf);
    double **Cfdifvecin2; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red extendida
    hor.espacio(Cfdifvecin2, muestra, Nf);
    int Cfdvecin2;
//     int Cfddori2; //Redes con duplicado que tienen fenotipo diferente al original
    
    //Mínimos, máximos y promedios (para condensar lo que va por número de factores)
    //proporciones
    double **pdv; 
    hor.espacio(pdv, muestra, 3);
    double **pdvin; 
    hor.espacio(pdvin, muestra, 3);
    double **pc1v; 
    hor.espacio(pc1v, muestra, 3);
    double **pc1vin; 
    hor.espacio(pc1vin, muestra, 3);
    double **pc2v; 
    hor.espacio(pc2v, muestra, 3);
    double **pc2vin; 
    hor.espacio(pc2vin, muestra, 3);
    //medias
    double **mdv; 
    hor.espacio(mdv, muestra, 3);
    double **mdvin; 
    hor.espacio(mdvin, muestra, 3);
    double **mc1v; 
    hor.espacio(mc1v, muestra, 3);
    double **mc1vin; 
    hor.espacio(mc1vin, muestra, 3);
    double **mc2v; 
    hor.espacio(mc2v, muestra, 3);
    double **mc2vin; 
    hor.espacio(mc2vin, muestra, 3);
    //desviaciones
    double **ddv; 
    hor.espacio(ddv, muestra, 3);
    double **ddvin; 
    hor.espacio(ddvin, muestra, 3);
    double **dc1v; 
    hor.espacio(dc1v, muestra, 3);
    double **dc1vin; 
    hor.espacio(dc1vin, muestra, 3);
    double **dc2v; 
    hor.espacio(dc2v, muestra, 3);
    double **dc2vin; 
    hor.espacio(dc2vin, muestra, 3);
    
    //////////Accesibilidad//////////

    int Fnuevom, /*Fnuevod, Fnuevo1, Fnuevo2,*/ Fnuevovecdup, Fnuevovecdupin, CFnuevovec, CFnuevovecin, CFnuevovec2, CFnuevovecin2, regresa1, regresa2, regresa3;
    int *Accesm; //Conteo de fenotipos nuevos en muestras debido a mutaciones. Por fenotipo nuevo excluye al fenotipo original y repetidos.
    hor.espacio(Accesm, muestra);
    
    int **Accesvecdup; //Conteo de f nuevos debido a mutaciones de duplicados que conservan f original.
    hor.espacio(Accesvecdup, muestra, Nf);
    int **Accesvecdupin; //Conteo de f nuevos debido a mutaciones de duplicados que no conservan f original.
    hor.espacio(Accesvecdupin, muestra, Nf);
    int **CAccesvec; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAccesvec, muestra, Nf);
    int **CAccesvecin; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAccesvecin, muestra, Nf);
    int **CAccesvec2; //Conteo de fenotipos nuevos debido a gen de más;
    hor.espacio(CAccesvec2, muestra, Nf);
    int **CAccesvecin2; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAccesvecin2, muestra, Nf);
    int **Regreso1;
    hor.espacio(Regreso1, muestra, Nf);
    int **Regreso2;
    hor.espacio(Regreso2, muestra, Nf);
    int **Regreso3;
    hor.espacio(Regreso3, muestra, Nf);
    
    int *Avecdup; //Conteo de f nuevos debido a mutaciones de duplicados que conservan f original.
    hor.espacio(Avecdup, muestra);
    int *Avecdupin; //Conteo de f nuevos debido a mutaciones de duplicados que no conservan f original.
    hor.espacio(Avecdupin, muestra);
    int *CAvec; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAvec, muestra);
    int *CAvecin; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAvecin, muestra);
    int *CAvec2; //Conteo de fenotipos nuevos debido a gen de más;
    hor.espacio(CAvec2, muestra);
    int *CAvecin2; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
    hor.espacio(CAvecin2, muestra);
    
    //Para el total de accesibilidad por evento por espécimen
    int **Accesdup, **CAcces, **CAcces2;
    hor.espacio(Accesdup, muestra, Nf);
    hor.espacio(CAcces, muestra, Nf);
    hor.espacio(CAcces2, muestra, Nf);
    
    ofstream proporciones, medias, desviaciones, comparaciones, duppasos, c1pasos, c2pasos, matrices, accesibilidad, accesibilidad1, accesibilidad2, accesibilidad3, accesibilidad4, accesibilidad5, accesibilidad6, accesibilidad7, accesibilidad8, accesibilidad9, accesibilidad10, accesibilidad11, accesibilidad12, accesibilidad13, especimenes/*, chequeo*/;
//     ofstream duppasos1, duppasos2, c1pasos1, c1pasos2, c2pasos1, c2pasos2;
    
    //Índice de Shannon
    ofstream dataesp, datadup, datac1, datac2, datadup2, datac12, datac22;
    int *valordekesp;
    hor.espacio(valordekesp, explorar);
    int *valordekdup;
    hor.espacio(valordekdup, explorar);
    int *valordekc1;
    hor.espacio(valordekc1, explorar);
    int *valordekc2;
    hor.espacio(valordekc2, explorar);
    double *shannonesp;
    hor.espacio(shannonesp, muestra);
    double **shannondup;
    hor.espacio(shannondup, muestra, Nf);
    double **shannonc1;
    hor.espacio(shannonc1, muestra, Nf);
    double **shannonc2;
    hor.espacio(shannonc2, muestra, Nf);
    double **shannondup2;
    hor.espacio(shannondup2, muestra, Nf);
    double **shannonc12;
    hor.espacio(shannonc12, muestra, Nf);
    double **shannonc22;
    hor.espacio(shannonc22, muestra, Nf);
    
//     ofstream nan;
        gat.condin(condinicial, Nf);
        
        for( int i = 0; i < Nf; i++){
            vectorinicio[i] = condinicial[i];
            trayectoria[0][i] = condinicial[i];
            vectoriniciod[i] = condinicial[i];
            Cvectorinicio[i] = condinicial[i];
            Cvectorinicio2[i] = condinicial[i];
            cout << condinicial[i];
        }
        cout << endl;
        
        cout << "////////////////////////////ELEGIR MATRIZ DE FACTORES CON PUNTO FIJO/////////////////////////////////////\n";
        
        //--OJO!, El que trayectoria sea finito puede provocar error a Nfs grandes--//
        
        do{
            n = 1;            
            for(int i = 0; i < Nf; i++)
                vectorentrada[i] = condinicial[i];
            
            lin.make(Matrizf, Nf, Nf, intnw);
            find.atractor(Matrizf, trayectoria, vectorentrada, Nf, n, p);
        }while(p!=1);
        
        for(int i = 0; i < Nf; i++)
            atractor[i] = trayectoria[n-1][i];
        
        cout << "//////////////////////////////////CREAMOS MATRIZ DE MARCADORES///////////////////////////////////////////\n";
        
        lin.make(Matrizp, Np, Nf, intnw);
        
        //Producto Matriz de fenotipo con vector(es) del atractor
        hor.espacio(Fenotipo, p, Np);
        
        lin.dot(Matrizp, trayectoria, Fenotipo, Nf, Np, n, p);
        
        cout << "//////////////////////////SELECCIONAR FENOTIPO///////////////////////////////////\n";
        
        hor.espacio(Wish, Np);       
        gat.escalon(Wish, Np);
        gat.cambiafenotipo(Matrizp, Wish, Fenotipo, Nf, Np, p);        
        lin.dot(Matrizp, trayectoria, Fenotipo, Nf, Np, n, p); //Corroboramos el nuevo fenotipo coincida con el desado
        
       cout << "///////////////////////////CAMINATA////////////////////////////////////\n";
        
        comp1[0] = 0;
        comp2[0] = 0;

        
        for(int m = 0; m < muestra; m++){
            cout << "Muestra: " << m;
            pasos = 0;
            standby = 0;
            
            do{
                //Trayectoria mantiene intacto trayectoria[0][0] = condición incial
                gat.walk(Matrizf, Matrizp, trayectoria, Fenotipo, Nf, Np, tolerancia, intnw, vectorinicio, caminata, n, p, pasos, standby, m);
                //Vector inicio como vector entrada para que sea condición inicial. n = 1.
                //Se muta aleatoriamente un elemento de la red, y la mutación permanece si la nueva red da como resultado el mismo fenotipo.
                //n cambia dentro de walk, pero no afecta a n a este nivel, ya que la n que introducimos aquí no cambia
            }while(pasos < caminata);
            
            cout << "Colectamos espécimen #: " << m + 1 << endl;
            especimenes.open("especimenes.txt", ios::app);
            for(int i = 0; i < Nf; i++){
                for(int j = 0; j < Nf; j++){
                    especimenes << Matrizf[i][j] << " ";
                }
                especimenes << endl;
            }
            for(int i = 0; i < Np; i++){
                for(int j = 0; j < Nf; j++){
                    especimenes << Matrizp[i][j] << " ";
                }
                especimenes << endl;
            }
            especimenes << endl;
            especimenes.close();
            
            pasosfuera[m] = (double)standby / caminata;
            
            if(m == 0){
                hor.rellenar(Matrizf, Matriz1f, Nf, Nf);
                hor.rellenar(Matrizp, Matriz1p, Np, Nf);
            }
            
            if(m > 0){
                comp1[1] = find.compara1(Matriz1f, Matriz1p, Matrizf, Matrizp, Nf, Np);
                if(comp1[1] > comp1[0]){comp1[0] = comp1[1];}
                comp2[1] = find.compara2(Matriz1f, Matriz1p, Matrizf, Matrizp, Nf, Np);
                if(comp2[1] > comp2[0]){comp2[0] = comp2[1];}
            }
            
            fdif = 0;
            hor.fillv0(valordekesp, explorar);
            gat.exploracion(Matrizf, Matrizp, trayectoria, Fenotipo, Distvecinos, Nf, Np, tolerancia, intnw, vectorinicio, n, p, m, explorar, fdif, Fnuevom, valordekesp);//Muta red de muestra, saca distancia fenotípica y la guarda en una matriz
            //               Matriz1: Factores
            //               Matriz2: Marcadores
            //               Matriz3: Trayectoria
            //               Matriz4: Fenotipo
            //               Matriz5: Distancias
            //               tam1 = hils1 = cols1 = cols2 = Nf
            //               tam2 = hils2 = Np
            //               num1: Número de muestra para saber en qué hilera va la distancia
            //               num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
            
            cout << ": Accesibilidad - ";
            //////////Accesibilidad después mutar, fmut diferentes al fori//////////
            
            Accesm[m] = Fnuevom;
            
            cout << "Robustez -";
            fdiferente[m] = (double)fdif / explorar;
            find.media(Distvecinos, Mediavecinos, explorar, m);
            find.desv(Distvecinos, Mediavecinos, Desvst, explorar, m);
            
            cout << "Shannon\n";
//             cout << "Antes\n";
            shannonesp[m] = gat.Shannon(valordekesp, explorar);
            
            fddup = 0;
            Cfd = 0;
            Cfd2 = 0;
            
            for(int gen = 0; gen < Nf; gen++){
                cout << "Muestra: " << m << " Gen: " << gen;
                ///////////////////////////DUPLICACIÓN////////////////////////////////////
                cout << " DUPLICACIÓN - ";
                gat.inducirduplicacion(Matrizfd, Matrizpd, Matrizf, Matrizp, gen, Nf, Np);
                //Matriz1: Factores duplicadas
                //Matriz2: Marcadores duplicadasFnuevo2
                //Matriz3: Factores
                //Matriz4: Marcadores
                //gen: gen a duplicar
                //Tam1: Nf
                //Tam2: Np
                
                //               Fenotipo de la red duplicada
                for( int i = 0; i < Nf; i++){
                    vectoriniciod[i] = condinicial[i];
                }
                vectoriniciod[Nf] = condinicial[gen]; //La condición inicial de ese gen se copia.
                nd = 1;
                
                for(int i = 0; i < (Nf + 1); i++){
                    vectorentradad[i] = vectoriniciod[i];
                    trayectoriad[0][i] = vectoriniciod[i];
                }
                
                find.atractor(Matrizfd, trayectoriad, vectorentradad, Nf + 1, nd, perd);
                
                //Producto Matriz de fenotipo con vector(es) del atractor
                
                hor.espacio(Fenotipod, perd, Np);
                lin.dot(Matrizpd, trayectoriad, Fenotipod, Nf + 1, Np, nd, perd);
                Distdup[gen] = 1 - find.distfenotipica(Fenotipo, p, Fenotipod, perd, Np);
                if(Distdup[gen] != 1){//Si fdup no es igual a fori
                    duppasos.open ("outs/Pasos/duppasosdi.txt", ios::app);
                    duppasos << "0 ";
                    duppasos.close ();
                    fddup++;
                }
                
                find.media(Distdup, Mediadup, Nf, m);
                find.desv(Distdup, Mediadup, Desvstdup, Nf, m);
                
                //Ya que tenemos la duplicada, la mutamos varias veces y tomamos las distancias
                fdvecdupin = 0;
                fdvecdup = 0;
                regresa1 = 0;
                hor.fillv0(valordekdup, explorar);
                gat.exploraciondup(Matrizfd, Matrizpd, trayectoriad, Fenotipod, Fenotipo, Distvecdupin, Distvecdup, Nf + 1, Np, tolerancia, intnw, vectoriniciod, nd, perd, p, gen, explorar, fdvecdup, fdvecdupin, regresa1, Fnuevovecdup, Fnuevovecdupin, 1, valordekdup);
                //     Matriz1: Factores
                //     Matriz2: Marcadores
                //     Matriz3: Trayectoria
                //     Matriz4: Fenotipo del duplicado
                //     Matriz5: Fenotipo original
                //     Vector1: Distancia con el fenotipo del duplicado
                //     Vector2: Distancia con el fenotipo original
                //     num: Máximo de steps
                //     tam1 = hils1 = cols1 = cols2 = Nf
                //     tam2 = hils2 = Np
                //     num1: Número de muestra para saber en qué hilera va la distancia
                //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
                //     cuenta1: Veces que no dio mismo fenotipo que el duplicado
                //     cuenta2: Veces que no dio mismo fenotipo que el original
                //                 hor.printmat(trayectoriad, nd, Nf + 1);
                //Vecinos de las redes duplicadas vs el f original
                
                
                cout << "Accesibilidad - ";
                Accesvecdup[m][gen] = Fnuevovecdup;
                Accesvecdupin[m][gen] = Fnuevovecdupin;
                Accesdup[m][gen] = Fnuevovecdup + Fnuevovecdupin;
                Regreso1[m][gen] = regresa1;
                
                cout << "Robustez - ";
                fdifvecdup[m][gen] = 1 - ((double)fdvecdup / explorar);
                find.media(Distvecdup, Mediavecdup, explorar, m, gen);
                find.desv(Distvecdup, Mediavecdup, Desvstvecdup, explorar, m, gen);
                

                //Vecinos de las redes duplicadas vs el f resultante de duplicar
                
                if(Distdup[gen] == 1){//Si el fenotipo es igual al original
                    duppasos.open ("outs/Pasos/duppasosdi.txt", ios::app);
                    duppasos << "1 ";
                    duppasos.close ();
                    fdifvecdupin[m][gen] = 2;
                    Mediavecdupin[m][gen] = 2;
                    Desvstvecdupin[m][gen] = 2;
                }
                else{
                    fdifvecdupin[m][gen] = 1 - ((double)(fdvecdupin + regresa1) / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                    find.media(Distvecdupin, Mediavecdupin, explorar, m, gen);
                    find.desv(Distvecdupin, Mediavecdupin, Desvstvecdupin, explorar, m, gen);
                }
                
                cout << "Shannon\n";
                shannondup[m][gen] = gat.Shannon(valordekdup, explorar);
                shannondup2[m][gen] = gat.Shannon2(valordekdup, explorar);
                
                /*
                if((m == 500) && (gen == 3)){
                    nan << "m: " << m << " evento: " << gen << endl;hábito
                    nan.open("nan.txt", ios::app);
                    for(int i = 0; i < explorar; i++){
                        nan << valordekdup[i] << endl;
                    }
                    nan.close();
                    gat.Shannon3(valordekdup, explorar);
                    
                }
                
                if((m == 947) && (gen == 5)){
                    nan << "m: " << m << " evento: " << gen << endl;
                    nan.open("nan.txt", ios::app);
                    for(int i = 0; i < explorar; i++){
                        nan << valordekdup[i] << endl;
                    }
                    nan.close();
                    gat.Shannon3(valordekdup, explorar);
                    
                }
                
                if(shannondup[m][gen] == 0) {
                    cout << "m: " << m << " evento: " << gen << endl;
                    nan.open("nan.txt", ios::app);
                    nan << "m: " << m << " evento: " << gen << endl;
                    for(int i = 0; i < explorar; i++){
                        nan << valordekdup[i] << endl;
                    }
                    nan.close();
                }
                
                if(shannondup2[m][gen] == 0) {
                    cout << "m2: " << m << " evento: " << gen << endl;
                    nan.open("nan.txt", ios::app);
                    nan << "m2: " << m << " evento: " << gen << endl;
                    for(int i = 0; i < explorar; i++){
                        nan << valordekdup[i] << endl;
                    }
                    nan.close();
                }*/
                
                ///////////////////////////CONTROL 1////////////////////////////////////
                cout << "CONTROL 1 - ";
                //Gen nuevo independiente

                gat.control1(CMatrizf, CMatrizp, Matrizf, Matrizp, Nf, Np, intnw);
                //Matriz1: Factores extendida
                //Matriz2: Marcadores extendida
                //Matriz3: Factores
                //Matriz4: Marcadores
                //Tam1: Nf
                //Tam2: Np
                //intnw
                //               Fenotipo de la red duplicada
                for( int i = 0; i < Nf; i++){
                    Cvectorinicio[i] = condinicial[i];
                }
                Cvectorinicio[Nf] = condinicial[gen]; //La condición inicial de ese gen se copia.
                
                Cn = 1;
                for(int i = 0; i < (Nf + 1); i++){
                    Cvectorentrada[i] = Cvectorinicio[i];
                    Ctrayectoria[0][i] = Cvectorinicio[i];
                }
                
                find.atractor(CMatrizf, Ctrayectoria, Cvectorentrada, Nf + 1, Cn, Cp);
               
                //Producto Matriz de fenotipo con vector(es) del atractor
                
                hor.espacio(CFenotipo, Cp, Np);             
                lin.dot(CMatrizp, Ctrayectoria, CFenotipo, Nf + 1, Np, Cn, Cp);             
                CDist[gen] = 1 - find.distfenotipica(Fenotipo, p, CFenotipo, Cp, Np);
                if(CDist[gen] != 1){//Si Fc1 diferente a fori
                    c1pasos.open ("outs/Pasos/c1pasosdi.txt", ios::app);
                    c1pasos << "0 ";
                    c1pasos.close ();
                    Cfd++;

                }            
                find.media(CDist, CMedia, Nf, m);             
                find.desv(CDist, CMedia, CDesvst, Nf, m);
                
                //Ya que tenemos la duplicada, la mutamos varias veces y tomamos las distancias
                Cfdvecin = 0;
                Cfdvec = 0;
                regresa2 = 0;
                hor.fillv0(valordekc1, explorar);
                gat.exploraciondup(CMatrizf, CMatrizp, Ctrayectoria, CFenotipo, Fenotipo, CDistvecin, CDistvec, Nf + 1, Np, tolerancia, intnw, Cvectorinicio, Cn, Cp, p, gen, explorar, Cfdvec, Cfdvecin, regresa2, CFnuevovec, CFnuevovecin, 2, valordekc1);
                //     Matriz1: Factores
                //     Matriz2: Marcadores
                //     Matriz3: Trayectoria
                //     Matriz4: Fenotipo del duplicado
                //     Matriz5: Fenotipo original
                //     Vector1: Distancia con el fenotipo del duplicado
                //     Vector2: Distancia con el fenotipo original
                //     num: Máximo de steps
                //     tam1 = hils1 = cols1 = cols2 = Nf
                //     tam2 = hils2 = Np
                //     num1: Número de muestra para saber en qué hilera va la distancia
                //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
                //     cuenta1: Veces que no dio mismo fenotipo que el duplicado
                //     cuenta2: Veces que no dio mismo fenotipo que el original
                
                //Vecinos de las redes duplicadas vs el f original

                cout << " Accesibilidad - ";
                CAccesvec[m][gen] = CFnuevovec;
                CAccesvecin[m][gen] = CFnuevovecin;
                CAcces[m][gen] = CFnuevovec + CFnuevovecin;
                Regreso2[m][gen] = regresa2;
                
                cout << "Robustez - ";
                Cfdifvec[m][gen] = 1 - ((double)Cfdvec / explorar);
                find.media(CDistvec, CMediavec, explorar, m, gen);
                find.desv(CDistvec, CMediavec, CDesvstvec, explorar, m, gen);
                
                if(CDist[gen] == 1){
                    c1pasos.open ("outs/Pasos/c1pasosdi.txt", ios::app);
                    c1pasos << "1 ";
                    c1pasos.close ();
                    Cfdifvecin[m][gen] = 2;
                    CMediavecin[m][gen] = 2;
                    CDesvstvecin[m][gen] = 2;
                }
                else{
                    Cfdifvecin[m][gen] = 1 - ((double) (Cfdvecin + regresa2) / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                    find.media(CDistvecin, CMediavecin, explorar, m, gen);
                    find.desv(CDistvecin, CMediavecin, CDesvstvecin, explorar, m, gen);
                }
                
                cout << "Shannon\n";
//                 cout << "CI\n";
                shannonc1[m][gen] = gat.Shannon(valordekc1, explorar);
                shannonc12[m][gen] = gat.Shannon2(valordekc1, explorar);
                
                ///////////////////////////CONTROL 2////////////////////////////////////
                cout << "CONTROL 2 - ";
                //Mismo número de conexiones del gen duplicado pero en orden diferente en el gen nuevo.

                gat.control2(CMatrizf2, CMatrizp2, Matrizfd, Matrizpd, Nf, Np);
                //Matriz1: Factores extendida
                //Matriz2: Marcadores extendida
                //Matriz3: Factores
                //Matriz4: Marcadores
                //Tam1: Nf
                //Tam2: Np
                //intnw
                //               Fenotipo de la red duplicada

                for( int i = 0; i < Nf; i++){
                    Cvectorinicio2[i] = condinicial[i];
                }
                Cvectorinicio2[Nf] = condinicial[gen]; //La condición inicial de ese gen se copia.
                
                Cn2 = 1;
                
                for(int i = 0; i < (Nf + 1); i++){
                    Cvectorentrada2[i] = Cvectorinicio2[i];
                    Ctrayectoria2[0][i] = Cvectorinicio2[i];
                }
                
                find.atractor(CMatrizf2, Ctrayectoria2, Cvectorentrada2, Nf + 1, Cn2, Cp2);
                
                //Producto Matriz de fenotipo con vector(es) del atractor
                
                hor.espacio(CFenotipo2, Cp2, Np);
                lin.dot(CMatrizp2, Ctrayectoria2, CFenotipo2, Nf + 1, Np, Cn2, Cp2);
                CDist2[gen] = 1 - find.distfenotipica(Fenotipo, p, CFenotipo2, Cp2, Np);
                if(CDist2[gen] != 1){//Si Fc2 diferente a fori
                    c2pasos.open ("outs/Pasos/c2pasosdi.txt", ios::app);
                    c2pasos << "0 ";
                    c2pasos.close ();
                    Cfd2++;
                }
                find.media(CDist2, CMedia2, Nf, m);
                find.desv(CDist2, CMedia2, CDesvst2, Nf, m);                
                //Ya que tenemos la duplicada, la mutamos varias veces y tomamos las distancias
                Cfdvecin2 = 0;
                Cfdvec2 = 0;
                regresa3 = 0;
                hor.fillv0(valordekc2, explorar);
                gat.exploraciondup(CMatrizf2, CMatrizp2, Ctrayectoria2, CFenotipo2, Fenotipo, CDistvecin2, CDistvec2, Nf + 1, Np, tolerancia, intnw, Cvectorinicio2, Cn2, Cp2, p, gen, explorar, Cfdvec2, Cfdvecin2, regresa3, CFnuevovec2, CFnuevovecin2, 3, valordekc2);
                //     Matriz1: Factores
                //     Matriz2: Marcadores
                //     Matriz3: Trayectoria
                //     Matriz4: Fenotipo del duplicado
                //     Matriz5: Fenotipo original
                //     Vector1: Distancia con el fenotipo del duplicado
                //     Vector2: Distancia con el fenotipo original
                //     num: Máximo de steps
                //     tam1 = hils1 = cols1 = cols2 = Nf
                //     tam2 = hils2 = Np
                //     num1: Número de muestra para saber en qué hilera va la distancia
                //     num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
                //     cuenta1: Veces que no dio mismo fenotipo que el duplicado
                //     cuenta2: Veces que no dio mismo fenotipo que el original
                
                //Vecinos de las redes duplicadas vs el f original
                
                cout << "Accesibilidad - ";
                CAccesvec2[m][gen] = CFnuevovec2;
                CAccesvecin2[m][gen] = CFnuevovecin2;
                CAcces2[m][gen] = CFnuevovec2 + CFnuevovecin2;
                Regreso3[m][gen] = regresa3;
                
                cout << "Robustez - ";
                Cfdifvec2[m][gen] = 1 - ((double)Cfdvec2 / explorar);
                find.media(CDistvec2, CMediavec2, explorar, m, gen);                          
                find.desv(CDistvec2, CMediavec2, CDesvstvec2, explorar, m, gen);
                
                if(CDist2[gen] == 1){
                    c2pasos.open ("outs/Pasos/c2pasosdi.txt", ios::app);
                    c2pasos << "1 ";
                    c2pasos.close ();
                    Cfdifvecin2[m][gen] = 2;
                    CMediavecin2[m][gen] = 2;
                    CDesvstvecin2[m][gen] = 2;
                }
                else{
                    Cfdifvecin2[m][gen] = 1 - ((double) (Cfdvecin2 + regresa3) / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                    find.media(CDistvecin2, CMediavecin2, explorar, m, gen);
                    find.desv(CDistvecin2, CMediavecin2, CDesvstvecin2, explorar, m, gen);
                }
                
                cout << " Shannon\n";
                shannonc2[m][gen] = gat.Shannon(valordekc2, explorar);
                shannonc22[m][gen] = gat.Shannon2(valordekc2, explorar);
                
            }//Cierra ciclo de Genes

            duppasos.open("outs/Pasos/duppasosdi.txt", ios::app);
            duppasos << endl;
            duppasos.close();
            c1pasos.open("outs/Pasos/c1pasosdi.txt", ios::app);
            c1pasos << endl;
            c1pasos.close();
            c2pasos.open("outs/Pasos/c2pasosdi.txt", ios::app);
            c2pasos << endl;
            c2pasos.close();
            
            comparaciones.open ("outs/comparaciones.txt", ios::app);
            comparaciones << comp1[0] << " " << comp2[0] << "\n";
            comparaciones.close ();
            fdifdup[m] = (double)fddup / Nf;
            Cfdif[m] = (double)Cfd / Nf; 
            Cfdif2[m] = (double)Cfd2 / Nf;

            hor.sumarmatriz(Accesvecdup, Avecdup, muestra, Nf);
            hor.sumarmatriz(Accesvecdupin, Avecdupin, muestra, Nf);
            hor.sumarmatriz(CAccesvec, CAvec, muestra, Nf);
            hor.sumarmatriz(CAccesvecin, CAvecin, muestra, Nf);
            hor.sumarmatriz(CAccesvec2, CAvec2, muestra, Nf);
            hor.sumarmatriz(CAccesvecin2, CAvecin2, muestra, Nf);

        }//Cierra ciclo de Muestras

        cout << "------------Condensamos matrices que van por factor a mín, máx y promedio------------\n";
        find.minmaxmed(fdifvecdup, pdv, muestra, Nf);
        find.minmaxmedin(fdifvecdupin, pdvin, muestra, Nf);
        find.minmaxmed(Cfdifvec, pc1v, muestra, Nf);
        find.minmaxmedin(Cfdifvecin, pc1vin, muestra, Nf);
        find.minmaxmed(Cfdifvec2, pc2v, muestra, Nf);
        find.minmaxmedin(Cfdifvecin2, pc2vin, muestra, Nf);
        
        find.minmaxmed(Mediavecdup, mdv, muestra, Nf);
        find.minmaxmedin(Mediavecdupin, mdvin, muestra, Nf);
        find.minmaxmed(CMediavec, mc1v, muestra, Nf);
        find.minmaxmedin(CMediavecin, mc1vin, muestra, Nf);
        find.minmaxmed(CMediavec2, mc2v, muestra, Nf);
        find.minmaxmedin(CMediavecin2, mc2vin, muestra, Nf);
        
        find.minmaxmed(Desvstvecdup, ddv, muestra, Nf);
        find.minmaxmedin(Desvstvecdupin, ddvin, muestra, Nf);
        find.minmaxmed(CDesvstvec, dc1v, muestra, Nf);
        find.minmaxmedin(CDesvstvecin, dc1vin, muestra, Nf);
        find.minmaxmed(CDesvstvec2, dc2v, muestra, Nf);
        find.minmaxmedin(CDesvstvecin2, dc2vin, muestra, Nf);
        
        proporciones.open ("outs/proporciones.txt", ios::app);
        medias.open ("outs/medias.txt", ios::app);
        desviaciones.open ("outs/desviaciones.txt", ios::app);
        accesibilidad.open ("outs/acces/accesibilidad.txt", ios::app);
        accesibilidad1.open ("outs/acces/accesesp.txt", ios::app);
        accesibilidad2.open ("outs/acces/accesdupori.txt", ios::app);
        accesibilidad3.open ("outs/acces/accesdupdif.txt", ios::app);
        accesibilidad4.open ("outs/acces/accesc1ori.txt", ios::app);
        accesibilidad5.open ("outs/acces/accesc1dif.txt", ios::app);
        accesibilidad6.open ("outs/acces/accesc2ori.txt", ios::app);
        accesibilidad7.open ("outs/acces/accesc2dif.txt", ios::app);
        accesibilidad8.open ("outs/acces/accesdup.txt", ios::app);
        accesibilidad9.open ("outs/acces/accesc1.txt", ios::app);
        accesibilidad10.open ("outs/acces/accesc2.txt", ios::app);
        accesibilidad11.open ("outs/acces/regdup.txt", ios::app);
        accesibilidad12.open ("outs/acces/regc1.txt", ios::app);
        accesibilidad13.open ("outs/acces/regc2.txt", ios::app);
        
        cout << "------------Archivos generándose------------\n";
        
        //Tabla
        for(int i = 0; i < muestra; i++){
            proporciones << 1 - pasosfuera[i] << " " << 1 - fdiferente[i] << " " << 1 - fdifdup[i] << " " << 1 - Cfdif[i] << " " <<  1 - Cfdif2[i] << " ";
            for(int j = 0; j < 3; j++){
                proporciones << pdv[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                proporciones << pdvin[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                proporciones << pc1v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                proporciones << pc1vin[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                proporciones << pc2v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                proporciones << pc2vin[i][j] << " ";
            }
            proporciones << "\n";
        }
        
        for(int i = 0; i < muestra; i++){
            medias << Mediavecinos[i] << " " << Mediadup[i] << " " << CMedia[i] << " " << CMedia2[i] << " ";
            for(int j = 0; j < 3; j++){
                medias << mdv[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                medias << mdvin[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                medias << mc1v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                medias << mc1vin[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                medias << mc2v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                medias << mc2vin[i][j] << " ";
            }
            medias << "\n";
        }
        
        for(int i = 0; i < muestra; i++){
            desviaciones << Desvst[i] << " " << Desvstdup[i] << " " << CDesvst[i] << " " << CDesvst2[i] << " ";
            for(int j = 0; j < 3; j++){
                desviaciones << ddv[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                desviaciones << ddvin[i][j] << " ";
            }      
            for(int j = 0; j < 3; j++){
                desviaciones << dc1v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                desviaciones << dc1vin[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                desviaciones << dc2v[i][j] << " ";
            }
            for(int j = 0; j < 3; j++){
                desviaciones << dc2vin[i][j] << " ";
            }
            desviaciones << "\n";
        }
        
        //Accesibilidad
        //Accesibilidad de los especímenes
        for(int i = 0; i < muestra; i++){
            accesibilidad1 << " " << Accesm[i] << endl;
        }
        //Accesibilidad de los duplicados que obtuvieron fenotipo original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad2 << Accesvecdup[i][j] << " ";
            }
            accesibilidad2 << endl;
        }
        //Accesibilidad de los duplicados que obtuvieron fenotipo diferente al original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad3 << Accesvecdupin[i][j] << " ";
            }
            accesibilidad3 << endl;
        }
        //Accesibilidad de los CI que obtuvieron fenotipo original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad4 << CAccesvec[i][j] << " ";
            }
            accesibilidad4 << endl;
        }
        //Accesibilidad de los CI que obtuvieron fenotipo diferente al original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad5 << CAccesvecin[i][j] << " ";
            }
            accesibilidad5 << endl;
        }
        //Accesibilidad de los CI que obtuvieron fenotipo original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad6 << CAccesvec2[i][j] << " ";
            }
            accesibilidad6 << endl;
        }
        //Accesibilidad de los CI que obtuvieron fenotipo diferente al original
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad7 << CAccesvecin2[i][j] << " ";
            }
            accesibilidad7 << endl;
        }
        //Totales
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad8 << Accesdup[i][j] << " ";
                accesibilidad9 << CAcces[i][j] << " ";
                accesibilidad10 << CAcces2[i][j] << " ";
            }
            accesibilidad8 << endl;
            accesibilidad9 << endl;
            accesibilidad10 << endl;
        }
        //Regresaron al fenotipo original después de agregar gen y mutar. Dup, C1 y C2 respectivamente
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad11 << Regreso1[i][j] << " ";
            }
            accesibilidad11 << endl;
        }
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad12 << Regreso2[i][j] << " ";
            }
            accesibilidad12 << endl;
        }
        for(int i = 0; i < muestra; i++){
            for(int j = 0; j < Nf; j++){
                accesibilidad13 << Regreso3[i][j] << " ";
            }
            accesibilidad13 << endl;
        }
        
        //Shannon
        dataesp.open("outs/Shannon/dataesp.txt", ios::app); 
        datadup.open("outs/Shannon/datadup.txt", ios::app);  
        datac1.open("outs/Shannon/datac1.txt", ios::app); 
        datac2.open("outs/Shannon/datac2.txt", ios::app); 
        datadup2.open("outs/Shannon/datadup2.txt", ios::app);  
        datac12.open("outs/Shannon/datac12.txt", ios::app); 
        datac22.open("outs/Shannon/datac22.txt", ios::app);

        for(int i = 0; i < muestra; i++){
            dataesp << shannonesp[i] << endl;
            for(int j = 0; j < Nf; j++){
                datadup << shannondup[i][j] << " ";
                datac1 << shannonc1[i][j] << " ";
                datac2 << shannonc2[i][j] << " ";
                datadup2 << shannondup2[i][j] << " ";
                datac12 << shannonc12[i][j] << " ";
                datac22 << shannonc22[i][j] << " ";
            }
            datadup << endl;
            datac1 << endl;
            datac2 << endl;
            datadup2 << endl;
            datac12 << endl;
            datac22 << endl;
        }
        
        proporciones.close ();
        medias.close ();
        desviaciones.close ();
        accesibilidad.close ();
        accesibilidad1.close ();
        accesibilidad2.close ();
        accesibilidad3.close ();
        accesibilidad4.close ();
        accesibilidad5.close ();
        accesibilidad6.close ();
        accesibilidad7.close ();
        accesibilidad8.close ();
        accesibilidad9.close ();
        accesibilidad10.close ();
        dataesp.close ();
        datadup.close ();
        datac1.close ();
        datac2.close ();
        datadup2.close ();
        datac12.close ();
        datac22.close ();
        ////////////////////////////////////////////////////////////////////////////
        cout << "------------Limpieza------------\n";
        hor.borrar(Matrizf, Nf);
        hor.borrar(Matrizp, Np);
        hor.borrar(Matriz1f, Nf);
        hor.borrar(Matriz1p, Np);
        hor.borrar(Matrizfd, Nf + 1);
        hor.borrar(Matrizpd, Np);
        hor.borrar(trayectoria, 1000);
        hor.borrar(trayectoriad, 1000);
        hor.borrar(Fenotipo, p);
        hor.borrar(Fenotipod, perd);
        hor.borrar(Wish);
        hor.borrar(condinicial);
        hor.borrar(vectorinicio);
        hor.borrar(vectorentrada);
        hor.borrar(vectoriniciod);
        hor.borrar(vectorentradad);
        hor.borrar(atractor);
        hor.borrar(pasosfuera);
        hor.borrar(fdiferente);
        hor.borrar(Distvecinos);
        hor.borrar(Mediavecinos);
        hor.borrar(Desvst);
        hor.borrar(Distdup);
        hor.borrar(Mediadup);
        hor.borrar(Desvstdup);
        hor.borrar(fdifdup);
        hor.borrar(Distvecdup);
        hor.borrar(Mediavecdup, muestra);
        hor.borrar(Desvstvecdup, muestra);
        hor.borrar(fdifvecdup, muestra);
        hor.borrar(Desvstvecdupin, muestra);
        hor.borrar(Distvecdupin);
        hor.borrar(Mediavecdupin, muestra);
        hor.borrar(fdifvecdupin, muestra);
        hor.borrar(CDistvec);
        hor.borrar(CMediavec, muestra);
        hor.borrar(CDesvstvec, muestra);
        hor.borrar(Cfdifvec, muestra);
        hor.borrar(CDistvecin);
        hor.borrar(CMediavecin, muestra);
        hor.borrar(CDesvstvecin, muestra);
        hor.borrar(Cfdifvecin, muestra);
        hor.borrar(CMatrizf, Nf + 1);
        hor.borrar(CMatrizp, Np);
        hor.borrar(Ctrayectoria, 1000);
        hor.borrar(Cvectorinicio);
        hor.borrar(Cvectorentrada);
        hor.borrar(CDistvec2);
        hor.borrar(CMediavec2, muestra);
        hor.borrar(CDesvstvec2, muestra);
        hor.borrar(Cfdifvec2, muestra);
        hor.borrar(CDistvecin2);
        hor.borrar(CMediavecin2, muestra);
        hor.borrar(CDesvstvecin2, muestra);
        hor.borrar(Cfdifvecin2, muestra);
        hor.borrar(CMatrizf2, Nf + 1);
        hor.borrar(CMatrizp2, Np);
        hor.borrar(Ctrayectoria2, 1000);
        hor.borrar(Cvectorinicio2);
        hor.borrar(Cvectorentrada2);
        hor.borrar(comp1);
        hor.borrar(comp2);
        hor.borrar(pdv, muestra);
        hor.borrar(pdvin, muestra);
        hor.borrar(pc1v, muestra);
        hor.borrar(pc1vin, muestra);
        hor.borrar(pc2v, muestra);
        hor.borrar(pc2vin, muestra);
        hor.borrar(mdv, muestra);
        hor.borrar(mdvin, muestra);
        hor.borrar(mc1v, muestra);
        hor.borrar(mc1vin, muestra);
        hor.borrar(mc2v, muestra);
        hor.borrar(mc2vin, muestra);
        hor.borrar(ddv, muestra);
        hor.borrar(ddvin, muestra);
        hor.borrar(dc1v, muestra);
        hor.borrar(dc1vin, muestra);
        hor.borrar(dc2v, muestra);
        hor.borrar(dc2vin, muestra);
        hor.borrar(Accesm);
        hor.borrar(Accesvecdup, muestra);
        hor.borrar(Accesvecdupin, muestra);
        hor.borrar(CAccesvec, muestra);
        hor.borrar(CAccesvecin, muestra);
        hor.borrar(CAccesvec2, muestra);
        hor.borrar(CAccesvecin2, muestra);
        hor.borrar(Avecdup);
        hor.borrar(Avecdupin);
        hor.borrar(CAvec2);
        hor.borrar(CAvecin2);
        hor.borrar(Accesdup, muestra);
        hor.borrar(CAcces, muestra);
        hor.borrar(CAcces2, muestra);
        hor.borrar(Regreso1, muestra);
        hor.borrar(Regreso2, muestra);
        hor.borrar(Regreso3, muestra);
        hor.borrar(valordekesp);
        hor.borrar(valordekdup);
        hor.borrar(valordekc1);
        hor.borrar(valordekc2);
        hor.borrar(shannonesp);
        hor.borrar(shannondup, muestra);
        hor.borrar(shannonc1, muestra);
        hor.borrar(shannonc2, muestra);
        cout << "------------Ejecución terminada!!-----------\n";

    hor.close_rng();
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    return 0;
}
