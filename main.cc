//Principal program. This program take the GRNs and performs general mutational experiments. It collects the distance between phenotypes after gene addition and mutations. It also takes the Shannon Index for other study.

// g++ -Wall main.cc libs/horse.o libs/linear.o libs/finder.o libs/gataca.o -lgsl -lgslcblas -lm
// ./a.out 10 1000 12 6

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
//     int caminata = 20 * ((Nf * Nf) + (Nf * Np)); // Número de redes a explorar antes de llegar al genotipo que será parte de la muestra
    int explorar = 2 * ((Nf * Nf) + (Nf * Np)); //Mutaciones para medir de la muestra
    //Lo anterior significa que se tomarán 3 muestras separadas por 10 pasos y se mutarán 3 veces, se les hará una duplicación y se mutarán 3 veces (para cada gen).
    double dist;
    int perd, Cp, Cp2; //Periodo del fenotipo
    int nd, Cn, Cn2; //Contador para lin.dot
    
    Horse hor(semilla);
    Finder find;
    Linear lin(semilla);
    GATACA gat(semilla);
    
    hor.start_rng(semilla);
    
    int /**pes, *eneses, */p, n;
//     hor.espacio(pes, muestra);
//     hor.espacio(eneses, muestra);
//     cout << "1\n";
    //Vectores y matrices para crear red de determiado fenotipo fenotipo
    
    double **Matrizf; //Matriz de red de genes factores de transcripción
    hor.espacio(Matrizf, Nf, Nf);
    double **Matrizp; //Matriz de red de genes reguladores del fenotipo
    hor.espacio(Matrizp, Np, Nf);
    int **trayectoria; //Trayectoria de condición inicial a atractor
    hor.espacio(trayectoria, 1000, Nf);
    int *condinicial;
    hor.espacio(condinicial, Nf);
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
    double *PromDist_up;
    hor.espacio(PromDist_up, muestra);
    double *PromDist_down;
    hor.espacio(PromDist_down, muestra);
    
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
    hor.espacio(Distvecdup, explorar);
    double **Mediavecdup; //Medias de vecinos por red duplicada
    hor.espacio(Mediavecdup, muestra, Nf);
    double **Desvstvecdup; //Desviación estándar por duplicado
    hor.espacio(Desvstvecdup, muestra, Nf);
    double **fdifvecdup; //Cuenta de fenotipos diferentes encontrados en los vecinos de una red duplicada
    hor.espacio(fdifvecdup, muestra, Nf);
    int fdvecdup;
    double **Distvecdup_gral;
    hor.espacio(Distvecdup_gral, muestra, Nf);
    
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
    hor.espacio(CDistvec, explorar);
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
    double **Distvecc1_gral;
    hor.espacio(Distvecc1_gral, muestra, Nf);
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
    double **Distvecc2_gral;
    hor.espacio(Distvecc2_gral, muestra, Nf);
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
    
//     int *Avecdup; //Conteo de f nuevos debido a mutaciones de duplicados que conservan f original.
//     hor.espacio(Avecdup, muestra);
//     int *Avecdupin; //Conteo de f nuevos debido a mutaciones de duplicados que no conservan f original.
//     hor.espacio(Avecdupin, muestra);
//     int *CAvec; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
//     hor.espacio(CAvec, muestra);
//     int *CAvecin; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
//     hor.espacio(CAvecin, muestra);
//     int *CAvec2; //Conteo de fenotipos nuevos debido a gen de más;
//     hor.espacio(CAvec2, muestra);
//     int *CAvecin2; //Conteo de fenotipos nuevos debido a mutar redes con un gen de más;
//     hor.espacio(CAvecin2, muestra);
    
    //Para el total de accesibilidad por evento por espécimen
    int **Accesdup, **CAcces, **CAcces2;
    hor.espacio(Accesdup, muestra, Nf);
    hor.espacio(CAcces, muestra, Nf);
    hor.espacio(CAcces2, muestra, Nf);
    
    ofstream proporciones, distancias, distanciasdup, distanciasc1, distanciasc2, desviaciones, comparaciones, duppasos, c1pasos, c2pasos, matrices, accesibilidad, accesibilidad1, accesibilidad2, accesibilidad3, accesibilidad4, accesibilidad5, accesibilidad6, accesibilidad7, accesibilidad8, accesibilidad9, accesibilidad10, accesibilidad11, accesibilidad12, accesibilidad13, regresan1, regresan2, regresan3, robmutdup, robmutc1, robmutc2, robmutupdup, robmutupc1, robmutupc2, robmutdowndup, robmutdownc1, robmutdownc2 ;
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

    //Densidad dos redes, red de arriba y red de abajo (Sólo para especímenes)
    ofstream densidades;
    double *densidadup, *densidaddown, *densidad;
    hor.espacio(densidadup, muestra);
    hor.espacio(densidaddown, muestra);
    hor.espacio(densidad, muestra);

     //Especímenes;
    ifstream especimenes, periodos, enes;
    double **bichos;
    hor.espacio(bichos, (muestra*(Nf + Np)), Nf);

    //Rmut dos redes, red de arriba y red de abajo (para especímenes)
    double *rmutespup;
    hor.espacio(rmutespup, muestra);
    double *rmutespdown;
    hor.espacio(rmutespdown, muestra);
    //Rmut general dos redes (+ 1 gen)
    double **rmutdupgeneral, **rmutc1general, **rmutc2general;
    hor.espacio(rmutdupgeneral, muestra, Nf);
    hor.espacio(rmutc1general, muestra, Nf);
    hor.espacio(rmutc2general, muestra, Nf);
    double **dvdg, **dvc1g, **dvc2g;
    hor.espacio(dvdg, muestra, 3);
    hor.espacio(dvc1g, muestra, 3);
    hor.espacio(dvc2g, muestra, 3);
    double **rmdg, **rmc1g, **rmc2g;
    hor.espacio(rmdg, muestra, 3);
    hor.espacio(rmc1g, muestra, 3);
    hor.espacio(rmc2g, muestra, 3);
    //Rmut general red de arriba (+ 1 gen)
    double **rmutupdup, **rmutupc1, **rmutupc2;
    hor.espacio(rmutupdup, muestra, Nf);
    hor.espacio(rmutupc1, muestra, Nf);
    hor.espacio(rmutupc2, muestra, Nf);
    double **rmupdup, **rmupc1, **rmupc2;
    hor.espacio(rmupdup, muestra, 3);
    hor.espacio(rmupc1, muestra, 3);
    hor.espacio(rmupc2, muestra, 3);
    //Rmut general red de abajo (+ 1 gen)
    double **rmutdowndup, **rmutdownc1, **rmutdownc2;
    hor.espacio(rmutdowndup, muestra, Nf);
    hor.espacio(rmutdownc1, muestra, Nf);
    hor.espacio(rmutdownc2, muestra, Nf);
    double **rmdowndup, **rmdownc1, **rmdownc2;
    hor.espacio(rmdowndup, muestra, 3);
    hor.espacio(rmdownc1, muestra, 3);
    hor.espacio(rmdownc2, muestra, 3);
    //Rmut fori red de arriba (+ 1 gen)
    double **rmutupdupori, **rmutupc1ori, **rmutupc2ori;
    hor.espacio(rmutupdupori, muestra, Nf);
    hor.espacio(rmutupc1ori, muestra, Nf);
    hor.espacio(rmutupc2ori, muestra, Nf);
    double **rmupdupori, **rmupc1ori, **rmupc2ori;
    hor.espacio(rmupdupori, muestra, 3);
    hor.espacio(rmupc1ori, muestra, 3);
    hor.espacio(rmupc2ori, muestra, 3);
    //Rmut fdif red de arriba (+ 1 gen)
    double **rmutupdupdif, **rmutupc1dif, **rmutupc2dif;
    hor.espacio(rmutupdupdif, muestra, Nf);
    hor.espacio(rmutupc1dif, muestra, Nf);
    hor.espacio(rmutupc2dif, muestra, Nf);
    double **rmupdupdif, **rmupc1dif, **rmupc2dif;
    hor.espacio(rmupdupdif, muestra, 3);
    hor.espacio(rmupc1dif, muestra, 3);
    hor.espacio(rmupc2dif, muestra, 3);

    //Rmut fori red de abajo (+ 1 gen)
    double **rmutdowndupori, **rmutdownc1ori, **rmutdownc2ori;
    hor.espacio(rmutdowndupori, muestra, Nf);
    hor.espacio(rmutdownc1ori, muestra, Nf);
    hor.espacio(rmutdownc2ori, muestra, Nf);
    double **rmdowndupori, **rmdownc1ori, **rmdownc2ori;
    hor.espacio(rmdowndupori, muestra, 3);
    hor.espacio(rmdownc1ori, muestra, 3);
    hor.espacio(rmdownc2ori, muestra, 3);
    //Rmut fdif red de abajo (+ 1 gen)
    double **rmutdowndupdif, **rmutdownc1dif, **rmutdownc2dif;
    hor.espacio(rmutdowndupdif, muestra, Nf);
    hor.espacio(rmutdownc1dif, muestra, Nf);
    hor.espacio(rmutdownc2dif, muestra, Nf);
    double **rmdowndupdif, **rmdownc1dif, **rmdownc2dif;
    hor.espacio(rmdowndupdif, muestra, 3);
    hor.espacio(rmdownc1dif, muestra, 3);
    hor.espacio(rmdownc2dif, muestra, 3);
    
    double **PromDistdup_up, **PromDistc1_up, **PromDistc2_up;
    hor.espacio(PromDistdup_up, muestra, Nf);
    hor.espacio(PromDistc1_up, muestra, Nf);
    hor.espacio(PromDistc2_up, muestra, Nf);
    double **PromDistdup_down, **PromDistc1_down, **PromDistc2_down;
    hor.espacio(PromDistdup_down, muestra, Nf);
    hor.espacio(PromDistc1_down, muestra, Nf);
    hor.espacio(PromDistc2_down, muestra, Nf);
    double **PDdup_up, **PDc1_up, **PDc2_up;
    hor.espacio(PDdup_up, muestra, 3);
    hor.espacio(PDc1_up, muestra, 3);
    hor.espacio(PDc2_up, muestra, 3);
    double **PDdup_down, **PDc1_down, **PDc2_down;
    hor.espacio(PDdup_down, muestra, 3);
    hor.espacio(PDc1_down, muestra, 3);
    hor.espacio(PDc2_down, muestra, 3);
    
    double **PromDistdup_up_tol, **PromDistdup_up_notol;
    hor.espacio(PromDistdup_up_tol, muestra, Nf);
    hor.espacio(PromDistdup_up_notol, muestra, Nf);
    double **PromDistc1_up_tol, **PromDistc1_up_notol;
    hor.espacio(PromDistc1_up_tol, muestra, Nf);
    hor.espacio(PromDistc1_up_notol, muestra, Nf);
    double **PromDistc2_up_tol, **PromDistc2_up_notol;
    hor.espacio(PromDistc2_up_tol, muestra, Nf);
    hor.espacio(PromDistc2_up_notol, muestra, Nf);
    
    double **PromDistdup_down_tol, **PromDistdup_down_notol;
    hor.espacio(PromDistdup_down_tol, muestra, Nf);
    hor.espacio(PromDistdup_down_notol, muestra, Nf);
    double **PromDistc1_down_tol, **PromDistc1_down_notol;
    hor.espacio(PromDistc1_down_tol, muestra, Nf);
    hor.espacio(PromDistc1_down_notol, muestra, Nf);
    double **PromDistc2_down_tol, **PromDistc2_down_notol;
    hor.espacio(PromDistc2_down_tol, muestra, Nf);
    hor.espacio(PromDistc2_down_notol, muestra, Nf);
    
    double **PDdup_up_tol, **PDc1_up_tol, **PDc2_up_tol;
    hor.espacio(PDdup_up_tol, muestra, 3);
    hor.espacio(PDc1_up_tol, muestra, 3);
    hor.espacio(PDc2_up_tol, muestra, 3);
    double **PDdup_down_tol, **PDc1_down_tol, **PDc2_down_tol;
    hor.espacio(PDdup_down_tol, muestra, 3);
    hor.espacio(PDc1_down_tol, muestra, 3);
    hor.espacio(PDc2_down_tol, muestra, 3);
    double **PDdup_up_notol, **PDc1_up_notol, **PDc2_up_notol;
    hor.espacio(PDdup_up_notol, muestra, 3);
    hor.espacio(PDc1_up_notol, muestra, 3);
    hor.espacio(PDc2_up_notol, muestra, 3);
    double **PDdup_down_notol, **PDc1_down_notol, **PDc2_down_notol;
    hor.espacio(PDdup_down_notol, muestra, 3);
    hor.espacio(PDc1_down_notol, muestra, 3);
    hor.espacio(PDc2_down_notol, muestra, 3);
    
    int up, down;
    ofstream redes, fo;
    ifstream fi;
    
    double alldists;
    
    int **Fn_esp, **Fn_dup, **Fn_c1, **Fn_c2;
    hor.espacio(Fn_esp, explorar, muestra);
    hor.espacio(Fn_dup, explorar, Nf);
    hor.espacio(Fn_c1, explorar, Nf);
    hor.espacio(Fn_c2, explorar, Nf);
    int *Fdif_esp, **Fdif_dup, **Fdif_c1, **Fdif_c2;
    hor.espacio(Fdif_esp, muestra);
    hor.espacio(Fdif_dup, muestra, Nf);
    hor.espacio(Fdif_c1, muestra, Nf);
    hor.espacio(Fdif_c2, muestra, Nf);
    
    double notoldists_noori;
    int notolrm_noori;
    
    int *evento_dup;
    hor.espacio(evento_dup, muestra);
    hor.open_ifstream(fi, "/outs/acces/eventos_dup.txt");
    for(int i = 0; i < muestra; i++){
        fi >> evento_dup[i];
    }
    fi.close();
    
    int *evento_c1;
    hor.espacio(evento_c1, muestra);
    hor.open_ifstream(fi, "/outs/acces/eventos_c1.txt");
    for(int i = 0; i < muestra; i++){
        fi >> evento_c1[i];
    }
    fi.close();
    
    int *evento_c2;
    hor.espacio(evento_c2, muestra);
    hor.open_ifstream(fi, "/outs/acces/eventos_c2.txt");
    for(int i = 0; i < muestra; i++){
        fi >> evento_c2[i];
    }
    fi.close();
        
    /////////////////////////////////CONDICIÓN INICIAL/////////////////////////////////
    //CI CERO
    for(int i = 0; i < Nf; i++){
        condinicial[i] = -1;
    }
    //OTRAS CI
// //     gat.condin(condinicial, Nf);
    
    for(int i = 0; i < Nf; i++){
        vectorentrada[i] = condinicial[i];
        trayectoria[0][i] = condinicial[i];
        vectoriniciod[i] = condinicial[i];
        Cvectorinicio[i] = condinicial[i];
        Cvectorinicio2[i] = condinicial[i];
    }
    
    hor.espacio(Wish, Np);
    
    Wish[0] = 1;
    Wish[1] = -1;
    Wish[2] = -1;
    Wish[3] = -1;
    Wish[4] = -1;
    Wish[5] = 1;
    
    //         cout << "////////////////////////////ELEGIR MATRIZ DE FACTORES CON PUNTO FIJO/////////////////////////////////////\n";
    
    //--OJO!, El que trayectoria sea finito puede provocar error a Nfs grandes--//
    
    especimenes.open("especimenes.txt");
    for(int i = 0; i < (muestra*(Nf + Np)); i++){
        for(int j = 0; j < Nf; j++){
            especimenes >> bichos[i][j];
        }
    }
    especimenes.close();
    
    //        cout << "///////////////////////////CAMINATA////////////////////////////////////\n";
    
    comp1[0] = 0;
    comp2[0] = 0;
    
    for(int m = 0; m < muestra; m++){
        cout << "Muestra: " << m << " evento " << evento_dup[m] << endl;
//         pasos = 0;
//         standby = 0;
        
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
        if(dist != 0){
            cout << "No coincide fenotipo\n";
        }
        
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
        gat.exploracion(Matrizf, Matrizp, trayectoria, Fenotipo, Distvecinos, Nf, Np, tolerancia, intnw, n, p, m, explorar, fdif, Fnuevom, valordekesp);//Muta red de muestra, saca distancia fenotípica y la guarda en una matriz           
        //               Matriz1: Factores
        //               Matriz2: Marcadores
        //               Matriz3: Trayectoria
        //               Matriz4: Fenotipo
        //               Matriz5: Distancias
        //               tam1 = hils1 = cols1 = cols2 = Nf
        //               tam2 = hils2 = Np
        //               num1: Número de muestra para saber en qué hilera va la distancia
        //               num2: Número de exploraciones, veces que se va a mutar la muestra, número de distancias a obtener
        
        //Fns 161019
        for(int hil = 0; hil < explorar; hil++){Fn_esp[hil][m] = valordekesp[hil];}
        Fdif_esp[m] = fdif;
        
        //Densidad
        densidadup[m] = hor.densidadparcial(Matrizf, Nf, Nf);
        densidaddown[m] = hor.densidadparcial(Matrizp, Np, Nf);
        densidad[m] = hor.densidadtotal(Matrizf, Nf, Nf, Matrizp, Np, Nf);

        //             cout << "//////////Accesibilidad después mutar, fmut diferentes al fori//////////\n";
        
        Accesm[m] = Fnuevom;
//        cout << Fnuevom << endl;
        
        //             cout << "///Robustez///\n";
        fdiferente[m] = (double)fdif / explorar;
        find.media(Distvecinos, Mediavecinos, explorar, m);
        find.desv(Distvecinos, Mediavecinos, Desvst, explorar, m);
        //General, arriba y abajo
        up = gat.rmutarriba(Matrizf, Nf, Nf, trayectoria, n, p, tolerancia, intnw, PromDist_up[m]);
        //Matriz1: Matriz a mutar
        //Matriz2: Trayectoria de la Matriz1 [n][nf]
        //Vector1: Condición inicial [Nf]
        //Nf: hils1, cols1
        //n original: hils2
        down = gat.rmutabajo(Matrizp, Np, Nf, Fenotipo, p, trayectoria, n, tolerancia, intnw, PromDist_down[m]);
        //Matriz1: Matriz a mutar
        //Matriz2: Fenotipo [p][np]
        //Matriz3: Trayectoria de la Matriz1 [n][nf]
        //Nf: cols1
        //Np: hils1
        //p: hils2
        //n: hils3
        
        rmutespup[m] = 1 - ((double)up / (2 * (Nf * Nf)));
        rmutespdown[m] = 1 - ((double)down / (2 * (Np * Nf)));
        
        //             cout << "//////////Shannon//////////\n";
        //             cout << "Antes\n";
        shannonesp[m] = gat.Shannon(valordekesp, explorar);
        
        fddup = 0;
        Cfd = 0;
        Cfd2 = 0;
        
        for(int gen = 0; gen < Nf; gen++){
//         for(int gen = 0; gen < 1; gen++){
//             cout << "Muestra: " << m << "Gen: " << gen << endl;
            //                 cout << "///////////////////////////DUPLICACIÓN////////////////////////////////////\n";
            gat.inducirduplicacion(Matrizfd, Matrizpd, Matrizf, Matrizp, gen, Nf, Np);
            //Matriz1: Factores duplicadas
            //Matriz2: Marcadores duplicadasFnuevo2
            //Matriz3: Factores
            //Matriz4: Marcadores
            //gen: gen a duplicar
            //Tam1: Nf
            //Tam2: Np
            
//             if(gen == evento_dup[m]){
//                 hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/especimenes_dup.txt");
//                 for(int hils_dup = 0; hils_dup < (Nf + 1); hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << Matrizfd[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 for(int hils_dup = 0; hils_dup < Np; hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << Matrizpd[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 fo << endl;
//                 fo.close();
//             }
            
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
                duppasos.open ("outs/duppasosdi.txt", ios::app);
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
            alldists = 0;
            notoldists_noori = 0;
            notolrm_noori = 0;
            hor.fillv0(valordekdup, explorar);
            gat.exploraciondup(Matrizfd, Matrizpd, trayectoriad, Fenotipod, Fenotipo, Distvecdupin, Distvecdup, alldists, notoldists_noori, notolrm_noori, Nf + 1, Np, tolerancia, intnw, nd, perd, p, gen, explorar, fdvecdup, fdvecdupin, regresa1, Fnuevovecdup, Fnuevovecdupin, 1, valordekdup, m);
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
            
//            if(gen == evento_dup[m]){
//                //01/01/21
//                hor.open_ofstream(fo, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/Rm/alldists_dup.txt");
//                fo << alldists << endl;
//                fo.close();
//            }
            //Fns 161019
            for(int hil = 0; hil < explorar; hil++){Fn_dup[hil][gen] = valordekdup[hil];}
            Fdif_dup[m][gen] = fdvecdup + fdvecdupin;
            
            //                 cout << "//////////Accesibilidad de duplicados después de mutación//////////\n";
            Accesvecdup[m][gen] = Fnuevovecdup;
            Accesvecdupin[m][gen] = Fnuevovecdupin;
            Accesdup[m][gen] = Fnuevovecdup + Fnuevovecdupin;
            Regreso1[m][gen] = regresa1;
            
            //                 cout << "///Robustez///\n";
            //General
            //             rmutdupgeneral[m][gen] = 1 - ((double)(fdvecdup - regresa1) / explorar); //Quité + fdvecdupin en numerador 250619
            //Parcial
            up = gat.rmutarriba(Matrizfd, Nf + 1, Nf + 1, trayectoriad, nd, perd, tolerancia, intnw, PromDistdup_up[m][gen]);
            down = gat.rmutabajo(Matrizpd, Np, Nf + 1, Fenotipod, perd, trayectoriad, nd, tolerancia, intnw, PromDistdup_down[m][gen]);
            rmutupdup[m][gen] = 1 - ((double)up / (2 * ((Nf + 1) * (Nf + 1))));
            rmutdowndup[m][gen] = 1- ((double)down / (2 * (Np * (Nf + 1))));
            
            //Vecinos de las redes duplicadas vs el f resultante de duplicar
            
            if(Distdup[gen] == 1){//Si el fenotipo es igual al original
                fdifvecdup[m][gen] = 1 - ((double)fdvecdup / explorar); //Fdif si es el mismo fenotipo
                rmutdupgeneral[m][gen] = fdifvecdup[m][gen]; //Quité + fdvecdupin en numerador 250619
                find.media(Distvecdup, Mediavecdup, explorar, m, gen);
                find.desv(Distvecdup, Mediavecdup, Desvstvecdup, explorar, m, gen);
                
                duppasos.open ("outs/duppasosdi.txt", ios::app);
                duppasos << "1 ";
                duppasos.close ();
                fdifvecdupin[m][gen] = 2;
                Mediavecdupin[m][gen] = 2;
                Desvstvecdupin[m][gen] = 2;
                
                rmutupdupori[m][gen] = rmutupdup[m][gen];
                rmutupdupdif[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdowndupori[m][gen] = rmutdowndup[m][gen];
                rmutdowndupdif[m][gen] = 2; //Para después poder descartar al sacar promedio
                PromDistdup_up_tol[m][gen] = PromDistdup_up[m][gen];
                PromDistdup_up_notol[m][gen] = 2;
                
                PromDistdup_down_tol[m][gen] = PromDistdup_down[m][gen];
                PromDistdup_down_notol[m][gen] = 2;
            }
            else{
                fdifvecdupin[m][gen] = 1 - ((double)fdvecdupin / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                rmutdupgeneral[m][gen] = (double) regresa1 / explorar;
                find.media(Distvecdupin, Mediavecdupin, explorar, m, gen);
                find.desv(Distvecdupin, Mediavecdupin, Desvstvecdupin, explorar, m, gen);
                
                rmutupdupori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutupdupdif[m][gen] = rmutupdup[m][gen];
                rmutdowndupori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdowndupdif[m][gen] = rmutdowndup[m][gen];
                
                fdifvecdup[m][gen] = 2;
                Mediavecdup[m][gen] = 2;
                Desvstvecdup[m][gen] = 2;
                
                PromDistdup_up_tol[m][gen] = 2;
                PromDistdup_up_notol[m][gen] = PromDistdup_up[m][gen];
                
                PromDistdup_down_tol[m][gen] = 2;
                PromDistdup_down_notol[m][gen] = PromDistdup_down[m][gen];
            }
            
            //General
            Distvecdup_gral[m][gen] = 1 - (alldists / explorar);
            hor.open_ofstream(fo, "/outs/notoldists_noori/sim_dup_noori.txt");
            fo << 1 - (notoldists_noori/explorar) << " ";
            fo.close();
            
            hor.open_ofstream(fo, "/outs/notoldists_noori/rm_dup_noori.txt");
            fo << notolrm_noori << " ";
            fo.close();
            
            
            //                 cout << "//////////Shannon//////////\n";
            shannondup[m][gen] = gat.Shannon(valordekdup, explorar);
            shannondup2[m][gen] = gat.Shannon2(valordekdup, explorar);
            
            //                 cout << "///////////////////////////CONTROL 1////////////////////////////////////\n";
            //Gen nuevo independiente
            
            gat.control1(CMatrizf, CMatrizp, Matrizf, Matrizp, Nf, Np, intnw);
            //Matriz1: Factores extendida
            //Matriz2: Marcadores extendida
            //Matriz3: Factores
            //Matriz4: Marcadores
            //Tam1: Nf
            //Tam2: Np
            //intnw
            
//             if(gen == evento_c1[m]){
//                 hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/especimenes_c1.txt");
//                 for(int hils_dup = 0; hils_dup < (Nf + 1); hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << CMatrizf[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 for(int hils_dup = 0; hils_dup < Np; hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << CMatrizp[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 fo << endl;
//                 fo.close();
//             }
            
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
                c1pasos.open ("outs/c1pasosdi.txt", ios::app);
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
            alldists = 0;
            notoldists_noori = 0;
            notolrm_noori = 0;
            hor.fillv0(valordekc1, explorar);
            gat.exploraciondup(CMatrizf, CMatrizp, Ctrayectoria, CFenotipo, Fenotipo, CDistvecin, CDistvec, alldists, notoldists_noori, notolrm_noori, Nf + 1, Np, tolerancia, intnw, Cn, Cp, p, gen, explorar, Cfdvec, Cfdvecin, regresa2, CFnuevovec, CFnuevovecin, 2, valordekc1, m);
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
            
//            if(gen == evento_c1[m]){
//                //01/01/21
//                hor.open_ofstream(fo, "/Users/yuridiaposadas/Documents/Doctorado/Resultados/outs070120/Rm/alldists_c1.txt");
//                fo << alldists << endl;
//                fo.close();
//            }
            
            //Fns 161019
            for(int hil = 0; hil < explorar; hil++){Fn_c1[hil][gen] = valordekc1[hil];}
            Fdif_c1[m][gen] = Cfdvec + Cfdvecin;
            
            //Vecinos de las redes duplicadas vs el f original
            
            //                 cout << "//////////Accesibilidad de c1 después de mutación//////////\n";
            CAccesvec[m][gen] = CFnuevovec;
            CAccesvecin[m][gen] = CFnuevovecin;
            CAcces[m][gen] = CFnuevovec + CFnuevovecin;
            Regreso2[m][gen] = regresa2;
            
            //                 cout << "///Robustez///\n";
            //General
//             rmutc1general[m][gen] = 1 - ((double)(Cfdvec + Cfdvecin - regresa2) / explorar);
            //Parcial
            up = gat.rmutarriba(CMatrizf, Nf + 1, Nf + 1, Ctrayectoria, Cn, Cp, tolerancia, intnw, PromDistc1_up[m][gen]);
            down = gat.rmutabajo(CMatrizp, Np, Nf + 1, CFenotipo, Cp, Ctrayectoria, Cn, tolerancia, intnw, PromDistc1_down[m][gen]);
            rmutupc1[m][gen] = 1 - ((double)up / (2 * ((Nf + 1) * (Nf + 1))));
            rmutdownc1[m][gen] = 1 - ((double)down / (2 * (Np * (Nf + 1))));
            
            if(CDist[gen] == 1){
                Cfdifvec[m][gen] = 1 - ((double)Cfdvec / explorar);
                rmutc1general[m][gen] = Cfdifvec[m][gen]; //Quité + fdvecdupin en numerador 250619
                find.media(CDistvec, CMediavec, explorar, m, gen);
                find.desv(CDistvec, CMediavec, CDesvstvec, explorar, m, gen);
                
                c1pasos.open ("outs/c1pasosdi.txt", ios::app);
                c1pasos << "1 ";
                c1pasos.close ();
                Cfdifvecin[m][gen] = 2;
                CMediavecin[m][gen] = 2;
                CDesvstvecin[m][gen] = 2;
                
                rmutupc1ori[m][gen] = rmutupc1[m][gen];
                rmutupc1dif[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdownc1ori[m][gen] = rmutdownc1[m][gen];
                rmutdownc1dif[m][gen] = 2; //Para después poder descartar al sacar promedio 
                PromDistc1_up_tol[m][gen] = PromDistc1_up[m][gen];
                PromDistc1_up_notol[m][gen] = 2;
                
                PromDistc1_down_tol[m][gen] = PromDistc1_down[m][gen];
                PromDistc1_down_notol[m][gen] = 2;
            }
            else{
                Cfdifvecin[m][gen] = 1 - ((double) Cfdvecin / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                rmutc1general[m][gen] = (double) regresa2 / explorar;
                find.media(CDistvecin, CMediavecin, explorar, m, gen);
                find.desv(CDistvecin, CMediavecin, CDesvstvecin, explorar, m, gen);
                
                rmutupc1ori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutupc1dif[m][gen] = rmutupc1[m][gen];
                rmutdownc1ori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdownc1dif[m][gen] = rmutdownc1[m][gen];
                
                Cfdifvec[m][gen] = 2;
                CMediavec[m][gen] = 2;
                CDesvstvec[m][gen] = 2;
                
                PromDistc1_up_tol[m][gen] = 2;
                PromDistc1_up_notol[m][gen] = PromDistc1_up[m][gen];
                
                PromDistc1_down_tol[m][gen] = 2;
                PromDistc1_down_notol[m][gen] = PromDistc1_down[m][gen];
            }
            
            //general
            Distvecc1_gral[m][gen] = 1 - (alldists / explorar);
            hor.open_ofstream(fo, "/outs/notoldists_noori/sim_c1_noori.txt");
            fo << 1 - (notoldists_noori/explorar) << " ";
            fo.close();
            hor.open_ofstream(fo, "/outs/notoldists_noori/rm_c1_noori.txt");
            fo << notolrm_noori << " ";
            fo.close();
            
            //                 cout << " //////////Shannon//////////\n";
            //                 cout << "CI\n";
            shannonc1[m][gen] = gat.Shannon(valordekc1, explorar);
            shannonc12[m][gen] = gat.Shannon2(valordekc1, explorar);
            
            //                 cout << "///////////////////////////CONTROL 2////////////////////////////////////\n";
            //Mismo número de conexiones del gen duplicado pero en orden diferente en el gen nuevo.
            
            gat.control2(CMatrizf2, CMatrizp2, Matrizfd, Matrizpd, Nf, Np);
            //Matriz1: Factores extendida
            //Matriz2: Marcadores extendida
            //Matriz3: Factores
            //Matriz4: Marcadores
            //Tam1: Nf
            //Tam2: Np
            //intnw
            
//             if(gen == evento_c2[m]){
//                 hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/especimenes_c2.txt");
//                 for(int hils_dup = 0; hils_dup < (Nf + 1); hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << CMatrizf2[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 for(int hils_dup = 0; hils_dup < Np; hils_dup++){
//                     for(int cols_dup = 0; cols_dup < (Nf + 1); cols_dup++){
//                         fo << CMatrizp2[hils_dup][cols_dup] << " ";
//                     }
//                     fo << endl;
//                 }
//                 fo << endl;
//                 fo.close();
//             }
            
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
                c2pasos.open ("outs/c2pasosdi.txt", ios::app);
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
            alldists = 0;
            notoldists_noori = 0;
            notolrm_noori = 0;
            hor.fillv0(valordekc2, explorar);
            gat.exploraciondup(CMatrizf2, CMatrizp2, Ctrayectoria2, CFenotipo2, Fenotipo, CDistvecin2, CDistvec2, alldists, notoldists_noori, notolrm_noori, Nf + 1, Np, tolerancia, intnw, Cn2, Cp2, p, gen, explorar, Cfdvec2, Cfdvecin2, regresa3, CFnuevovec2, CFnuevovecin2, 3, valordekc2, m);
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
            
//            if(gen == evento_c2[m]){
//                //01/01/21
//                hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs070120/Rm/alldists_c2.txt");
//                fo << alldists << endl;
//                fo.close();
//            }
            
            //Fns 161019
            for(int hil = 0; hil < explorar; hil++){Fn_c2[hil][gen] = valordekc2[hil];}
            Fdif_c2[m][gen] = Cfdvec2 + Cfdvecin2;
            
            //Vecinos de las redes duplicadas vs el f original
            
            //                 cout << "//////////Accesibilidad de c2 después de mutación//////////\n";
            CAccesvec2[m][gen] = CFnuevovec2;
            CAccesvecin2[m][gen] = CFnuevovecin2;
            CAcces2[m][gen] = CFnuevovec2 + CFnuevovecin2;
            Regreso3[m][gen] = regresa3;
            
            //                 cout << "///Robustez///\n";
            //General
//             rmutc2general[m][gen] = 1 - ((double)(Cfdvec2 + Cfdvecin2 - regresa3) / explorar);
            //Parcial
            up = gat.rmutarriba(CMatrizf2, Nf + 1, Nf + 1, Ctrayectoria2, Cn2, Cp2, tolerancia, intnw, PromDistc2_up[m][gen]);
            down = gat.rmutabajo(CMatrizp2, Np, Nf + 1, CFenotipo2, Cp2, Ctrayectoria2, Cn2, tolerancia, intnw, PromDistc2_down[m][gen]);
            rmutupc2[m][gen] = 1 - ((double)up / (2 * ((Nf + 1) * (Nf + 1))));
            rmutdownc2[m][gen] = 1 - ((double)down / (2 * (Np * (Nf + 1))));
            
            if(CDist2[gen] == 1){
                Cfdifvec2[m][gen] = 1 - ((double)Cfdvec2 / explorar);
                rmutc2general[m][gen] = Cfdifvec2[m][gen]; //Quité + fdvecdupin en numerador 250619
                find.media(CDistvec2, CMediavec2, explorar, m, gen);                          
                find.desv(CDistvec2, CMediavec2, CDesvstvec2, explorar, m, gen);
                
                c2pasos.open ("outs/c2pasosdi.txt", ios::app);
                c2pasos << "1 ";
                c2pasos.close ();
                Cfdifvecin2[m][gen] = 2;
                CMediavecin2[m][gen] = 2;
                CDesvstvecin2[m][gen] = 2;
                
                rmutupc2ori[m][gen] = rmutupc1[m][gen];
                rmutupc2dif[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdownc2ori[m][gen] = rmutdownc1[m][gen];
                rmutdownc2dif[m][gen] = 2; //Para después poder descartar al sacar promedio
                PromDistc2_up_tol[m][gen] = PromDistc2_up[m][gen];
                PromDistc2_up_notol[m][gen] = 2;
                
                PromDistc2_down_tol[m][gen] = PromDistc2_down[m][gen];
                PromDistc2_down_notol[m][gen] = 2;
            }
            else{
                Cfdifvecin2[m][gen] = 1 - ((double) Cfdvecin2 / explorar);//Esto es la proporción de redes duplicadas que al duplicar perdieron el fenotipo original pero que mantenieron en fenotipo nuevo después de mutar.
                rmutc2general[m][gen] = (double) regresa3 / explorar;
                find.media(CDistvecin2, CMediavecin2, explorar, m, gen);
                find.desv(CDistvecin2, CMediavecin2, CDesvstvecin2, explorar, m, gen);
                
                rmutupc2ori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutupc2dif[m][gen] = rmutupc2[m][gen];
                rmutdownc2ori[m][gen] = 2; //Para después poder descartar al sacar promedio
                rmutdownc2dif[m][gen] = rmutdownc2[m][gen];
                
                Cfdifvec2[m][gen] = 2;
                CMediavec2[m][gen] = 2;
                CDesvstvec2[m][gen] = 2;
                
                PromDistc2_up_tol[m][gen] = 2;
                PromDistc2_up_notol[m][gen] = PromDistc2_up[m][gen];
                
                PromDistc2_down_tol[m][gen] = 2;
                PromDistc2_down_notol[m][gen] = PromDistc2_down[m][gen];
            }
            
            //general
            Distvecc2_gral[m][gen] = 1 - (alldists / explorar);
            hor.open_ofstream(fo, "/outs/notoldists_noori/sim_c2_noori.txt");
            fo << 1 - (notoldists_noori/explorar) << " ";
            fo.close();
            hor.open_ofstream(fo, "/outs/notoldists_noori/rm_c2_noori.txt");
            fo << notolrm_noori << " ";
            fo.close();
            //                 cout << "//////////Shannon//////////\n";
            shannonc2[m][gen] = gat.Shannon(valordekc2, explorar);
            shannonc22[m][gen] = gat.Shannon2(valordekc2, explorar);
            
        }//Cierra ciclo de Genes
        
        hor.open_ofstream(fo, "/outs/notoldists_noori/sim_dup_noori.txt");
        fo << endl;
        fo.close();
        hor.open_ofstream(fo, "/outs/notoldists_noori/sim_c1_noori.txt");
        fo << endl;
        fo.close();
        hor.open_ofstream(fo, "/outs/notoldists_noori/sim_c2_noori.txt");
        fo << endl;
        fo.close();
        
        hor.open_ofstream(fo, "/outs/notoldists_noori/rm_dup_noori.txt");
        fo << endl;
        fo.close();
        hor.open_ofstream(fo, "/outs/notoldists_noori/rm_c1_noori.txt");
        fo << endl;
        fo.close();
        hor.open_ofstream(fo, "/outs/notoldists_noori/rm_c2_noori.txt");
        fo << endl;
        fo.close();
        
        
        //070120 <<< Para tener la distancia después de agregar un gen con original
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/dist/dist_tolerancia_dup.txt");
//        for(int hils = 0; hils < Nf; hils++){
//            fo << Distdup[hils] << " ";
//        }
//        fo << endl;
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/dist/dist_tolerancia_c1.txt");
//        for(int hils = 0; hils < Nf; hils++){
//            fo << CDist[hils] << " ";
//        }
//        fo << endl;
//        fo.close();
//
//        hor.open_ofstream(fo, "/home/rdaneel/Documentos/Doctorado/Resultados/outs/dist/dist_tolerancia_c2.txt");
//        for(int hils = 0; hils < Nf; hils++){
//            fo << CDist2[hils] << " ";
//        }
//        fo << endl;
//        fo.close();
        
        
        hor.open_ofstream(fo, "/outs/acces/valordek/valordek_dup_"+ hor.inttostring(m) +".txt");
        for(int hils = 0; hils < explorar; hils++){
            for(int cols = 0; cols < Nf; cols++){
                fo << Fn_dup[hils][cols] << " ";
            }
            fo << endl;
        }
        fo.close();

        hor.open_ofstream(fo, "/outs/acces/valordek/valordek_c1_"+ hor.inttostring(m) +".txt");
        for(int hils = 0; hils < explorar; hils++){
            for(int cols = 0; cols < Nf; cols++){
                fo << Fn_c1[hils][cols] << " ";
            }
            fo << endl;
        }
        fo.close();

        hor.open_ofstream(fo, "/outs/acces/valordek/valordek_c2_"+ hor.inttostring(m) +".txt");
        for(int hils = 0; hils < explorar; hils++){
            for(int cols = 0; cols < Nf; cols++){
                fo << Fn_c2[hils][cols] << " ";
            }
            fo << endl;
        }
        fo.close();
        
        duppasos.open("outs/duppasosdi.txt", ios::app);
        duppasos << endl;
        duppasos.close();
        c1pasos.open("outs/c1pasosdi.txt", ios::app);
        c1pasos << endl;
        c1pasos.close();
        c2pasos.open("outs/c2pasosdi.txt", ios::app);
        c2pasos << endl;
        c2pasos.close();
        
        comparaciones.open ("outs/comparaciones.txt", ios::app);
        comparaciones << comp1[0] << " " << comp2[0] << "\n";
        comparaciones.close ();
        fdifdup[m] = (double)fddup / Nf;
        Cfdif[m] = (double)Cfd / Nf; 
        Cfdif2[m] = (double)Cfd2 / Nf;
        
//         hor.sumarmatriz(Accesvecdup, Avecdup, muestra, Nf);
//         hor.sumarmatriz(Accesvecdupin, Avecdupin, muestra, Nf);
//         hor.sumarmatriz(CAccesvec, CAvec, muestra, Nf);
//         hor.sumarmatriz(CAccesvecin, CAvecin, muestra, Nf);
//         hor.sumarmatriz(CAccesvec2, CAvec2, muestra, Nf);
//         hor.sumarmatriz(CAccesvecin2, CAvecin2, muestra, Nf);
        
    }//Cierra ciclo de Muestras
    
    hor.open_ofstream(fo, "/outs/acces/valordek/valordek_esp.txt");
    for(int hils = 0; hils < explorar; hils++){
        for(int cols = 0; cols < muestra; cols++){
            fo << Fn_esp[hils][cols] << " ";
        }
        fo << endl;
    }
    fo.close();

    hor.open_ofstream(fo, "/outs/acces/valordek/Fdif_esp.txt");
    for(int hils = 0; hils < muestra; hils++){
        fo << Fdif_esp[hils]<< endl;
    }
    fo.close();


    hor.open_ofstream(fo, "/outs/acces/valordek/Fdif_dup.txt");
    for(int hils = 0; hils < muestra; hils++){
        for(int cols = 0; cols < Nf; cols++){
            fo << Fdif_dup[hils][cols] << " ";
        }
        fo << endl;
    }
    fo.close();


    hor.open_ofstream(fo, "/outs/acces/valordek/Fdif_c1.txt");
    for(int hils = 0; hils < muestra; hils++){
        for(int cols = 0; cols < Nf; cols++){
            fo << Fdif_c1[hils][cols] << " ";
        }
        fo << endl;
    }
    fo.close();


    hor.open_ofstream(fo, "/outs/acces/valordek/Fdif_c2.txt");
    for(int hils = 0; hils < muestra; hils++){
        for(int cols = 0; cols < Nf; cols++){
            fo << Fdif_c2[hils][cols] << " ";
        }
        fo << endl;
    }
    fo.close();

             cout << "------------Condensamos matrices que van por factor a mín, máx y promedio------------\n";
    
    find.minmaxmedin(rmutdupgeneral, rmdg, muestra, Nf);
    find.minmaxmedin(rmutc1general, rmc1g, muestra, Nf);
    find.minmaxmedin(rmutc2general, rmc2g, muestra, Nf);
    
    find.minmaxmedin(Distvecdup_gral, dvdg, muestra, Nf);
    find.minmaxmedin(Distvecc1_gral, dvc1g, muestra, Nf);
    find.minmaxmedin(Distvecc2_gral, dvc2g, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_up, PDdup_up, muestra, Nf);
    find.minmaxmedin(PromDistc1_up, PDc1_up, muestra, Nf);
    find.minmaxmedin(PromDistc2_up, PDc2_up, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_down, PDdup_down, muestra, Nf);
    find.minmaxmedin(PromDistc1_down, PDc1_down, muestra, Nf);
    find.minmaxmedin(PromDistc2_down, PDc2_down, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_up_tol, PDdup_up_tol, muestra, Nf);
    find.minmaxmedin(PromDistc1_up_tol, PDc1_up_tol, muestra, Nf);
    find.minmaxmedin(PromDistc2_up_tol, PDc2_up_tol, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_up_notol, PDdup_up_notol, muestra, Nf);
    find.minmaxmedin(PromDistc1_up_notol, PDc1_up_notol, muestra, Nf);
    find.minmaxmedin(PromDistc2_up_notol, PDc2_up_notol, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_down_tol, PDdup_down_tol, muestra, Nf);
    find.minmaxmedin(PromDistc1_down_tol, PDc1_down_tol, muestra, Nf);
    find.minmaxmedin(PromDistc2_down_tol, PDc2_down_tol, muestra, Nf);
    
    find.minmaxmedin(PromDistdup_down_notol, PDdup_down_notol, muestra, Nf);
    find.minmaxmedin(PromDistc1_down_notol, PDc1_down_notol, muestra, Nf);
    find.minmaxmedin(PromDistc2_down_notol, PDc2_down_notol, muestra, Nf);
    
    find.minmaxmedin(rmutupdup, rmupdup, muestra, Nf);
    find.minmaxmedin(rmutupc1, rmupc1, muestra, Nf);
    find.minmaxmedin(rmutupc2, rmupc2, muestra, Nf);
    
    find.minmaxmedin(rmutdowndup, rmdowndup, muestra, Nf);
    find.minmaxmedin(rmutdownc1, rmdownc1, muestra, Nf);
    find.minmaxmedin(rmutdownc2, rmdownc2, muestra, Nf);
    
    find.minmaxmedin(rmutupdupori, rmupdupori, muestra, Nf); //270219
    find.minmaxmedin(rmutupc1ori, rmupc1ori, muestra, Nf);
    find.minmaxmedin(rmutupc2ori, rmupc2ori, muestra, Nf);
    
    find.minmaxmedin(rmutupdupdif, rmupdupdif, muestra, Nf); //270219
    find.minmaxmedin(rmutupc1dif, rmupc1dif, muestra, Nf);
    find.minmaxmedin(rmutupc2dif, rmupc2dif, muestra, Nf);
    
    find.minmaxmedin(rmutdowndupori, rmdowndupori, muestra, Nf); //270219
    find.minmaxmedin(rmutdownc1ori, rmdownc1ori, muestra, Nf);
    find.minmaxmedin(rmutdownc2ori, rmdownc2ori, muestra, Nf);
    
    find.minmaxmedin(rmutdowndupdif, rmdowndupdif, muestra, Nf); //270219
    find.minmaxmedin(rmutdownc1dif, rmdownc1dif, muestra, Nf);
    find.minmaxmedin(rmutdownc2dif, rmdownc2dif, muestra, Nf);
    
    find.minmaxmedin(fdifvecdup, pdv, muestra, Nf);
    find.minmaxmedin(fdifvecdupin, pdvin, muestra, Nf);
    find.minmaxmedin(Cfdifvec, pc1v, muestra, Nf);
    find.minmaxmedin(Cfdifvecin, pc1vin, muestra, Nf);
    find.minmaxmedin(Cfdifvec2, pc2v, muestra, Nf);
    find.minmaxmedin(Cfdifvecin2, pc2vin, muestra, Nf);
    
    find.minmaxmedin(Mediavecdup, mdv, muestra, Nf);
    find.minmaxmedin(Mediavecdupin, mdvin, muestra, Nf);
    find.minmaxmedin(CMediavec, mc1v, muestra, Nf);
    find.minmaxmedin(CMediavecin, mc1vin, muestra, Nf);
    find.minmaxmedin(CMediavec2, mc2v, muestra, Nf);
    find.minmaxmedin(CMediavecin2, mc2vin, muestra, Nf);
    
    find.minmaxmedin(Desvstvecdup, ddv, muestra, Nf);
    find.minmaxmedin(Desvstvecdupin, ddvin, muestra, Nf);
    find.minmaxmedin(CDesvstvec, dc1v, muestra, Nf);
    find.minmaxmedin(CDesvstvecin, dc1vin, muestra, Nf);
    find.minmaxmedin(CDesvstvec2, dc2v, muestra, Nf);
    find.minmaxmedin(CDesvstvecin2, dc2vin, muestra, Nf);
    
    proporciones.open ("outs/proporciones.txt", ios::app);
    robmutdup.open("outs/Rm/robmutdup.txt", ios::app);
    robmutc1.open("outs/Rm/robmutc1.txt", ios::app);
    robmutc2.open("outs/Rm/robmutc2.txt", ios::app);
    robmutupdup.open("outs/Rm/robmutupdup.txt", ios::app);
    robmutupc1.open("outs/Rm/robmutupc1.txt", ios::app);
    robmutupc2.open("outs/Rm/robmutupc2.txt", ios::app);
    robmutdowndup.open("outs/Rm/robmutdowndup.txt", ios::app);
    robmutdownc1.open("outs/Rm/robmutdownc1.txt", ios::app);
    robmutdownc2.open("outs/Rm/robmutdownc2.txt", ios::app);
    regresan1.open("outs/Rm/regresan1.txt", ios::app);
    regresan2.open("outs/Rm/regresan2.txt", ios::app);
    regresan3.open("outs/Rm/regresan3.txt", ios::app);
    densidades.open("outs/densidad.txt", ios::app);
    distancias.open("outs/distancias.txt", ios::app);
    distanciasdup.open("outs/dist/distdup.txt", ios::app);
    distanciasc1.open("outs/dist/distc1.txt", ios::app);
    distanciasc2.open("outs/dist/distc2.txt", ios::app);
    desviaciones.open("outs/desviaciones.txt", ios::app);
    accesibilidad.open("outs/acces/accesibilidad.txt", ios::app);
    accesibilidad1.open("outs/acces/accesesp.txt", ios::app);
    accesibilidad2.open("outs/acces/accesdupori.txt", ios::app);
    accesibilidad3.open("outs/acces/accesdupdif.txt", ios::app);
    accesibilidad4.open("outs/acces/accesc1ori.txt", ios::app);
    accesibilidad5.open("outs/acces/accesc1dif.txt", ios::app);
    accesibilidad6.open("outs/acces/accesc2ori.txt", ios::app);
    accesibilidad7.open("outs/acces/accesc2dif.txt", ios::app);
    accesibilidad8.open("outs/acces/accesdup.txt", ios::app);
    accesibilidad9.open("outs/acces/accesc1.txt", ios::app);
    accesibilidad10.open("outs/acces/accesc2.txt", ios::app);
    accesibilidad11.open("outs/acces/regdup.txt", ios::app);
    accesibilidad12.open("outs/acces/regc1.txt", ios::app);
    accesibilidad13.open("outs/acces/regc2.txt", ios::app);
    
    cout << "------------Archivos generándose------------\n";
    
    //Tabla
    //         proporciones << "muestra Caminata Robmut Robdup RobC1 RobC2 Robmutdupmín Robmutdupmáx Robmutdupmed Robinmutdupmín Robinmutdupmáx Robinmutdupmed Robmutc1mín Robmutc1máx Robmutc1med Robinmutc1mín Robinmutc1máx Robinmutc1med Robmutc2mín Robmutc2máx Robmutc2med Robinmutc2mín Robinmutc2máx Robinmutc2med\n";
    
    for(int i = 0; i < muestra; i++){
        proporciones << 1 - fdiferente[i] << " " << 1 - fdifdup[i] << " " << 1 - Cfdif[i] << " " <<  1 - Cfdif2[i] << " ";
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
        proporciones << rmutespup[i] << " " << rmutespdown[i] << " ";
        for(int j = 0; j < 3; j++){
            proporciones << rmdg[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmc1g[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmc2g[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupdup[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc1[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc2[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdowndup[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc1[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc2[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupdupori[i][j] << " "; //270219
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc1ori[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc2ori[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdowndupori[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc1ori[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc2ori[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupdupdif[i][j] << " "; //270219
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc1dif[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmupc2dif[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdowndupdif[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc1dif[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            proporciones << rmdownc2dif[i][j] << " ";
        }
        proporciones << "\n";
    }
    
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            robmutdup << rmutdupgeneral[i][j] << " ";
            robmutc1 << rmutc1general[i][j] << " ";
            robmutc2 << rmutc2general[i][j] << " ";
            robmutupdup << rmutupdup[i][j] << " ";
            robmutupc1 << rmutupc1[i][j] << " ";
            robmutupc2 << rmutupc2[i][j] << " ";
            robmutdowndup << rmutdowndup[i][j] << " ";
            robmutdownc1 << rmutdownc1[i][j] << " ";
            robmutdownc2 << rmutdownc2[i][j] << " ";
            regresan1 << Regreso1[i][j] << " ";
            regresan2 << Regreso2[i][j] << " ";
            regresan3 << Regreso3[i][j] << " ";
        }
        robmutdup << "\n";
        robmutc1 << "\n";
        robmutc2 << "\n";
        robmutupdup << "\n";
        robmutupc1 << "\n";
        robmutupc2 << "\n";
        robmutdowndup << "\n";
        robmutdownc1 << "\n";
        robmutdownc2 << "\n";
        regresan1 << "\n";
        regresan2 << "\n";
        regresan3 << "\n";
    }
    
    for(int i = 0; i < muestra; i++){
        densidades << densidad[i] << " " << densidadup[i] << " " << densidaddown[i] << endl;
    }
    
    for(int i = 0; i < muestra; i++){
        distancias << Mediavecinos[i] << " " << Mediadup[i] << " " << CMedia[i] << " " << CMedia2[i] << " ";
        for(int j = 0; j < 3; j++){
            distancias << mdv[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << mdvin[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << mc1v[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << mc1vin[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << mc2v[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << mc2vin[i][j] << " ";
        }
        distancias << PromDist_up[i] << " " << PromDist_down[i] << " ";
        for(int j = 0; j < 3; j++){
            distancias << dvdg[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << dvc1g[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << dvc2g[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_up[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_up[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_up[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_down[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_down[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_down[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_up_tol[i][j] << " "; //270219
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_up_tol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_up_tol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_down_tol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_down_tol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_down_tol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_up_notol[i][j] << " "; //270219
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_up_notol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_up_notol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDdup_down_notol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc1_down_notol[i][j] << " ";
        }
        for(int j = 0; j < 3; j++){
            distancias << PDc2_down_notol[i][j] << " ";
        }
        distancias << "\n";
    }
    
    for(int i = 0; i < muestra; i++){
        for(int j = 0; j < Nf; j++){
            if(Mediavecdup[i][j] == 2){distanciasdup << Mediavecdupin[i][j] << " ";}
            else{distanciasdup << Mediavecdup[i][j] << " ";}
            if(CMediavec[i][j] == 2){distanciasc1 << CMediavecin[i][j] << " ";}
            else{distanciasc1 << CMediavec[i][j] << " ";}
            if(CMediavec2[i][j] == 2){distanciasc2 << CMediavecin2[i][j] << " ";}
            else{distanciasc2 << CMediavec2[i][j] << " ";}
        }
        distanciasdup << "\n";
        distanciasc1 << "\n";
        distanciasc2 << "\n";
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
    robmutdup.close ();
    robmutc1.close ();
    robmutc2.close ();
    robmutupdup.close ();
    robmutupc1.close ();
    robmutupc2.close ();
    robmutdowndup.close ();
    robmutdownc1.close ();
    robmutdownc2.close ();
    regresan1.close ();
    regresan2.close ();
    regresan3.close ();
    densidades.close ();
    distancias.close ();
    distanciasdup.close ();
    distanciasc1.close ();
    distanciasc2.close ();
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
    hor.borrar(Fenotipo, 1); //Recuerda que en condin esto puede variar
    hor.borrar(Fenotipod, perd);
    hor.borrar(Wish);
    hor.borrar(condinicial);
    hor.borrar(vectorentrada);
    hor.borrar(vectoriniciod);
    hor.borrar(vectorentradad);
    hor.borrar(atractor);
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
//     hor.borrar(Avecdup);
//     hor.borrar(Avecdupin);
//     hor.borrar(CAvec2);
//     hor.borrar(CAvecin2);
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
    hor.borrar(densidadup);
    hor.borrar(densidaddown);
    hor.borrar(densidad);
    hor.borrar(rmutespup);
    hor.borrar(rmutespdown);
    hor.borrar(bichos, (muestra*(Nf + Np)));
    hor.borrar(rmutdupgeneral, muestra);
    hor.borrar(rmutc1general, muestra);
    hor.borrar(rmutc2general, muestra);
    hor.borrar(rmdg, muestra);
    hor.borrar(rmc1g, muestra);
    hor.borrar(rmc2g, muestra);
    hor.borrar(rmutupdup, muestra);
    hor.borrar(rmutupc1, muestra);
    hor.borrar(rmutupc2, muestra);
    hor.borrar(rmutdowndup, muestra);
    hor.borrar(rmutdownc1, muestra);
    hor.borrar(rmutdownc2, muestra);
    hor.borrar(rmupdup, muestra);
    hor.borrar(rmupc1, muestra);
    hor.borrar(rmupc2, muestra);
    hor.borrar(rmdowndup, muestra);
    hor.borrar(rmdownc1, muestra);
    hor.borrar(rmdownc2, muestra);
    hor.borrar(rmutdowndupori, muestra);//270219
    hor.borrar(rmutdownc1ori, muestra);
    hor.borrar(rmutdownc2ori, muestra);
    hor.borrar(rmdowndupori, muestra);
    hor.borrar(rmdownc1ori, muestra);
    hor.borrar(rmdownc2ori, muestra);
    hor.borrar(rmutdowndupdif, muestra);
    hor.borrar(rmutdownc1dif, muestra);
    hor.borrar(rmutdownc2dif, muestra);
    hor.borrar(rmdowndupdif, muestra);
    hor.borrar(rmdownc1dif, muestra);
    hor.borrar(rmdownc2dif, muestra);
    hor.borrar(PromDist_up);
    hor.borrar(PromDist_down);
    hor.borrar(Distvecdup_gral, muestra);
    hor.borrar(Distvecc1_gral, muestra);
    hor.borrar(Distvecc2_gral, muestra);
    hor.borrar(dvdg, muestra);
    hor.borrar(dvc1g, muestra);
    hor.borrar(dvc2g, muestra);
    hor.borrar(PromDistdup_up, muestra);
    hor.borrar(PromDistc1_up, muestra);
    hor.borrar(PromDistc2_up, muestra);
    hor.borrar(PromDistdup_down, muestra);
    hor.borrar(PromDistc1_down, muestra);
    hor.borrar(PromDistc2_down, muestra);
    hor.borrar(PDdup_up, muestra);
    hor.borrar(PDc1_up, muestra);
    hor.borrar(PDc2_up, muestra);
    hor.borrar(PDdup_down, muestra);
    hor.borrar(PDc1_down, muestra);
    hor.borrar(PDc2_down, muestra);
    hor.borrar(PromDistdup_up_tol, muestra);
    hor.borrar(PromDistdup_up_notol, muestra);
    hor.borrar(PromDistc1_up_tol, muestra);
    hor.borrar(PromDistc1_up_notol, muestra);
    hor.borrar(PromDistc2_up_tol, muestra);
    hor.borrar(PromDistc2_up_notol, muestra);
    hor.borrar(PromDistdup_down_tol, muestra);
    hor.borrar(PromDistdup_down_notol, muestra);
    hor.borrar(PromDistc1_down_tol, muestra);
    hor.borrar(PromDistc1_down_notol, muestra);
    hor.borrar(PromDistc2_down_tol, muestra);
    hor.borrar(PromDistc2_down_notol, muestra);
    hor.borrar(PDdup_up_tol, muestra);
    hor.borrar(PDc1_up_tol, muestra);
    hor.borrar(PDc2_up_tol, muestra);
    hor.borrar(PDdup_down_tol, muestra);
    hor.borrar(PDc1_down_tol, muestra);
    hor.borrar(PDc2_down_tol, muestra);
    hor.borrar(PDdup_up_notol, muestra);
    hor.borrar(PDc1_up_notol, muestra);
    hor.borrar(PDc2_up_notol, muestra);
    hor.borrar(PDdup_down_notol, muestra);
    hor.borrar(PDc1_down_notol, muestra);
    hor.borrar(PDc2_down_notol, muestra);
    hor.borrar(Fn_esp, explorar);
    hor.borrar(Fn_c1, explorar);
    hor.borrar(Fn_c2, explorar);
    hor.borrar(Fdif_esp);
    hor.borrar(Fdif_dup, muestra);
    hor.borrar(Fdif_c1, muestra);
    hor.borrar(Fdif_c2, muestra);
    hor.borrar(evento_dup);
    hor.borrar(evento_c1);
    hor.borrar(evento_c2);
    //#########################################################################################################//
    
    cout << "------------Ejecución terminada!!-----------\n";
    
    hor.close_rng();
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    return 0;
}
