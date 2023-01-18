/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#ifndef __SISLIN_H__
#define __SISLIN_H__

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

#define ALOCACAO 64
#define TAM_PADRAO  10
#define STRING_SIZE 100
#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
     // coeficientes nos sistemas lineares.

typedef struct {
int n;
double **A;
double **L;
double **AI;
double **AC;
double **I;
} SL;

SL* criaNovoSl (int n); // Aloca espaço na memória para a struct SL.
void liberaSisLin(SL *SL);
int ehInvertivel(SL *sl);
void geraIdentidade(double **I, int tam);
int refinamento (SL *SL, double *tTotal, int iteracoes, FILE* file);
void criaL(double **L, double *m, int tam, double *tTotal);
void fatoracaoLU(double **restrict I, double **restrict AI, double **restrict L, double **restrict A, int tam, double *tTotal);
int eliminacaoGauss (SL *sl, double *tTotal);
void iniciarSlEntrada(FILE* input, SL *sl);
void geraSlAleatoria (SL *SL, real_t coef_max);
void imprimeSL(double **A, int tam);
void imprimeSlSaida(double **A, int tam, FILE * output);
void pegarcoluna(double *restrict y, double **restrict R, int j, int tam);

#endif
