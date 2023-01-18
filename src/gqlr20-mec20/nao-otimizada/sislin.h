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

#define ERRO 0.000000008
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

// Tipos de matrizes de coeficientes usados pela função 'inicializaSistLinear()'
typedef enum {
    generico = 0,
    hilbert,
    diagDominante,
    eqNula,
    eqProporcional,
    eqCombLinear
} tipoSL;

// Tipo de alocação para matrizes
typedef enum {
  pontPont=0, // Matriz como vetor de N ponteiros para vetores de tamanho N
  pontVet     // Matriz como vetor de N ponteiros para um único vetor de tamanho N*N
} tipoAloc_t;
void liberaSisLin(SL *SL);
int ehInvertivel(SL *sl);
void geraIdentidade(SL *SL);
int refinamento (SL *SL, double *tTotal, int iteracoes, FILE* file);
void fatoracaoLU(double **AI, double **L, double **A, int tam, double *tTotal);
int eliminacaoGauss (double **A, double **L, double **AC, double **I, int tam, double *tTotal);
char *nomesArquivos(int argc, char *argv[]);
void criaNovoSl(SL *sl, int tamanho); // Aloca espaço na memória para a struct SL.
FILE *abrirEntrada(int argc, char *argv[]);
void iniciarSlEntrada(FILE* input, SL *sl);
void geraSlAleatoria (SL *SL, tipoSL tipo, real_t coef_max);
void imprimeSL(double **A, int tam);
void imprimeSlSaida(double **A, int tam, FILE * output);

#endif
