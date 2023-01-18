/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#include <likwid.h>
#include "utils.h"
#include "sislin.h"
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

void criaNovoSl(SL *sl, int tamanho) { // Aloca espaço na memória para a struct SL.
        int i;
	sl->n=tamanho;
	if (!(sl->A=malloc(sizeof(void **) * tamanho))) {
                perror ("Erro alocacao Sl->A: ");
                exit (1);
        }
	if (!(sl->AI=malloc(sizeof(void **) * tamanho))) {
                perror ("Erro alocacao Sl->AI: ");
                exit (1);
        }
	if (!(sl->L=(double **)malloc(sizeof(double *)*tamanho))) {
                perror ("Erro alocacao Sl->L: ");
                exit (1);
        }
	if (!(sl->AC=malloc(sizeof(double **)*tamanho))) {
                perror ("Erro alocacao Sl->AC: ");
                exit (1);
        }
	if (!(sl->I=malloc(sizeof(double **)*tamanho))) {
                perror ("Erro alocacao Sl->I: ");
                exit (1);
        }
  	for (i=0;i<tamanho;i++) {
    		if (!(sl->A[i] = malloc(sizeof(void *) * tamanho))) {
                        perror ("Erro alocacao Sl->A[i]: ");
                        exit (1);
                }
		if (!(sl->AI[i] = malloc(sizeof(void *) * tamanho))) {
                        perror ("Erro alocacao Sl->AI[i]: ");
                        exit (1);
                }
		if (!(sl->L[i] = malloc(sizeof(void *) * tamanho))) {
                        perror ("Erro alocacao Sl->L[i]: ");
                        exit (1);
                }
  		if (!(sl->AC[i] = malloc(sizeof(void *) * tamanho))) {
                        perror ("Erro alocacao Sl->AC[i]: ");
                        exit (1);
                }
		if (!(sl->I[i] = malloc(sizeof(void *) * tamanho))) {
                        perror ("Erro alocacao Sl->I[i]: ");
                        exit (1);
                }
	}
}

void geraSlAleatoria (SL *SL, tipoSL tipo, real_t coef_max) { // Gera uma SL com valores aleatórios
  unsigned int n = SL->n;
  // para gerar valores no intervalo [0,coef_max]
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);
    
  if (tipo == hilbert) {
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
        SL->A[i][j] = 1.0 / (real_t)(i+j+1);
	SL->AC[i][j] = SL->A[i][j];
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
             SL->A[i][j] = (real_t)rand() * invRandMax;
	     SL->AC[i][j] = SL->A[i][j];
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % n;
      for (unsigned int j=0; j<n; ++j) {
             SL->A[nula][j] = 0.0;
      }
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % n;
      unsigned int propSrc = (propDst + 1) % n;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[propDst][j] = SL->A[propSrc][j] * mult;
      }
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % n;
      unsigned int combSrc1 = (combDst + 1) % n;
      unsigned int combSrc2 = (combDst + 2) % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[combDst][j] = SL->A[combSrc1][j] + SL->A[combSrc2][j];
      }
    }
    else if (tipo == diagDominante) {
      // aumenta o valor dos termos da diagonal principal
      for (unsigned int i=0; i<n; ++i) {
        real_t soma = 0.0;
        for (unsigned int j=0; j < i; ++j) soma += SL->A[i][j];
        for (unsigned int j=i+1; j < n; ++j) soma += SL->A[i][j];
        SL->A[i][i] += soma;
      }
    }
  }
}

void iniciarSlEntrada(FILE* input, SL *sl) {
  	int i, j;

  	for(i=0;i<sl->n;i++) {
		for(j=0;j<sl->n;j++)
    			fscanf(input,"%le", &sl->A[i][j]);
  	}
	for(i=0;i<sl->n;i++) {
                for(j=0;j<sl->n;j++)
                        sl->AC[i][j]=sl->A[i][j];
        }
}

int encontraMax(double **A, int i, int tam){
        int iMaior=i;
        for(int j=i+1;j<tam;j++){
                if(A[j][i] > A[iMaior][i]){
                        iMaior = j;
                }
        }
        return iMaior;
}

void trocaLinha(double **A, int i, int iPivo, int tam){
        double *tmp = A[i];
        A[i] = A[iPivo];
        A[iPivo] = tmp;
}

int eliminacaoGauss (double **A, double **L, double **AC, double **I, int tam, double *tTotal){
        *tTotal = timestamp();

        for (int i = 0; i < tam; i++)
                for (int j = 0; j < tam; j++)
                        L[i][j] = 0;

        for(int i=0;i<tam-1;i++){
                int iPivo = encontraMax(A, i, tam);
                if(i != iPivo){
                        trocaLinha(A, i, iPivo, tam);
                        trocaLinha(AC, i, iPivo, tam);
                        trocaLinha(I, i, iPivo, tam);
                        trocaLinha(L, i, iPivo, tam);
		}
                for(int k=i+1;k<tam;k++){
                        double m = A[k][i] / A[i][i];
                        L[k][i]=m;
                        if (isnan (m) || isinf (m)) {
                                perror ("Erro, divisao por 0 na ElinacaoGauss: ");
                                exit (1);
                        }
                        A[k][i] = 0.0;
                        for(int j=i+1;j<tam;++j)
                                A[k][j] -= A[i][j] * m;
                }
        }

        for (int i = 0; i < tam; i++)
                for (int j = 0; j < tam; j++)
                        if (i == j)
                                L[i][j] = 1;

        *tTotal = timestamp() - *tTotal;
        return 0;
}   

void retrossubsU(double **A, double *b, double *x, int tam){
        for(int i=tam-1;i>=0;--i){
                x[i]=b[i];
                for(int j=i+1;j<tam;++j){
                        x[i] = x[i] - (A[i][j] * x[j]);
                }
                x[i] = x[i]/A[i][i];
                if (isnan (x[i]) || isinf (x[i])) {
                        perror ("Erro, divisao por 0 na retrossubs: ");
                        exit (1);
                }
        }
}

void retrossubsL(double **A, double *b, double *x, int tam){
        for(int i=0;i<tam;i++){
                x[i]=b[i];
                for(int j=0;j<i;j++){
                        x[i] = x[i] - (A[i][j] * x[j]);
                }
                x[i] = x[i]/A[i][i];
                if (isnan (x[i]) || isinf (x[i])) {
                        perror ("Erro, divisao por 0 na retrossubs: ");
                        exit (1);
                }
        }
}

void fatoracaoLU(double **AI, double **L, double **A, int tam, double *tTotal){
        *tTotal = timestamp();
        char nome[20]="teste_fatoracaoLU_";
	char text[20];
	sprintf(text, "%d", tam);
	strcat(nome, text);
	LIKWID_MARKER_INIT;
        LIKWID_MARKER_START (nome);
	double x[tam], y[tam];
        for(int i=0;i<tam;i++)
                y[i]=0;
        for(int i=0;i<tam;i++){
                for(int k=0;k<tam;k++)
                        y[k]=0;
                y[i]=1;
                retrossubsL(L, y, x, tam);
                retrossubsU(A, x, y, tam);
                for(int j=0;j<tam;j++){
                        AI[j][i]=y[j];
                }
        }
	LIKWID_MARKER_STOP (nome);
        LIKWID_MARKER_CLOSE;
        *tTotal = timestamp() - *tTotal;
}

void somaMatrizes(double **A, double **B, int tam){
    for(int i=0;i<tam;i++){
        for(int j=0;j<tam;j++){
            A[i][j] = A[i][i] + B[i][j];
        }
    }
}

void pegarColuna(double *y, double **R, int j, int tam){
	for(int i=0;i<tam;i++){
		if(i==j)
			y[i]=R[i][j];
		else
			y[i]=R[i][j];
	}
}

double normaL2Residuo(double **R, int tam){
        double soma=0;
        for(int i=0;i<tam;i++){
                for(int j=0;j<tam;j++){
                        soma = soma + (R[i][j] * R[i][j]);  
                }
        }
        return sqrt(soma);
}

void residuo(double **R, double **A, double **AI, double **I, int tam, double *tTotal){
	char nome[20]="teste_residuo_";
	char text[20];
	sprintf(text, "%d", tam);  
	strcat(nome, text);
        LIKWID_MARKER_INIT;
        LIKWID_MARKER_START(nome);
        *tTotal = timestamp();
        for (int i=0;i<tam;i++){
                for (int k=0;k<tam;k++){
                        R[i][k] = 0;
                        for (int j=0;j<tam;j++) {
                                R[i][k] += A[i][j] * AI[j][k];
                        }
                }
        }
        for(int i=0;i<tam;i++){
                for(int j=0;j<tam;j++){
                        R[i][j]=I[i][j]-R[i][j];
                }
        }
        *tTotal = timestamp() - *tTotal;

       	LIKWID_MARKER_STOP(nome);
        LIKWID_MARKER_CLOSE;
}

int refinamento (SL *SL, double *tTotal, int iteracoes, FILE* file){
	double **R, r;
	if (!(R = malloc(sizeof(void **) * SL->n))) {
                perror ("Erro de alocacao no refinamento: ");
                exit (1);
        }
	for (int i=0;i<SL->n;i++) {
                if (!(R[i] = malloc(sizeof(void *) * SL->n))) {
                        perror ("Erro de alocaco no refinamento: ");
                        exit (1);
                }
	}
	int i=1;
	int tam=SL->n;
	residuo(R, SL->AC, SL->AI, SL->I, SL->n, tTotal);
        double tempo = *tTotal;
	r = normaL2Residuo(R, tam);
	while(i <= iteracoes){
                if (file != NULL) // imprime iter na saida ou na tela
                        fprintf(file, "iter %d: %.15g\n", i, r);
                else
		      printf("iter %d: %.15g\n", i, r);
		residuo(R, SL->AC, SL->AI, SL->I, SL->n, tTotal);
                tempo += *tTotal;
		double x[tam], y[tam];
		int k=0;
		for(int j=0;j<tam;j++){
			pegarColuna(y, R, j, tam);
			retrossubsL(SL->L, y, x, tam);
                	retrossubsU(SL->A, x, y, tam);
			for(int c=0;c<tam;c++){	
			     SL->AI[c][k]=SL->AI[c][k] + y[c];
			}
			k++;
               	}
		i++;
		r = normaL2Residuo(R, tam);
	}
        if (file != NULL) // pula linha
                fprintf(file, "\n");
         else
                printf("\n");
        *tTotal = tempo;
	return 1;
}

void imprimeSL(double **A, int tam){
        for(int i=0;i<tam;i++){
                for(int j=0;j<tam;j++)
                        printf("%.15g ", A[i][j]);
                printf("\n");
        }
        printf("\n\n\n");
}

void imprimeSlSaida(double **A, int tam, FILE * output){
    for(int i=0;i<tam;i++){
        for(int j=0;j<tam;j++)
            fprintf(output, "%.15g ", A[i][j]);
                fputs ("\n", output); /* Pula linha */
    }
    fputs ("\n\n\n", output);
}
void geraIdentidade(SL *SL){
	for(int i=0;i<SL->n;i++){
        	for(int j=0;j<SL->n;j++)
                	if(i==j)
                        	SL->I[i][j]=1;
                        else
                                SL->I[i][j]=0;
	}
}
int ehInvertivel(SL *sl){
	double soma1=0, soma2=0;
	for(int i=0;i<sl->n;i++){
		for(int j=0;j<sl->n;j++){
			if(i==j)
				soma1 = soma1 + sl->A[i][j];
			if(i+j==sl->n-1)
				soma2 = soma2 + sl->A[i][j];
		}
	}
	if(soma2-soma1 == 0.0000000000000)
		return 0;
	return 1;
}
void liberaSisLin(SL *SL){
	if(SL){
  			for(int i=0; i < SL->n; ++i){
    				free(SL->A[i]);
				free(SL->L[i]);
				free(SL->AI[i]);
				free(SL->AC[i]);
				free(SL->I[i]);
			}
      			free(SL->A);
			free(SL->L);
			free(SL->AI);
			free(SL->AC);
			free(SL->I);
  	}
}
