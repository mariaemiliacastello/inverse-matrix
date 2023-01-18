/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#include "utils.h"
#include "sislin.h"
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <likwid.h>

// Aloca espaço na memória para a struct SL
SL* criaNovoSl (int n) {
        SL *sl = (SL *) malloc(sizeof(SL));
 
        if ( sl ) {
                sl->n = n;
                sl->A = malloc(n * sizeof(double *));
                sl->AI = malloc(n * sizeof(double *));
                sl->AC = malloc(n * sizeof(double *));
                sl->I = malloc(n * sizeof(double *));
                sl->L = malloc(n * sizeof(double *));

                if (!(sl->A) || !(sl->AI) || !(sl->AC) || !(sl->I) || !(sl->L)) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao 1");
                        exit (1);
                }

                // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
                sl->A[0] = (double *) aligned_alloc(ALOCACAO, n * n * sizeof(double));
                if (!(sl->A[0])) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao A");
                        exit (1);
                }
                for (int i=1; i < n; ++i)
                        sl->A[i] = sl->A[i-1]+n;

                // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
                sl->AI[0] = (double *) aligned_alloc(ALOCACAO, n * n * sizeof(double));
                if (!(sl->AI[0])) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao AI");
                        exit (1);
                }
                for (int i=1; i < n; ++i)
                        sl->AI[i] = sl->AI[i-1]+n;

                // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
                sl->AC[0] = (double *) aligned_alloc(ALOCACAO, n * n * sizeof(double));
                if (!(sl->AC[0])) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao AC");
                        exit (1);
                }
                for (int i=1; i < n; ++i)
                        sl->AC[i] = sl->AC[i-1]+n;

                // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
                sl->I[0] = (double *) aligned_alloc(ALOCACAO, n * n * sizeof(double));
                if (!(sl->I[0])) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao I");
                        exit (1);
                }
                for (int i=1; i < n; ++i)
                        sl->I[i] = sl->I[i-1]+n;

                // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
                sl->L[0] = (double *) aligned_alloc(ALOCACAO, n * n * sizeof(double));
                if (!(sl->L[0])) {
                        liberaSisLin(sl);
                        perror ("Erro alocacao L");
                        exit (1);
                }
                for (int i=1; i < n; ++i)
                        sl->L[i] = sl->L[i-1]+n;
        }
        return (sl);
}

// Liberacao de memória
void liberaSisLin (SL *sl) {
        if (sl) {
                if (sl->A) {
                        if (sl->A[0])
                                free (sl->A[0]);
                        free(sl->A);
                }

                if (sl->AI) {
                        if (sl->AI[0])
                                free (sl->AI[0]);
                        free(sl->AI);
                }

                if (sl->AC) {
                        if (sl->AC[0])
                                free (sl->AC[0]);
                        free(sl->AC);
                }

                if (sl->I) {
                        if (sl->I[0])
                                free (sl->I[0]);
                        free(sl->I);
                }

                if (sl->L) {
                        if (sl->L[0])
                                free (sl->L[0]);
                        free(sl->L);
                }

                free(sl);
        }
}

void geraSlAleatoria (SL *SL, double coef_max) { // Gera uma SL com valores aleatórios
        unsigned int n = SL->n;
        // para gerar valores no intervalo [0,coef_max]
        double invRandMax = ((double)coef_max / (double)RAND_MAX);
        // inicializa a matriz A
        for (unsigned int i=0; i<n; ++i) {
                for (unsigned int j=0; j<n; ++j)  {
                        SL->A[i][j] = (double)rand() * invRandMax;
                        SL->AC[i][j] = SL->A[i][j];
                }
        }
}

void iniciarSlEntrada(FILE* input, SL *sl) {
        int i, j;

        if (input != NULL) {
                for(i=0;i<sl->n;i++) {
                        for(j=0;j<sl->n;j++)
                                fscanf(input,"%le", &sl->A[i][j]);
                }
                for(i=0;i<sl->n;i++) {
                        for(j=0;j<sl->n;j++)
                                sl->AC[i][j]=sl->A[i][j];
                }
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

int eliminacaoGauss (SL *sl, double *tTotal){
        *tTotal = timestamp();
        int tam = sl->n;
        memset(sl->L[0], 0, tam * sizeof(double));
        for(int i = 0; i < tam - 1; i++) {
                int iPivo = encontraMax(sl->A, i, tam);
                if(i != iPivo) {
                        trocaLinha(sl->A, i, iPivo, tam);
                        trocaLinha(sl->AC, i, iPivo, tam);
                        trocaLinha(sl->I, i, iPivo, tam);
                        trocaLinha(sl->L, i, iPivo, tam);
                }

                for(int k = i + 1; k < tam; k++){
                        double m = sl->A[k][i] / sl->A[i][i];
                        sl->L[k][i] = m;
                        if (isnan (m) || isinf (m)) {
                                perror ("Erro, divisao por 0 na ElinacaoGauss: ");
                                exit (1);
                        }
                        sl->A[k][i] = 0.0;
                        for(int j = i + 1; j < tam; ++j)
                                sl->A[k][j] -= sl->A[i][j] * m;
                }
        }
        for(int i = 0; i < tam; i++){
                sl->L[i][i]=1;
        }
        *tTotal = timestamp() - *tTotal;
        return 0;
}

void retrossubsU (double **A, double *b, double *x, int tam){
for(int i=tam-1;i>=0;i--){
                x[i]=b[i];
                for(int j=i+1;j<tam;j++){
                        x[i] = x[i] - (A[i][j] * x[j]);
                }
                x[i] = x[i]/A[i][i];
                if (isnan (x[i]) || isinf (x[i])) {
                        perror ("Erro, divisao por 0 na retrossubs: ");
                        exit (1);
                }
        }
}

void retrossubsL (double **A, double *b, double *x, int tam){
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

void fatoracaoLU(double **restrict I, double **restrict AI, double **restrict L, double **restrict A, int tam, double *tTotal){
        *tTotal=timestamp();
	char nome[20]="teste_fatoracaoLU_";
	char text[20];
	sprintf(text, "%d", tam);  
	strcat(nome, text);
	LIKWID_MARKER_INIT;
        LIKWID_MARKER_START (nome);
        double x[tam], y[tam];
        for(int i=0;i<tam;i++){
                memset (&y,0,sizeof(y)); 
                y[i]=1;
                retrossubsL(L, y, x, tam);
                retrossubsU(A, x, y, tam);
                for (int j = 0; j < tam; j++)
                        AI[j][i]=y[j];
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

void pegarcoluna(double *restrict y, double **restrict R, int j, int tam) {
        for(int i=0;i<tam;i++)
                y[i]=R[i][j];
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

void residuo(double **restrict R, double **restrict A, double **restrict AI, double **restrict I, int tam, double *tTotal){
        *tTotal = timestamp();
	char nome[20]="teste_residuo_";
	char text[20];
	sprintf(text, "%d", tam);  
	strcat(nome, text);
		
        LIKWID_MARKER_INIT;
        LIKWID_MARKER_START (nome);

        for (int i = 0; i < tam; i++) {
                for (int k = 0; k < tam - (tam % 4); k+= 4) {
                        R[i][k] = R[i][k+1] = R[i][k+2] = R[i][k+3] = 0;
                        for (int j = 0; j < tam; j++) {
                                R[i][k] += A[i][j] * AI[j][k];
                                R[i][k+1] += A[i][j] * AI[j][k+1];
                                R[i][k+2] += A[i][j] * AI[j][k+2];
                                R[i][k+3] += A[i][j] * AI[j][k+3];
                        }
                }
                for (int k = tam - (tam % 4); k < tam; k++) {
                        R[i][k] = 0;
                        for (int j = 0; j < tam; j++)
                                R[i][k] += A[i][j] * AI[j][k];
                }
        }

        for(int i=0;i<tam;i++) {
                for(int j = 0; j < tam - (tam % 4); j += 4) { // loop unroll
                        R[i][j] = I[i][j] - R[i][j];
                        R[i][j+2] = I[i][j+2] - R[i][j+2];
                        R[i][j+1] = I[i][j+1] - R[i][j+1];
                        R[i][j+3] = I[i][j+3] - R[i][j+3];
                }        
        }
        for(int i=0;i<tam;i++) { // Resíduo do laço
                for(int j = tam - (tam % 4); j < tam; j++) {
                        R[i][j] = I[i][j]-R[i][j];
                }
        }
	*tTotal = timestamp() - *tTotal;
        LIKWID_MARKER_STOP (nome);
        LIKWID_MARKER_CLOSE;
}

int refinamento (SL *SL, double *tTotal, int iteracoes, FILE* file){ // Vai calculando o resíduo
        double **R, r;
        // Alocação
        R = malloc (SL->n * sizeof(double *));
        if (!R) {
                perror ("Erro alocacao Residuo");
                exit (1);
        }
        // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
        R[0] = (double *) malloc(SL->n * SL->n * sizeof(double));
        if (!(R[0])) {
                perror ("Erro alocacao Residuo");
                exit (1);
        }
        for (int i=1; i < SL->n; ++i)
                R[i] = R[i-1] + SL->n;

        int i=1;
        int tam = SL->n;
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
                        pegarcoluna(y, R, j, tam);
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

        free (R[0]);
        free(R);

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

void geraIdentidade(double **I, int tam) {

        memset (I[0], 0, tam * tam * sizeof (double));

        for(int i = 0; i < (tam - (tam % 4)); i+= 4) {
                I[i][i] = 1;
                I[i+1][i+1] = 1;
                I[i+2][i+2] = 1;
                I[i+3][i+3] = 1;
        }
        for (int i = (tam - (tam % 4)); i < tam; i++) // Resíduo do laço
                I[i][i] = 1;
}

int ehInvertivel(SL *sl){ // Se determinante for diferente de 0 é invertível, já se não for não é invertível
        double soma1=1, soma2=1;
        int i, j;

        for (i = 0; i < sl->n; i++)
                soma1 = soma1 * sl->A[i][i];

        for (i = 0; i < sl->n; i++)
                for (j = 0; j < i; j++)
                       if(i+j==sl->n-1)
                                soma2 = soma2 * sl->A[i][j];
                               
        for (i = 0; i < sl->n; i++)
                for(j = i + 1; j < sl->n; j++)
                        if(i+j==sl->n-1)
                                soma2 = soma2 * sl->A[i][j];

        if(fabs(soma2-soma1) <= FLT_EPSILON)
                return 0;
        return 1;
}
