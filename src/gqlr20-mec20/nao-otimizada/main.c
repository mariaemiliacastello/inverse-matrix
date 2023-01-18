/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#include "sislin.h"
#include "linhaComando.h"
#include "utils.h"
#include <stdio.h>

int main(int argc, char *argv[]){
	srand(202201);

	real_t LU_tempo, iter_tempo, residuo_tempo;

	int cont = 0; // Conta quantas matrizes já foram calculadas

	FILE *input = NULL;
	FILE *output = NULL;

	SL sl;

	int tam = le_tamanho (argc, argv);

	int iteracoes = le_iteracoes (argc, argv);

	if (tam == -1) { // Se nao especificou tamanho da SL na linha de comando tenta ler da entrada
		input = abre_input(argc, argv,"r");
	}
	output = abre_output(argc, argv,"w"); /* tenta achar saida na linha de comando */

	if (tam == -1 && input == NULL) // Se não indicou nem o tam e nem arquivo de entrada cria uma SL com um tam padrão
		tam = TAM_PADRAO;

  	if (input != NULL) { // Se for ler SL do arquivo de entrada
  		while (!feof(input) && fscanf(input,"%d\n", &tam) == 1) {
        	criaNovoSl(&sl, tam); // Aloca SL
  			iniciarSlEntrada(input, &sl);

  			cont++;

			if(ehInvertivel(&sl) == 0){
                        	if (output != NULL) // Imprime na saida ou na tela
                			fprintf(output, "matriz %d nao invertivel\n\n", cont);
               			else 
                        		printf("matriz %d nao invertivel\n\n", cont);      
			}
			else{
				if (output != NULL) { // Imprime na saida ou na tela
  					fprintf(output, "matriz %d\n\ntam: %d\niteracoes: %d\n\n", cont, tam, iteracoes);
  				}
  				else {
  					printf("matriz %d\n\ntam : %d\niteracoes: %d\n\n", cont, tam, iteracoes);
  				}

				geraIdentidade(&sl); // Gera matriz identidade
				eliminacaoGauss(sl.A, sl.L, sl.AC, sl.I, sl.n, &LU_tempo);		
				fatoracaoLU(sl.AI, sl.L, sl.A, sl.n, &iter_tempo);

				//imprimeSL (sl.AI, sl.n);
				//printf("\n\n");

				refinamento(&sl, &residuo_tempo, iteracoes,output);
			
				if (output != NULL) { // Imprime tempos na saida ou na tela
  					fprintf(output, "LU_tempo: %.15g\niter_tempo: %.15g\nresiduo_tempo: %.15g\n\n", LU_tempo, iter_tempo, residuo_tempo);
				}
  				else {
  					printf("LU_tempo: %.15g\niter_tempo: %.15g\nresiduo_tempo: %.15g\n\n", LU_tempo, iter_tempo, residuo_tempo);
  				}
				
				/*if (output != NULL) // Imprime matriz inversa na saida ou na tela
                                	imprimeSlSaida (sl.AI, tam, output);
                        	else
                                	imprimeSL(sl.AI, tam);*/
			}
			liberaSisLin(&sl);
  		}
  	}
  	else { // Cria SL aleatória de tamanho fixo
  		criaNovoSl (&sl, tam); // Aloca SL
  		geraSlAleatoria (&sl, generico, COEF_MAX);

  		cont++;
		
  		if(ehInvertivel(&sl) == 0){
                        if (output != NULL) // Imprime na saida ou na tela
                		fprintf(output, "matriz %d nao invertivel\n\n", cont);
               		else 
                        	printf("matriz %d nao invertivel\n\n", cont);      
		}
		else{
			if (output != NULL) { // Imprime na saida ou na tela
                		fprintf(output, "matriz %d\n\ntam: %d\niteracoes: %d\n\n", cont, tam, iteracoes);
                	}
               		else {
                        	printf("matriz %d\n\ntam : %d\niteracoes: %d\n\n", cont, tam, iteracoes);
                	}

			geraIdentidade(&sl);
			eliminacaoGauss(sl.A, sl.L, sl.AC, sl.I, sl.n, &LU_tempo);
			fatoracaoLU(sl.AI, sl.L, sl.A, sl.n, &iter_tempo);

			refinamento(&sl, &residuo_tempo, iteracoes,output);

			if (output != NULL) { // Imprime tempos na saida ou na tela
  				fprintf(output, "LU_tempo: %.15g\niter_tempo: %.15g\nresiduo_tempo: %.15g\n\n", LU_tempo, iter_tempo, residuo_tempo);
  			}
  			else {
  				printf("LU_tempo: %.15g\niter_tempo: %.15g\nresiduo_tempo: %.15g\n\n", LU_tempo, iter_tempo, residuo_tempo);
  			}

/*			if (output != NULL) // Imprime matriz inversa na saida ou na tela
                        	imprimeSlSaida (sl.AI, tam, output);
                	else
                        	imprimeSL(sl.AI, tam);*/
  		}
		liberaSisLin(&sl);
	}
  	// Caso entrada ou saida tenham sido abertos, fecha
  	if (input != NULL)
  		fclose (input);
	if (output != NULL)
		fclose (output);
}
