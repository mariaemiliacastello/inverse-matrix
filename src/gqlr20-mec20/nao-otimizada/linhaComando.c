/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#include "linhaComando.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define ITERACAO_PADRAO 50

int le_tamanho (int argc, char **argv) {
	for (int i = 0; i < argc - 1; i++)
		if (strcmp (argv[i], "-r") == 0)
			return atoi(argv[i+1]);
	return -1; // Se não achar o tamanho da matriz retorna -1
}

int le_iteracoes (int argc, char **argv) {
	for (int i = 0; i < argc - 1; i++)
		if (strcmp (argv[i], "-i") == 0)
			return atoi(argv[i+1]);
	return ITERACAO_PADRAO; // Se não achar o num de iterações na linha de comando retorna um num padrão 
}

FILE* abre_input (int argc, char **argv, char* modo) {
	for (int i = 0; i < argc - 1; i++)
		if (strcmp (argv[i], "-e") == 0) {
			FILE *file;
			file = fopen(argv[i+1], modo);
			if (!file) { /* Se não conseguir abrir arquivo emite mensagem de erro e fecha o programa */
				perror ("Erro ao abrir arquivo de entrada: ");
				exit (1);
			}
			return file;
		}
	return NULL;
}

FILE* abre_output (int argc, char **argv, char* modo) {
	for (int i = 0; i < argc - 1; i++)
		if (strcmp (argv[i], "-s") == 0) {
			FILE *file;
			file = fopen(argv[i+1], modo);
			if (!file) { /* Se não conseguir abrir arquivo emite mensagem de erro e fecha o programa */
				perror ("Erro ao abrir arquivo de saida: ");
				exit (1);
			}
			return file;
		}
	return NULL;
}
