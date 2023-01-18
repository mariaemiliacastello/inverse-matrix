/* Guilherme de Queiroz Lima Roth - GRR20206149
Maria Emilia Castello - GRR20203921
Introdução à computação Científica - 2o semestre de 2022
*/

#ifndef __linhaComando_H__
#define __linhaComando_H__

#include <stdio.h>

FILE* abre_input (int argc, char **argv, char* modo);
FILE* abre_output (int argc, char **argv, char* modo);
int le_tamanho (int argc, char **argv);
int le_iteracoes (int argc, char **argv);

#endif
