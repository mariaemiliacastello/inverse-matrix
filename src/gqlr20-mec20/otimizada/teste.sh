#!/bin/bash

METRICA="L3 L2CACHE FLOPS_DP"
TAMANHO="32 33 64 65 128 129 256"
LIKWID_HOME=/home/soft/likwid
CFLAGS="-I${LIKWID_HOME}/include -DLIKWID_PERFMON"
LFLAGS="-L${LIKWID_HOME}/lib -llikwid"

gcc ${CFLAGS} -c main.c -o main.o ${LFLAGS}
gcc ${CFLAGS} -c sislin.c utils.c ${LFLAGS}
gcc ${CFLAGS} -c linhaComando.c ${LFLAGS}
gcc ${CFLAGS} -c utils.c ${LFLAGS}
for k in $METRICA
do
    		likwid-perfctr -C 1 -g ${k} -m ./invmat -i 10 -e ../entradaSimples.txt >${k}_teste.log
done

