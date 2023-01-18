#!/bin/bash

METRICA="L3 L2CACHE FLOPS_DP"
TAMANHO="6000"
LIKWID_HOME=/home/soft/likwid
CFLAGS="-I${LIKWID_HOME}/include -DLIKWID_PERFMON"
LFLAGS="-L${LIKWID_HOME}/lib -llikwid"

cd ./otimizada
gcc ${CFLAGS} -c main.c -o main.o ${LFLAGS}
gcc ${CFLAGS} -c sislin.c utils.c ${LFLAGS}
gcc ${CFLAGS} -c linhaComando.c ${LFLAGS}
gcc ${CFLAGS} -c utils.c ${LFLAGS}
for r in $TAMANHO
do
     for k in $METRICA
     do
         echo ${r}
                echo ${k}
                likwid-perfctr -C 1 -g ${k} -m ./invmat -i 10 -r ${r} >>o${k}_otimizada.log
     done
done
cd ../nao-otimizada
gcc ${CFLAGS} -c main.c -o main.o ${LFLAGS}
gcc ${CFLAGS} -c sislin.c utils.c ${LFLAGS}
gcc ${CFLAGS} -c linhaComando.c ${LFLAGS}
gcc ${CFLAGS} -c utils.c ${LFLAGS}
for r in $TAMANHO
do
     for k in $METRICA
     do
         echo ${r}
                echo ${k}
                likwid-perfctr -C 1 -g ${k} -m ./invmat -i 10 -r ${r} >>${k}_naootimizada.log
     done
done
