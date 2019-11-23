#!/bin/bash

#gcc practica3_1b.c -o practica3_1b -lm

#L=16
TRANSIT=1000000
DECORR=1000000
Ti="1.80"
Tf="3.00"
DELTA_T="0.05"

for L in 16 32 64 128; do

    OUTPUT="L${L}_Ti${Ti}_Tf${Tf}_deltaT${DELTA_T}_transit${TRANSIT}_decorr${DECORR}.dat"

    echo $L
    echo $TRANSIT
    echo $DECORR
    echo $Ti
    echo $Tf
    echo $DELTA_T
    echo $OUTPUT

    ./practica3_1b ${OUTPUT} ${L} ${TRANSIT} ${DECORR} ${Ti} ${Tf} ${DELTA_T}

done