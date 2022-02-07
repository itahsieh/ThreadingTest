#!/bin/bash
# Prerequisites: gcc

SRC=ThreadingTest.c
EXE=ThreadingTest

clear
rm -f $EXE

gcc $SRC -o $EXE \
-Wall -lm -O3 \
-fopenmp -lpthread 
./$EXE
