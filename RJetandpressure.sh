#!/bin/bash
qcc -O2 -Wall -D_FORTIFY_SOURCE=0 Jetandpressure.c -o Jetandpressure -lm
./Jetandpressure


