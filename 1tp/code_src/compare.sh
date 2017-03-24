#!/bin/bash

sigma="0.5 1 5 10"
Wlist=`seq 1 150`
# Spatial
for s in $sigma
do
    for w in $Wlist
    do
        ./test.sh "$1" $s $w
    done
done

## Frequentiel
for s in $sigma
do
    ./test.sh "$1" $s
done