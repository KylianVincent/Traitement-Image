#!/bin/bash

for s in `ls *.pgm` 
do
	convert $s "`basename $s '.pgm'`.jpg" 
done

