#!/bin/bash

for s in `ls p2*` 
do
	convert $s "`basename $s '.pgm'`.jpg" 
done

