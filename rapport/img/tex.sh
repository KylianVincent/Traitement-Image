#!/bin/bash

for s in `ls p2*` 
do
echo "\\begin{minipage}[c]{0.20\\linewidth}
	\\begin{center}
		\\includegraphics[width = 33mm]{./img/$s}
		\\textit{origine}
	\\end{center}
\\end{minipage}" 
done

