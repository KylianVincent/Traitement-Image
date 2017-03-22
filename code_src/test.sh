#!/bin/bash

IMAGES_PATH=../imagestp

./filtreGaussien "$IMAGES_PATH/formes$1.pgm" "../rapport/img/$1.pgm" $2 $3
##eog "../rapport/img/$1.pgm" &
##echo "PSNR pour $IMAGES_PATH/formes$1.pgm  (sigma=$2, W=$3):"
./testPsnr "../imagestp/formes2.pgm" "../rapport/img/$1.pgm"
echo ""
