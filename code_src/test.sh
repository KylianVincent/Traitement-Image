#/bin/sh

IMAGES_PATH=../imagestp

for f in $(find "$IMAGES_PATH" -type f -name "*.pgm" ); do
    ./filtreGaussien ../imagestp/formes1pets10.pgm test.pgm 0.5
    echo "PSNR pour $f :"
    ./testPsnr ../imagestp/formes1pets10.pgm test.pgm
    echo ""
    echo ""
done;
