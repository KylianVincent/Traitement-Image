#/bin/sh

IMAGES_PATH=../imagestp

for f in $(find "$IMAGES_PATH" -type f -name "formes1*.pgm" ); do
    ./filtreGaussien "$f" test.pgm 0.5 10
    echo "PSNR pour $f :"
    ./testPsnr ../imagestp/formes1.pgm test.pgm
    echo ""
    echo ""
done;

for f in $(find "$IMAGES_PATH" -type f -name "formes2*.pgm" ); do
    ./filtreGaussien "$f" test.pgm 0.5 10
    echo "PSNR pour $f :"
    ./testPsnr ../imagestp/formes2.pgm test.pgm
    echo ""
    echo ""
done;
