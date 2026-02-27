#!/bin/bash
for f in *.svg; do
    t=$(echo $f | sed 's/\.svg$/.png/')
    echo "Move $f to $t"
    # Convert stopped working after installing inkscape
    # convert $f -resize 48x48 $t
    inkscape -z -e $t -w 48 -h 48 $f
    rm $f
done
identify -verbose *.png | grep -E 'Image:|Geometry' > IM.txt
