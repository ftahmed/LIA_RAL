#!/bin/bash

wav=$1
if [ ! -e $wav ]; then
        echo "can't find the sphere file ($wav)"
fi

show=`basename $wav .wav`
mfcc=$show/$show.mfcc

mkdir $show
sfbcep -v -F WAVE -p 12 -m -e  -f 16000 $wav $mfcc
scopy -D $mfcc $mfcc.delta
nbf=`scopy $mfcc -o ascii - | wc -l`

echo "$show 1 0 $nbf U U U S1" > $show/$show.uem.seg
