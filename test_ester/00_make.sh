#!/bin/bash

wav=$1

./01_make_mfcc.sh $wav
./02_make_segCLR.sh $wav

