#!/bin/bash

PATH=$PATH:.

wav=$1
if [ ! -e $wav ]; then
        echo "can't find the sphere file ($wav)"
fi

show=`basename $wav .wav`
datadir=${show}
mkdir ./$datadir >& /dev/null

featureExt=.mfcc
fileFormat=SPRO4
mask=0-12
clrMask=0-11,13-24
size=13
segPath=${datadir}/
inSegType=.seg
segType=.seg
outSegType=.seg
saveSegType=XML
extSegInitOut=.i
extSegInitIn=.uem
extSegOut=.s
extSegIn=.i
extSegClustOut=.l
extSegClustIn=.s
extSegClustHIn=.l
extSegClustHOut=.h
extSegTrainInitIn=.h
extSegTrainEMIn=.h
extSegDecodeOut=.d
extSegClustCLRIn=.d
extSegClustCLROut=.c
mergeMax=2147483647

uem=sph/$show.uem.seg

echo "#####################################################"
echo "#   $show"
echo "#####################################################"

mkdir ./$datadir >& /dev/null

#first version with vectsize
../src/segInit --debug=true --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --loadSegServerFileFormat=$segType --extSegClustIn=$extSegInitIn --saveSegServerFileExtension=$outSegType --saveSegServerFileFormat=$saveSegType --extSegClustOut=$extSegInitOut --show=$show

#GLR based segmentation, make small segments
#Alize version
../src/seg --debug=true --minLLK=-200 --maxLLK=200 --segThr=-1000 --segMinWSize=250 --segWSize=250 --distribType=GF --segMethod=GLR --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --loadSegServerFileFormat=$segType --extSegClustIn=$extSegIn --saveSegServerFileExtension=$outSegType --saveSegServerFileFormat=$saveSegType --extSegClustOut=$extSegOut --show=$show

# linear clustering
../src/clust --minLLK=-200 --maxLLK=200 --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --loadSegServerFileFormat=$segType --extSegClustIn=$extSegClustIn --saveSegServerFileExtension=$outSegType --saveSegServerFileFormat=$saveSegType --extSegClustOut=$extSegClustOut --method=l --thr=2 --distribType=GF --mergeMax=$mergeMax --show=$show

# hierarchical clustering
../src/clust --minLLK=-200 --maxLLK=200 --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --loadSegServerFileFormat=$segType --extSegClustIn=$extSegClustHIn --saveSegServerFileExtension=$outSegType --saveSegServerFileFormat=$saveSegType --extSegClustOut=$extSegClustHOut --method=h --thr=4.5 --distribType=GF --mergeMax=$mergeMax --show=$show

# initialize GMM
#LIA Version
../src/trainInit --minLLK=-200 --maxLLK=200 --mixtureDistribCount=8 --distribType=GD --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --extSegClustIn=$extSegTrainInitIn --loadSegServerFileFormat=$segType --mixtureFilesPath=$datadir/ --saveMixtureFileExtension=.init.xml --dMethod=std --show=$show

# EM computation
#LIA Version
../src/trainEM --mixtureDistribCount=8 --minLLK=-200 --maxLLK=200 --distribType=GD --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --extSegClustIn=$extSegTrainEMIn --loadSegServerFileFormat=$segType --mixtureFilesPath=$datadir/ --saveMixtureFileExtension=.xml --loadMixtureFileExtension=.init.xml --gain=0.01 --maxIt=10 --minIt=3  --show=$show

#Viterbi decoding
../src/decode  --mixtureDistribCount=8 --minLLK=-200 --maxLLK=200 --distribType=GD --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt --loadFeatureFileFormat=$fileFormat --featureServerMask=$mask --loadFeatureFileVectSize=$size --vectSize=$size --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --extSegClustIn=$extSegClustHOut --loadSegServerFileFormat=$segType --mixtureFilesPath=$datadir/ --loadMixtureFileExtension=.xml --extSegDecodeOut=$extSegDecodeOut --fudge=1 --viterbiBufferLength=30 --HMMPenalityIJ=-250 --show=$show

#CLR clustering
# Features contain static and delta and are centered and reduced (--fdesc)
../src/clust --computeLLKWithTopDistribs=COMPLETE --meanAdapt=true --MAPAlgo=MAPModelBased --MAPRegFactorMean=15 --initVarianceFlooring=0 --initVarianceCeiling=10 --nbTrainIt=5 --normalizeModel=true --mixtureDistribCount=512 --topDistribsCount=5 --minLLK=-200 --maxLLK=200 --featureFilesPath=$datadir/ --loadFeatureFileExtension=$featureExt.delta --loadFeatureFileFormat=$fileFormat --featureServerMask=$clrMask --loadFeatureFileVectSize=24 --vectSize=24 --segServerFilesPath=$segPath --loadSegServerFileExtension=$inSegType --loadSegServerFileFormat=$segType --extSegClustIn=$extSegClustCLRIn --saveSegServerFileExtension=$outSegType --saveSegServerFileFormat=$saveSegType --extSegClustOut=$extSegClustCLROut --loadMixtureFileExtension=.xml --mixtureFilesPath=./gmm_ubm/ --world=ubm --method=c --thr=1.75 --distribType=GD --mergeMax=2147483647 --minSpk=0 --featureServerMode=FEATURE_WRITABLE --featureServerMemAlloc=50000000 --show=$show

