#include <iostream>
#include "liatools.h"
#include "ScorpioGenerator.h"

int main(int argc, char* argv[]){

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    //cc.addIntegerParam("verboseLevel",false,true,"level of the berose information 0=no verbose, 1=normal, 2=more");
    //cc.addStringParam("inputFeatureFilename",false, true,"feature filename or filename of a text file with the list of feature filenames");
    //cc.addStringParam("targetIdList",true,true,"The file with the list of models to train. A line is composed by client_id file1 file2 ...");
    //cc.addStringParam("inputWorldFilename",false,true,"if set, the init is based on a model get from this file, else frrom scratch");
    //cc.addStringParam("mixtureServer",false,true,"If set save the complete mixture server in the filename (FUTURE USED, TODO)");
    //cc.addBooleanParam("initByClient",false,true,"For by lael option. Modify the initial model for statistic estimation (EM), default world, if set client");
    //cc.addBooleanParam("saveEmptyModel",false,true,"If no data is available for a model (or a lable model), save the not adapted model (world or global client)");
    //cc.addBooleanParam("useIdForSelectedFrame",false,true,"If set, the segments with label ID are used for training the client model ID");
    //cc.addStringParam("labelSelectedFrames",false,true,"The segments with this label are used for training the worldmodel (if UseIdForSelectedFrame is not used)"); 
    //cc.addFloatParam("baggedFrameProbability",false,true,"Defines the % of frames taken for each iterations (default 1)");
    //cc.addIntegerParam("nbTrainIt",false,true,"number of it (default=1)"); 
    //cc.addBooleanParam("normalizeModel",false,true,"if set to true,  normalize the world (at each iteration)");
    //cc.addBooleanParam("normalizeModelMeanOnly",false,true,"Used only if normalizeModel is On, says if only mean parameters should be normalized"); 
    //cc.addIntegerParam("normalizeModelNbIt",false,true,"Used only if noramlizeModelMeanOnly is set, nb of normalization it");
    //cc.addBooleanParam("meanAdapt",false,true,"Mean adaptation (default false)");
    //cc.addBooleanParam("varAdapt",false,true,"Variance adaptation (default false)");
    //cc.addBooleanParam("weightAdapt",false,true,"Weight adaptation (default false)");
    //cc.addStringParam("MAPAlgo",true,true,"Adaptation method (MAPConst,MAPConst2,MAPOccDep,MLLR)");
    //cc.addFloatParam("MAPAlphaMean",false,true,"a priori proba for world");	    
    //cc.addFloatParam("MAPAlphaVar",false,true,"a priori proba for world");
    //cc.addFloatParam("MAPAlphaWeight",false,true,"a priori proba for world");
    //cc.addFloatParam("MAPRegFactorMean",false,true,"Reg factor");			
    //cc.addFloatParam("MAPRegFactorVar",false,true,"Reg factor");	
    //cc.addFloatParam("MAPRegFactorWeight",false,true,"Reg factor");	
    //cc.addBooleanParam("info",false,false,"If info is requested, just info on the train set is outputed");
    //cc.addBooleanParam("useModelData",false,true,"New MAP algo based on ML estimate of the training data");
    //cc.addStringParam("initModel",false,true,"With model based, use a specific model for initialize the EM estimate (default=world");
    //cc.addBooleanParam("outputAdaptParam",false,true,"Saving a vector (matrix if MLLR, weights if MAP) instead of a mixture");
 try{
    CmdLine cmdLine(argc, argv);
    if (cmdLine.displayHelpRequired()){
      cout <<"ScorpioGenerator.exe"<<endl<<"This program is used to generate Header file for Scorpio"
	   <<endl<<cc.getParamList()<<endl;
      return 0;  
    }
    if (cmdLine.displayVersionRequired()){
    cout <<"Version 1.0"<<endl;
    } 
    Config tmp;
    cmdLine.copyIntoConfig(tmp);
    Config config;
    if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
    cmdLine.copyIntoConfig(config);
    cc.check(config);
    debug=config.getParam_debug();
    
	ScorpioGenerator(config);  
 
    
    return 0;
  }
  catch(alize::Exception & e){ cout <<"ScorpioGenerator "<< e.toString() << endl << cc.getParamList()<<endl;}
}
