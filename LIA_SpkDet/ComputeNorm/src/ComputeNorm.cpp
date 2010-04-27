#if !defined(ALIZE_ComputeNorm_cpp)
#define ALIZE_ComputeNorm_cpp

#include <iostream>
#include <fstream>  
#include <cstdio> 
#include <cassert> 
#include <cmath>
#include <cstdlib>
#include <liatools.h>
#include "ComputeNorm.h"


// JF Bonastre 4/11/2009 - cleaning and add dichotomic search in Norm
// JF Bonastre 27/4/2010 -cleaning and add 
//                          * impostor (target independent) selection
//                          * High score discard
//                          * Median computation instead of aryth mean
                          

using namespace alize;
using namespace std;


//TODO move mean() and std() to getMean() and getStd()
class DistribNorm{
	LKVector tabScore; // table of norm scores
	/* Use of LKVector ALIZE object, it provides sort functions */
	String *tabIdImp; // table of norm scores
//	double meanVal;   // storage of mean value to avoid recomputation
//	double stdVal;    // storage of std value to avoid recomputation
	Config conf;
//	static int compare(const void*, const void*);
	public:
	void init();
	void addScore(double, const String&);
	/*void selectScore(XLine *, String);*/
	//long findImp(String);
	bool computeMeanStd(double &, double&, char,double);
//	double std();
//	double sum();
//	double sum2(); 
	unsigned long size();
	void print();
	DistribNorm(Config &);  
	~DistribNorm();
};

void DistribNorm::init(){
 	tabScore.clear();
//	meanVal = 0.0;
//	stdVal = 1.0;	 
}

unsigned long DistribNorm::size(){
	return tabScore.size();
}

void DistribNorm::addScore(double d, const String& idImp){
	LKVector::type t;
 	unsigned long ind;
	/* idx value allows to keep in memory the impostor name related to each score, even after applying the sort function */
        /* After sorting, the name of an impostor is retrieved thanks to the tabIdImp[idx] */
        
	if(tabScore.size() == 0)
		ind = (unsigned long)0;
	else
		ind = (tabScore.size());
	t.lk = d;
	t.idx = ind;
	tabIdImp[ind] = idImp;
	tabScore.addValue(t);
}

/*
void DistribNorm::selectScore(XLine *l, String selectType){
	
	if(selectType == "noSelect"){		
	 	return;	 
	}	 
	
	if(selectType == "selectNBestByTestSegment"){
		// Sort of impostor scores	
		if(debug){
			cout << "Sort the score" << endl;		 
		}
	 	sortScore();		 
		// Size of impostor cohort
		String sMax = conf.getParam("cohortNb");
		unsigned long maxCohortNb = (unsigned long)sMax.toLong();	  	
	 	
		if(debug){
			cout << "Size of the cohort " << maxCohortNb <<endl;		 
		}
		tabScore.setSize(maxCohortNb);
		
		if(debug){
			cout << "Score distribution after selection" << endl;
		 	print();
		}
		return;
	}
	
       if(selectType == "selectTargetDependentCohortInFile"){
		// Scores are selected according to an impostor name list, which is dependent on target model 
	
	 	if(debug){
			cout << "Retrieve client dependent cohort score" << endl;		 
		}
		String impFile = conf.getParam("cohortFilePath")+l->getElement(1)+conf.getParam("cohortFilesExt");
		cout << "ouverture file : " << impFile << endl;
		XList impFileList(impFile, conf);
	 	XLine * linep;
		unsigned long impNb=0;
		while ((linep=impFileList.getLine()) != NULL){
			String impName = linep->getElement(0);
			// necessary to find where is stored the informationof this impostor 
			//
			long indImp = findImp(impName);			
			
		        if(indImp == -1){
				cout << "dependent Impostor name does not correspond to any impostor" << endl; 
				exit(0);
			 }
			else{
			 	//permutation of structures in tabScore 
				LKVector::type t = tabScore[impNb];
				tabScore[impNb] = tabScore[indImp];
				tabScore[indImp] = t;
				impNb++;				
			}
		}
	    tabScore.setSize(impNb);
		if(debug){
			cout << "Score distribution after Cohort selection" << endl;
		 	print();
		}
	       	
	}	
}

long DistribNorm::findImp(String impName){
	unsigned int i=0;
	while((i < tabScore.size()) && (tabIdImp[tabScore[i].idx] != impName)) i++;
	if(i == tabScore.size()) return -1;
	return i;	
}
*/
/*
double DistribNorm::sum(){
 	double sum=0.0;
	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}
	for(unsigned int i=0; i<tabScore.size(); i++)
		sum += tabScore[i].lk;   	
    return sum;
}

double DistribNorm::sum2(){
 	double sum2=0.0;	
	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}
	for(unsigned int i=0; i<tabScore.size(); i++)
		sum2 += tabScore[i].lk * tabScore[i].lk;   	
	return sum2;
}
*/
bool DistribNorm::computeMeanStd(double &mean,double &std,char mode, double percent){
	if(tabScore.size() == 0) return false;
	else{
		unsigned long size=tabScore.size();
		unsigned long begin=0;
		if (percent<1.0){
			tabScore.descendingSort();
			size=((double) size) * percent;
			begin=tabScore.size()-size;
		   }
		double sum=0;
		double sum2=0;
		switch (mode){
			case 0: //classical mean computation	
					for(unsigned int i=begin; i<tabScore.size(); i++){
						sum += tabScore[i].lk;  
						sum2 += tabScore[i].lk * tabScore[i].lk;
						}
					mean=sum/(double)(size);
					std=sqrt((sum2/(double)(size)-(mean*mean)));  
					break;
			case 1: // Median instead of mean
			 		mean=tabScore[begin+(size/2)].lk;
			 		for(unsigned int i=begin; i<tabScore.size(); i++)
						sum += abs(tabScore[i].lk-mean); 
			 		std=sum/(double)(size);
			 		break;
			default: cout <<"mean comptation mode["<<mode<<"] unknown"<<endl;
					 exit(0);
		}  				
	
	    return true;
	}
}

/*
double DistribNorm::std(){
	double m;	

	if(tabScore.size() == 0){
		cout << "No normalization score available for std computation: std=1.0" << endl;
		return 1.0;	 
	}
	m = mean();
	return 
}*/


void DistribNorm::print(){
	cout << "Distrib Score Print (nbScore=" << tabScore.size() << ")" << endl;
 	for(unsigned int i=0; i<tabScore.size(); i++){
		cout << tabIdImp[tabScore[i].idx] << " ";
	 	printf("%.8lE ", tabScore[i].lk);
	}
	cout << endl;
}

DistribNorm::DistribNorm(Config &config){
	conf = config;
	unsigned long nbScoreMax = (unsigned long)(conf.getParam("maxScoreDistribNb").toLong());
	tabIdImp = new String[nbScoreMax];	
//	meanVal = 0.0;
//	stdVal = 1.0;
}

DistribNorm::~DistribNorm() {
 	delete []tabIdImp;
}

class NormSeg{
	String name;
 	double mean;
	double std;
	bool computed;
//	double sum;
//	double sum2;
	unsigned long nb;
	DistribNorm *distribNorm; // A pointer on an imp score distrib used once for computing the norm parameters
	public:
	double getMean(char,double);
	double getStd(char,double);
//	double getDistribNormMean();
//	double getDistribNormStd();
//	double getSum();
//	double getSum2();
	unsigned long getNb();
	String getName();
	void setMean(double);
	void setStd(double);
//	void setSum(double);
//	void setSum2(double);
	void setNb(unsigned long);
	void setName(String); 
	void addScore(double,const String &imp);
        void newDistribNorm(Config &);
	void deleteDistribNorm();
	NormSeg();
	~NormSeg();
};
void NormSeg::addScore(double score,const String & imp){
	distribNorm->addScore(score,imp);
}
double NormSeg::getMean(char computeMode,double computePercent){
	if (!computed){
			if (!distribNorm->computeMeanStd(mean,std,computeMode,computePercent)){
			cout << "Problem: empty impostor cohort"<<endl;
			exit(0);
		    }
		computed=true;
		}
	return mean; 
}

double NormSeg::getStd(char computeMode,double computePercent){
	if (!computed){
		if (!distribNorm->computeMeanStd(mean,std,computeMode,computePercent)){
			cout << "Problem: empty impostor cohort"<<endl;
			exit(0);
		    }
		computed=true;
		}
	return std;	
}
/*
double NormSeg::getDistribNormMean(){
	return distribNorm->mean(); 
}

double NormSeg::getDistribNormStd(){
	return distribNorm->std(); 
}
*/
/*double NormSeg::getSum(){
	return sum; 
}*/

/*double NormSeg::getSum2(){
	return sum2; 
}
*/
unsigned long NormSeg::getNb(){
	return nb; 
}

String NormSeg::getName(){
	return name; 
}

void NormSeg::setMean(double m){
	mean = m; 
}

void NormSeg::setStd(double s){
	std = s; 
}
/*
void NormSeg::setSum(double m){
	sum = m; 
}
*/
/*void NormSeg::setSum2(double m){
	sum2 = m; 
}
*/
void NormSeg::setNb(unsigned long m){
	nb = m; 
}

void NormSeg::setName(String n){
	name = n; 
}

void NormSeg::newDistribNorm(Config &config){
	distribNorm=new DistribNorm(config);
}
void NormSeg::deleteDistribNorm(){
	delete distribNorm;
	distribNorm=NULL;
}
NormSeg::NormSeg(){
	 distribNorm=NULL; // initialize at NULL the distribNorm, initialized independently of NormSeg but freezed automatically
}
NormSeg::~NormSeg(){ // Clean, if needed, the DistribNorm
	if (distribNorm) deleteDistribNorm();
}

class Norm {
 	NormSeg *normTab;	 // table of normSeg => seg (or Id for znorm) name associated with mean and std values
	unsigned long nbNorm;    // number of normSeg in tabScore
	unsigned long nbNormMax; // Max number of normSeg in tabScore
	unsigned long findEntityIdxInNorm(String);             // Find where to add the seg info
        void moveEntity(unsigned long,unsigned long); // Make the room for the new Entity
	unsigned long idxLastFind;             // Save the idx of the last segment search
	char computeMode;             // The mode for mean/std computation
								  // 0=classical
								  // 1=MedianBased
								  
    double computePercent;        // % of scores used for computation.
    							  // 1-computePercent high scores are discarded								  
	public:
	bool findEntityInNorm(String,unsigned long &);   
	double getMean(unsigned long);
	double getStd(unsigned long);
	//double getSum(unsigned long);
	//double getSum2(unsigned long);
	unsigned long getNb(unsigned long);
	void print();
	unsigned long addSeg(String, double, double);
//	unsigned long addSeg(String, double, double, unsigned long);
	unsigned long findAddSeg(String,Config&); // TAKE CARE, DO NOT INITIALIZE MEAN AND COV
	void addScore(unsigned long,double,const String &);
	void deleteAllDistribNorm(); // Clean the memory (done automatically - in ~NormSeg() - but could be asked manually)
	void setNormAndFreeDistrib();// Compute Mean and Cov and free the DistribNorm (per seg)
//	void setNorm();              // Same, but not clean the DistribNorm
	Norm(unsigned long);         // Create a Norm object, with x NormSeg and Classical mean/std computation
	Norm(unsigned long, char,double); // Id but with a defined type of computation for mean/std
	~Norm();
};
void Norm::setNormAndFreeDistrib(){
	for(unsigned long idx=0; idx<nbNorm; idx++){
		normTab[idx].getMean(computeMode,computePercent);
  		normTab[idx].deleteDistribNorm();	
		}
}
/*
void Norm::setNorm(){
	for(unsigned long idx=0; idx<nbNorm; idx++){
		normTab[idx].setMean(normTab[idx].getDistribNormMean());
		normTab[idx].setStd(normTab[idx].getDistribNormStd()); 
		}
}
*/
void Norm::addScore(unsigned long idx,double score,const String & imp)
{		
	normTab[idx].addScore(score,imp);
}
void Norm::deleteAllDistribNorm(){		// clean the distribNorm to save space
	for(unsigned long i=0; i<nbNorm; i++)
		normTab[i].deleteDistribNorm();	
}
// Make the room for the new element
// Could be time consumming. good to sort the score file (.res) before
void Norm::moveEntity(unsigned long beg,unsigned long end){ 
	
        if (debug || (verboseLevel>2)) 
		cout << "moveEntity beg["<<beg<<"] end["<<end<<"]"<<endl;
	for (unsigned idx=end;idx>beg;idx--)
		normTab[idx]=normTab[idx-1];
	//memmove((void*)&(normTab[beg+1]),(void*)&(normTab[beg]),sizeof(NormSeg)*((beg-end)+1));
}
unsigned long Norm::findEntityIdxInNorm(String name){
	unsigned long beg=0;	
	unsigned long end=nbNorm-1;
	unsigned long idx=0;	
  	if (nbNorm==0) return 0; // empty...
	while ((end-beg)>1){
		idx=(beg+end)/2;
		if (name > normTab[idx].getName())
			beg=idx;
			else end=idx;
		}
	if (name <=normTab[beg].getName()) return beg;
		else if (name>normTab[end].getName()) return end+1;
			else return end; 
}
bool Norm::findEntityInNorm(String name,unsigned long & idx){
	if(nbNorm == 0) {idx=0;return false;}
	// last one searched?
	if (normTab[idxLastFind].getName() == name){
		idx=idxLastFind;
	 	return true;		
	}
        idx=findEntityIdxInNorm(name);
	if(normTab[idx].getName() == name){
		idxLastFind=idx; 
		return true;
		}
	else return false; 
}

void Norm::print(){
	cout << "nbnorm: " << nbNorm << endl; 
 	for(unsigned long i=0; i<nbNorm; i++){
		cout << normTab[i].getName() << " " << normTab[i].getMean(computeMode,computePercent) << " " << normTab[i].getStd(computeMode,computePercent) << endl; 
	} 
}

double Norm::getMean(unsigned long ind){
 	return normTab[ind].getMean(computeMode,computePercent);
}

double Norm::getStd(unsigned long ind){
 	return normTab[ind].getStd(computeMode,computePercent);
}
/*
double Norm::getSum(unsigned long ind){
 	return normTab[ind].getSum();
}

double Norm::getSum2(unsigned long ind){
 	return normTab[ind].getSum2();
}
*/
unsigned long Norm::getNb(unsigned long ind){
 	return normTab[ind].getNb();
}

unsigned long Norm::findAddSeg(String name,Config & config){ // TAKE CARE, JUST ADD THE SEG, DO NOT INITIALIZE MEAN/COV
        if(nbNorm == (nbNormMax)){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
	unsigned long idx;
	if (!findEntityInNorm(name,idx)){       // the seg is missing
		if (debug) cout <<"addSeg:seg["<<name<<"] idx["<<idx<<"]"<<endl;
		moveEntity(idx,nbNorm);                      // Make the room for the new Entity
		normTab[idx].setName(name);
		nbNorm++;
		normTab[idx].newDistribNorm(config);
	}
    return idx;
}

unsigned long Norm::addSeg(String name, double mean, double std){
        if(nbNorm == (nbNormMax)){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
	unsigned long idx=findEntityIdxInNorm(name);   // Find where to add the seg info
	if (debug) cout <<"addSeg:seg["<<name<<"] idx["<<idx<<"]"<<endl;
	moveEntity(idx,nbNorm);                      // Make the room for the new Entity
	normTab[idx].setName(name);
	normTab[idx].setMean(mean);
	normTab[idx].setStd(std); 
	nbNorm++;	
    return (idx);
}
/*
unsigned long Norm::addSeg(String name, double sum, double sum2, unsigned long nb){
	if(nbNorm == (nbNormMax)){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
	unsigned long idx=findEntityIdxInNorm(name); // Find where to add the seg info
	moveEntity(idx,nbNorm);                      // Make the room for the new Entity
	normTab[idx].setName(name);
	normTab[idx].setSum(sum);
	normTab[idx].setSum2(sum2); 
	normTab[idx].setNb(nb); 
	nbNorm++;
    return (idx);
}
*/
Norm::Norm(unsigned long nbMax){
	 computeMode=0;
	 computePercent=1.0;
	 nbNorm = 0;
	 nbNormMax = nbMax;
	 normTab = new NormSeg[nbNormMax];
	 idxLastFind=0;
}
Norm::Norm(unsigned long nbMax, char compute,double percent){
	 computeMode=compute;
	 computePercent=percent;
	 nbNorm = 0;
	 nbNormMax = nbMax;
	 normTab = new NormSeg[nbNormMax];
	 idxLastFind=0;
}
Norm::~Norm(){
  	delete [] normTab; 
}

// Get the imp scores, for one seg or once for all
//-------------------------------------------------------------------------------------------------------
bool selectImp(String &fieldOne,String& fieldTwo,XList &impList,char selectMode)
{
	bool select;
	switch (selectMode){
		case 0:select=true;break; // No selection
		case 1:select=(impList.findLine(fieldTwo) !=NULL);break; // target independent selection 
		default:select=true;cout <<"selectMode unknown. No impostor scores selection is applied"<<endl;
		}     
	return select;		       	     
}																
void getAllScores(XList &list, unsigned long indFieldOne, unsigned long indFieldTwo, unsigned long indScore, Norm &norm,XList &impList,char selectMode,Config &config){
	XLine *linep;
	// ind* are used to locate information in NIST File in the case of format change
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		String fieldOne=linep->getElement(indFieldOne);
		String fieldTwo=linep->getElement(indFieldTwo);
		double score=linep->getElement(indScore).toDouble();		      	 
		if (selectImp(fieldOne,fieldTwo,impList,selectMode)){ 
		   // If no selection of impostor scores is needed or if there is a selection and the impostor is in the list 
		   // Find the index in norm or add a seg/id in norm if needed, return the idx of the seg
		   // And initialize (memory ans others) the corresponding DistribNorm if needed
           unsigned long idx=norm.findAddSeg(fieldOne,config);
		   norm.addScore(idx,score, fieldTwo);
		   }
		}	
}	
		
// getAllScoresFirstNormed computes the score distributions for t or znorm AFTER applying firstNorm	normalization
// used for ztnorm and tznorm												
void getAllScoresFirstNormed(XList &list, unsigned long indFieldOne, unsigned long indFieldTwo,unsigned long indScore, Norm &firstNorm, 
 			Norm &outNorm,XList &impList,char selectMode,Config &config){
	XLine *linep;
	// TODO Replace 0 and 3 by the indElt, indIdImp
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		String fieldOne=linep->getElement(indFieldOne);
		String fieldTwo=linep->getElement(indFieldTwo);
		double score=linep->getElement(indScore).toDouble();	 	 
		unsigned long ind; 		
		if (selectImp(fieldOne,fieldTwo,impList,selectMode)){ 
		   	// If no selection of impostor scores is needed or if there is a selection and the impostor is in the list 
			if(firstNorm.findEntityInNorm(fieldTwo,ind)){
			    unsigned long idx=outNorm.findAddSeg(fieldOne,config);
				outNorm.addScore(idx,(score - firstNorm.getMean(ind)) / firstNorm.getStd(ind), fieldTwo);
		        }
			else{
				cout << "distribution for["<<fieldTwo<<"] not found!" << endl;
				cout <<"line:"<<linep<<endl;
				exit(-1); 
				}
		}
	}
}		
//-------------------------------------------------------------------------------------------------------
int ComputeNorm(Config& config)
{

	using namespace alize;
	using namespace std;

	int decision=0; // TODO change ASAP
    String outputNISTFileName = config.getParam("outputFileBaseName");                        // Result file BAsenamein NIST style (.nist) result file format
    String znormFilesExtension = config.getParam("znormFilesExtension");                      // Result file znorm extension
    String tnormFilesExtension = config.getParam("tnormFilesExtension");                      // Result file tnorm extension
    String tznormFilesExtension = config.getParam("tznormFilesExtension");                    // Result file tztnorm extension
    String ztnormFilesExtension = config.getParam("ztnormFilesExtension");                    // Result file tztnorm extension
  
	// Norm Type = znorm, tnorm ou ztnorm
	unsigned long maxIdNb = (unsigned long)(config.getParam("maxIdNb").toLong());
	unsigned long maxSegNb = (unsigned long)(config.getParam("maxSegNb").toLong());
	String normType = config.getParam("normType");
	String testNistFile = config.getParam("testNistFile"); 
    char selectMode=0;       // 0 if no selection, 1 if target independent selection
    XList impList;
    if (config.existsParam("impostorIDList")){
    	selectMode=1; // Target independent Impostor score selection mode  	
    	impList.load(config.getParam("impostorIDList"),config);
    }  
	char computeMode = (char)(config.getParam("meanMode").toLong()); // 0 classical, 1 Median
	double computePercent=(double)(config.getParam("percentScoreUsed").toDouble());
	// Define the field of the score files
	char fieldGender=(char) config.getParam("fieldGender").toLong();
	char fieldName=(char)config.getParam("fieldName").toLong();
	//char fieldDecision=(char)config.getParam("fieldDecision").toLong(); // not used
	char fieldSeg=(char)config.getParam("fieldSeg").toLong();
	char fieldLLR=(char)config.getParam("fieldLLR").toLong();
	
	XList testList(testNistFile, config);
	testList.getLine(0); // goto first line still to debug
	
	try{		
	if(normType == "tnorm"){	
		String outputFilename=outputNISTFileName+tnormFilesExtension;	   // Create the complete output file filename
		ofstream outFile(outputFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
	 	Norm tnorm(maxSegNb,computeMode,computePercent);// storage of mean and std for all segments
		String tnormNistFile = config.getParam("tnormNistFile"); 
		if (verbose) cout << "Tnorm, reading tnormList"<<endl;
		XList tnormList(tnormNistFile, config);
		getAllScores(tnormList, fieldSeg, fieldName, fieldLLR, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)

		if (verbose) cout << "Tnorm, begin the test list"<<endl;
	 	XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);
			unsigned long ind;
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm impostor score not found for seg["<<seg<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 										
			double score = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);
	        outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFile);	
		}
		outFile.close();    			
	    if(debug||verbose){
			cout << "tnorm distrib" << endl;
	       	tnorm.print();       	 
		}
		      
	}
    else if(normType == "znorm"){	
    	String outputFilename=outputNISTFileName+znormFilesExtension;	   // Create the complete output file filename
		ofstream outFile(outputFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute Znorm" << endl; 
		// storage of mean and std for all segments
	 	Norm znorm(maxIdNb,computeMode,computePercent);
		String znormNistFile = config.getParam("znormNistFile"); 	
		XList znormList(znormNistFile, config);
	 	XLine * linep;		
	 	if (verbose) cout << "Computing once the Znorm distributions"<<endl<<
							 "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(znormList, fieldName, fieldSeg, fieldLLR, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String id=linep->getElement(fieldName);
			unsigned long ind; 	
			if(verboseLevel>2) cout << endl << "Compute Znorm [" << id << "]"<< endl; 							
			if(debug) cout << endl << "findEntityInNorm Id[" << id << "]"; 				
			if (!znorm.findEntityInNorm(id,ind)){// Problem, the seg parameters are not present...
				cout << "znorm impostor score not found for id["<<id<<endl;
				exit(-1); 
			    }
			if(debug) cout << " idx["<< ind<< "]"<<endl; 						
			double score = ((linep->getElement(fieldLLR)).toDouble() - znorm.getMean(ind)) / znorm.getStd(ind);
	        outputResultLine(score, id, linep->getElement(fieldSeg), linep->getElement(fieldGender), decision, outFile);
		    }	 
		outFile.close();		    
	    if(debug||verbose){
			cout << "znorm distrib" << endl;
	       	znorm.print();       	 
		    }	
	    }
     else if(normType == "ztnorm"){ // Ztnorm outputs ztnorm AND tnorm results
     	String outputTnormFilename=outputNISTFileName+tnormFilesExtension;	    // Create the complete output file filename
     	String outputZtnormFilename=outputNISTFileName+ztnormFilesExtension;	// Create the complete output file filename
     	
		ofstream outFileTnorm(outputTnormFilename.c_str(),ios::out | ios::trunc);      // Initialise the output file
		ofstream outFileZtnorm(outputZtnormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute ZTnorm (and tnorm)" << endl; 
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb,computeMode,computePercent);
	 	Norm znorm(maxIdNb,computeMode,computePercent);
		Norm ztnorm(maxSegNb,computeMode,computePercent);
		
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String znormNistFile = config.getParam("znormNistFile"); 
		String ztnormNistFile = config.getParam("ztnormNistFile"); 
		XList tnormList(tnormNistFile, config);
		XList znormList(znormNistFile, config);
		XList ztnormList(ztnormNistFile, config);
				
		// tnorm distribution for the impostor segments (from imp_imp) in ztnorm		
		if (verbose) cout << "Computing once the Tnorm distribution for impostor segments"<<endl<<
							 "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(ztnormList, fieldSeg, fieldName, fieldLLR, ztnorm,impList,selectMode,config);
		ztnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		// compute the tnorm distributions for the (test) seg in tnorm
	    if (verbose) cout << "Computing once all the tnorm distributions for the test segments "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(tnormList, fieldSeg, fieldName, fieldLLR, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		
		// compute the znorm distrib for the tnormed scores in znorm
		if (verbose) cout << "Computing once all the znorm (tnormed) distributions "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScoresFirstNormed(znormList, fieldName, fieldSeg, fieldLLR, ztnorm, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib();// finalize the mean/cov computation and free the score (the per seg DistribNorm)
	
		// for all the test lines (one test seg against several target speakers)
		XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);    // the test seg 
			String id= linep->getElement(fieldName);   // the Id (target spk)
			if (verboseLevel>2) cout << "ztnorm seg["<<seg<<"] Id["<<id<<"]"<<endl;
			unsigned long ind; 		
			// segment tnorm distribution already computed ?
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm impostor distribution not found for seg["<<seg<<"]"<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 	
			double scoreTNorm = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);

			// ID mean and std already computed ?
			if(!znorm.findEntityInNorm(id,ind)){// Problem, the seg parameters are not present...
				cout << "znorm(tnormed) distribution not found for id["<<id<<"]"<<endl;
				exit(-1); 
				}
			double score = (scoreTNorm - znorm.getMean(ind)) / znorm.getStd(ind);
			
		    outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileZtnorm);
		    outputResultLine(scoreTNorm, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileTnorm);
		}	      
		outFileZtnorm.close();
		outFileTnorm.close();
	    if(debug||verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();       	 
		}
	}     
	   else if(normType == "tznorm"){ // tznorm outputs tznorm and znorm scores
     	String outputTznormFilename=outputNISTFileName+tznormFilesExtension;	// Create the complete output file filename  	
     	String outputZnormFilename=outputNISTFileName+znormFilesExtension;	    // Create the complete output file filename  	
		ofstream outFileTznorm(outputTznormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		ofstream outFileZnorm(outputZnormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute TZnorm" << endl; 
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb,computeMode,computePercent);
	 	Norm znorm(maxIdNb,computeMode,computePercent);
		Norm tznorm(maxSegNb,computeMode,computePercent);
		
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String znormNistFile = config.getParam("znormNistFile"); 
		String tznormNistFile = config.getParam("ztnormNistFile"); 
		XList tnormList(tnormNistFile, config);
		XList znormList(znormNistFile, config);
		XList tznormList(tznormNistFile, config);
		
		// compute znorm distribs for the target speakers ID in znorm 
		if (verbose) cout << "Compute Znorm distribs for the target speakers (using znorm scores)"<<endl
						  <<"Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(znormList, fieldName, fieldSeg, fieldLLR, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)		
		
		// compute znorm distribs for the impostor ID in tznorm 
		if (verbose) cout << "Compute znorm distribs for the impostor speakers (using ztnorm scores)"
						  << "Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(tznormList, fieldName, fieldSeg, fieldLLR, tznorm,impList,selectMode,config);
		tznorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)		
	
		// compute tnorm distribs for the znormed scores
		if (verbose) cout << "Compute tnorm distribs of the znormed scores "<<endl
						  << "Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScoresFirstNormed(tnormList, fieldSeg, fieldName, fieldLLR, tznorm, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib();// finalize the mean/cov computation and free the score (the per seg DistribNorm)
		if(debug|| verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();    
		    cout << "tznorm distrib" << endl;
		   	tznorm.print();       	   	 
		}		
		// for all the test lines (one test seg against several target speakers)
		if (verbose) cout << "begin the test list"<<endl;
		XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);    // the test seg 
			String id= linep->getElement(fieldName);   // the Id (target spk)
			if (verboseLevel>2) cout << "tznorm seg["<<seg<<"] Id["<<id<<"]"<<endl;
			unsigned long ind; 		
			// id mean and std already computed ?
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!znorm.findEntityInNorm(id,ind)){// Problem, the id parameters are not present...
				cout << "znorm distribution not found for id["<<id<<"]"<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 	
			double scoreZNorm = ((linep->getElement(fieldLLR)).toDouble() - znorm.getMean(ind)) / znorm.getStd(ind);
			// segment mean and std already computed ?
			if(!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm(znormed) distribution not found for seg["<<seg<<"]"<<endl;
				exit(-1); 
				}
			double score = (scoreZNorm - tnorm.getMean(ind)) / tnorm.getStd(ind);
			
		    outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileTznorm);
		    outputResultLine(scoreZNorm, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileZnorm);
		    if (verboseLevel>2) outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, cout);

		}	      
		outFileTznorm.close();
		outFileZnorm.close();
	    if(debug || verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();       	 
		}
	}   
	else {
		  cout << "unknown normalization mode:"<<normType<<endl;
		  exit(-1); 
		  }   
	} //try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}

              
    return 0;
}

#endif //!defined(ALIZE_ComputeNorm_cpp)