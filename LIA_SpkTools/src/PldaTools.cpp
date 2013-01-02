/*
This file is part of LIA_RAL which is a set of software based on ALIZE
toolkit for speaker recognition. ALIZE toolkit is required to use LIA_RAL.

LIA_RAL project is a development project was initiated by the computer
science laboratory of Avignon / France (Laboratoire Informatique d'Avignon -
LIA) [http://lia.univ-avignon.fr <http://lia.univ-avignon.fr/>]. Then it
was supported by two national projects of the French Research Ministry:
	- TECHNOLANGUE program [http://www.technolangue.net]
	- MISTRAL program [http://mistral.univ-avignon.fr]

LIA_RAL is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of
the License, or any later version.

LIA_RAL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with LIA_RAL.
If not, see [http://www.gnu.org/licenses/].

The LIA team as well as the LIA_RAL project team wants to highlight the
limits of voice authentication in a forensic context.
The "Person Authentification by Voice: A Need of Caution" paper
proposes a good overview of this point (cf. "Person
Authentification by Voice: A Need of Caution", Bonastre J.F.,
Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-
chagnolleau I., Eurospeech 2003, Genova].
The conclusion of the paper of the paper is proposed bellow:
[Currently, it is not possible to completely determine whether the
similarity between two recordings is due to the speaker or to other
factors, especially when: (a) the speaker does not cooperate, (b) there
is no control over recording equipment, (c) recording conditions are not
known, (d) one does not know whether the voice was disguised and, to a
lesser extent, (e) the linguistic content of the message is not
controlled. Caution and judgment must be exercised when applying speaker
recognition techniques, whether human or automatic, to account for these
uncontrolled factors. Under more constrained or calibrated situations,
or as an aid for investigative purposes, judicious application of these
techniques may be suitable, provided they are not considered as infallible.
At the present time, there is no scientific process that enables one to
uniquely characterize a persones voice or to identify with absolute
certainty an individual from his or her voice.]

Copyright (C) 2004-2010
Laboratoire d'informatique d'Avignon [http://lia.univ-avignon.fr]
LIA_RAL admin [alize@univ-avignon.fr]
Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
*/

#if !defined(ALIZE_PldaTools_cpp)
#define ALIZE_PldaTools_cpp

#if defined(_WIN32)
  #include <cfloat> // for _isnan()
  #define ISNAN(x) _isnan(x)
  #define ISINF(x) (!_finite(x))
#elif defined(linux) || defined(__linux) || defined(__CYGWIN__) || defined(__APPLE__)
  #define ISNAN(x) isnan(x)
  #define ISINF(x) isinf(x)
#else
  #error "Unsupported OS\n"
#endif

#include<PldaTools.h>
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
#include "RealVector.h"
#include <liatools.h>
#include <limits>
#ifdef THREAD
#include <pthread.h>
#endif

using namespace alize;
using namespace std;



//***************************************************************//
//			Constructeurs et Destructeurs				 //
//***************************************************************//

PldaDev::PldaDev(){}
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev::PldaDev(String & ndxFilename,Config & config){

	// Create the XList and sort speakers by decreasing number of sessions
	_fileList.load(ndxFilename,config);
	_fileList.sortByElementNumber("descend");

	_n_speakers = _fileList.getLineCount();
	_session_per_speaker.setSize(_n_speakers);
	_n_sessions = 0;

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	linep=_fileList.getLine();
	String fileName = linep->getElement(0);

	String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
	Matrix<double> clientV(fName,config);
	_vectSize = clientV.cols();

	// 	the total number of sessions
	_fileList.getLine(0);
	unsigned long c = 0;
	while ((linep=_fileList.getLine()) != NULL){
		_session_per_speaker[c] = linep->getElementCount();
		_n_sessions += linep->getElementCount();
		c++;
	}

	// Initialize data, mean vector and matrix
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);

	_speaker_means.setDimensions(_vectSize,_n_speakers);

	_speaker_means.setAllValues(0.0);

	_class.setSize(_n_sessions);
	_style.setSize(_n_sessions);

	// Read vectors and fill the _data Matrix
	_data.setDimensions(_vectSize,_n_sessions);

	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	
	double *data, *speaker_mean, *tv;
	data = _data.getArray();

	speaker_mean = _speaker_means.getArray();

	while ((linep=_fileList.getLine()) != NULL){

		if(verboseLevel>0) cout<<"	Load class [ "<<LCounter+1<<" / "<<_n_speakers<<" ]"<<endl; 

		for(unsigned long s=0;s<linep->getElementCount();s++){

			// Read the vector from file
			fileName = linep->getElement(s);
			String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");

			if(verboseLevel>1) cout<<"		Load file "<<fName<<endl;

			Matrix<double> tmpVect(fName,config);
			tv = tmpVect.getArray();

			if((tmpVect.rows() !=1)||(tmpVect.cols()!=_vectSize))
				throw Exception("Incorrect dimension of vector to load",__FILE__,__LINE__);

			_class[SCounter]	= LCounter;
			_style[SCounter]	= 0;	// to modify for future implementation of Tied_FA and Tied-PLDA

			for(unsigned long k=0;k<_vectSize;k++){
				data[k*_n_sessions+SCounter] = tv[k];
				speaker_mean[k*_n_speakers+LCounter] += tv[k];
				_mean[k] += tmpVect(0,k);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev::PldaDev(XList & fileList,Config & config){

	_fileList.sortByElementNumber("descend");

	_n_speakers = _fileList.getLineCount();
	_session_per_speaker.setSize(_n_speakers);
	_n_sessions = 0;

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	linep=_fileList.getLine();
	String fileName = linep->getElement(0);

	String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
	Matrix<double> clientV(fName,config);
	_vectSize = clientV.cols();

	// Compute the total number of sessions
	_fileList.getLine(0);
	unsigned long c = 0;
	while ((linep=_fileList.getLine()) != NULL){
		_session_per_speaker[c] = linep->getElementCount();
		_n_sessions += linep->getElementCount();
		c++;
	}

	// Initialize data, mean vector and matrix
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);

	_speaker_means.setDimensions(_vectSize,_n_speakers);
	_speaker_means.setAllValues(0.0);

	_class.setSize(_n_sessions);
	_style.setSize(_n_sessions);

	// Read vectors and fill the _data Matrix
	_data.setDimensions(_vectSize,_n_sessions);

	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){

			// Read the vector from file
			fileName = linep->getElement(s);
			String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
			Matrix<double> tmpVect(fName,config);

			// verifier la dimension de la matrice lue...
			if((tmpVect.rows() !=1)||(tmpVect.cols()!=_vectSize))
				throw Exception("Incorrect dimension of vector to load",__FILE__,__LINE__);

			_class[SCounter]	= LCounter;
			_style[SCounter]	= 0;	// to modify for future implementation of Tied_FA and Tied-PLDA

			for(unsigned long k=0;k<_vectSize;k++){
				_data(k,SCounter) = tmpVect(0,k);
				_speaker_means(k,LCounter) += tmpVect(0,k);
				_mean[k] += tmpVect(0,k);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev::~PldaDev(){}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaDev::getClassName() const{	return "PldaDev";}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::load(String & ndxFilename,Config & config){

	// Create the XList and sort speakers by decreasing number of sessions
	_fileList.load(ndxFilename,config);
	_fileList.sortByElementNumber("descend");

	_n_speakers = _fileList.getLineCount();
	_session_per_speaker.setSize(_n_speakers);
	_n_sessions = 0;

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	linep=_fileList.getLine();
	String fileName = linep->getElement(0);

	String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
	Matrix<double> clientV(fName,config);
	_vectSize = clientV.cols();

	// 	the total number of sessions
	_fileList.getLine(0);
	unsigned long c = 0;
	while ((linep=_fileList.getLine()) != NULL){
		_session_per_speaker[c] = linep->getElementCount();
		_n_sessions += linep->getElementCount();
		c++;
	}

	// Initialize data, mean vector and matrix
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);

	_speaker_means.setDimensions(_vectSize,_n_speakers);
	_speaker_means.setAllValues(0.0);

	_class.setSize(_n_sessions);
	_style.setSize(_n_sessions);

	// Read vectors and fill the _data Matrix
	_data.setDimensions(_vectSize,_n_sessions);

	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){

			// Read the vector from file
			fileName = linep->getElement(s);
			String fName = config.getParam("loadVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
			Matrix<double> tmpVect(fName,config);

			if((tmpVect.rows() !=1)||(tmpVect.cols()!=_vectSize))
				throw Exception("Incorrect dimension of vector to load",__FILE__,__LINE__);

			_class[SCounter]	= LCounter;
			_style[SCounter]	= 0;	// to modify for future implementation of Tied_FA and Tied-PLDA

			for(unsigned long k=0;k<_vectSize;k++){
				_data(k,SCounter) = tmpVect(0,k);
				_speaker_means(k,LCounter) += tmpVect(0,k);
				_mean[k] += tmpVect(0,k);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeAll(){

	//Update the size of vectors
	_vectSize = _data.rows();
	
	// Reset all mean values
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);
	_speaker_means.setDimensions(_vectSize,_n_speakers);
	_speaker_means.setAllValues(0.0);

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){
			for(unsigned long k=0;k<_vectSize;k++){
				_speaker_means(k,LCounter) += _data(k,SCounter);
				_mean[k] += _data(k,SCounter);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getVectSize(){
	return _vectSize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getSpeakerNumber(){
	return  _n_speakers;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getSessionNumber(){
	return _n_sessions;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getSpeakerSessionNumber(unsigned long spkIdx){
	return _session_per_speaker[spkIdx];
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaDev::getData(){
	return _data;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
double PldaDev::getData(unsigned long i, unsigned long j){
	return _data(i,j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
RealVector<double>& PldaDev::getMean(){
	return _mean;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::getSpeakerMean(unsigned long spk,RealVector<double>& spkMean){
	spkMean.setSize(_vectSize);
	for(unsigned long i=0;i<_vectSize;i++)
		spkMean[i] = _speaker_means(i,spk);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::lengthNorm(){

	if(verboseLevel>0) cout<<"	DEV: Normalize length of ... ";

	// Compute the norm of all vectors in _data
	RealVector<double> vecNorm(0);
	vecNorm.setSize(_n_sessions);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_sessions;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _data(k,s)*_data(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_sessions;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_data(k,s) /=vecNorm[s];
		}
	}
	// Recompute the means
	this->computeAll();

	if(verboseLevel>0) cout<<"done"<<endl;;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::center(RealVector<double> &mu){
	for(unsigned long s =0;s<_n_sessions;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_data(k,s) -= mu[k];
		}
	}
	// Recompute the means
	this->computeAll();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::center(Eigen::VectorXd &mu){
	for(unsigned long s =0;s<_n_sessions;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_data(k,s) -= mu(k);
		}
	}
	// Recompute the means
	this->computeAll();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::rotateLeft(Matrix<double> &M){

	if(M.cols() != _data.rows())	cerr<<"Rotation dimension mismatch !"<<endl<<"	Rootation matrix: "<<M.rows()<<" x "<<M.cols()<<endl<<"	Vector:	"<<_data.rows()<<endl;

	Matrix<double> tmpData(_data);
	_data.setDimensions(M.rows(),_n_sessions);
	_data.setAllValues(0.0);

	// Rotate _data
	_data = M*tmpData;

	_vectSize = M.rows();

	// Recompute the means
	this->computeAll();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMat(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeCovMatThreaded(Sigma,W,B,config.getParam("numThread").toULong());
	else computeCovMatUnThreaded(Sigma,W,B);
	#else
	computeCovMatUnThreaded(Sigma,W,B);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatUnThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B){

	if(verboseLevel>0) cout<<"	DEV: Compute Covariance Matrices unthreaded ... ";

	// Initialize matrices
	Sigma.setSize(_vectSize);
	Sigma.setAllValues(0.0);
	W.setSize(_vectSize);
	W.setAllValues(0.0);
	B.setSize(_vectSize);
	B.setAllValues(0.0);

	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			for(unsigned long s=0;s<_n_sessions;s++){
				Sigma(i,j)	+= (_data(i,s)-_mean[i])*(_data(j,s)-_mean[j]);
				// Sigma(i,j)	+= (_data(i,s)-_mean[i])*(_data(j,s)-_mean[j]);
				W(i,j)		+= (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
			}

			for(unsigned long c=0;c<_n_speakers;c++){
				B(i,j)		+= _session_per_speaker[c]* (_speaker_means(i,c)-_mean[i]) * (_speaker_means(j,c)-_mean[j]);	
			}
		}
	}

	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			Sigma(i,j)	/= _n_sessions;
			W(i,j)		/= _n_sessions;
			B(i,j)		/= _n_sessions;
		}
	}
	if(verboseLevel>0) cout<<"done"<<endl;
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct Covthread_data{

	double *data;
	double *speaker_means;
	double *mean;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *sigma;
	RefVector <DoubleSquareMatrix> *w;
	RefVector <DoubleSquareMatrix> *b;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *Covthread(void *threadarg){
	struct Covthread_data *my_data;
	my_data = (struct Covthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	double *mean = my_data->mean;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &_sigma=(*(my_data->sigma))[threadNb];
	DoubleSquareMatrix &_w=(*(my_data->w))[threadNb];
	DoubleSquareMatrix &_b=(*(my_data->b))[threadNb];

	double *sigma = _sigma.getArray();
	double *w = _w.getArray();
	double *b = _b.getArray();

	for(unsigned long i=0;i<vectSize;i++){
		for(unsigned long j=0;j<vectSize;j++){
			for(unsigned long s=startSession;s<stopSession;s++){
				_sigma(i,j)	+= (data[i*n_sessions+s]-mean[i])*(data[j*n_sessions+s]-mean[j]);
				_w(i,j)		+= (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
			}
			for(unsigned long c=_class[startSession];c<_class[stopSession-1]+1;c++){
				_b(i,j)		+= session_per_speaker[c]* (speaker_means[i*n_speakers+c]-mean[i]) * (speaker_means[j*n_speakers+c]-mean[j]);
			}
		}
	}
	
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	Sigma.setSize(_vectSize); W.setSize(_vectSize); B.setSize(_vectSize);
	Sigma.setAllValues(0.0); W.setAllValues(0.0); B.setAllValues(0.0);

	// Compute the covariance matrices
	if(verboseLevel>0) cout<<"	DEV: Compute Covariance Matrices threaded ... ";

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the 3 RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> sigMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		sigMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		sigMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> wMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		wMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		wMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> bMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		bMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		bMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct Covthread_data *thread_data_array = new Covthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data, *mean;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); mean = _mean.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t].mean = mean;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].sigma=&(sigMatrices);
		thread_data_array[t].w=&(wMatrices);
		thread_data_array[t].b=&(bMatrices);

		if (verboseLevel > 0) cout<<"(computeCovMat) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, Covthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeCovMat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				Sigma(i,j) += sigMatrices[t](i,j);
				W(i,j) += wMatrices[t](i,j);
				B(i,j) += bMatrices[t](i,j);
			}
		}
	}
	sigMatrices.deleteAllObjects();
	wMatrices.deleteAllObjects();
	bMatrices.deleteAllObjects();

	// Normalize 
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			Sigma(i,j)	/= _n_sessions;
			W(i,j)		/= _n_sessions;
			B(i,j)		/= _n_sessions;
		}
	}

	if(verboseLevel>0) cout<<"done"<<endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatEigen(Eigen::MatrixXd &Sigma, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeCovMatEigenThreaded(Sigma,config.getParam("numThread").toULong());
	else computeCovMatEigenUnThreaded(Sigma);
	#else
	computeCovMatEigenUnThreaded(Sigma);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatEigenUnThreaded(Eigen::MatrixXd &Sigma){

	if(verboseLevel>0) cout<<"	DEV: Compute Total covariance matrix unthreaded ... ";

	// Initialize matrices
	Sigma = Eigen::MatrixXd::Zero(_vectSize,_vectSize);

	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			for(unsigned long s=0;s<_n_sessions;s++){
				Sigma(i,j)	+= (_data(i,s)-_mean[i])*(_data(j,s)-_mean[j]);
			}
		}
	}
	if(verboseLevel>0) cout<<"done"<<endl;
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct CovEigenthread_data{

	double *data;
	double *speaker_means;
	double *mean;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *sigma;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *CovEigenthread(void *threadarg){
	struct CovEigenthread_data *my_data;
	my_data = (struct CovEigenthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	double *mean = my_data->mean;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &_sigma=(*(my_data->sigma))[threadNb];

	double *sigma = _sigma.getArray();

	for(unsigned long i=0;i<vectSize;i++){
		for(unsigned long j=0;j<vectSize;j++){
			for(unsigned long s=startSession;s<stopSession;s++){
				_sigma(i,j)	+= (data[i*n_sessions+s]-mean[i])*(data[j*n_sessions+s]-mean[j]);
			}
		}
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatEigenThreaded(Eigen::MatrixXd &Sigma, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	DoubleSquareMatrix tmpSigma(_vectSize);
	tmpSigma.setAllValues(0.0);
	Sigma.resize(_vectSize,_vectSize);

	// Compute the covariance matrices
	if(verboseLevel>0) cout<<"	DEV: Compute Covariance Matrices threaded ... ";

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the 3 RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> sigMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		sigMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		sigMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct CovEigenthread_data *thread_data_array = new CovEigenthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data, *mean;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); mean = _mean.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t].mean = mean;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].sigma=&(sigMatrices);

		if (verboseLevel > 0) cout<<"(computeCovMat) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, CovEigenthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeCovMat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				Sigma(i,j) += sigMatrices[t](i,j);
			}
		}
	}
	sigMatrices.deleteAllObjects();

	if(verboseLevel>0) cout<<"done"<<endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnChol(DoubleSquareMatrix &WCCN, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeWccnCholThreaded(WCCN,config.getParam("numThread").toULong());
	else computeWccnCholUnThreaded(WCCN);
	#else
	computeWccnCholUnThreaded(WCCN);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnCholUnThreaded(DoubleSquareMatrix &WCCN){

	if(verbose) cout<<"Compute WCCN Matrix unThreaded"<<endl;

	WCCN.setSize(_vectSize);
	WCCN.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize);
	W.setAllValues(0.0);

	// Compute the covariance matrix per class
	if(verboseLevel > 0) cout<<"	Compute WCCN matrix"<<endl;
	unsigned long spk = 0;
	unsigned long s = 0;

	while(s < _n_sessions){
		DoubleSquareMatrix covSpk(_vectSize);
		covSpk.setAllValues(0.0);
		unsigned long sessionNumber = _session_per_speaker[spk];
		while((s<_n_sessions)&&(_class[s] == spk)){
			for(unsigned long i=0;i<_vectSize;i++){
				for(unsigned long j=0;j<_vectSize;j++){
					covSpk(i,j) += (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
				}
			}
			s++;
		}
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				W(i,j) += covSpk(i,j)/(double)sessionNumber;
			}
		}
		if(s<_n_sessions)	spk = _class[s];
	}

	// Normalize the WCCN Matrix by the number of speakers
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			W(i,j) /= _n_speakers;
		}
	}

    // Invert the WCCN matrix
	if(verboseLevel > 0) cout<<"	Invert WCCN matrix"<<endl;
	DoubleSquareMatrix invW(_vectSize);
	W.invert(invW);

	// Choleski decompostion of WCCN inverse
	if(verboseLevel > 0) cout<<"	Cholesky decomposition of inv(WCCN)"<<endl;
	invW.upperCholesky(WCCN);

}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct WCCNthread_data{

	double *data;
	double *speaker_means;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *cov;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *WCCNthread(void *threadarg){
	struct WCCNthread_data *my_data;
	my_data = (struct WCCNthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &cov=(*(my_data->cov))[threadNb];
	double *c = cov.getArray();
	unsigned long spk = _class[startSession];
	unsigned long s = startSession;

	while(s < stopSession){
		DoubleSquareMatrix covSpk(vectSize);
		covSpk.setAllValues(0.0);
		unsigned long sessionNumber = session_per_speaker[spk];
		while((s<stopSession)&&(_class[s] == spk)){
			for(unsigned long i=0;i<vectSize;i++){
				for(unsigned long j=0;j<vectSize;j++){
					covSpk(i,j) += (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
				}
			}
			s++;
		}
		for(unsigned long i=0;i<vectSize;i++){
			for(unsigned long j=0;j<vectSize;j++){
				c[i*vectSize+j] += covSpk(i,j)/(double)sessionNumber;
			}
		}
		if(s<n_sessions)	spk = _class[s];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnCholThreaded(DoubleSquareMatrix &WCCN, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	WCCN.setSize(_vectSize);
	WCCN.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize);
	W.setAllValues(0.0);

	// Compute the covariance matrix per class
	if(verboseLevel > 0) cout<<"	Compute WCCN matrix"<<endl;

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> covMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		covMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		covMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct WCCNthread_data *thread_data_array = new WCCNthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].cov=&(covMatrices);

		if (verboseLevel > 0) cout<<"(computeWccnChol) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, WCCNthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >0) cout <<"(computeWccnChol) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				W(i,j) += covMatrices[t](i,j);
			}
		}
	}
	covMatrices.deleteAllObjects();

	// Normalize the WCCN Matrix by the number of speakers
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			W(i,j) /= _n_speakers;
		}
	}

    // Invert the WCCN matrix
	if(verboseLevel > 0) cout<<"	Invert WCCN matrix"<<endl;
	DoubleSquareMatrix invW(_vectSize);
	W.invert(invW);

	// Choleski decompostion of WCCN inverse
	if(verboseLevel > 0) cout<<"	Cholesky decomposition of inv(WCCN)";
	invW.upperCholesky(WCCN);
	if(verboseLevel > 0) cout<<"	-> done"<<endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeMahalanobis(DoubleSquareMatrix &M, Config& config){

	// Initialize matrices
	M.setSize(_vectSize);
	M.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize), B(_vectSize), Sigma(_vectSize);
	W.setAllValues(0.0); B.setAllValues(0.0); Sigma.setAllValues(0.0);

	this->computeCovMat(Sigma,W,B,config);

	W.invert(M);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeLDA(Matrix<double> &ldaMat, long ldaRank, Config& config){

	// Initialize matrices
	if (verboseLevel >0)  cout << "LDA reduction from  " << _vectSize << " --> "  << ldaRank <<endl;

	DoubleSquareMatrix W(_vectSize), B(_vectSize), Sigma(_vectSize);
	W.setAllValues(0.0); B.setAllValues(0.0); Sigma.setAllValues(0.0);

	this->computeCovMat(Sigma,W,B,config);
	if (verboseLevel >0) cout << "- done " <<endl;

	DoubleSquareMatrix invW(_vectSize);
	invW.setAllValues(0.0);
	W.invert(invW);

	//  use Matrix object  ( invert and computeCovMat worked on DoubleSquareMatrix but now we work with   Matrix <double> )
	Matrix <double> MB(B);
	Matrix <double> MinvW(invW);

	Matrix <double> EP = MinvW*MB;
	Matrix<double> tmpLdaMat;
	tmpLdaMat.setDimensions (_vectSize,ldaRank)  ;
	tmpLdaMat.setAllValues(0.0); 
	
	Matrix<double> eigenVal;
	eigenVal.setDimensions (ldaRank,ldaRank);
	eigenVal.setAllValues(0.0);
	
	this->computeEigenProblem(EP,tmpLdaMat,eigenVal,ldaRank,config);
	ldaMat = tmpLdaMat.transpose();
  
	/*
	cout << "LDA done some EVal and Evect " <<endl;
	for(int i=0;i<10;i++){
		for(int j=0;j<3;j++){  cout << eigenVect(i,j) << "\t" ;}
	cout << endl;}
	*/
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect, long rank, Config& config)
{
	Matrix <double> eigenVal;
	eigenVal.setDimensions (rank,rank);
	eigenVal.setAllValues(0.0);
	this->computeEigenProblem(EP,eigenVect,eigenVal,rank,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect,Matrix<double> &eigenVal, long rank, Config& config)
{
	if(verboseLevel>0)	cout<<"(PldaDev) Compute Eigen Problem"<<endl;

	// convert ALIZE matrix into Eigen::MatrixXd
	Eigen::MatrixXd A(EP.rows(),EP.cols());
	for(unsigned long i=0;i<EP.rows();i++){
		for(unsigned long j=0;j<EP.cols();j++)
		A(i,j) = EP(i,j);
	}	
		
	// Compute Eigen Decomposition
	Eigen::EigenSolver<Eigen::MatrixXd> es(A);
	complex<double> lambda = es.eigenvalues()[0];
	Eigen::MatrixXcd V = es.eigenvectors();

	// get and check eigen values
	LKVector EV(0,0);
	unsigned long imagPart = 0;
	for(unsigned long i=0;i<_vectSize;i++){ 
		if(imag(es.eigenvalues()[i])!=0)	imagPart++;
		LKVector::type s;
		s.idx = i;
		s.lk = real(es.eigenvalues()[i]);
		EV.addValue(s);
	}
	if(imagPart>0) cout << "WARNING "<<imagPart<<" eigenvalues have an imaginary part" << endl;
	eigenVal.setAllValues(0.0);

	// Order the EigenValues
	EV.descendingSort();
	for(unsigned long k=0; k<_vectSize;k++){
		for(int j=0;j<rank;j++){
			eigenVect(k,j)= real(V(k,j));
		}
	}
	for(int j=0;j<rank;j++){ 
		eigenVal(j,0) = EV[j].lk;
	}
	
	if(verboseLevel>1){
		cerr<<"EigenValues"<<endl;
		for(int i=0;i<rank;i++){ 
			cerr<<eigenVal(i,0)<<endl;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getClass(unsigned long s){
	return(_class[s]);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex){

//	unsigned long spkPerThread = floor((double)_n_speakers/(double)nbThread);
	unsigned long spkPerThread = _n_speakers/nbThread;
	startIndex.setSize(nbThread);
	startIndex[0] = 0;
	
	unsigned long firstSpeakerNextList = spkPerThread;
	unsigned long spkCounter = 1;
	unsigned long currentSession = 0;

	while((currentSession<_n_sessions)&&(spkCounter<nbThread)){
		while(_class[currentSession]<spkCounter*spkPerThread){
			currentSession++;
		}
		startIndex[spkCounter] = currentSession;
		spkCounter++;
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMat(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeScatterMatThreaded(SB,SW,config.getParam("numThread").toULong());
	else computeScatterMatUnThreaded(SB,SW);
	#else
	computeScatterMatUnThreaded(SB,SW);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMatUnThreaded(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW){

	// Initialize matrices
	SB.setSize(_vectSize);
	SB.setAllValues(0.0);
	SW.setSize(_vectSize);
	SW.setAllValues(0.0);

	unsigned long sessionCounter = 0;
	for(unsigned long c=0;c<_n_speakers;c++){

		DoubleSquareMatrix tmpSB;
		tmpSB.setSize(_vectSize);
		tmpSB.setAllValues(0.0);

		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SB(i,j)		+= (_speaker_means(i,c)-_mean[i]) * (_speaker_means(j,c)-_mean[j]);	

				for(unsigned long s=0;s<_session_per_speaker[c];s++){
					tmpSB(i,j)		+= (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
					sessionCounter++;
				}
			}
		}
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SW(i,j) = tmpSB(i,j) / _session_per_speaker[c];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct Scatterthread_data{

	double *data;
	double *speaker_means;
	double *mean;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *sb;
	RefVector <DoubleSquareMatrix> *sw;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *Scatterthread(void *threadarg){
	struct Scatterthread_data *my_data;
	my_data = (struct Scatterthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	double *mean = my_data->mean;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &_sb=(*(my_data->sb))[threadNb];
	DoubleSquareMatrix &_sw=(*(my_data->sw))[threadNb];

	double *sb = _sb.getArray();
	double *sw = _sw.getArray();

	unsigned long spk = _class[startSession];
	unsigned long s = startSession;

	DoubleSquareMatrix tmpSW;
	tmpSW.setSize(vectSize);
	tmpSW.setAllValues(0.0);

	unsigned long currentClass = _class[startSession];
	while(s<stopSession){
		while(currentClass == _class[s]){
			for(unsigned long i=0;i<vectSize;i++){
				for(unsigned long j=0;j<vectSize;j++){
					tmpSW(i,j)		+= (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
				}
			}
			s++;
		}
		for(unsigned long i=0;i<vectSize;i++){
			for(unsigned long j=0;j<vectSize;j++){
				sw[i*vectSize+j]		+= tmpSW(i,j) / session_per_speaker[currentClass];
				sb[i*vectSize+j]		+= (speaker_means[i*n_speakers+currentClass] - mean[i]) * (speaker_means[j*n_speakers+currentClass] - mean[j]);
			}
		}
		tmpSW.setAllValues(0.0);
		if(s<stopSession) currentClass = _class[s];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMatThreaded(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	SB.setSize(_vectSize); SB.setAllValues(0.0);
	SW.setSize(_vectSize); SW.setAllValues(0.0);

	// Compute the covariance matrices
	if(verboseLevel > 0) cout<<"	Compute Scatter matrices"<<endl;

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the 2 RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> sbMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		sbMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		sbMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> swMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		swMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		swMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct Scatterthread_data *thread_data_array = new Scatterthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data, *mean;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); mean = _mean.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t].mean = mean;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].sb=&(sbMatrices);
		thread_data_array[t].sw=&(swMatrices);

		if (verboseLevel > 0) cout<<"(computeCovMat) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, Scatterthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeCovMat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SB(i,j) += sbMatrices[t](i,j);
				SW(i,j) += swMatrices[t](i,j);
			}
		}
	}
	sbMatrices.deleteAllObjects();
	swMatrices.deleteAllObjects();
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::sphericalNuisanceNormalization(Config& config){
	
	// Chose the number of iterations
	unsigned long nb_iterations = 1;
	if(config.existsParam("ivNormIterationNb"))	nb_iterations = config.getParam("ivNormIterationNb").toULong();

	// Chose normalization mode (EFR or SphNorm)
	String mode = "EFR";
	if(config.existsParam("ivNormEfrMode"))	mode = config.getParam("ivNormEfrMode");
	
	//Prepare filename to save normalization matrices and means
	String sphNormMatrixBaseName = "ivNormEfrMatrix_it";
	if(config.existsParam("ivNormEfrMatrixBaseName"))	sphNormMatrixBaseName = config.getParam("ivNormEfrMatrixBaseName");
	
	String sphNormMeanBaseName = "ivNormEfrMean_it";
	if(config.existsParam("ivNormEfrMeanBaseName"))	sphNormMeanBaseName = config.getParam("ivNormEfrMeanBaseName");
	
	
	for(unsigned long it=0;it<nb_iterations;it++){
		
		if(verboseLevel>0) cout<<"	Eigen Factor Radial Normalization - start iteration "<<it<<endl;
	
		//Compute covariance matrices
		DoubleSquareMatrix W(_vectSize), B(_vectSize), Sigma(_vectSize);
		W.setAllValues(0.0); B.setAllValues(0.0); Sigma.setAllValues(0.0);
		this->computeCovMat(Sigma,W,B,config);
		
		Matrix<double> sphNormMat;
		sphNormMat.setDimensions (_vectSize,_vectSize)  ;
		sphNormMat.setAllValues(0.0); 
		
// AJOUTER UNE OPTION POUR RECUPERER LES VALEURS PROPRES ET POUVOIR TRACER
		
		if(mode == "EFR"){
			if(verboseLevel>0) cout<<"	Eigen Factor Radial "<<it<<endl;
			if(verboseLevel>0) cout<<"		compute whitening matrix"<<endl;
			
			//Compute EigenVectors and values of Sigma
			Matrix<double> EP(Sigma);

			Matrix<double> eigenVect, eigenVal,tmp;
			tmp.setDimensions(_vectSize,_vectSize);
			tmp.setAllValues(0.0);
			eigenVect.setDimensions(_vectSize,_vectSize);
			eigenVect.setAllValues(0.0);
			eigenVal.setDimensions(_vectSize,_vectSize);
			eigenVal.setAllValues(0.0);
			this->computeEigenProblem(EP,eigenVect,eigenVal,_vectSize,config);

			//Compute inverse square root of Sigma
			Matrix<double> invSqrtEigenValues(_vectSize,_vectSize);
			invSqrtEigenValues.setAllValues(0.0);
			for(unsigned long ii=0;ii<_vectSize;ii++){			//				sqSigma = bsxfun(@power, diag(Eval_sigma), 0.5);
				invSqrtEigenValues(ii,ii) = 1/sqrt(eigenVal(ii,0));			//				sqrInv_Eval_sigma = 1./sqSigma;

			}
			tmp =  eigenVect*invSqrtEigenValues;						//				sqrInvSigma = Evec_sigma*diag(sqrInv_Eval_sigma);	

			sphNormMat = tmp.transpose();
		}
		else{
			if(verboseLevel>0) cout<<"	Spherical Nuisance Normalization "<<it<<endl;
			if(verboseLevel>0) cout<<"		compute normalization matrix"<<endl;
			
			//Compute EigenVectors and values of W
			Matrix<double> EP(W);
			Matrix<double> eigenVect, eigenVal,tmp;
			tmp.setDimensions(_vectSize,_vectSize);
			tmp.setAllValues(0.0);
			eigenVect.setDimensions(_vectSize,_vectSize);
			eigenVect.setAllValues(0.0);
			eigenVal.setDimensions(_vectSize,_vectSize);
			eigenVal.setAllValues(0.0);
			this->computeEigenProblem(EP,eigenVect,eigenVal,_vectSize,config);

			//Compute inverse square root of Sigma
			Matrix<double> invSqrtEigenValues(_vectSize,_vectSize);
			invSqrtEigenValues.setAllValues(0.0);
			for(unsigned long ii=0;ii<_vectSize;ii++){			//				sqSigma = bsxfun(@power, diag(Eval_sigma), 0.5);
				invSqrtEigenValues(ii,ii) = 1/sqrt(eigenVal(ii,0));			//				sqrInv_Eval_sigma = 1./sqSigma;
			}

			tmp =  eigenVect*invSqrtEigenValues;						//				sqrInvSigma = Evec_sigma*diag(sqrInv_Eval_sigma);	
			sphNormMat = tmp.transpose();

		}

		//Save the  inverse square root of Sigma
		String s;
		String sphNormMatrixFilename = config.getParam("matrixFilesPath") + sphNormMatrixBaseName + s.valueOf(it) + config.getParam("saveMatrixFilesExtension");
		sphNormMat.save(sphNormMatrixFilename,config);
		
		//Save the mean of the iVector development set
		String sphNormMeanFilename = config.getParam("matrixFilesPath") + sphNormMeanBaseName + s.valueOf(it) + config.getParam("saveMatrixFilesExtension");
		Matrix<double> sphNormMean(_mean);
		sphNormMean.save(sphNormMeanFilename,config);
		
		//Center iVectors
		if(verboseLevel>0) cout<<"		center the data"<<endl;
		this->center(_mean);
	
		//Rotate using the inverse square root of Sigma
		if(verboseLevel>0) cout<<"		rotate the data"<<endl;
		this->rotateLeft(sphNormMat);
		
		//Normalize the length of iVectors
		if(verboseLevel>0) cout<<"		normalize the length of the data"<<endl;
		this->lengthNorm();
	}
}































//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//		PLDA TRAINING
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaModel::PldaModel(){
	_rankF		= 0;
	_rankG		= 0;
	_vectSize	= 0;

	// Set to 0 all other matrices used only for training
	_F.resize(0,0);
	_G.resize(0,0);
	_Sigma.resize(0,0);
	_originalMean.resize(0);
	_Delta.resize(0);

	_sigmaObs.resize(0,0);
	_Eh.resize(0,0);
	_EhhSum.resize(0,0);
	_xhSum.resize(0,0);
	_U.resize(0,0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaModel::PldaModel(String mode, Config &config){

	if(mode == "train"){
		//Read the NDX file
		String ndxFilename = config.getParam("pldaNdxFilename");

		PldaDev dev(ndxFilename, config);
		initTrain(dev, config);
	}
	else{
		initTest(config);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::initTrain(PldaDev dev,Config &config){

	 try {
		_rankF		= config.getParam("pldaEigenVoiceNumber").toULong();
		_rankG		= config.getParam("pldaEigenChannelNumber").toULong();

		//Initialize accumulators
		_Dev = dev;
		_vectSize	= _Dev.getVectSize();

		_Eh			= Eigen::MatrixXd::Zero(1,_Dev.getSpeakerNumber());
		_EhhSum		= Eigen::MatrixXd::Zero(_rankF+_rankG,_rankF+_rankG);
		_xhSum		= Eigen::MatrixXd::Zero(_vectSize,_rankF+_rankG);
		_U			= Eigen::MatrixXd::Zero(_rankF+_rankG,1);


		// Store Mean of the original data
		_originalMean	= Eigen::VectorXd::Zero(_vectSize);
		_sigmaObs		= Eigen::MatrixXd::Zero(_vectSize,_vectSize);
		_Delta			= Eigen::VectorXd::Zero(_vectSize);

		// Initialize temporary accumulators
		_invSigma.resize(_vectSize,_vectSize);
		_Ftweight.resize(_rankF, _vectSize);
		_Gtweight.resize(_rankG,_vectSize);
		_GtweightG.resize(_rankG,_rankG);
		_FtweightG.resize(_rankF, _rankG);
		_invGtweightGplusEye.resize(_rankG,_rankG);

		// Load matrices
		String tmpF,tmpG,tmpSigma,tmpMeanVec;

	// MODIFIER LES OPTIONS POUR SEPARER  eigenVoiceMatrix ET loadPldaMatrices (PAREIL POUR eigenChannel)

		//if EigenVoice matrix exists, load
		if(config.getParam("pldaLoadInitMatrices").toBool() && config.existsParam("pldaEigenVoiceMatrixInit")){
			Matrix<double> F;
			tmpF		= config.getParam("matrixFilesPath")+config.getParam("pldaEigenVoiceMatrixInit")+config.getParam("loadMatrixFilesExtension");
			if(verbose) cout<<"PldaModel load EigenVoiceMatrix from	[ "<<tmpF<<" ]"<<endl;
 			F.load(tmpF,config);
			_F.resize(_vectSize,_rankF);
			for(unsigned long i=0;i<F.rows();i++)
				for(unsigned long j=0;j<F.cols();j++)
					_F(i,j) = F(i,j);
		}
		else{		// if this matrix doesn't exists, then initialize randomly
			initF(config);
		}

		if(config.getParam("pldaLoadInitMatrices").toBool() && config.existsParam("eigenChannelMatrixInit")){
			Matrix<double> G;
			tmpG		= config.getParam("matrixFilesPath")+config.getParam("eigenChannelMatrixInit")+config.getParam("loadMatrixFilesExtension");
			if(verbose) cout<<"PldaModel load EigenChannelMatrix from	[ "<<tmpG<<" ]"<<endl;
			G.load(tmpG,config);
			_G.resize(_vectSize,_rankG);
			for(unsigned long i=0;i<G.rows();i++)
				for(unsigned long j=0;j<G.cols();j++)
					_G(i,j) = G(i,j);
		}
		else{		// if this matrix doesn't exists, then initialize randomly
			initG(config);
		}

		if(config.getParam("pldaLoadInitMatrices").toBool() && config.existsParam("sigmaMatrixInit")){
			tmpSigma	= config.getParam("matrixFilesPath")+config.getParam("sigmaMatrixInit")+config.getParam("loadMatrixFilesExtension");
			Matrix<double> TS(tmpSigma,config);
			_Sigma.resize(_vectSize,_vectSize);
			if(verbose) cout<<"PldaModel load Precision Matrix from	[ "<<tmpSigma<<" ]"<<endl;
			if((TS.rows() == _vectSize)&&(TS.cols()==_vectSize)){
				for(unsigned long i=0;i<_vectSize;i++)
					for(unsigned long j=0;j<_vectSize;j++)
						_Sigma(i,j) = TS(i,j);
			}
		}
		else{		// Initialize _Sigma to identity matrix
			_Sigma = Eigen::MatrixXd::Identity(_vectSize,_vectSize);
		}

		if(config.getParam("pldaLoadInitMatrices").toBool() && config.existsParam("meanVecInit")){
			tmpMeanVec	= config.getParam("matrixFilesPath")+config.getParam("meanVecInit")+config.getParam("loadMatrixFilesExtension");
			Matrix<double> tM(tmpMeanVec,config);
			_originalMean.resize(_vectSize);
			if(verbose) cout<<"PldaModel load Mean Vector from	[ "<<tmpMeanVec<<" ]"<<endl;
			for(unsigned long i=0;i<_vectSize;i++)
				_originalMean(i) = tM(i,0);
		}
		else{		// Initialize _originalMean to mean of pldaDev data
			RealVector<double> tM;
			tM = _Dev.getMean();
			_originalMean.resize(_vectSize);
			for(unsigned long i=0;i<_vectSize;i++)
				_originalMean(i) = tM[i];
		}
	}	// end try
	catch (Exception& e) {cout << e.toString().c_str() << endl;}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::initTest(Config &config){

	// Get dimensions from the config
	_rankF		= config.getParam("pldaEigenVoiceNumber").toULong();
	_rankG		= config.getParam("pldaEigenChannelNumber").toULong();
	_vectSize	= config.getParam("iVectSize").toULong();

	String pldaOriginalMeanMatrix = "pldaOriginalMeanMatrix";
	String pldaEigenVoiceMatrix	= "pldaEigenVoiceMatrix";
	String pldaEigenChannelMatrix = "pldaEigenChannelMatrix";
	String pldaPrecisionMatrix	= "pldaPrecisionMatrix";
	String pldaMinDivMeanMatrix	= "pldaMinDivMeanMatrix";
	if(config.existsParam("pldaOriginalMeanMatrix")) pldaOriginalMeanMatrix	= config.getParam("pldaOriginalMeanMatrix");
	if(config.existsParam("pldaEigenVoiceMatrix")) pldaEigenVoiceMatrix		= config.getParam("pldaEigenVoiceMatrix");
	if(config.existsParam("pldaEigenChannelMatrix")) pldaEigenChannelMatrix = config.getParam("pldaEigenChannelMatrix");
	if(config.existsParam("pldaPrecisionMatrix")) pldaPrecisionMatrix		= config.getParam("pldaPrecisionMatrix");
	if(config.existsParam("pldaMinDivMeanMatrix")) pldaMinDivMeanMatrix		= config.getParam("pldaMinDivMeanMatrix");

	// Load matrices
	loadF(pldaEigenVoiceMatrix,config);
	if(config.getParam("pldaEigenChannelNumber").toULong() == 0 ){
		_G.resize(_vectSize,0);
		if(verboseLevel>0) cout<<" (PldaModel) No EigenChannel";
	}
	else{	loadG(pldaEigenChannelMatrix,config);}
	loadNoise(pldaPrecisionMatrix,config);
	loadOriginalMean(pldaOriginalMeanMatrix,config);
	loadDelta(pldaMinDivMeanMatrix,config);

	// Set to 0 all other matrices used only for training
	_sigmaObs.resize(0,0);
	_Eh.resize(0,0);
	_EhhSum.resize(0,0);
	_xhSum.resize(0,0);
	_U.resize(0,0);

	// Initialize temporary accumulators
	_invSigma.resize(_vectSize,_vectSize);
	_Ftweight.resize(_rankF, _vectSize);
	_Gtweight.resize(_rankG,_vectSize);
	_GtweightG.resize(_rankG,_rankG);
	_FtweightG.resize(_rankF, _rankG);
	_invGtweightGplusEye.resize(_rankG,_rankG);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaModel::~PldaModel(){

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaModel::getClassName() const{	return "PldaModel";}

////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//void PldaModel::splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex){
//
//	unsigned long spkPerThread = _Dev.getSpeakerNumber()/nbThread;
//	startIndex.setSize(nbThread);
//	startIndex[0] = 0;
//	
//	unsigned long firstSpeakerNextList = spkPerThread;
//	unsigned long spkCounter = 1;
//	unsigned long currentSession = 0;
//
//	while((currentSession<getSessionNumber)&&(spkCounter<nbThread)){
//		while(_class[currentSession]<spkCounter*spkPerThread){
//			currentSession++;
//		}
//		startIndex[spkCounter] = currentSession;
//		spkCounter++;
//	}
//}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::initF(Config& config){		///random initialisation of the EigenVoices Matrix
		_rankF=config.getParam("pldaEigenVoiceNumber").toULong();
		_F.resize(_vectSize,_rankF);

		//Two type of random initializations
		//	normal  : random with a normal law by using a Box-Muller Generator
		//	uniform : random with a uniform distribution
try{
		String randomLaw = "normal";
		if(config.existsParam("pldaRandomInitLaw"))	randomLaw = config.getParam("pldaRandomInitLaw");

		//Initialize the matrix by generating random values following a uniform law
		if(randomLaw == "uniform"){
			_F.Random();
		}

		//Initialize the matrix by generating random values following a normal law whose mean is 0 and variance is 1
		else if (randomLaw == "normal"){

			boxMullerGeneratorInit();
			for(unsigned long i=0; i<_vectSize; i++){
				for(unsigned long j=0; j<_rankF; j++){
					double val = boxMullerGenerator(0.0, 1.0);
					while((ISNAN(val)) || (ISINF(val))){
						val = boxMullerGenerator(0.0, 1.0);
					}
					_F(i,j) = val;
				}
			}
		}
		else{
			throw Exception("Selected random initialization law does not exist",__FILE__,__LINE__);
		}
		if (verboseLevel >=1) cout << "(PldaModel) Random Init for EigenVoices Matrix with "<<randomLaw<<" law "<<", rank: ["<<_F.rows() << "] sv size: [" << _F.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::initG(Config& config){		///random initialisation of the EigenVoices Matrix
		_rankG=config.getParam("pldaEigenChannelNumber").toULong();
		_G.resize(_vectSize,_rankG);

		//Two type of random initializations
		//	normal  : random with a normal law by using a Box-Muller Generator
		//	uniform : random with a uniform distribution
try{
		String randomLaw = "normal";
		if(config.existsParam("pldaRandomInitLaw"))	randomLaw = config.getParam("pldaRandomInitLaw");

		//Initialize the matrix by generating random values following a uniform law
		if(randomLaw == "uniform"){

			srand48(_vectSize*_rankG);
			_G.Random();
		}

		//Initialize the matrix by generating random values following a normal law whose mean is 0 and variance is 1
		else if (randomLaw == "normal"){

			boxMullerGeneratorInit();
			for(int i=0; i<_G.rows(); i++){
				for(int j=0; j<_G.cols(); j++){
					double val = boxMullerGenerator(0.0, 1.0);
					while((ISNAN(val)) || (ISINF(val))){
						val = boxMullerGenerator(0.0, 1.0);
					}
					_G(i,j) = val;
				}
			}
		}
		else{
			throw Exception("Selected random initialization law does not exist",__FILE__,__LINE__);
		}
		if (verboseLevel >=1) cout << "(PldaModel) Random Init for EigenChannel Matrix with "<<randomLaw<<" law "<<", rank: ["<<_G.rows() << "] sv size: [" << _G.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::updateMean(){

	for(unsigned long i=0;i<_vectSize;i++)
		_originalMean(i) = _Dev.getMean()[i];
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::centerData(){
	_Dev.center(_originalMean);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::updateModel(Config& config){

	//Initialize accumulators
	_vectSize	= _Dev.getVectSize();
	_xhSum		= Eigen::MatrixXd::Zero(_vectSize,_rankF+_rankG);

	// Store Mean of the original data
	_originalMean	= Eigen::VectorXd::Zero(_vectSize);
	_sigmaObs		= Eigen::MatrixXd::Zero(_vectSize,_vectSize);
	_Delta			= Eigen::VectorXd::Zero(_vectSize);

	// Initialize temporary accumulators
	_invSigma.resize(_vectSize,_vectSize);
	_Ftweight.resize(_rankF, _vectSize);
	_Gtweight.resize(_rankG,_vectSize);

	// Load matrices
	String tmpF,tmpG,tmpSigma,tmpMeanVec;

	//if EigenVoice matrix exists, load
	if(config.existsParam("eigenVoiceMatrixInit")){
		Matrix<double> F;
		tmpF		= config.getParam("matrixFilesPath")+config.getParam("eigenVoiceMatrixInit")+config.getParam("loadMatrixFilesExtension");
		if(verbose) cout<<"PldaModel load EigenVoiceMatrix from	[ "<<tmpF<<" ]"<<endl;
 		F.load(tmpF,config);
		_F.resize(_vectSize,_rankF);
		for(unsigned long i=0;i<F.rows();i++)
			for(unsigned long j=0;j<F.cols();j++)
				_F(i,j) = F(i,j);
	}
	else{		// if this matrix doesn't exists, then initialize randomly
		initF(config);
	}

	if(config.existsParam("eigenChannelMatrixInit")){
		Matrix<double> G;
		tmpG		= config.getParam("matrixFilesPath")+config.getParam("eigenChannelMatrixInit")+config.getParam("loadMatrixFilesExtension");
		if(verbose) cout<<"PldaModel load EigenChannelMatrix from	[ "<<tmpG<<" ]"<<endl;
		G.load(tmpG,config);
		_G.resize(_vectSize,_rankG);
		for(unsigned long i=0;i<G.rows();i++)
			for(unsigned long j=0;j<G.cols();j++)
				_G(i,j) = G(i,j);
	}
	else{		// if this matrix doesn't exists, then initialize randomly
		initG(config);
	}

	if(config.existsParam("sigmaMatrixInit")){
		tmpSigma	= config.getParam("matrixFilesPath")+config.getParam("sigmaMatrixInit")+config.getParam("loadMatrixFilesExtension");
		Matrix<double> TS(tmpSigma,config);
		_Sigma.resize(_vectSize,_vectSize);
		if(verbose) cout<<"PldaModel load Precision Matrix from	[ "<<tmpSigma<<" ]"<<endl;
		if((TS.rows() == _vectSize)&&(TS.cols()==_vectSize)){
			for(unsigned long i=0;i<_vectSize;i++)
				for(unsigned long j=0;j<_vectSize;j++)
					_Sigma(i,j) = TS(i,j);
		}
	}
	else{		// Initialize _Sigma to identity matrix
		_Sigma = Eigen::MatrixXd::Identity(_vectSize,_vectSize);
	}

	if(config.existsParam("meanVecInit")){
		tmpMeanVec	= config.getParam("matrixFilesPath")+config.getParam("meanVecInit")+config.getParam("loadMatrixFilesExtension");
		Matrix<double> tM(tmpMeanVec,config);
		_originalMean.resize(_vectSize);
		if(verbose) cout<<"PldaModel load Mean Vector from	[ "<<tmpMeanVec<<" ]"<<endl;
		for(unsigned long i=0;i<_vectSize;i++)
			_originalMean(i) = tM(i,0);
	}
	else{		// Initialize _originalMean to mean of pldaDev data
		RealVector<double> tM;
		tM = _Dev.getMean();
		_originalMean.resize(_vectSize);
		for(unsigned long i=0;i<_vectSize;i++)
			_originalMean(i) = tM[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::em_iteration(Config &config){

	// Remove the mean re-estimated by minimum divergence
	if(verboseLevel>1) cout<<" (PldaModel) Center data using minimum divergence estimate"<<endl;
	_Dev.center(_Delta);

	// Compute Co-Variance matrix on the data
	_Dev.computeCovMatEigen(_sigmaObs,config);

	// Get statistics (E-step)
	this->getExpectedValues(config);

	// Re-estimate _F, _G, _Sigma and _Delta (M-step)
	this->mStep();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::getExpectedValues(Config &config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0){
		getExpectedValuesThreaded(config.getParam("numThread").toULong(),config);
	}
	else	getExpectedValuesUnThreaded(config); 			//unthreaded version
	#else
		getExpectedValuesUnThreaded(config);			//accumute stats
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::getExpectedValuesUnThreaded(Config& config){

	// Create temporary matrices
	Eigen::MatrixXd FtweightF(_rankF,_rankF);
	Eigen::MatrixXd GtweightG(_rankG,_rankG);
	Eigen::MatrixXd inGtweightG(_rankG,_rankG);
	Eigen::MatrixXd S(_rankG,_rankF);
	Eigen::MatrixXd A(_rankF,_rankF);

	// Compute temporary variables
	this->preComputation();
	FtweightF = _Ftweight * _F;
	S = _invGtweightGplusEye * _FtweightG.transpose();
	A = FtweightF - _FtweightG*_invGtweightGplusEye*_FtweightG.transpose();

	// Compute Eigen Decomposition
	Eigen::EigenSolver<Eigen::MatrixXd> es(A);
	Eigen::MatrixXd V = es.eigenvectors().real();
	Eigen::MatrixXd D = es.eigenvalues().real().asDiagonal();

	//Init accumulators
	_EhhSum = Eigen::MatrixXd::Zero(_rankF+_rankG,_rankF+_rankG);
	_xhSum = Eigen::MatrixXd::Zero(_vectSize,_rankF+_rankG);
	_U = Eigen::MatrixXd::Zero(_rankF+_rankG,1);	// TO DO changer pour un vecteur

	Eigen::MatrixXd M, invJDIVt,MsT;
	unsigned long currentSessionNb = 0;

// IL y a surement des choses qu'on peut sortir de la boucle pour accelerer ?
// multi-threader cette boucle
	unsigned long sessionCounter = 0;		// counter used to know to which speaker belongs the current session
	for(unsigned long spk=0;spk<_Dev.getSpeakerNumber();spk++){

		// Initialize_Eh for the current speaker
		unsigned long spkSessionNumber = _Dev.getSpeakerSessionNumber(spk);
		Eigen::MatrixXd _Eh =  Eigen::MatrixXd::Zero(_rankF+_rankG,spkSessionNumber);		//initialize for each speaker

		// Compute M, the first part of ATISigA^-1 which depends on the number of session per speaker
		if(spkSessionNumber != currentSessionNb){

			currentSessionNb = spkSessionNumber;
			M = Eigen::MatrixXd::Zero(_rankF,_rankF);

			// JDIVt = inv(J*D+eye(N_F))*V';
			Eigen::MatrixXd JDI = (currentSessionNb*D+Eigen::MatrixXd::Identity(_rankF,_rankF));
			invJDIVt = JDI.inverse() * V.transpose();	// optimize: don't use inverse for a diagonal matrix

			// Compute  M = V * Inv_JD_I_VT;
			M = V*invJDIVt;
	
			//Compute MsT = M*s';
			MsT = M*S.transpose();
		}

		//Compute Sigma_x = inv_Sigma * x{cInd};
				// get data of current speaker: session start from sessionCounter until 
		unsigned long spkFirstSession = sessionCounter;
		unsigned long spkClass = _Dev.getClass(sessionCounter); 
		while((sessionCounter < _Dev.getSessionNumber())&&(_Dev.getClass(sessionCounter) == spkClass)){
			sessionCounter++;
		}
		unsigned long spkLastSession = sessionCounter-1;

		Eigen::MatrixXd Sigma_x = Eigen::MatrixXd::Zero(_vectSize,_Dev.getSpeakerSessionNumber(spk));
		for(unsigned long i=0;i<_vectSize;i++)
			for(unsigned long j=0;j<_Dev.getSpeakerSessionNumber(spk);j++)
				for(unsigned long k=0;k<_vectSize;k++)
					Sigma_x(i,j) += _invSigma(i,k) * _Dev.getData(k,spkFirstSession+j);

		//Compute f = sum(F'* Sigma_x,2);
		Eigen::MatrixXd fi	= _F.transpose()* Sigma_x;
		Eigen::VectorXd f	= Eigen::VectorXd::Zero(_rankF); 
		for(int r=0;r<fi.cols();r++)
			f += fi.col(r);

		//Compute gi = G' * Sigma_x;
		//Compute g = sum(gi,2);
		Eigen::MatrixXd gi	= _G.transpose()* Sigma_x;
		Eigen::VectorXd g	= Eigen::VectorXd::Zero(_rankG); 
		for(int r=0;r<gi.cols();r++)
			g += gi.col(r);

		//Compute thisEh = M *(f - s'*g);
		Eigen::MatrixXd thisEh = M * (f - S.transpose()*g);

		//Compute thisEw = bsxfun(@minus,c*gi, s*thisEh);
		Eigen::MatrixXd invGtweightGplusEyegi = _invGtweightGplusEye*gi;
		Eigen::VectorXd SthisEh = S*thisEh;
		for(int c=0;c<invGtweightGplusEyegi.cols();c++)
			invGtweightGplusEyegi.col(c) -= SthisEh;

		//Create Eh{cInd} = [repmat(thisEh,1,J); thisEw];
		Eigen::MatrixXd Eh(thisEh.rows()+invGtweightGplusEyegi.rows(),invGtweightGplusEyegi.cols());
		for(int c=0;c<invGtweightGplusEyegi.cols();c++)
			Eh.block(0,c,thisEh.rows(),1)=thisEh;
		for(int r = 0;r<invGtweightGplusEyegi.rows();r++)
			Eh.row(thisEh.rows()+r) = invGtweightGplusEyegi.row(r);


		//Compute EhhSum = EhhSum + J*[M, -MsT; -MsT', c + s*MsT] + Eh{cInd}*Eh{cInd}'; 
			// Create [M, -MsT; -MsT', c + s*MsT]
// TO DO modifier les dimensions de la matrice tmpM et Eh (ci-dessus) pour utiliser les valeurs constantes ou plus parlantes
		Eigen::MatrixXd tmpM(thisEh.rows()+invGtweightGplusEyegi.rows(),thisEh.rows()+invGtweightGplusEyegi.rows());
		tmpM.block(0,0,_rankF,_rankF)			= M;
		tmpM.block(0,_rankF,_rankF,_rankG)		= -MsT;
		tmpM.block(_rankF,0,_rankG,_rankF)		= -MsT.transpose();
		tmpM.block(_rankF,_rankF,_rankG,_rankG)	= _invGtweightGplusEye + S*MsT;

		_EhhSum +=  spkSessionNumber*tmpM + Eh*Eh.transpose();

		//Compute xhSum = xhSum + x{cInd}* Eh{cInd}';   % Calculate terms for updates
		for(unsigned long i=0;i<_vectSize;i++)
			for(unsigned long j=0;j<_rankF+_rankG;j++)
				for(unsigned long k=0;k<_Dev.getSpeakerSessionNumber(spk);k++)
					_xhSum(i,j) += _Dev.getData(i,spkFirstSession+k) * Eh(j,k);

		//Compute u = u + sum(Eh{cInd},2);
		for(int c=0;c<gi.cols();c++)
			_U += Eh.col(c);
	}	// end speaker loop
}

#ifdef THREAD
pthread_mutex_t mutexPLDA=PTHREAD_MUTEX_INITIALIZER;				// Mutex for PLDA
	
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct getExpectedValuesThread_data{

	double *v;
	double *d;
	double *s;
	double *invsigma;
	double *F;
	double *G;
	double * invGtweightGplusEye;

	unsigned long rankF;
	unsigned long rankG;
	unsigned long vectSize;

	PldaDev *dev;

	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long numThread;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *getExpectedValuesThread(void *threadarg) {

// Routine en cours--------------------------------------

	pthread_mutex_init(&mutexPLDA, NULL);		// Mutex for ??? tout ou juste un des accumulateurs ?
	
	unsigned long currentSessionNb = 0;		// a inclure dans chaque thread
	unsigned long sessionCounter = 0;		// counter used to know to which speaker belongs the current session

	struct getExpectedValuesThread_data *my_data;
	my_data = (struct getExpectedValuesThread_data *) threadarg;

	unsigned long rankF = my_data->rankF;
	unsigned long rankG = my_data->rankG;
	unsigned long vectSize = my_data->vectSize;

	Eigen::Map<Eigen::MatrixXd> v(my_data->v,rankF,rankF);
	Eigen::Map<Eigen::MatrixXd> d(my_data->d,rankF,rankF);
	Eigen::Map<Eigen::MatrixXd> s(my_data->s,rankG,rankF);
	Eigen::Map<Eigen::MatrixXd> invsigma(my_data->invsigma,vectSize,vectSize);
	Eigen::Map<Eigen::MatrixXd> F(my_data->F,vectSize,rankF);
	Eigen::Map<Eigen::MatrixXd> G(my_data->G,vectSize,rankG);
	Eigen::Map<Eigen::MatrixXd> invGtweightGplusEye(my_data->invGtweightGplusEye,rankG,rankG);

	Eigen::Map<Eigen::MatrixXd> EhhSum(my_data->invGtweightGplusEye,rankF+rankG,rankF+rankG);
	Eigen::Map<Eigen::MatrixXd> xhSum(my_data->invGtweightGplusEye,vectSize,rankF+rankG);
	Eigen::Map<Eigen::MatrixXd> U(my_data->invGtweightGplusEye,rankF+rankG,1);

	PldaDev *dev = my_data->dev;

	unsigned long spkBottom = my_data->spkBottom;
	unsigned long spkUp = my_data->spkUp;	
	unsigned long numThread = my_data->numThread;
	
	Eigen::MatrixXd M, invJDIVt,MsT;

	for(unsigned long spk=spkBottom;spk<spkUp;spk++){

			// Initialize_Eh for the current speaker
			unsigned long spkSessionNumber = dev->getSpeakerSessionNumber(spk);
			Eigen::MatrixXd _Eh =  Eigen::MatrixXd::Zero(rankF+rankG,spkSessionNumber);		//initialize for each speaker

			// Compute M, the first part of ATISigA^-1 which depends on the number of session per speaker
			if(spkSessionNumber != currentSessionNb){

				currentSessionNb = spkSessionNumber;
				M = Eigen::MatrixXd::Zero(rankF,rankF);

				// JDIVt = inv(J*D+eye(N_F))*V';
				Eigen::MatrixXd JDI = (currentSessionNb*d+Eigen::MatrixXd::Identity(rankF,rankF));
				invJDIVt = JDI.inverse() * v.transpose();	// optimize: don't use inverse for a diagonal matrix

				// Compute  M = V * Inv_JD_I_VT;
				M = v*invJDIVt;
	
				//Compute MsT = M*s';
				MsT = M*s.transpose();
			}

			//Compute Sigma_x = inv_Sigma * x{cInd}; 			// get data of current speaker: session start from sessionCounter until 
			unsigned long spkFirstSession = sessionCounter;
			unsigned long spkClass = dev->getClass(sessionCounter); 
			while((sessionCounter < dev->getSessionNumber())&&(dev->getClass(sessionCounter) == spkClass)){
				sessionCounter++;
			}
			unsigned long spkLastSession = sessionCounter-1;

			Eigen::MatrixXd Sigma_x = Eigen::MatrixXd::Zero(vectSize,dev->getSpeakerSessionNumber(spk));
			for(unsigned long i=0;i<vectSize;i++)
				for(unsigned long j=0;j<dev->getSpeakerSessionNumber(spk);j++)
					for(unsigned long k=0;k<vectSize;k++)
						Sigma_x(i,j) += invsigma(i,k) * dev->getData(k,spkFirstSession+j);

			//Compute f = sum(F'* Sigma_x,2);
			Eigen::MatrixXd fi	= F.transpose()* Sigma_x;
			Eigen::VectorXd f	= Eigen::VectorXd::Zero(rankF); 
			for(int r=0;r<fi.cols();r++)
				f += fi.col(r);

			//Compute gi = G' * Sigma_x;
			//Compute g = sum(gi,2);
			Eigen::MatrixXd gi	= G.transpose()* Sigma_x;
			Eigen::VectorXd g	= Eigen::VectorXd::Zero(rankG); 
			for(int r=0;r<gi.cols();r++)
				g += gi.col(r);

			//Compute thisEh = M *(f - s'*g);
			Eigen::MatrixXd thisEh = M * (f - s.transpose()*g);

			//Compute thisEw = bsxfun(@minus,c*gi, s*thisEh);
			Eigen::MatrixXd invGtweightGplusEyegi = invGtweightGplusEye*gi;
			Eigen::VectorXd SthisEh = s*thisEh;
			for(int c=0;c<invGtweightGplusEyegi.cols();c++)
				invGtweightGplusEyegi.col(c) -= SthisEh;

			//Create Eh{cInd} = [repmat(thisEh,1,J); thisEw];
			Eigen::MatrixXd Eh(thisEh.rows()+invGtweightGplusEyegi.rows(),invGtweightGplusEyegi.cols());
			for(int c=0;c<invGtweightGplusEyegi.cols();c++)
				Eh.block(0,c,thisEh.rows(),1)=thisEh;
			for(int r = 0;r<invGtweightGplusEyegi.rows();r++)
				Eh.row(thisEh.rows()+r) = invGtweightGplusEyegi.row(r);


			//Compute EhhSum = EhhSum + J*[M, -MsT; -MsT', c + s*MsT] + Eh{cInd}*Eh{cInd}'; 
				// Create [M, -MsT; -MsT', c + s*MsT]
	// TO DO modifier les dimensions de la matrice tmpM et Eh (ci-dessus) pour utiliser les valeurs constantes ou plus parlantes
			Eigen::MatrixXd tmpM(thisEh.rows()+invGtweightGplusEyegi.rows(),thisEh.rows()+invGtweightGplusEyegi.rows());
			tmpM.block(0,0,rankF,rankF)			= M;
			tmpM.block(0,rankF,rankF,rankG)		= -MsT;
			tmpM.block(rankF,0,rankG,rankF)		= -MsT.transpose();
			tmpM.block(rankF,rankF,rankG,rankG)	= invGtweightGplusEye + s*MsT;

			// Update Accumulators
			pthread_mutex_lock(&mutexPLDA);		// Lock Mutex

			EhhSum +=  spkSessionNumber*tmpM + Eh*Eh.transpose();

			//Compute xhSum = xhSum + x{cInd}* Eh{cInd}';   % Calculate terms for updates
			for(unsigned long i=0;i<vectSize;i++)
				for(unsigned long j=0;j<rankF+rankG;j++)
					for(unsigned long k=0;k<dev->getSpeakerSessionNumber(spk);k++)
						xhSum(i,j) += dev->getData(i,spkFirstSession+k) * Eh(j,k);

			//Compute u = u + sum(Eh{cInd},2);
			for(int c=0;c<gi.cols();c++)
				U += Eh.col(c);

			pthread_mutex_unlock(&mutexPLDA);	// Unlock Mutex
	}	// end speaker loop

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void PldaModel::getExpectedValuesThreaded(unsigned long NUM_THREADS, Config& config){

	try{
		if (verboseLevel >= 1) cout << "(PldaModel) Estimate Compute Expected Values for PLDA Threaded"<<endl;
		if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

		// Create temporary matrices
		Eigen::MatrixXd FtweightF(_rankF,_rankF);
		Eigen::MatrixXd GtweightG(_rankG,_rankG);
		Eigen::MatrixXd inGtweightG(_rankG,_rankG);
		Eigen::MatrixXd S(_rankG,_rankF);
		Eigen::MatrixXd A(_rankF,_rankF);

		// Compute temporary variables
		this->preComputation();
		FtweightF = _Ftweight * _F;
		S = _invGtweightGplusEye * _FtweightG.transpose();
		A = FtweightF - _FtweightG*_invGtweightGplusEye*_FtweightG.transpose();

		// Compute Eigen Decomposition
		Eigen::EigenSolver<Eigen::MatrixXd> es(A);
		Eigen::MatrixXd V = es.eigenvectors().real();
		Eigen::MatrixXd D = es.eigenvalues().real().asDiagonal();

		//Init accumulators
		_EhhSum = Eigen::MatrixXd::Zero(_rankF+_rankG,_rankF+_rankG);
		_xhSum = Eigen::MatrixXd::Zero(_vectSize,_rankF+_rankG);
		_U = Eigen::MatrixXd::Zero(_rankF+_rankG,1);	// TO DO changer pour un vecteur
		
		int rc, status;
		if (NUM_THREADS > _Dev.getSpeakerNumber()) NUM_THREADS=_Dev.getSpeakerNumber();

		// split the list
		RealVector<unsigned long> startIndex;
		_Dev.splitPerSpeaker(NUM_THREADS, startIndex);

		struct getExpectedValuesThread_data *thread_data_array = new getExpectedValuesThread_data[NUM_THREADS];
		pthread_t *threads = new pthread_t[NUM_THREADS];

		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		//Create threads
		for(unsigned long t=0; t<NUM_THREADS; t++){

			double *v = V.data();	// a remplacer par l'equivalent dans Eigen
			double *d = D.data();
			double *s = S.data();
			double *invsigma = _invSigma.data();
			double *F = _F.data();
			double *G = _G.data();
			double * invGtweightGplusEye = _invGtweightGplusEye.data();
			
			//PldaDev &dev = _Desv;

			unsigned long spkBottom		= startIndex[t];
			unsigned long spkUp = 0;
			if(t<NUM_THREADS-1)	spkUp	= startIndex[t+1];
			else{				spkUp	= _Dev.getSpeakerNumber();}
			unsigned long numThread = t;

			thread_data_array[t].v = v;
			thread_data_array[t].d = d;
			thread_data_array[t].s = s;
			thread_data_array[t].invsigma = invsigma;
			thread_data_array[t].F = F;
			thread_data_array[t].G = G;
			thread_data_array[t].invGtweightGplusEye = invGtweightGplusEye;
			thread_data_array[t].rankF = _rankF;
			thread_data_array[t].rankG = _rankG;
			thread_data_array[t].vectSize = _vectSize;
			thread_data_array[t].dev = &_Dev;
			thread_data_array[t].spkBottom = spkBottom;
			thread_data_array[t].spkUp = spkUp;	
			thread_data_array[t].numThread = numThread;

			if (verboseLevel>1) cout<<"(PldaModel) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
			rc = pthread_create(&threads[t], &attr, getExpectedValuesThread, (void *)&thread_data_array[t]);
			if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		}
	
		pthread_attr_destroy(&attr);
		for(unsigned long t=0; t<NUM_THREADS; t++) {
			rc = pthread_join(threads[t], (void **)&status);
			if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
			if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
		}
		
		free(thread_data_array);
		free(threads);
	
		if (verboseLevel >= 1) cout << "(PldaModel) Done " << endl;
	}
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::mStep(){

	// Update F and G 
	Eigen::MatrixXd FGEst = _xhSum * _EhhSum.inverse();
	_F = FGEst.block(0,0,_vectSize,_rankF);
	_G = FGEst.block(0,_rankF,_vectSize,_rankG);

	// Update Sigma
	Eigen::MatrixXd SigmaLat = FGEst * _xhSum.transpose();
	_Sigma = (_sigmaObs - SigmaLat)/_Dev.getSessionNumber();

	// Minimum Divergence
	_U = _U/_Dev.getSessionNumber();
	Eigen::MatrixXd c = _EhhSum/_Dev.getSessionNumber() - _U*_U.transpose();
	
	Eigen::MatrixXd Rh = c.block(0,0,_rankF,_rankF).llt().matrixL().transpose();	//upper triangular matrix from the Cholesky decomposition
	Eigen::MatrixXd Rw = c.block(_rankF,_rankF,_rankG,_rankG).llt().matrixL().transpose();
	_F *= Rh.transpose();
	_G *= Rw.transpose();
	_Delta += FGEst * _U;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::saveModel(Config &config){

	// Convert matrices from Eigen to ALIZE format
	Matrix<double> originalMean(_vectSize,1);
	Matrix<double> F(_vectSize,_rankF);
	Matrix<double> G(_vectSize,_rankG);
	Matrix<double> Sigma(_vectSize,_vectSize);
	Matrix<double> minDivMean(_vectSize,1);

	for(unsigned long i=0;i<_vectSize;i++){

		originalMean(i,0)	= _originalMean(i);
		minDivMean(i,0)		= _Delta(i);

		for(unsigned long j=0;j<_rankF;j++){
			F(i,j) = _F(i,j);
		}

		for(unsigned long j=0;j<_rankG;j++){
			G(i,j) = _G(i,j);
		}

		for(unsigned long j=0;j<_vectSize;j++){
			Sigma(i,j) = _Sigma(i,j);
		}
	}

	String pldaOriginalMeanMatrix = "pldaOriginalMeanMatrix";
	String pldaEigenVoiceMatrix	= "pldaEigenVoiceMatrix";
	String pldaEigenChannelMatrix = "pldaEigenChannelMatrix";
	String pldaPrecisionMatrix	= "pldaPrecisionMatrix";
	String pldaMinDivMeanMatrix	= "pldaMinDivMeanMatrix";
	if(config.existsParam("pldaOriginalMeanMatrix")) pldaOriginalMeanMatrix	= config.getParam("pldaOriginalMeanMatrix");
	if(config.existsParam("pldaEigenVoiceMatrix")) pldaEigenVoiceMatrix		= config.getParam("pldaEigenVoiceMatrix");
	if(config.existsParam("pldaEigenChannelMatrix")) pldaEigenChannelMatrix = config.getParam("pldaEigenChannelMatrix");
	if(config.existsParam("pldaPrecisionMatrix")) pldaPrecisionMatrix		= config.getParam("pldaPrecisionMatrix");
	if(config.existsParam("pldaMinDivMeanMatrix")) pldaMinDivMeanMatrix		= config.getParam("pldaMinDivMeanMatrix");

	// Save Original Mean _originalMean
	String originalMeanFilename = config.getParam("matrixFilesPath") + pldaOriginalMeanMatrix + config.getParam("saveMatrixFilesExtension");
	originalMean.save(originalMeanFilename, config);

	// Save EigenVoice matrix _F
	String eigenVoiceFilename = config.getParam("matrixFilesPath") + pldaEigenVoiceMatrix + config.getParam("saveMatrixFilesExtension");
	F.save(eigenVoiceFilename, config);

	// Save EigenChannel matrix _G
	String eigenChannelFilename = config.getParam("matrixFilesPath") + pldaEigenChannelMatrix + config.getParam("saveMatrixFilesExtension");
	G.save(eigenChannelFilename, config);

	// Save Precision matrix _Sigma
	String precisionFilename = config.getParam("matrixFilesPath") + pldaPrecisionMatrix + config.getParam("saveMatrixFilesExtension");
	Sigma.save(precisionFilename, config);

	// Save minimum divergence mean estimation _Delta
	String minDivMeanFilename = config.getParam("matrixFilesPath") + pldaMinDivMeanMatrix + config.getParam("saveMatrixFilesExtension");
	minDivMean.save(minDivMeanFilename, config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::loadF(String filename, Config& config){
//TO DO: exception for dimension mismatch

	String f = config.getParam("matrixFilesPath") + filename + config.getParam("loadMatrixFilesExtension");
	Matrix<double> F; F.load (f,config);
	_rankF = F.cols();
	_F.resize(F.rows(),F.cols());

	for(unsigned long i=0;i<F.rows();i++)
		for(unsigned long j=0;j<F.cols();j++)
			_F(i,j) = F(i,j);

	if(verboseLevel>0) cout<<" (PldaModel) Load EigenVoice Matrix, rank = "<<_rankF<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::loadG(String filename, Config& config){
//TO DO: exception for dimension mismatch, if G is 0 then set dimension to _vectSize x 0

	String f = config.getParam("matrixFilesPath") + filename + config.getParam("loadMatrixFilesExtension");

	Matrix<double> G; G.load (f,config);
	_rankG = G.cols();
	_G.resize(G.rows(),G.cols());
	if(G.rows() ==0){
		_G.resize(_vectSize,G.cols());
	}


	for(unsigned long i=0;i<G.rows();i++)
		for(unsigned long j=0;j<G.cols();j++)
			_G(i,j) = G(i,j);

	if(verboseLevel>0) cout<<" (PldaModel) Load EigenChannel Matrix, rank = "<<_rankG<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::loadNoise(String filename, Config& config){

	String f = config.getParam("matrixFilesPath") + filename + config.getParam("loadMatrixFilesExtension");

	Matrix<double> tmpNoise;
	tmpNoise.load (f,config);
	_Sigma.resize(tmpNoise.rows(),tmpNoise.rows());
	for(unsigned long i=0;i<tmpNoise.rows();i++)
		for(unsigned long j=0;j<tmpNoise.rows();j++)
			_Sigma(i,j) = tmpNoise(i,j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::loadOriginalMean(String filename, Config& config){

	String f = config.getParam("matrixFilesPath") + filename + config.getParam("loadMatrixFilesExtension");

	Matrix<double> tmp;
	tmp.load (f,config);
	_originalMean.resize(tmp.rows());
	for(unsigned long i=0;i<tmp.rows();i++)
		_originalMean[i] = tmp(i,0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::loadDelta(String filename, Config& config){

	String f = config.getParam("matrixFilesPath") + filename + config.getParam("loadMatrixFilesExtension");

	Matrix<double> tmp;
	tmp.load (f,config);
	_Delta.resize(tmp.rows());
	for(unsigned long i=0;i<tmp.rows();i++)
		_Delta[i] = tmp(i,0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaModel::preComputation(){

	// Inverse Sigma
	_invSigma = _Sigma.inverse();

	// Compute _Ftweight(_rankF, _vectSize)
	_Ftweight = _F.transpose()*_invSigma;

	// Compute _Gtweight(_rankG,_vectSize)
	_Gtweight = _G.transpose()*_invSigma;

	// Compute _Gtweight(_rankG,_rankG)
	_GtweightG = _Gtweight*_G;

	// Compute _FtweightG(_rankF, _rankG)
	_FtweightG = _Ftweight * _G;

	// Compute invGtweightGplusEye
	Eigen::MatrixXd GtweightGplusEye(_rankG,_rankG);
	GtweightGplusEye = _GtweightG + Eigen::MatrixXd::Identity(_rankG,_rankG);
	_invGtweightGplusEye = GtweightGplusEye.inverse();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd PldaModel::getFtweight(){
	return(_Ftweight);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd PldaModel::getFtweightG(){
	return(_FtweightG);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd PldaModel::getInvGtweightGplusEye(){
	return(_invGtweightGplusEye);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd PldaModel::getGtweight(){
	return(_Gtweight);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd PldaModel::getF(){
	return(_F);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaModel::getRankF(){
	return(_rankF);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaModel::getRankG(){
	return(_rankG);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev& PldaModel::getDev(){
	return(_Dev);
}
















//***********************************************
//	PldaTest methods
//***********************************************

// Create the PldaTest object and load model and test segments.
// Case where each user get several enrollment segments or that the model name is the different from the segment name
//
PldaTest::PldaTest(String &ndxTrialsFilename, String &ndxModelId ,Config &config){

	_segLine.reset();
	_modelIDLine.reset();
	_modelIndexLine.reset();
	_modelSessionslLine.reset();

	// Create the XList for model definition
	XList enrollList;
	enrollList.load(ndxModelId,config);
	enrollList.sortByElementNumber("descend");

	// Read targetIdList and fill _modelIDLine and _modelSessionslLine
	enrollList.getLine(0);
	String model, enrollSession;
	XLine *linep;
	linep = enrollList.getLine();

	// Get the maximum number of sessions per speaker
	_n_enrollSessions_max = linep->getElementCount();

	while ((linep=enrollList.getLine()) != NULL){
		
		model = linep->getElement(0);
		if(_modelIDLine.getIndex(model) == -1){	//add the model ID to the list if it doesn't exists
			_modelIDLine.addElement(model);
		}

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			enrollSession = linep->getElement(e);
			_modelIndexLine.addElement(model);
			_modelSessionslLine.addElement(enrollSession);
			e++;
		}
	}

	// Create the XList for trials
	XList testList;
	testList.load(ndxTrialsFilename,config);

	// Read the XList and fill the XLines when required
	testList.getLine(0);
	String vect;

	while ((linep=testList.getLine()) != NULL){
		
		vect = linep->getElement(0);
		if(_segLine.getIndex(vect) == -1){
			_segLine.addElement(vect);
		}

		// Verify that all models involved in the testing have at least one enrollment session
		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			if(_modelIDLine.getIndex(vect) == -1){
				cout<<"	WARNING: model "<<vect<<"	has no enrollment session"<<endl;
				//_modelIDLine.addElement(vect);
			}
			e++;
		}
	}
	_n_models				= _modelIDLine.getElementCount();
	_n_enrollment_segments	= _modelSessionslLine.getElementCount();
	_n_test_segments		= _segLine.getElementCount();

	// Get _vectSize from the first model
	String tmp = _modelIDLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("vectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_enrollment_segments);
	_segments.setDimensions(_vectSize,_n_test_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<" (PldaTest) found: "<<_n_models<<" models and "<<_n_test_segments<<" test segments"<<endl;
		cout<<" (PldaTest) Maximum number of sessions per speaker = "<<_n_enrollSessions_max<<endl;
	}

	// Load model enrollment vectors
	String fileName;
	for(unsigned long m=0;m<_n_models;m++){

		// Read the vector from file
		//fileName = _modelIDLine.getElement(m);
		fileName = _modelSessionslLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_test_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++)
			_segments(k,s) = tmpVect(0,k);
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_test_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	testList.getLine(0);
	while ((linep=testList.getLine()) != NULL){
	
		vect = linep->getElement(0);
		unsigned long segIdx = (unsigned long)_segLine.getIndex(vect);

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			unsigned long modelIdx = (unsigned long)_modelIDLine.getIndex(vect);
			_trials(modelIdx,segIdx) = true; 
			e++;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Create the PldaTest object and load model and test segments.
// Case where each user get one only enrollment segment and that the model name is the same as the segment name
//
PldaTest::PldaTest(String &ndxTrialsFilename ,Config &config){

	bool existIdList = false;
	String targetIdList;
	if(config.existsParam("targetIdList")){
		existIdList = true;
		targetIdList = config.getParam("targetIdList");
	}

	_segLine.reset();
	_modelIDLine.reset();
	_modelIndexLine.reset();
	_modelSessionslLine.reset();

	XLine *linep;

	if(existIdList){

		// Create the XList for model definition
		XList enrollList;
		enrollList.load(targetIdList,config);
		enrollList.sortByElementNumber("descend");

		// Read targetIdList and fill _modelIDLine and _modelSessionslLine
		enrollList.getLine(0);
		linep=enrollList.getLine();
		String model, enrollSession;

		// Get the maximum number of sessions per speaker
		_n_enrollSessions_max = linep->getElementCount()-1;
		enrollList.getLine(0);
		while ((linep=enrollList.getLine()) != NULL){
		
			model = linep->getElement(0);
			if(_modelIDLine.getIndex(model) == -1){	//add the model ID to the list if it doesn't exists
				_modelIDLine.addElement(model);
			}

			unsigned long e = 1;
			while(e<linep->getElementCount()){
				enrollSession = linep->getElement(e);
				_modelIndexLine.addElement(model);
				_modelSessionslLine.addElement(enrollSession);
				e++;
			}
		}
	}
	else{
		// Get the maximum number of sessions per speaker
		_n_enrollSessions_max = 1;
	}

	// Create the XList for trials
	XList testList;
	testList.load(ndxTrialsFilename,config);

	// Read the XList and fill the XLines when required
	testList.getLine(0);
	String vect;
	while ((linep=testList.getLine()) != NULL){

		vect = linep->getElement(0);
		if(_segLine.getIndex(vect) == -1){
			_segLine.addElement(vect);
		}

		// Verify that all models involved in the testing have at least one enrollment session
		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			if(_modelIDLine.getIndex(vect) == -1){
				if(existIdList){
					cout<<"	WARNING: model "<<vect<<"	has no enrollment session"<<endl;
				}
				else{
					_modelIDLine.addElement(vect);
					_modelIndexLine.addElement(vect);
					_modelSessionslLine.addElement(vect);
				}
			}
			e++;
		}
	}
	_n_models				= _modelIDLine.getElementCount();
	_n_enrollment_segments	= _modelSessionslLine.getElementCount();
	_n_test_segments		= _segLine.getElementCount();

	// Get _vectSize from the first model
	String tmp = _modelSessionslLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("vectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_enrollment_segments);
	_segments.setDimensions(_vectSize,_n_test_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<" (PldaTest) found: "<<_n_models<<" models and "<<_n_test_segments<<" test segments"<<endl;
		cout<<" (PldaTest) Maximum number of sessions per speaker = "<<_n_enrollSessions_max<<endl;
	}

	// Load model enrollment vectors
	String fileName;
	for(unsigned long m=0;m<_n_enrollment_segments;m++){
		// Read the vector from file
		fileName = _modelSessionslLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_test_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) = tmpVect(0,k);
		}
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_test_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	testList.getLine(0);
	while ((linep=testList.getLine()) != NULL){
	
		vect = linep->getElement(0);
		unsigned long segIdx = (unsigned long)_segLine.getIndex(vect);

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			unsigned long modelIdx = (unsigned long)_modelIDLine.getIndex(vect);
			_trials(modelIdx,segIdx) = true; 
			e++;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::PldaTest(XList &testList,Config & config){

	_segLine.reset();
	_modelIDLine.reset();
	_modelIndexLine.reset();
	_modelSessionslLine.reset();

	// Read the XList and fill the XLines when required
	XLine *linep;
	testList.getLine(0);
	String vect;

	while ((linep=testList.getLine()) != NULL){
		
		vect = linep->getElement(0);
		if(_segLine.getIndex(vect) == -1){
			_segLine.addElement(vect);
		}

		unsigned long e = 1;
		while(e<linep->getElementCount()){

			vect = linep->getElement(e);
			if(_modelIDLine.getIndex(vect) == -1){
				_modelIDLine.addElement(vect);
			}
			e++;
		}
	}
	_n_models	= _modelIDLine.getElementCount();
	_n_test_segments = _segLine.getElementCount();

	// Get _vectSize from the first model
	String tmp = _modelIDLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("vectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_models);
	_segments.setDimensions(_vectSize,_n_test_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<"PldaTest found: "<<_n_models<<" models and "<<_n_test_segments<<" test segments"<<endl;
	}

	// Load model vectors
	String fileName;
	for(unsigned long m=0;m<_n_models;m++){

		// Read the vector from file
		fileName = _modelIDLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_test_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) = tmpVect(0,k);
		}
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_test_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_scores.setDimensions(_n_models,_n_test_segments);
	_scores.setAllValues(0);

	testList.getLine(0);
	while ((linep=testList.getLine()) != NULL){
	
		vect = linep->getElement(0);
		unsigned long seg = (unsigned long)_segLine.getIndex(vect);

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			unsigned long model = (unsigned long)_modelIDLine.getIndex(vect);
			_trials(model,seg) = 1; 
			e++;
		}
	}
} // TO MODIFY TO DEAL WITH MULTIPLE SEGMENTS

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::PldaTest(Config& config){
	this->load(config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::PldaTest(){}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::~PldaTest(){}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getClassName() const{	return "PldaTest";}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::load(Config &config){

	bool existIdList = false;

	// Create the XList for trials
	XList testList;
	XList enrollList;
	_segLine.reset();
	_modelIDLine.reset();
	_modelIndexLine.reset();
	_modelSessionslLine.reset();

	//Load from one only file list
	if(config.existsParam("inputVectorFilename")){

		config.setParam("testVectorFilesPath",config.getParam("loadVectorFilesPath"));

		String ndxTrialsFilename = config.getParam("inputVectorFilename");
		testList.load(ndxTrialsFilename,config);

		//Initialize enrollList with all files (by copying testList data)
		_modelIDLine = testList.getAllElements();
		_modelIndexLine = testList.getAllElements();
		_modelSessionslLine = testList.getAllElements();
		_n_enrollSessions_max = 1;

		//Initialize the _segLine with all segments from testList
		_segLine = testList.getAllElements();

	}
	else{ //load existing testList 
		String ndxTrialsFilename = config.getParam("ndxFilename");
		testList.load(ndxTrialsFilename,config);
	}

	String targetIdList;
	if(config.existsParam("targetIdList")){
		existIdList = true;
		targetIdList = config.getParam("targetIdList");
	}

	XLine *linep;

	if(existIdList){

		// si le parametre targetIdList exist, on re-initialise la enrollList avec les vrais donnees des modeles
		// Create the XList for model definition
//		if(config.existsParam("targetIdList")){
			enrollList.reset();
			enrollList.load(targetIdList,config);
//		}

		// Sort XList
		enrollList.sortByElementNumber("descend");

		// Read targetIdList and fill _modelIDLine and _modelSessionslLine
		enrollList.getLine(0);
		linep=enrollList.getLine();
		String model, enrollSession;

		// Get the maximum number of sessions per speaker
		_n_enrollSessions_max = linep->getElementCount()-1;
		enrollList.getLine(0);
		while ((linep=enrollList.getLine()) != NULL){
		
			model = linep->getElement(0);
			if(_modelIDLine.getIndex(model) == -1){	//add the model ID to the list if it doesn't exists
				_modelIDLine.addElement(model);
			}

			unsigned long e = 1;
			while(e<linep->getElementCount()){
				enrollSession = linep->getElement(e);
				_modelIndexLine.addElement(model);
				_modelSessionslLine.addElement(enrollSession);
				e++;
			}
		}
	}
	else{
		// Get the maximum number of sessions per speaker
		_n_enrollSessions_max = 1;
	}

	//if(config.existsParam("inputVectorFilename")){
	//	// Fill the XLines from the vector list
	//	String vect;
	//	while ((linep=testList.getLine()) != NULL){
	//		vect = linep->getElement(0);
	//		_segLine.addElement(vect);
	//	}
	//}
	//else{
	if(!config.existsParam("inputVectorFilename")){
		// Read the XList and fill the XLines when required
		testList.getLine(0);
		String vect;
		while ((linep=testList.getLine()) != NULL){

			vect = linep->getElement(0);
			if(_segLine.getIndex(vect) == -1){
				_segLine.addElement(vect);
			}

			// Verify that all models involved in the testing have enrollment data
			unsigned long e = 1;
			while(e<linep->getElementCount()){
				vect = linep->getElement(e);
				if(_modelIDLine.getIndex(vect) == -1){
					if(existIdList){
						cout<<"	WARNING: model "<<vect<<"	has no enrollment session"<<endl;
					}
					else{
						_modelIDLine.addElement(vect);
						_modelIndexLine.addElement(vect);
						_modelSessionslLine.addElement(vect);
					}
				}
				e++;
			}
		}
	}
	_n_models				= _modelIDLine.getElementCount();
	_n_enrollment_segments	= _modelSessionslLine.getElementCount();
	_n_test_segments		= _segLine.getElementCount();


	// Get _vectSize from the first model
	String tmp = _modelSessionslLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("vectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_enrollment_segments);
	_segments.setDimensions(_vectSize,_n_test_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<" (PldaTest) found: "<<_n_models<<" models and "<<_n_test_segments<<" test segments"<<endl;
		cout<<" (PldaTest) Maximum number of sessions per speaker = "<<_n_enrollSessions_max<<endl;
	}

	// Load model enrollment vectors
	String fileName;
	for(unsigned long m=0;m<_n_enrollment_segments;m++){
		// Read the vector from file
		fileName = _modelSessionslLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_test_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) = tmpVect(0,k);
		}
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_test_segments);
	_trials.setAllValues(false);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_test_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_trials.setDimensions(_n_models,_n_test_segments);
	
	if(config.existsParam("inputVectorFilename")){
		_trials.setAllValues(true);
	}
	else{
		_trials.setAllValues(false);
		testList.getLine(0);
		String vect;
		while ((linep=testList.getLine()) != NULL){
			
			vect = linep->getElement(0);
			unsigned long segIdx = (unsigned long)_segLine.getIndex(vect);

			unsigned long e = 1;
			while(e<linep->getElementCount()){
				vect = linep->getElement(e);
				unsigned long modelIdx = (unsigned long)_modelIDLine.getIndex(vect);
				_trials(modelIdx,segIdx) = true; 
				e++;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getVectSize(){
	return _vectSize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getModelsNumber(){
	return _n_models;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getMaxEnrollmentSession(){
	return _n_enrollSessions_max;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getSegmentsNumber(){
	return _n_test_segments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getModels(){
	return _models;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getSegments(){
	return _segments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getScores(){
	return _scores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
BoolMatrix PldaTest::getTrials(){
	return _trials;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getModelName(unsigned long index){
	return(_modelIDLine.getElement(index));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getSegmentName(unsigned long index){
	return(_segLine.getElement(index));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::lengthNorm(){

	if(verboseLevel>0) cout<<"	TEST: Normalize length ... ";

	// Compute the norm of all models
	RealVector<double> vecNorm(0);
	vecNorm.setSize(_n_models);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_models;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _models(k,s)*_models(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_models;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,s) /=vecNorm[s];
		}
	}

	// Compute the norm of all segments
	vecNorm.setSize(_n_test_segments);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_test_segments;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _segments(k,s)*_segments(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_test_segments;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) /=vecNorm[s];
		}
	}
	if(verboseLevel>0) cout<<"done"<<endl;;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::center(RealVector<double> &mu){
	// center all models and test segments
	for(unsigned long k=0;k<_vectSize;k++){
		for(unsigned long s =0;s<_n_models;s++){
			_models(k,s) -= mu[k];
		}
		for(unsigned long s =0;s<_n_test_segments;s++){
			_segments(k,s) -= mu[k];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::rotateLeft(Matrix<double> &M){

	Matrix<double> tmpModels(_models);
	Matrix<double> tmpSegments(_segments);

	_models.setDimensions(M.rows(),_n_models);
	_models.setAllValues(0.0);
	
	_segments.setDimensions(M.rows(),_n_test_segments);
	_segments.setAllValues(0.0);

	// Rotate _models
	_models = M*tmpModels;

	// Rotate _segments
	_segments = M * tmpSegments;
	
	_vectSize = M.rows();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::sphericalNuisanceNormalization(Config& config){
	
	// Chose the number of iterations
	unsigned long nb_iterations = 1;
	if(config.existsParam("ivNormIterationNb"))	nb_iterations = config.getParam("ivNormIterationNb").toULong();

	//Prepare filename to load normalization matrices and means
	String sphNormMatrixBaseName = "ivNormEfrMatrix_it";
	if(config.existsParam("ivNormEfrMatrixBaseName"))	sphNormMatrixBaseName = config.getParam("ivNormEfrMatrixBaseName");
	
	String sphNormMeanBaseName = "ivNormEfrMean_it";
	if(config.existsParam("ivNormEfrMeanBaseName"))	sphNormMeanBaseName = config.getParam("ivNormEfrMeanBaseName");

	for(unsigned long it=0;it<nb_iterations;it++){
		
		if(verbose) cout<<"Eigen Factor Radial Normalization - start iteration "<<it<<endl;

		//Load matrix and mean for current iteration
		String s;
		String sphNormMatrixFilename = config.getParam("matrixFilesPath") + sphNormMatrixBaseName + s.valueOf(it) + config.getParam("loadMatrixFilesExtension");
		Matrix<double>  sphNormMat(sphNormMatrixFilename,config);
			
		//Load the mean of the iVector development set
		String sphNormMeanFilename = config.getParam("matrixFilesPath") + sphNormMeanBaseName + s.valueOf(it) + config.getParam("loadMatrixFilesExtension");
		Matrix<double> sphNormMean(sphNormMeanFilename,config);	// add an exception if the size of the mean is not the vectSize

		RealVector<double> mean(sphNormMean.cols(),sphNormMean.cols());
		for(unsigned long i=0;i<sphNormMean.cols();i++)
			mean[i] = sphNormMean(0,i);

		//Center iVectors
		if(verboseLevel>0) cout<<"		center the test data"<<endl;
		this->center(mean);
		
		//Rotate using the inverse square root of Sigma
		if(verboseLevel>0) cout<<"		rotate the test data"<<endl;
		this->rotateLeft(sphNormMat);
		
		//Normalize the length of iVectors
		if(verboseLevel>0) cout<<"		normalize the length of the test data"<<endl;
		this->lengthNorm();
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::cosineDistance(Config &config){

	RealVector<double> normM(0),normS(0);
	normM.setSize(_n_models);
	normS.setSize(_n_test_segments);
	normM.setAllValues(0.0);
	normS.setAllValues(0.0);

	// Compute Norm of models and segments
	for(unsigned long k=0;k<_vectSize;k++){
		for(unsigned long m=0;m<_n_models;m++){
			normM[m] += _models(k,m)*_models(k,m);
		}
		for(unsigned long s=0;s<_n_test_segments;s++){
			normS[s] += _segments(k,s)*_segments(k,s);
		}
	}

	for(unsigned long m=0;m<_n_models;m++){
		normM[m] = sqrt(normM[m]);
	}
	for(unsigned long s=0;s<_n_test_segments;s++){
		normS[s] = sqrt(normS[s]);
	}

	for(unsigned long m=0;m<_n_models;m++){
		for(unsigned long s=0;s<_n_test_segments;s++){
			if(_trials(m,s)){
				for(unsigned long k=0;k<_vectSize;k++){
					_scores(m,s) += _models(k,m)*_segments(k,s);
				}
				_scores(m,s) /= (normM[m]*normS[s]);
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::mahalanobisDistance(Matrix<double>& Mah, Config &config){

	for(unsigned long m=0;m<_n_models;m++){
		for(unsigned long s=0;s<_n_test_segments;s++){
	
			if(_trials(m,s)){
				DoubleVector tmp(_vectSize,_vectSize);
				for(unsigned long k=0;k<_vectSize;k++){
					tmp[k] = _models(k,m)-_segments(k,s);
				}

				DoubleVector t(_vectSize,_vectSize);
				t.setAllValues(0.0);
				for(unsigned long k=0;k<_vectSize;k++){
					for(unsigned long i=0;i<_vectSize;i++){
						t[k] += -0.5*tmp[i]*Mah(i,k);
					}
				}
	
				for(unsigned long i=0;i<_vectSize;i++){
					_scores(m,s) += t[i]*tmp[i];
				}	
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::twoCovScoring(DoubleSquareMatrix& W, DoubleSquareMatrix& B, Config &config){

	// Compute G and H matrices for two covariance scoring
	DoubleSquareMatrix invW(_vectSize);
	DoubleSquareMatrix invB(_vectSize);
	W.invert(invW);
	B.invert(invB);
	
	DoubleSquareMatrix tmp1, tmp2, tmp3, tmp4, sum1, sum2, G, H;
	sum1.setSize(_vectSize);sum2.setSize(_vectSize);tmp1.setSize(_vectSize);tmp2.setSize(_vectSize);tmp3.setSize(_vectSize);tmp4.setSize(_vectSize);G.setSize(_vectSize);H.setSize(_vectSize);
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			sum1(i,j) = invB(i,j)+2*invW(i,j);
			sum2(i,j) = invB(i,j)+invW(i,j);
		}
	}
	sum1.invert(tmp1);
	sum2.invert(tmp2);
		
	tmp3.setAllValues(0.0);
	tmp4.setAllValues(0.0);
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			for(unsigned long k=0;k<_vectSize;k++){
				tmp3(i,j) += invW(i,k)*tmp2(k,j);
				tmp4(i,j) += invW(i,k)*tmp1(k,j);
			}
		}
	}

	G.setAllValues(0.0);
	H.setAllValues(0.0);
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			for(unsigned long k=0;k<_vectSize;k++){
				G(i,j) += tmp3(i,k)*invW(k,j);
				H(i,j) += tmp4(i,k)*invW(k,j);
			}
		}
	}

	for(unsigned long m=0;m<_n_models;m++){

		// test.model_vec(ii,:)*H*test.model_vec(ii,:)'
		double scorePart1 = 0;
		DoubleVector a(_vectSize,_vectSize);
		a.setAllValues(0.0);
		for(unsigned long k=0;k<_vectSize;k++){
			for(unsigned long i=0;i<_vectSize;i++){
						a[k] += _models(i,m)*H(i,k);
			}
		}
		for(unsigned long j=0;j<_vectSize;j++){
			scorePart1 += a[j]*_models(j,m);
		}		

		for(unsigned long s=0;s<_n_test_segments;s++){
	
			// test.seg_vec(jj,:)*H*test.seg_vec(jj,:)'
			double scorePart2 = 0;
			DoubleVector b(_vectSize,_vectSize);
			b.setAllValues(0.0);
			for(unsigned long k=0;k<_vectSize;k++){
				for(unsigned long i=0;i<_vectSize;i++){
							b[k] += _segments(i,s)*H(i,k);
				}
			}
			for(unsigned long j=0;j<_vectSize;j++){
				scorePart1 += a[j]*_models(j,m);
			}
	
			if(_trials(m,s)){
				
				double scorePart3 = 0.0;
				
				// (test.model_vec(ii,:) + test.seg_vec(jj,:))*G*(test.model_vec(ii,:)+test.seg_vec(jj,:))'
				DoubleVector diff(_vectSize,_vectSize), c(_vectSize,_vectSize);
				c.setAllValues(0.0);
				diff.setAllValues(0.0);
				for(unsigned long j=0;j<_vectSize;j++){
					diff[j] = _models(j,m) + _segments(j,s);
				}
				for(unsigned long k=0;k<_vectSize;k++){
					for(unsigned long i=0;i<_vectSize;i++){
								c[k] += diff[i]*G(i,k);
					}
				}
				for(unsigned long j=0;j<_vectSize;j++){
					scorePart3 += c[j]*diff[j];
				}
				
				//S(ii,jj)    = (test.model_vec(ii,:) + test.seg_vec(jj,:))*G*(test.model_vec(ii,:)+test.seg_vec(jj,:))' - test.model_vec(ii,:)*H*test.model_vec(ii,:)' - test.seg_vec(jj,:)*H*test.seg_vec(jj,:)';
				_scores(m,s) = scorePart3 - scorePart1 - scorePart2;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::pldaNativeScoring(PldaModel& pldaModel, Config &config){

	//Precomputation for scoring
	pldaModel.preComputation();

	Eigen::MatrixXd FTJ		= pldaModel.getFtweight() - pldaModel.getFtweightG() * pldaModel.getInvGtweightGplusEye() * pldaModel.getGtweight();
	Eigen::MatrixXd FTJF	= FTJ * pldaModel.getF();

	//	Apply transformation to the models and test segments
	Matrix<double> tmpFTJ(pldaModel.getRankF(),_vectSize);
	for(unsigned long ii=0;ii<pldaModel.getRankF();ii++)
		for(unsigned long jj=0;jj<_vectSize;jj++)
			tmpFTJ(ii,jj) = FTJ(ii,jj);
	this->rotateLeft(tmpFTJ);

	// Compute K1
	Eigen::MatrixXd tmpK	= FTJF + Eigen::MatrixXd::Identity(pldaModel.getRankF(),pldaModel.getRankF());
	Eigen::MatrixXd K_one	= tmpK.inverse();

	// Compute alpha1
	Eigen::MatrixXd a = K_one.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
	Eigen::VectorXd b = a.diagonal();
	double alpha_one = 0.0;
	for(unsigned long i=0;i<pldaModel.getRankF();i++)
		alpha_one += log(b(i));
	alpha_one *= 2.0;

	// Define variables for loop on speakers
	Eigen::MatrixXd K_L, K_L_1;
	double alpha_L, alpha_L_1, Constant ;
	Eigen::VectorXd model,testSeg;

	Eigen::MatrixXd tmpTest_L(_segments.rows(),_segments.cols());
	for(unsigned long i=0;i<_segments.rows();i++)
		for(unsigned long j=0;j<_segments.cols();j++)
			tmpTest_L(i,j) = _segments(i,j);

	unsigned long sessionCounter = 0;		// counter used to know to which speaker belongs the current session
	unsigned long currentSessionNb = 0;
	for(unsigned long mod=0; mod < this->getModelsNumber(); mod++){

		String previousModel = _modelIndexLine.getElement(sessionCounter);

		// Initialize the model for the current speaker
		model	= Eigen::VectorXd::Zero(_vectSize);
		unsigned long spkSessionNumber = 0;

		// Get number of sessions, first and last session indexes of the current speaker
		while((sessionCounter < _n_enrollment_segments) && (_modelIndexLine.getElement(sessionCounter) == previousModel)){

			// Accumulate i-vectors for the current speaker
			for(unsigned long d=0;d<_vectSize;d++)
				model(d) += _models(d,sessionCounter);

			spkSessionNumber++;
			previousModel = _modelIndexLine.getElement(sessionCounter);
			sessionCounter++;
		}

		// If the number if session is different from the previous one, compute alpha and K for this configuration
		if(spkSessionNumber != currentSessionNb){

			currentSessionNb = spkSessionNumber;
			
			// Compute K(L+1)
			tmpK	= currentSessionNb*FTJF + Eigen::MatrixXd::Identity(pldaModel.getRankF(),pldaModel.getRankF());
			K_L	= tmpK.inverse();

			// Compute K(L)
			tmpK	= (currentSessionNb+1)*FTJF + Eigen::MatrixXd::Identity(pldaModel.getRankF(),pldaModel.getRankF());
			K_L_1	= tmpK.inverse();

			// Compute alpha(L)
			a = K_L.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
			b = a.diagonal();
			alpha_L = 0.0;
			for(unsigned long i=0;i<pldaModel.getRankF();i++)
				alpha_L += log(b(i));
			alpha_L *= 2.0;

			// Compute alpha(L+1)
			a = K_L_1.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
			b = a.diagonal();
			alpha_L_1 = 0.0;
			for(unsigned long i=0;i<pldaModel.getRankF();i++)
				alpha_L_1 += log(b(i));
			alpha_L_1 *= 2.0;

			Constant = (alpha_L_1 - alpha_L - alpha_one)/2.0;
		}

		double s2 = model.transpose()*K_L*model;

		// add the model to all test segments (create a new matrix)
		Eigen::MatrixXd tmpTest_L_1(_segments.rows(),_segments.cols());
		for(unsigned long c=0;c<_segments.cols();c++)
			tmpTest_L_1.col(c) = tmpTest_L.col(c) + model;

		Eigen::MatrixXd tmp = tmpTest_L.transpose() * K_one;
		Eigen::VectorXd s1(_segments.cols());

		Eigen::MatrixXd tmp_1 = tmpTest_L_1.transpose() * K_L_1;
		Eigen::VectorXd s3(_segments.cols());

		for(unsigned long i=0;i<_segments.cols();i++){
			s1(i) = tmp.row(i).dot(tmpTest_L.col(i));
			s3(i) = tmp_1.row(i).dot(tmpTest_L_1.col(i));
			_scores(mod,i) = (s3(i)-s2-s1(i))/2 + Constant;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::pldaMeanScoring(PldaModel& pldaModel, Config &config){

	//Precomputation for scoring
	pldaModel.preComputation();

	Eigen::MatrixXd FTJ		= pldaModel.getFtweight() - pldaModel.getFtweightG() * pldaModel.getInvGtweightGplusEye() * pldaModel.getGtweight();
	Eigen::MatrixXd FTJF	= FTJ * pldaModel.getF();

	//	Apply transformation to the models and test segments
	Matrix<double> tmpFTJ(pldaModel.getRankF(),_vectSize);
	for(unsigned long ii=0;ii<pldaModel.getRankF();ii++)
		for(unsigned long jj=0;jj<_vectSize;jj++)
			tmpFTJ(ii,jj) = FTJ(ii,jj);
	this->rotateLeft(tmpFTJ);

	// Compute K1
	Eigen::MatrixXd tmpK	= FTJF + Eigen::MatrixXd::Identity(pldaModel.getRankF(),pldaModel.getRankF());
	Eigen::MatrixXd K_one	= tmpK.inverse();

	// Compute alpha1
	Eigen::MatrixXd a = K_one.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
	Eigen::VectorXd b = a.diagonal();
	double alpha_one = 0.0;
	for(unsigned long i=0;i<pldaModel.getRankF();i++)
		alpha_one += log(b(i));
	alpha_one *= 2.0;

	// Compute K2
	tmpK	= 2*FTJF + Eigen::MatrixXd::Identity(pldaModel.getRankF(),pldaModel.getRankF());
	Eigen::MatrixXd K_two	= tmpK.inverse();

	// Compute alpha2
	a = K_two.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
	b = a.diagonal();
	double alpha_two = 0.0;
	for(unsigned long i=0;i<pldaModel.getRankF();i++)
		alpha_two += log(b(i));
	alpha_two *= 2.0;

	//Compute the constant component
	double Constant = (alpha_two - 2*alpha_one)/2.0;

	// Define variables for loop on speakers
	Eigen::VectorXd model,testSeg;

	// Convert test segment in Eigen format
	Eigen::MatrixXd tmpTest_L(_segments.rows(),_segments.cols());
	for(unsigned long i=0;i<_segments.rows();i++)
		for(unsigned long j=0;j<_segments.cols();j++)
			tmpTest_L(i,j) = _segments(i,j);

	unsigned long sessionCounter = 0;		// counter used to know to which speaker belongs the current session
	unsigned long currentSessionNb = 0;
	for(unsigned long mod=0; mod < this->getModelsNumber(); mod++){

		String previousModel = _modelIndexLine.getElement(sessionCounter);

		// Initialize the model for the current speaker
		model	= Eigen::VectorXd::Zero(_vectSize);
		unsigned long spkSessionNumber = 0;

		// Get number of sessions, first and last session indexes of the current speaker
		while((sessionCounter < _n_enrollment_segments) && (_modelIndexLine.getElement(sessionCounter) == previousModel)){

			// Accumulate i-vectors for the current speaker
			for(unsigned long d=0;d<_vectSize;d++)
				model(d) += _models(d,sessionCounter);

			spkSessionNumber++;
			previousModel = _modelIndexLine.getElement(sessionCounter);
			sessionCounter++;
		}

		//Compute the mean of the speaker i-vectors
		model /= spkSessionNumber;

		double s2 = model.transpose()*K_one*model;

		// add the model to all test segments (create a new matrix)
		Eigen::MatrixXd tmpTest_L_1(_segments.rows(),_segments.cols());
		for(unsigned long c=0;c<_segments.cols();c++)
			tmpTest_L_1.col(c) = tmpTest_L.col(c) + model;

		Eigen::MatrixXd tmp = tmpTest_L.transpose() * K_one;
		Eigen::VectorXd s1(_segments.cols());

		Eigen::MatrixXd tmp_1 = tmpTest_L_1.transpose() * K_two;
		Eigen::VectorXd s3(_segments.cols());

		for(unsigned long i=0;i<_segments.cols();i++){
			s1(i) = tmp.row(i).dot(tmpTest_L.col(i));
			s3(i) = tmp_1.row(i).dot(tmpTest_L_1.col(i));
			_scores(mod,i) = (s3(i)-s2-s1(i))/2 + Constant;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::saveSegments(String outputDir, Config& config){

	if(verboseLevel>0) cout<<"(PldaTest) Save test segments"<<endl;

	String vectExtension = config.getParam("vectorFilesExtension");
	String vectName;

	Matrix<double> tmp(1,_vectSize);

	for(unsigned long s=0;s<_n_test_segments;s++){
		//Get output vector filename
		vectName =  outputDir + _segLine.getElement(s) + vectExtension;

		//Get the vector to save from _segments
		for(unsigned long i=0;i<_vectSize;i++){
			tmp(0,i) = _segments(i,s);
		}
		tmp.save(vectName,config);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::saveVectors(String outputDir, Config& config){

	if(verboseLevel>0) cout<<"(PldaTest) Save enrollment and test segments"<<endl;

	String vectExtension = config.getParam("vectorFilesExtension");
	String vectName;

	Matrix<double> tmp(1,_vectSize);

	//Save enrollment data
	for(unsigned long s=0;s<_n_enrollment_segments;s++){
		
		//Get output vector filename
		vectName =  outputDir + _modelSessionslLine.getElement(s) + vectExtension;

		//Get the vector to save from _segments
		for(unsigned long i=0;i<_vectSize;i++)
			tmp(0,i) = _models(i,s);

		tmp.save(vectName,config);
	}

	//Save test segments
	for(unsigned long s=0;s<_n_test_segments;s++){
		
		//Get output vector filename
		vectName =  outputDir + _segLine.getElement(s) + vectExtension;

		//Get the vector to save from _segments
		for(unsigned long i=0;i<_vectSize;i++)
			tmp(0,i) = _segments(i,s);

		tmp.save(vectName,config);
	}
}



















#endif
