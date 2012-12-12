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


#if !defined(ALIZE_ScorpioTest_cpp)
#define ALIZE_ScorpioTest_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "ScorpioTest.h"

#include "IvExtractor.h"

using namespace std;
using namespace alize;

using namespace Eigen;


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int ScorpioTest(Config &config){

	if(verboseLevel>0) cout<<"(ScorpioTest) ScorpioTest"<<endl;

	// Load statistics into VectorXd

	String statZero = "system32/mat/stat_N_32.matx";
	String statOne  = "system32/mat/stat_F_X_32.matx";
	alize::Matrix<double >stat0(statZero,config);
	alize::Matrix<double >stat1(statOne,config);

	cerr<<"size of statZero : "<<stat0.rows()<<"	"<<stat0.cols()<<endl;
	cerr<<"size of statOne : "<<stat1.rows()<<"	"<<stat1.cols()<<endl;

	VectorXd S0(stat0.cols());
	VectorXd S1(stat1.cols());

	for(unsigned long ii=0;ii<stat0.cols();ii++)
		S0(ii) = stat0(0,ii);
	for(unsigned long ii=0;ii<stat1.cols();ii++)
		S1(ii) = stat1(0,ii);

	real_t *s_0 = S0.data();
	real_t *s_1 = S1.data();

	unsigned long enrollSessionNb = 1;

	// Load test i-vector
	String testSegFilename  = "system32/ivec32/NIST05/jaac_A.y";
	alize::Matrix<double >testSeg(testSegFilename,config);
	VectorXd testIV(RANKT);
	for(unsigned long i=0;i<RANKT;i++)
		testIV(i) = testSeg(0,i);

	// Initialize temporary variables
	MatrixXd L = MatrixXd::Identity(RANKT,RANKT);
	Map<VectorXd> vL(L.data(),RANKT*RANKT);

//-------------------------------------------------------------------------------------------	
	// TRAIN FUNCTION TO ENCAPSULATE

		// Load meta-parameters for I-Vector extraction
		Map<VectorXd> mean_min_div(MEAN_MIN_DIV[0],DISTRIB_SIZE*VECT_SIZE);
		Map<MatrixXd> tett(TETt[0],DISTRIB_SIZE,RANKT*RANKT);
		Map<MatrixXd> T(TE[0],DISTRIB_SIZE*VECT_SIZE,RANKT);

		// Load meta-parameters for I-Vector normalization
		Map<VectorXd> sphnorm_mean(SN_MEAN[0],RANKT);
		Map<MatrixXd> sphnorm_mat(SN_MAT[0],RANKT,RANKT);

		// Compute i-vector from statistics

			// Substract mean_min_div
			double *mu = mean_min_div.data();
			for(unsigned long i=0; i<DISTRIB_SIZE; i++){
				for(unsigned long j = 0; j< VECT_SIZE;j++){
					s_1[i*VECT_SIZE+j] -= mu[i*VECT_SIZE+j]*s_0[i];
				}
			}

			// Accumulate vL
			for(unsigned long dis=0; dis<DISTRIB_SIZE;dis++){
				vL += S0(dis)* tett.row(dis);
			}

			// Compute auxiliary accumulator
			VectorXd aux = VectorXd::Zero(RANKT);

			aux = T.transpose()*S1;

			// Estimate w
			VectorXd modelIV = aux.transpose() * L.inverse();

		// Normalise i-vector 
			//Center
			modelIV -= sphnorm_mean;
			//Whitening			
			modelIV = sphnorm_mat * modelIV;
			//LengthNorm
			modelIV /= modelIV.norm();

		// Add to the existing SpeakerModel and incremente the number of sessions


//-------------------------------------------------------------------------------------------
	// Normalise test i-vector	(done during the i-vector extraction in TRAIN FUNCTION)
			//Center
			testIV -= sphnorm_mean;
			//Whitening			
			testIV = sphnorm_mat * modelIV;
			//LengthNorm
			testIV /= modelIV.norm();

//-------------------------------------------------------------------------------------------

	// TEST FUNCTION TO ENCAPSULATE

		// Load meta-parameters
		Map<MatrixXd> FtJ(FTJ[0],RANKF,RANKT);
		Map<MatrixXd> K_one(K_ONE[0],RANKF,RANKF);
		Map<MatrixXd> K_two(K_TWO[0],RANKF,RANKF);
		double pldaCst = PLDA_CST[0];

		// Compute Test using PLDA

			// Divide model i-vector by number of enrollment session
			modelIV /= enrollSessionNb;

			// Compute score
			VectorXd model, test;
			model = FtJ * modelIV;
			test = FtJ * testIV;

			VectorXd ModelPlusTest = model + test;

			double s3 = ModelPlusTest.transpose() * K_two * ModelPlusTest;
			double s2 = model.transpose()* K_one * model;
			double s1 = test.transpose() * K_one * test;

			double score = (s3-s2-s1)/2 + pldaCst;

cerr<<"SCORE = "<<score<<endl;

return 0;
}


#endif 
