#if !defined(ALIZE_ScorpioGenerator_cpp)
#define ALIZE_ScorpioGenerator_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "ScorpioGenerator.h"

using namespace alize;
using namespace std;




//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int ScorpioGenerator(Config& config)
{

// TO DO: Add UBM writing

try{

	config.setParam("iVectSize", config.getParam("totalVariabilityNumber"));

	//----------------------------------------
	// LOAD WORLD MODELS (MALE & FEMALE)
	//----------------------------------------
	String inputMaleWorldFilename = config.getParam("maleInputWorldFilename");
	String inputFemaleWorldFilename = config.getParam("femaleInputWorldFilename");
 
	MixtureServer ms(config);
	if (verbose) cout << "(ScorpioGenerator) Load male world model [" << inputMaleWorldFilename<<"]"<<endl;
	MixtureGD& worldMale = ms.loadMixtureGD(inputMaleWorldFilename);
	if (verbose) cout << "(ScorpioGenerator) Load female world model [" << inputMaleWorldFilename<<"]"<<endl;
	MixtureGD& worldFemale = ms.loadMixtureGD(inputFemaleWorldFilename);  


	//----------------------------------------
	// LOAD MALE META PARAMETERS
	//----------------------------------------
	config.setParam("inputWorldFilename", config.getParam("maleInputWorldFilename"));
	config.setParam("totalVariabilityMatrix", config.getParam("maleTotalVariabilityMatrix"));
	config.setParam("meanEstimate", config.getParam("maleMeanEstimate"));
	config.setParam("sphNormMat", config.getParam("maleSphNormMat"));
	config.setParam("sphNormMean", config.getParam("maleSphNormMean"));
	config.setParam("pldaEigenVoiceMatrix", config.getParam("malePldaEigenVoiceMatrix"));
	//config.setParam("pldaEigenChannelMatrix", config.getParam("malePldaEigenChannelMatrix"));
	config.setParam("pldaSigmaMatrix", config.getParam("malePldaSigmaMatrix"));
	config.setParam("pldaOriginalMean", config.getParam("malePldaOriginalMean"));
	config.setParam("pldaDelta", config.getParam("malePldaDelta"));

	// Load I-Vector extractor parameters
		TVAcc maleTvAcc(config);

		//Load TotalVariability matrix
		maleTvAcc.loadT(config.getParam("totalVariabilityMatrix"), config);

		// Then load the meanEstimate computed by minDiv if required
		DoubleVector meanEstimate = maleTvAcc.getUbmMeans();
		if(config.existsParam("minDivergence")&& config.getParam("minDivergence").toBool()){
			String minDivName = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("loadMatrixFilesExtension");
			Matrix<double> tmpMean(minDivName,config);
			for(unsigned long i=0;i<meanEstimate.size();i++){
				meanEstimate[i] = tmpMean(0,i);
			}
		}
		//Update the mean Estimate
		cout<<"	(TrainTargetiVector) Load Mean Estimate"<<endl;
		maleTvAcc.loadMeanEstimate(meanEstimate);

		// Compute T*invSigma
		unsigned long _svSize = maleTvAcc.getSvSize();
		Matrix<double> maleTE(_svSize,maleTvAcc.getRankT());
		Matrix<double> mT = maleTvAcc.getT();
		double* tp= mT.getArray();
		DoubleVector invVar = maleTvAcc.getUbmInvVar();
		maleTE.setAllValues(0.0);
		double *te = maleTE.getArray();
		for(unsigned long i=0;i<maleTvAcc.getRankT();i++){
			for(unsigned long k=0;k<_svSize;k++) {
				te[i*_svSize+k] = invVar[k] * tp[i*_svSize+k];
			}
		}

		//Compute vEvT for each session
		maleTvAcc.estimateTETt(config);

	// Load Male Normalization parameters
		// Load Whitening matrix and mean
		String SphNormMat_filename = config.getParam("matrixFilesPath") + config.getParam("sphNormMat") + config.getParam("loadMatrixFilesExtension");
		String SphNormMean_filename = config.getParam("matrixFilesPath") + config.getParam("sphNormMean") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> maleSphNormMat(SphNormMat_filename,config);
		Matrix<double> maleSphNormMean(SphNormMean_filename,config);


	// Load Male PLDA scoring meta-parameters
		PldaModel malePldaModel("test",config);	// Load in test mode

		//Precomputation for scoring
		malePldaModel.preComputation();

		Eigen::MatrixXd maleFTJ		= malePldaModel.getFtweight() - malePldaModel.getFtweightG() * malePldaModel.getInvGtweightGplusEye() * malePldaModel.getGtweight();
		Eigen::MatrixXd FTJF		= maleFTJ * malePldaModel.getF();

		// Compute K1
		Eigen::MatrixXd tmpK		= FTJF + Eigen::MatrixXd::Identity(malePldaModel.getRankF(),malePldaModel.getRankF());
		Eigen::MatrixXd maleK_one	= tmpK.inverse();

		// Compute alpha1
		Eigen::MatrixXd a = maleK_one.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
		Eigen::VectorXd b = a.diagonal();
		double alpha_one = 0.0;
		for(unsigned long i=0;i<malePldaModel.getRankF();i++)
			alpha_one += log(b(i));
		alpha_one *= 2.0;

		// Compute K2
		Eigen::MatrixXd tmpK2	= 2*FTJF + Eigen::MatrixXd::Identity(malePldaModel.getRankF(),malePldaModel.getRankF());
		Eigen::MatrixXd maleK_two	= tmpK2.inverse();

		// Compute alpha1
		Eigen::MatrixXd a2 = maleK_two.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
		Eigen::VectorXd b2 = a2.diagonal();
		double alpha_two = 0.0;
		for(unsigned long i=0;i<malePldaModel.getRankF();i++)
			alpha_two += log(b2(i));
		alpha_two *= 2.0;

		double maleConstant = (alpha_two - 2*alpha_one)/2.0;



	//----------------------------------------
	// LOAD FEMALE META PARAMETERS
	//----------------------------------------
	config.setParam("inputWorldFilename", config.getParam("femaleInputWorldFilename"));
	config.setParam("totalVariabilityMatrix", config.getParam("femaleTotalVariabilityMatrix"));
	config.setParam("meanEstimate", config.getParam("femaleMeanEstimate"));
	config.setParam("sphNormMat", config.getParam("femaleSphNormMat"));
	config.setParam("sphNormMean", config.getParam("femaleSphNormMean"));
	config.setParam("pldaEigenVoiceMatrix", config.getParam("femalePldaEigenVoiceMatrix"));
	//config.setParam("pldaEigenChannelMatrix", config.getParam("femalePldaEigenChannelMatrix"));
	config.setParam("pldaSigmaMatrix", config.getParam("femalePldaSigmaMatrix"));
	config.setParam("pldaOriginalMean", config.getParam("femalePldaOriginalMean"));
	config.setParam("pldaDelta", config.getParam("femalePldaDelta"));

	// Load I-Vector extractor parameters
		TVAcc femaleTvAcc(config);

		//Load TotalVariability matrix
		femaleTvAcc.loadT(config.getParam("totalVariabilityMatrix"), config);

		// Then load the meanEstimate computed by minDiv if required
		meanEstimate = femaleTvAcc.getUbmMeans();
		if(config.existsParam("minDivergence")&& config.getParam("minDivergence").toBool()){
			String minDivName = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("loadMatrixFilesExtension");
			Matrix<double> tmpMean(minDivName,config);
			for(unsigned long i=0;i<meanEstimate.size();i++){
				meanEstimate[i] = tmpMean(0,i);
			}
		}
		//Update the mean Estimate
		cout<<"	(TrainTargetiVector) Load Mean Estimate"<<endl;
		femaleTvAcc.loadMeanEstimate(meanEstimate);

		// Compute T*invSigma
		_svSize = femaleTvAcc.getSvSize();
		Matrix<double> femaleTE(_svSize,femaleTvAcc.getRankT());
		Matrix<double> fT = femaleTvAcc.getT();
		tp= fT.getArray();
		invVar = femaleTvAcc.getUbmInvVar();
		femaleTE.setAllValues(0.0);
		te = femaleTE.getArray();
		for(unsigned long i=0;i<femaleTvAcc.getRankT();i++){
			for(unsigned long k=0;k<_svSize;k++) {
				te[i*_svSize+k] = invVar[k] * tp[i*_svSize+k];
			}
		}

		//Compute vEvT for each session
		femaleTvAcc.estimateTETt(config);

	// Load Male Normalization parameters
		// Load Whitening matrix and mean
		SphNormMat_filename = config.getParam("matrixFilesPath") + config.getParam("sphNormMat") + config.getParam("loadMatrixFilesExtension");
		SphNormMean_filename = config.getParam("matrixFilesPath") + config.getParam("sphNormMean") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> femaleSphNormMat(SphNormMat_filename,config);
		Matrix<double> femaleSphNormMean(SphNormMean_filename,config);


	// Load Male PLDA scoring meta-parameters
		PldaModel femalePldaModel("test",config);	// Load in test mode

		//Precomputation for scoring
		femalePldaModel.preComputation();

		Eigen::MatrixXd femaleFTJ	= femalePldaModel.getFtweight() - femalePldaModel.getFtweightG() * femalePldaModel.getInvGtweightGplusEye() * femalePldaModel.getGtweight();
		FTJF		= femaleFTJ * femalePldaModel.getF();

		// Compute K1
		tmpK		= FTJF + Eigen::MatrixXd::Identity(femalePldaModel.getRankF(),femalePldaModel.getRankF());
		Eigen::MatrixXd femaleK_one	= tmpK.inverse();

		// Compute alpha1
		a = femaleK_one.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
		b = a.diagonal();
		alpha_one = 0.0;
		for(unsigned long i=0;i<femalePldaModel.getRankF();i++)
			alpha_one += log(b(i));
		alpha_one *= 2.0;

		// Compute K2
		tmpK2	= 2*FTJF + Eigen::MatrixXd::Identity(malePldaModel.getRankF(),malePldaModel.getRankF());
		Eigen::MatrixXd femaleK_two	= tmpK2.inverse();

		// Compute alpha1
		a2 = femaleK_two.llt().matrixL();	//lower triangular matrix from the Cholesky decomposition
		b2 = a2.diagonal();
		alpha_two = 0.0;
		for(unsigned long i=0;i<femalePldaModel.getRankF();i++)
			alpha_two += log(b2(i));
		alpha_two *= 2.0;

		double femaleConstant = (alpha_two - 2*alpha_one)/2.0;



	//----------------------------------------
	// WRITE OUTPUT FILE
	//----------------------------------------

	//Open output file
	String fname = config.getParam("outputFile");
	ofstream outfile(fname.c_str(), ios::out);

	// WRITE head of the file
	outfile<<"#define VECT_SIZE "<<config.getParam("featureVectSize").toULong()<<endl;
	outfile<<"#define DISTRIB_SIZE "<<config.getParam("distribSize").toULong()<<endl;
	outfile<<"#define RANKT "<<config.getParam("totalVariabilityNumber").toULong()<<endl;
	outfile<<"#define RANKF "<<config.getParam("pldaEigenVoiceNumber").toULong()<<endl<<endl;

	// SAVE the mean vector to use for mean substraction
	outfile<<"double MEAN_MIN_DIV[2][VECT_SIZE*DISTRIB_SIZE] = {"<<endl;

	// Write FEMALE MEAN_MIN_DIV
		DoubleVector mean = femaleTvAcc.getUbmMeans();
		outfile<<"	{ ";
		for(unsigned long ii=0;ii<mean.size()-1;ii++)
			outfile<<mean[ii]<<",";
		outfile<<mean[mean.size()-1]<<" },"<<endl;

	// Write MALE MEAN_MIN_DIV
		mean = maleTvAcc.getUbmMeans();
		outfile<<"	{ ";
		for(unsigned long ii=0;ii<mean.size()-1;ii++)
			outfile<<mean[ii]<<",";
		outfile<<mean[mean.size()-1]<<" }"<<endl;
	outfile<<"};"<<endl<<endl;


	// SAVE the TETt matrices in one big ,matrix
	outfile<<"double TETt[2][DISTRIB_SIZE*RANKT*RANKT] = {"<<endl;

		// Write FEMALE TETt
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleTvAcc.getRankT()*femaleTvAcc.getRankT()-1;ii++){
			for(unsigned long dd=0;dd<femaleTvAcc.getNDistrib();dd++){
				DoubleSquareMatrix TETt = femaleTvAcc.getTETt(dd);
				double *tett = TETt.getArray();
				outfile<<tett[ii]<<",";
			}
		}
		// Dernier coefficient
		for(unsigned long dd=0;dd<femaleTvAcc.getNDistrib()-1;dd++){
			DoubleSquareMatrix TETt = femaleTvAcc.getTETt(dd);
			double *tett = TETt.getArray();
			outfile<<tett[femaleTvAcc.getRankT()*femaleTvAcc.getRankT()-1]<<",";
		}
		//dernier coefficient, derniere matrice
		DoubleSquareMatrix TETt = femaleTvAcc.getTETt(femaleTvAcc.getNDistrib()-1);
		double *tett = TETt.getArray();
		outfile<<tett[femaleTvAcc.getRankT()*femaleTvAcc.getRankT()-1]<<" },"<<endl;


		// Write MALE TETt
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleTvAcc.getRankT()*maleTvAcc.getRankT()-1;ii++){
			for(unsigned long dd=0;dd<maleTvAcc.getNDistrib();dd++){
				DoubleSquareMatrix TETt = maleTvAcc.getTETt(dd);
				double *tett = TETt.getArray();
				outfile<<tett[ii]<<",";
			}
		}
		// Dernier coefficient
		for(unsigned long dd=0;dd<maleTvAcc.getNDistrib()-1;dd++){
			DoubleSquareMatrix TETt = maleTvAcc.getTETt(dd);
			double *tett = TETt.getArray();
			outfile<<tett[maleTvAcc.getRankT()*maleTvAcc.getRankT()-1]<<",";
		}
		//dernier coefficient, derniere matrice
		TETt = maleTvAcc.getTETt(maleTvAcc.getNDistrib()-1);
		tett = TETt.getArray();
		outfile<<tett[maleTvAcc.getRankT()*maleTvAcc.getRankT()-1]<<" }"<<endl;
	outfile<<"};"<<endl<<endl;


	// SAVE the TE matrix
	outfile<<"double TE[2][RANKT*DISTRIB_SIZE*VECT_SIZE] = {"<<endl;
		// Write FEMALE TE
		double *t = femaleTE.getArray();
	
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleTE.rows()*femaleTE.cols()-1;ii++)
			outfile<<t[ii]<<",";
		outfile<<t[femaleTE.rows()*femaleTE.cols()-1]<<"},"<<endl;

		// Write MALE TE
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleTE.rows()*maleTE.cols()-1;ii++)
			outfile<<t[ii]<<",";
		outfile<<t[maleTE.rows()*maleTE.cols()-1]<<"}"<<endl;
	outfile<<"};"<<endl<<endl;

	
	// SAVE FTJ
	outfile<<"double FTJ[2][RANKF*RANKT] = {"<<endl;
		// Write FEMALE FTJ
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleFTJ.cols()-1;ii++){
			for(unsigned long jj=0;jj<femaleFTJ.rows();jj++){
				outfile<<femaleFTJ(jj,ii)<<",";
			}
		}
			for(unsigned long jj=0;jj<femaleFTJ.rows()-1;jj++){		// last column
				outfile<<femaleFTJ(jj,femaleFTJ.cols()-1)<<",";
			}
		outfile<<femaleFTJ(femaleFTJ.rows()-1,femaleFTJ.cols()-1)<<" },"<<endl;

		// Write MALE FTJ
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleFTJ.cols()-1;ii++)
			for(unsigned long jj=0;jj<maleFTJ.rows();jj++)
				outfile<<maleFTJ(jj,ii)<<",";
			for(unsigned long jj=0;jj<maleFTJ.rows()-1;jj++)		// last column
				outfile<<maleFTJ(jj,maleFTJ.cols()-1)<<",";
		outfile<<maleFTJ(maleFTJ.rows()-1,maleFTJ.cols()-1)<<" }"<<endl;
	outfile<<"};"<<endl<<endl;

	// SAVE SphNorm parameters
	
	// WRITE SPHNORM Matrix
	outfile<<"double SN_MAT[2][RANKT*RANKT] = {"<<endl;
		// Write FEMALE SN_MAT
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleSphNormMat.cols()-1;ii++){
			for(unsigned long jj=0;jj<femaleSphNormMat.rows();jj++){
				outfile<<femaleSphNormMat(jj,ii)<<",";
			}
		}
			for(unsigned long jj=0;jj<femaleSphNormMat.rows()-1;jj++){		// last column
				outfile<<femaleSphNormMat(jj,femaleSphNormMat.cols()-1)<<",";
			}
		outfile<<femaleSphNormMat(femaleSphNormMat.rows()-1,femaleSphNormMat.cols()-1)<<" },"<<endl;

		// Write MALE SN_MAT
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleSphNormMat.cols()-1;ii++)
			for(unsigned long jj=0;jj<maleSphNormMat.rows();jj++)
				outfile<<maleSphNormMat(jj,ii)<<",";
			for(unsigned long jj=0;jj<maleSphNormMat.rows()-1;jj++)		// last column
				outfile<<maleSphNormMat(jj,maleSphNormMat.cols()-1)<<",";
		outfile<<maleSphNormMat(maleSphNormMat.rows()-1,maleSphNormMat.cols()-1)<<" }"<<endl;
	outfile<<"};"<<endl<<endl;	

	// WRITE SPHNORM Mean vectors
	outfile<<"double SN_MEAN[2][RANKT] = {"<<endl;
		// Write FEMALE SN_MEAN
		outfile<<"	{";
		double *snm = femaleSphNormMean.getArray();
		for(unsigned long ii=0;ii<femaleSphNormMean.cols()-1;ii++)
			outfile<<snm[ii]<<",";
		outfile<<snm[femaleSphNormMean.cols()-1]<<" },"<<endl;

		// Write MALE SN_MEAN
		outfile<<"	{";
		snm = maleSphNormMean.getArray();
		for(unsigned long ii=0;ii<maleSphNormMean.cols()-1;ii++)
			outfile<<snm[ii]<<",";
		outfile<<snm[maleSphNormMean.cols()-1]<<" }"<<endl;
	outfile<<"};"<<endl<<endl;

	// SAVE K_one
	outfile<<"double K_ONE[2][RANKF*RANKF] = {"<<endl;

		// Write FEMALE K_one
		double *k1 = femaleK_one.data();
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleK_one.rows()*femaleK_one.cols()-1;ii++)
			outfile<<k1[ii]<<",";
		outfile<<k1[femaleK_one.rows()*femaleK_one.cols()-1]<<" },"<<endl;

		// Write MALE K_one
		k1 = maleK_one.data();
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleK_one.rows()*maleK_one.cols()-1;ii++)
			outfile<<k1[ii]<<",";
		outfile<<k1[maleK_one.rows()*maleK_one.cols()-1]<<" }"<<endl;
	outfile<<"};"<<endl<<endl;


	// SAVE K_two
	outfile<<"double K_TWO[2][RANKF*RANKF] = {"<<endl;

		// Write FEMALE K_two
		double *k2 = femaleK_two.data();
		outfile<<"	{";
		for(unsigned long ii=0;ii<femaleK_two.rows()*femaleK_two.cols()-1;ii++)
			outfile<<k2[ii]<<",";
		outfile<<k2[femaleK_two.rows()*femaleK_two.cols()-1]<<" },"<<endl;

		// Write MALE K_two
		k2 = maleK_two.data();
		outfile<<"	{";
		for(unsigned long ii=0;ii<maleK_two.rows()*maleK_two.cols()-1;ii++)
			outfile<<k2[ii]<<",";
		outfile<<k2[maleK_two.rows()*maleK_two.cols()-1]<<" }"<<endl;
	outfile<<"};"<<endl<<endl;

	// SAVE Constant
	outfile<<"double PLDA_CST[] = {"<<femaleConstant<<","<<maleConstant<<"};"<<endl;

	outfile.close();

} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

#endif //!defined(ALIZE_TrainTarget_cpp)
