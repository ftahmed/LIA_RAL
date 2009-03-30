/*!
\file seg.cpp
\author LIUM
\author Meignier Sylvain
\author sylvain.meignier@lium.univ-lemans.fr

\author LIUM
\author Dufour Richard
\author richard.dufour@lium.univ-lemans.fr

\author LIUM
\author Jousse Vincent
\author vincent.jousse@lium.univ-lemans.fr

* This file is part of LIUM_Seg.
* Copyright (c) 2007 Université du Maine. All rights reserved.

* GNU General Public License, as published by the Free Software
* Foundation; either version 2 of the License.
* 
* this program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with main.cpp.  If not, write to:
* 	The Free Software Foundation, Inc.,
* 	51 Franklin Street, Fifth Floor
* 	Boston, MA  02110-1301, USA.
*/

#include "alize.h"
#include <math.h>
#include <vector>

#include "distance.h"

#include "clusterServer.h"
#include "segment.h"
#include "cluster.h"

#define MIN(x,y) ((x)<(y))?(x):(y)
#define MAX(x,y) ((x)>(y))?(x):(y)

using namespace alize;

/// \brief storage for similarity measures
typedef vector < double >Measures;
/// \brief storage for segment border measures 
typedef map < long, double >Borders;

/// \brief select and compute the similarity method 
/// \param g1 the first Gaussien
/// \param g2 the second Gaussien
/// \param param the parameter structure
/// \param BICCst the constant need in BIC similarity
/// \return the similarity
double getSimilarity(FrameAcc & accG1, FrameAcc & accG2, FrameAcc & accG12, String segMethod,
					 double BICCst)
{
	
	if (segMethod == "GLR") {
		return Distance::GLR(accG1, accG2, accG12);
	} 
		
	else {
		if (segMethod == "BIC") {
			long len = accG1.getCount() + accG2.getCount();
			return Distance::BIC(accG1, accG2, accG12, BICCst, len);
		} else {
			if (segMethod == "KL2") {
				return Distance::KL2(static_cast<FrameAccGD&>(accG1), static_cast<FrameAccGD&>(accG2));
			} else {
				if (segMethod == "GD") {
					return Distance::GD(static_cast<FrameAccGD&>(accG1), static_cast<FrameAccGD&>(accG2));
				} else {
					if (segMethod == "H2") {
						return Distance::H2(static_cast<FrameAccGD&>(accG1), static_cast<FrameAccGD&>(accG2));
					} else {
						throw Exception("mSeg unknown segMethod param ", __FILE__, __LINE__);
					}
				}
			}
		}
	}
	return 0.0;
}

/// \brief Compute all the similarity
/// \return a array of similarity
//Measures doMeasures(Features & features, const Seg & seg, Param & param)
Measures doMeasures(FeatureServer & fs, const Segment & seg, ClusterServer & clusterServer, Config & config)
{
	unsigned long vectDimension = fs.getVectSize();
	
	//A changer prendre le debut du segment
	unsigned long start=seg.getBegin();
	unsigned long nb = MIN(seg.getLength(), fs.getFeatureCount());
	unsigned long segWSize = config.getIntegerParam("segWSize");
	
	Measures measures(nb);
	
#ifdef DEBUG
	cerr << "trace[mSeg] \t start="<<start<<" len="<<nb<<endl;
	cerr << "trace[mSeg] \t dim="<<vectDimension<<endl;
	cerr << "trace[mSeg] \t segWSize="<<segWSize<<endl;
#endif
	 
		
	//If the segment is to small to contain the 2 windows
	//return an empty measures
	if(nb<(2*segWSize)) {
		for(unsigned long i=0;i<nb;i++) {
			measures[i] = -1.0;
		}
		return measures;
	}
	
	
	double cst = Distance::BICCst(config.getParam("distribType"),vectDimension,config.getFloatParam("segThr"));
	double s;
	
	FrameAcc* fAccG1;
	FrameAcc* fAccG2;
	FrameAcc* fAccG12;
	
	if(config.getParam("distribType") == "GF") {
		fAccG1 = new FrameAccGF();
		fAccG2 = new FrameAccGF();
		fAccG12 = new FrameAccGF();
	} else {
		fAccG1 = new FrameAccGD();
		fAccG2 = new FrameAccGD();
		fAccG12 = new FrameAccGD();
	}
	
	
	fAccG1->reset();
	fAccG2->reset();
	fAccG12->reset();
	
	//we check the name of the source
	String sourceNameTemp = clusterServer.getSourceName(seg.getIdxSrcName());

	//in order to manage multiple stream source, we add to the begin and the end
	//the position in feature of the source
	unsigned long position = fs.getFirstFeatureIndexOfASource(sourceNameTemp);	
	start+=position;
	
	for(unsigned int i=start;i<segWSize;i++) {
		fs.seekFeature(i);
		Feature f;
		fs.readFeature(f);
		fAccG1->accumulate(f);
		fAccG12->accumulate(f);
	}	
	
	for(unsigned int i=start+segWSize;i<start+(2*segWSize);i++) {
		fs.seekFeature(i);
		Feature f;
		Feature f2;
		fs.readFeature(f);
		fAccG2->accumulate(f);
		fAccG12->accumulate(f);
	}	
	
	s = getSimilarity(*fAccG1,*fAccG2,*fAccG12,config.getParam("segMethod"),cst);
		
	#ifdef DEBUG
	cerr << "trace[mSeg] \t similarity="<<s<<endl;
	#endif
	
	for (unsigned long i = 0; i < segWSize; i++) {
		measures[i] = s;
	}	
	
	
	for (unsigned long i = position + segWSize; i < nb - segWSize; i++) {

		fs.seekFeature(start+i - segWSize);
		Feature f;
		fs.readFeature(f);
		fAccG1->deaccumulate(f);
		fAccG12->deaccumulate(f);
		
		fs.seekFeature(start+i);
		fs.readFeature(f);
		fAccG1->accumulate(f);
			
		fs.seekFeature(start+i);
		fs.readFeature(f);
		fAccG2->deaccumulate(f);
		
		fs.seekFeature(start+i + segWSize);
		fs.readFeature(f);
		fAccG2->accumulate(f);
		fAccG12->accumulate(f);
			
		s = getSimilarity(*fAccG1,*fAccG2,*fAccG12,config.getParam("segMethod"),cst);	
		
		measures[i] = s;
		

	}	
	
	for (unsigned long i = nb - segWSize; i < nb; i++) {
		measures[i] = s;
	}

	delete(fAccG1);
	delete(fAccG2);
	delete(fAccG12);

	return measures;
}

/// \brief select the true border from the array of similarities.
/// \return the borders
Borders doBorders(Measures & measures, Config & config)
{
	unsigned long size = measures.size();
	Borders borders;
	borders[0] = 0.0;
	long j = 0;
	borders[measures.size() - 1] = 0.0;
	unsigned long segMinWSize = config.getIntegerParam("segMinWSize");
	double thr = config.getFloatParam("segThr");
	
	unsigned long i = segMinWSize - 1;
	while(i < size) {
		double curr = measures[i];
		unsigned long start = MAX(0, i - segMinWSize);
		unsigned long end = MIN(size, i + segMinWSize);
		double max = measures[start];
		for (unsigned long m = start + 1; m < end; m++) {
			double v = measures[m];
			if ((i != m) && (v > max)) {
				max = v;
			}
		}
		if ((curr > max) && (curr > thr)) {

			borders[i] = curr;
			i += segMinWSize;
			j++;
		} else {
			i++;
		}
	}
	return borders;
}

/// \brief add Borders to Clusters.
/// \return a copy of \e clusters.
void doClusters(Borders & borders, const Segment & seg,
					ClusterServer & clusters, Config & config, String show)
{

	Borders::iterator itPrev = borders.begin();
	Borders::iterator itCur = borders.begin();
	itCur++;
						
	clusters.setIdxSrcName(show);
	long begin = seg.getBegin();
	long idxShow = clusters.getIdxSrcName(show);
	long i=clusters.getMaxIdxClusterName();
				
						
	for (; itCur != borders.end(); itCur++) {
		
		String name;
		name="S"+String::valueOf(i);
		i++;
		
		long st = begin + itPrev->first;
		long ln = itCur->first - itPrev->first;
		
		//Create a new locutor name, and get the corresponding id (index of the clusterMap)
		long clusterId = clusters.setIdxClusterName(name);		
		Cluster c;
		
		Segment newSeg;
		newSeg.setIdxSrcName(idxShow);
		
		newSeg.setBegin(st);
		newSeg.setLength(ln);
		
		// we add it to the list of the cluster
		c.insertSeg(newSeg);		
		
		clusters.addCluster(clusterId,c);	
		
		itPrev = itCur;
	}

//	return clusters;
}

/* seg
 */
int main(int argc, char* argv[])
{
	
	try {
		//permit to check the configuration parameters
		#define MANDATORY true
		#define OPTIONAL false
		#define ARG_REQUIRED true
		#define ARG_IS_OPTIONAL false
		
		ConfigChecker cc;
		//description of the parameters list
		//--fInMask=path.extension (path.mfcc)
		cc.addStringParam("featureFilesPath", MANDATORY, ARG_IS_OPTIONAL);
		cc.addStringParam("loadFeatureFileExtension", MANDATORY, ARG_IS_OPTIONAL);
		
		//--fdesc="spro4,1:1:0:0:0:0,13,0:0:0
		cc.addStringParam("loadFeatureFileFormat", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("featureServerMask", MANDATORY, ARG_REQUIRED);
		//cc.addIntegerParam("loadFeatureFileVectSize", MANDATORY, ARG_REQUIRED);
		//cc.addIntegerParam("vectSize", MANDATORY, ARG_REQUIRED);
		
		//--sInMask=path.uem.seg (sph/myFile.uem.seg)
		cc.addStringParam("segServerFilesPath", MANDATORY, ARG_REQUIRED); //+output -> idem
		cc.addStringParam("loadSegServerFileExtension", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("extSegClustIn", MANDATORY, ARG_REQUIRED);
		
		//--sInFormat=0 (seg)
		cc.addStringParam("loadSegServerFileFormat", MANDATORY, ARG_REQUIRED);//RAW only ? --> TODO .seg
		
		//--sOutMask=path.uem.seg (sph/myFile.uem.seg)
		cc.addStringParam("saveSegServerFileExtension", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("extSegClustOut", MANDATORY, ARG_REQUIRED);
		
		//--sOutFormat=0 (seg)
		cc.addStringParam("saveSegServerFileFormat", MANDATORY, ARG_REQUIRED); //RAW or XML only ? --> TODO .seg
		
		//--kind (full/diag)
		cc.addStringParam("distribType", MANDATORY, ARG_REQUIRED); //GF or GD
		
		//--segWSize
		cc.addStringParam("segWSize", MANDATORY, ARG_REQUIRED); //Default : 250
		
		//--sMethod (GLR)
		cc.addStringParam("segMethod", MANDATORY, ARG_REQUIRED);
		
		//--show
		cc.addStringParam("show", MANDATORY, ARG_REQUIRED);
		
		//get the line command
		CmdLine cmdLine(argc, argv);
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "seg 1.0" << endl
			<< "Usage: seg [option]..." << endl << endl
			<< "       A measure based segmentation software." << endl 
			<< "       The segmentation process consists to find the instantaneous " << endl
			<< "       speaker change points corresponding to segment boundaries. " << endl
			<< "       The proposed technique is to detect the change point in a " << endl
			<< "       sliding segment over the whole signal. Given a two sliding " << endl
			<< "       segments, a measure is computed between segments. A change " << endl
			<< "       point is present when the measure exceed a threshold. " << endl
			<<  endl
			<<cc.getParamList();
		}
		else {
			if (cmdLine.displayVersionRequired()) {
				cout << "Version: 1.0" << endl;
			}

			else {
				//define the configuration
				Config config;
				cmdLine.copyIntoConfig(config);

				//check all the arguments : if one is missing, program stop
				cc.check(config);
				
				if (config.getParam_debug()) {
					cout << "mode debug" << endl;
				}
      
				// performs job
				//create segments and cluster
				ClusterServer inputClusterServer(config);
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				inputClusterServer.read(filename);
				
				cout << endl << "-------------- Start seg " << config.getParam("show") << " -----------" << endl << endl;

				#ifdef DEBUG
				cout << "debug mode" << endl << endl;
				#else
				cout << "release mode" << endl << endl;
				#endif

				FeatureServer fs(config, inputClusterServer.getSourceNameInXLine());
				ClusterServer outputClusterServer(config);
				
				ClusterMap& inputClusterMap = inputClusterServer.getContainer();
				ClusterMap::iterator inputIterCluster;
				
				for(inputIterCluster=inputClusterMap.begin();inputIterCluster != inputClusterMap.end();
					inputIterCluster++) {
					Cluster inputCluster = inputIterCluster->second;
						
					for(inputCluster.seekBeginSeg(); inputCluster.isEof() == false; inputCluster.seekNextSeg()){
						const Segment & inSeg = inputCluster.getCurrentSeg();
						 
						time_t rawtime;
						time ( &rawtime );
						printf ( "trace[time] \t Start: %s", ctime (&rawtime) );
						 
						Measures m = doMeasures(fs,inSeg,inputClusterServer,config);
						time ( &rawtime );
						printf ( "trace[time] \t After doMeasures: %s", ctime (&rawtime) );
						 
						Borders b = doBorders(m, config);
						time ( &rawtime );
						printf ( "trace[time] \t After doBorders: %s", ctime (&rawtime) );
						
						Segment seg = inSeg;
						
						String currentShow = inputClusterServer.getSourceName(seg.getIdxSrcName());
							
						doClusters(b, seg, outputClusterServer, config,currentShow);
						time ( &rawtime );
						printf ( "trace[time] \t After doClusters: %s", ctime (&rawtime) );						
					}
					
				}
					
				//save into a file in format .seg
				//we try to write
				filename = config.getParam("show")+config.getParam("extSegClustOut");
				outputClusterServer.write(filename);
				cout << endl << "-------------- End seg -----------" << endl << endl;
			
				


			}
		}
	}//end of try
	catch (ConfigCheckException& e) {
		cout << e.msg << endl
			<< "Try test --help for more informations" << endl;
	} 
	catch (Exception& e) {
		cout << e.toString() << endl;
	}

	return 0;
}
