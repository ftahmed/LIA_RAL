/*!
\file trainInit.cpp
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

#include <vector>

#include "alize.h"
#include "clusterServer.h"

using namespace lium_seg;
using namespace alize;

//Sera surement a deplacer dans Mixture ...
void initMixture(Mixture& mixture, Config& config, Cluster & cluster, ClusterServer & clusterServer,
				 FeatureServer & fs, String kind) {
	
	unsigned long distribCount = mixture.getDistribCount(); 
	//cerr << "Distrib count : " << distribCount << endl;
  	FrameAcc** frameAcc=new FrameAcc* [distribCount]; 
	
	for (unsigned long i=0;i<distribCount;i++) {        
    	frameAcc[i]= new FrameAccGD();                           // Create the accumulators
    	frameAcc[i]->reset();                                    // Reset it
    	
		//For each accumulator, add the features
		for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
			const Segment & seg = cluster.getCurrentSeg();
			//we check the name of the source
			String sourceNameTemp = clusterServer.getSourceName(seg.getIdxSrcName());

			//in order to manage multiple stream source, we add to the begin and the end
			//the position in feature of the source
			unsigned long position = fs.getFirstFeatureIndexOfASource(sourceNameTemp);
			
			long s = seg.getBegin();
			s+=position;
			
			long e = s + seg.getLength();

			for (long j = s; j < e; j++) {
				fs.seekFeature(j);
				Feature f;
				fs.readFeature(f);
				frameAcc[i]->accumulate(f);
			}
		}
  	}
	
	long nbComp = distribCount;
	long len = 5;
	double weight = 1.0 / nbComp;
					 
	long nb = 0;
	vector<long> pFreatures;
	vector<long> idxShow;
	for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
		const Segment & seg = cluster.getCurrentSeg();
		//we check the name of the source
		String sourceNameTemp = clusterServer.getSourceName(seg.getIdxSrcName());

		//in order to manage multiple stream source, we add to the begin and the end
		//the position in feature of the source
		unsigned long position = fs.getFirstFeatureIndexOfASource(sourceNameTemp);
		
		nb += seg.getLength();
		long s = seg.getBegin();
		s+=position;
		long e = s + seg.getLength();
		for(long i = s; i < e; i++){
			pFreatures.push_back(i);
		}
	}
					 
	//std::cerr << "trace[uniforme] \t pFreatures size : " << pFreatures.size() << std::endl;

	long step = (nb - 1) / nbComp;
	long dim = fs.getVectSize();
	for (long i = 0; i < nbComp; i++) {
		Distrib* g1 = mixture.getTabDistrib()[i];
		for (long k = 0; k < dim; k++) {
			double s = 0.0;
			for (long j = 0; j < len; j++) {
				long idx = pFreatures[i * step + j];
				Feature feature;
				fs.seekFeature(idx);
				fs.readFeature(feature);
				s += feature[k];
			}
			double v = s / (double) len;
			g1->setMean(v,k);
			
			if(kind == "GD") {
				DistribGD* g1D = (DistribGD*)g1;
				FrameAccGD* fAccGD = (FrameAccGD*)frameAcc[i];
				g1D->setCov(fAccGD->getCovVect()[k],k);
				g1D->computeAll();
			}
		}
		mixture.getTabWeight()[i]=weight;
	}
	
}

int main(int argc, char* argv[])
{
	try {
		//we get all the params
		//permit to check the configuration parameters
		#define MANDATORY true
		#define OPTIONAL false
		#define ARG_REQUIRED true
		#define ARG_IS_OPTIONAL false
		
		ConfigChecker cc;

		//kind (full or diagonal)
		cc.addStringParam("distribType", MANDATORY, ARG_REQUIRED);

		//nbComp / nb gaussiennes
		cc.addStringParam("mixtureDistribCount", MANDATORY, ARG_REQUIRED);

		//--fInMask=path.extension (path.mfcc)
		cc.addStringParam("featureFilesPath", MANDATORY, ARG_IS_OPTIONAL);
		cc.addStringParam("loadFeatureFileExtension", MANDATORY, ARG_IS_OPTIONAL);
		cc.addStringParam("show", MANDATORY, ARG_REQUIRED);

		//--fdesc="spro4,1:1:0:0:0:0,13,0:0:0
		cc.addStringParam("loadFeatureFileFormat", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("featureServerMask", MANDATORY, ARG_REQUIRED);
		cc.addIntegerParam("loadFeatureFileVectSize", MANDATORY, ARG_REQUIRED);
		cc.addIntegerParam("vectSize", MANDATORY, ARG_REQUIRED);

		//--sInMask=path.uem.seg (sph/myFile.uem.seg)
		cc.addStringParam("segServerFilesPath", MANDATORY, ARG_REQUIRED); //+output -> idem
		cc.addStringParam("loadSegServerFileExtension", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("extSegClustIn", MANDATORY, ARG_REQUIRED);
		
		//--sInFormat=0 (seg)
		cc.addStringParam("loadSegServerFileFormat", MANDATORY, ARG_REQUIRED);//RAW only ? --> TODO .seg

		//--tOutMask=./$datadir/%s.gmms
		cc.addStringParam("mixtureFilesPath", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("saveMixtureFileExtension", MANDATORY, ARG_REQUIRED);

		//method
		cc.addStringParam("dMethod", MANDATORY, ARG_REQUIRED);

		//get the line command
		CmdLine cmdLine(argc, argv);
		
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "trainInit 1.0" << endl
			<< "Usage: trainInit [option]..." << endl << endl
			<< "       This software is useful to initialize a set of GMM."  << endl << endl
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
				cout << endl << "-------------- Start trainInit -----------" << endl << endl;
				// clusters
				//create segments and cluster
				ClusterServer inputClusterServer(config);
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				inputClusterServer.read(filename);

				// Features
				FeatureServer features(config, inputClusterServer.getSourceNameInXLine());

				ClusterMap& inputClusterMap = inputClusterServer.getContainer();
				ClusterMap::iterator inputIterCluster;
				
				MixtureServer ms(config);
				
				long mixtureDistribCount=config.getIntegerParam("mixtureDistribCount");
				
				for(inputIterCluster=inputClusterMap.begin();inputIterCluster != inputClusterMap.end();
					inputIterCluster++) {
					Cluster inputCluster = inputIterCluster->second;
					
					DistribType distribType=DistribType_GD;
					if(config.getParam("distribType")=="GF")
						distribType=DistribType_GF;
						
					Mixture& gmm = ms.createMixture(mixtureDistribCount,distribType);
					
					initMixture(gmm,config,inputCluster,inputClusterServer,features,config.getParam("distribType"));
						
					String clusterId = String::valueOf(inputCluster.getId());
					gmm.setId(clusterId);
					//A voir si on a besoin de mettre le genre ?

				}
				
				ms.save(config.getParam("show"));
			}
		}
				cout << endl << "-------------- End trainInit -----------" << endl << endl;
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
