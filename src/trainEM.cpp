/*!
\file trainEM.cpp
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
#include "clusterServer.h"
#include "liaToolsWrapper.h"

using namespace lium_seg;
using namespace alize;

DoubleVector computeGlobalCov(ClusterServer & clusters, FeatureServer & features)
{
	FrameAccGD accu; 
	ClusterMap::iterator it = clusters.getContainer().begin();
	while(it != clusters.getContainer().end()) {
		Cluster & clust = it->second;
		//we get the segments
		for(clust.seekBeginSeg(); clust.isEof() == false; clust.seekNextSeg()) {
			const Segment & seg = clust.getCurrentSeg();
			String sourceName = clusters.getSourceName(seg.getIdxSrcName());
			unsigned long position = features.getFirstFeatureIndexOfASource(sourceName);
			unsigned long start = seg.getBegin() + position; // begin of the segment
			unsigned long end = start + seg.getLength();
			Feature f;
			for (unsigned long j = start; j<end; j++) {
				features.seekFeature(j);
				features.readFeature(f);
				accu.accumulate(f);
			}
		}
		it++;
	}
	return accu.getCovVect();
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

		//--tInMask=./$datadir/%s.init.gmms
		cc.addStringParam("loadMixtureFileExtension", MANDATORY, ARG_REQUIRED);

		//--minLLK
		cc.addStringParam("minLLK", MANDATORY, ARG_REQUIRED);

		//--maxLLK
		cc.addStringParam("maxLLK", MANDATORY, ARG_REQUIRED);

		//gain for EM
		cc.addFloatParam("gain", MANDATORY, ARG_REQUIRED);

		//max iteration for EM
		cc.addIntegerParam("maxIt", MANDATORY, ARG_REQUIRED);

		//min iteration for EM
		cc.addIntegerParam("minIt", MANDATORY, ARG_REQUIRED);

		//variance Flooring for EM
//		cc.addFloatParam("varFloor", MANDATORY, ARG_IS_OPTIONAL);
		
		//variance celling for EM
//		cc.addFloatParam("varCel", MANDATORY, ARG_IS_OPTIONAL);

		//normalize model for EM
//		cc.addStringParam("nomModel", MANDATORY, ARG_IS_OPTIONAL);

		//get the line command
		CmdLine cmdLine(argc, argv);
		
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "trainEM 1.0" << endl
			<< "Usage: trainEM [option]..." << endl << endl
			<< "       This software is used to train a set of GMM using EM " << endl
			<<"        algorithm (previously initialized with trainInit software)." << endl << endl
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
				cout << endl << "-------------- Start trainEM -----------" << endl << endl;
				// clusters
				//create segments and cluster
				ClusterServer clusters(config);
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				clusters.read(filename);

				// Features
				FeatureServer features(config, clusters.getSourceNameInXLine());
				DoubleVector globalCov = computeGlobalCov(clusters, features);
				//Compute Model
				//GMMVect initVect = MainTools::getInputGMMVect(param);;

				MixtureServer mxServerInit(clusters.getConfig());
				cerr<<"Filename : "<<filename <<endl;
				mxServerInit.reset();
				mxServerInit.load(config.getParam("show"));
				
				//we must have the same number of GMM than clusters
				if (mxServerInit.getMixtureCount() != clusters.getContainer().size()) {
					cerr<<"!!! error \t initial model number is not good "<<endl;
					cerr<<"Nb Mixtures : "<< mxServerInit.getMixtureCount() <<endl;
					cerr<<"Nb Clusters : "<< clusters.getContainer().size() <<endl;

					exit(0);
				}
				
				//Compute Model
				StatServer statServer(config);
				//training
				int maxIt = config.getIntegerParam("maxIt");
				int minIt = config.getIntegerParam("minIt");
				float gain = config.getFloatParam("gain");
				bool diag = (config.getParam("distribType") == "diagonal" ? true : false);
				ClusterMap::iterator it = clusters.getContainer().begin();
				while(it != clusters.getContainer().end()) {
					//we do it several times in order to gain in precision
					Cluster & cluster = it->second;
					long idMxInit = mxServerInit.getMixtureIndex(String::valueOf(it->first));
					Mixture& mixture = mxServerInit.getMixture(idMxInit);

					//the first mixture, to do the tests
					//create a MixtureStat in order to calculate the EM
					MixtureStat &emAcc=statServer.createAndStoreMixtureStat(mixture);
					emAcc.resetEM();
					mixture.computeAll();
					//we accumulate EM for every segments in the cluster
					double old = 0;
					double cur = 0;

					for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
						const Segment & seg = cluster.getCurrentSeg();

						//we check the name of the source
						String sourceName = clusters.getSourceName(seg.getIdxSrcName());
						//we accumulate the stats
						old += LiaToolsWrapper::accStatEM(statServer, features,emAcc,seg,sourceName,config);
					}
					mixture = emAcc.getEM();
					if (diag) 
						LiaToolsWrapper::checkModel(static_cast<MixtureGD&>(mixture), globalCov, config);						
					//now we try to be more accurate
					for (int i = 1;i<maxIt;i++)
					{
						//create a MixtureStat in order to calculate the EM
						MixtureStat &emAccNew=statServer.createAndStoreMixtureStat(mixture);
						emAccNew.resetEM();
						mixture.computeAll();

						cur=0;
						//we get the segments
						for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
							const Segment & seg = cluster.getCurrentSeg();

							//we check the name of the source
							String sourceName = clusters.getSourceName(seg.getIdxSrcName());
							//we accumulate the stats
							cur += LiaToolsWrapper::accStatEM(statServer, features,emAccNew,seg,sourceName,config);
						}
						double g = cur - old;
						mixture = emAccNew.getEM();
						if (diag) 
							LiaToolsWrapper::checkModel(static_cast<MixtureGD&>(mixture), globalCov, config);						
						//we the result satisfies us, we don't continue
						if ((i >=minIt) && (g < gain)) {
							break;
						}
						old = cur;
					}
					it++;
				}
				mxServerInit.save(config.getParam("show"));
			}
		}
		cout << endl << "-------------- End trainEM -----------" << endl << endl;
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
