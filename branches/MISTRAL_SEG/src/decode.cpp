/*! 
\file decode.cpp

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

\date 12 Oct 2006
\version 0.0.1

\brief EM trainer program for the GMMs

*/

#include <vector>

#include "alize.h"
#include "clusterServer.h"
#include "liaToolsWrapper.h"

using namespace lium_seg;
using namespace alize;


ClusterServer getClusters(const ClusterServer & clusters, MixtureServer& mxServer, ViterbiAccum& va, SegLst& segLst,Config& config)
{
cerr << "get cluster 1" << endl;
	ULongVector path = va.getPath();
	unsigned long pathSize = path.size();
	
	ClusterServer res(config,clusters.getSrcNameMap(), clusters.getMaxIdxSrcName(), 
	clusters.getClusterNameMap(), clusters.getMaxIdxClusterName());
	
	map<long, SegLst::iterator> m;
	SegLst::iterator its = segLst.begin();
	for(; its != segLst.end(); its++){
		long end = its->getBegin()+its->getLength();
		for(long i = its->getBegin(); i < end; i++){
			m[i] = its;
		}
	}
cerr << "get cluster 2" << endl;
	
	for(unsigned long i=0;i<pathSize;) {
//		cerr << "get cluster 2 i=" << i << endl;
		unsigned long mixtureIndex = path[i];
		
		Cluster & c = res.getContainer()[mixtureIndex];
		
		c.setId(mixtureIndex);
		Segment s = *(m[i]);
		s.setBegin(i);
		unsigned long len = 1;
		s.setLength(len);
		c.insertSeg(s);
		i += len;
	}
cerr << "colapse" << endl;
	res.collapse();
	return res;
}

// A function for create and initialize a Viterbi accumulator
ViterbiAccum& createAndInitializeViterbiAccum(StatServer &ss, hmm &cHmm, double penality){
  ViterbiAccum& va=ss.createViterbiAccum(); // CREATE 
  for(unsigned int i=0;i<cHmm.getNbState();i++) // Add the states
    va.addState(cHmm.getDensity(i));
  for(unsigned int i=0;i<cHmm.getNbState();i++) // Set the transitions
    for(unsigned int j=0;j<cHmm.getNbState();j++){
      double p = 0.0;
      if (i != j) {
      	p = penality;
      }
	  cHmm.setTransition(p, i,j);
      va.logTransition(i,j) = p;
    }
  return va;
}


// accumulateStatViterbi() is used for computing the viterbi path
// TAKE CARE: THE accumulator should be reseted before the call
void   accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,unsigned long beginIdx,unsigned long length,Config &config)
{	
	for (unsigned int i=beginIdx;i<beginIdx+length;i++) {
		Feature f;
		fs.seekFeature(i);
		fs.readFeature(f);
		va.computeAndAccumulate(f, 1.0, 0.0);
	}

}
// On a segment
void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,const Segment &seg,Config &config, 
						   ClusterServer& clusterServer)
{  
	String sourceNameTemp = clusterServer.getSourceName(seg.getIdxSrcName());

	//in order to manage multiple stream source, we add to the begin and the end
	//the position in feature of the source
	unsigned long position = fs.getFirstFeatureIndexOfASource(sourceNameTemp);	
  	unsigned long begin=seg.getBegin()+position;              // Find the index of the first frame of the file in the buffer
  								   
	accumulateStatViterbi(fs,va,begin,seg.getLength(),config);
}

/// END FROM LIA_SpkTools/SegTools

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
		
		cc.addStringParam("extSegDecodeOut", MANDATORY, ARG_REQUIRED);
		
		//--sInFormat=0 (seg)
		cc.addStringParam("loadSegServerFileFormat", MANDATORY, ARG_REQUIRED);//RAW only ? --> TODO .seg

		
		//--tInMask=./$datadir/%s.gmms
		cc.addStringParam("loadMixtureFileExtension", MANDATORY, ARG_REQUIRED);

		//HMM Penality I-->J
		cc.addStringParam("HMMPenalityIJ", MANDATORY, ARG_REQUIRED);

		//get the line command
		CmdLine cmdLine(argc, argv);
		
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "decode 1.0" << endl
			<< "Usage: decode [option]..." << endl << endl
			<< "       This software generates a new segmentation based on Viterbi " << endl
			<< "       algorithm using a set of GMMs. Each state of the HMM corresponds " << endl
			<< "       to a GMM (a speaker cluster). " << endl << endl
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
				cout << endl << "-------------- Start decode -----------" << endl << endl;
				// clusters
				//create segments and cluster
				ClusterServer inputClusterServer(config);
				ClusterServer outputClusterServer;
				
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				inputClusterServer.read(filename);
				// Features
				FeatureServer features(config, inputClusterServer.getSourceNameInXLine());
				
				MixtureServer mxServerInit(config);
				cerr<<"Filename : "<<filename <<endl;
				mxServerInit.reset();
				mxServerInit.load(config.getParam("show"));
				
				cerr << "Get the current HMM" << endl;
				// Get the current HMM
  				hmm hmmactuel(mxServerInit,config);
				StatServer ss(config, mxServerInit);
								
				int nbMixtureCount = mxServerInit.getMixtureCount();
				ClusterMap::iterator inputIterCluster;
				
				for(int i =0; i < nbMixtureCount;i++) {
				
					MixtureGD& mixtureGD = mxServerInit.getMixtureGD(i);
					mixtureGD.computeAll();
						
					hmmactuel.LoadState(mixtureGD,String::valueOf(i));
					
				}
				
				cerr << "Initialize the viterbi accumulator" << endl;
				
				//Initialize the viterbi accumulator
				ViterbiAccum& va=createAndInitializeViterbiAccum(ss, hmmactuel, config.getParam("HMMPenalityIJ").toDouble());
				
				SegLst segLstDecode = inputClusterServer.getSegLst(false);
				SegLst::iterator its = segLstDecode.begin();
				for(; its != segLstDecode.end(); its++){
					Segment s = *its;
					cerr << "--- new seg, start=" << s.getBegin() << " " << s.getLength() << endl;
					accumulateStatViterbi(features,va,s,config,inputClusterServer);
				}
				
				outputClusterServer = getClusters(inputClusterServer,mxServerInit,va,segLstDecode,config);
				filename = config.getParam("show")+config.getParam("extSegDecodeOut");
				outputClusterServer.write(filename);
				cout << endl << "-------------- End decode toto-----------" << endl << endl;
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
