/*! 
\file clust.cpp
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

/*! \mainpage &nbsp;
  \htmlinclude README

 */

#include "clusterServer.h"
#include "biclclust.h"
#include "bichclust.h"
#include "clrhclust.h"
#include <cstdlib>
#include "float.h"

using namespace lium_seg;

#define MIN(x,y) ((x)<(y))?(x):(y)
#define MAX(x,y) ((x)>(y))?(x):(y)

/// \brief BIC linear clustering
/// \return a Clusters
ClusterServer lclust(ClusterServer & clusters, FeatureServer & features)
{
	Config & config = clusters.getConfig();
	real_t threshold = config.getParam("thr").toDouble();
	String kind = config.getParam("distribType");
	BICLClust clust(clusters, features, threshold, kind);
	long ci = 0;
	long cj = 1;
	real_t score;
	clust.init();
	clust.next(score, ci, cj);
	while (score < DBL_MAX) {
		cerr << "trace[lclust] \t score = " << score;
		cerr << " ci = " << ci << "(" << clust.getClusterName(ci) << ")";
		cerr << " cj = " << cj << "(" << clust.getClusterName(cj) << ")" << endl;
		if (score < 0.0) {
			clust.update(ci, cj);
		} else {
			ci++;
			cj++;
		}
		clust.next(score, ci, cj);
	}
	return clust.getContainer();
}

/// \brief BIC Hierarchical clustering
/// \return a Clusters
ClusterServer hclust(ClusterServer & clusters, FeatureServer & features)
{
	cerr<<"[clust] \t begin function hclust"<<endl;
	Config & config = clusters.getConfig();
	real_t threshold = config.getParam("thr").toDouble();
	String kind = config.getParam("distribType");

	BICHClust clust(clusters, features, threshold, kind);
	long nb = 0;
	long ci, cj;
	real_t score;
	String strMergeMax = config.getParam("mergeMax");
	long mergeMax = strMergeMax.toLong();
	clust.init();
	clust.next(score, ci, cj);
	while ((score < 0.0) && (nb < mergeMax)) {
		cerr << "trace[hclust] \t merge = " << nb;
		cerr << " score = " << -1.0 * score;
		cerr << " ci = " << ci << "(" << clust.getClusterName(ci) << ")";
		cerr << " cj = " << cj << "(" << clust.getClusterName(cj) << ")" << endl;
		clust.update(ci, cj);
		clust.next(score, ci, cj);
		nb++;
	}
	return clust.getContainer();
}

/// \brief CLR Hierarchical clustering
/// \return a Clusters
ClusterServer cclust(ClusterServer & clusters, FeatureServer & features, bool nclr = false)
{
	Config & config = clusters.getConfig();
	real_t threshold = config.getParam("thr").toDouble();
	String kind = config.getParam("distribType");
	long mergeMax = config.getIntegerParam("mergeMax");
	long minSpk = config.getIntegerParam("minSpk");
	CLRHClust clust(clusters, features, threshold, kind, nclr);

	long nb = 0;
	long ci, cj;
	real_t score;
	clust.init();
	clust.next(score, ci, cj);
	double old_score = score;
	long nbSpk = clusters.getContainer().size();
	while ((score < -1.0 * threshold) && (nb < mergeMax) && (nbSpk > minSpk)) {
		cerr << "trace[cclust] \t merge = " << nb;
		cerr << " score = " << old_score << " thr=" << threshold;
		cerr << " ci = " << ci << "(" << clusters.getSourceName(ci) << ")";
		cerr << " cj = " << cj << "(" << clusters.getSourceName(cj) << ")" << endl;
		old_score = score;
		clust.update(ci, cj);
		clust.next(score, ci, cj);
		nb++;
		nbSpk = clust.getContainer().getContainer().size();
	}
	cerr << " score_old = " << old_score << " thr=" << threshold << endl;
	cerr << " score = " << score << " thr=" << threshold  << endl;
	cerr << " nbSpk = " << minSpk << endl;
	cerr << " nb = " << mergeMax << endl;
	clust.reset();
	return clust.getContainer();
}

//Surement a mettre dans une autre classe de normalisation
//On le laisse ici pour l'instant, voir si ca marche
void computeNorm(ClusterServer & clusterServer, const Segment & seg, FeatureServer & features) {

	unsigned long start = seg.getBegin();
	//we check the name of the source
	String sourceNameTemp = clusterServer.getSourceName(seg.getIdxSrcName());
	//in order to manage multiple stream source, we add to the begin and the end
	//the position in feature of the source
	unsigned long position = features.getFirstFeatureIndexOfASource(sourceNameTemp);	
	start += position;
	unsigned long end = MIN(seg.getBegin() + seg.getLength(), features.getFeatureCount());
	
	FrameAccGD frameAcc;
	frameAcc.reset();
	
	//On accumule de quoi normaliser
	for(unsigned int i = start; i < end; i++) {
		features.seekFeature(i);
		Feature f;
		features.readFeature(f);
		frameAcc.accumulate(f);
	}	
	
	bool reduce = true;
	//On normalise chaque trame
	for(unsigned int i = start; i < end; i++) {
		features.seekFeature(i);
		Feature f;
		features.readFeature(f);
		for (unsigned long j = 0; j < features.getVectSize(); j++) {
			f[j] -= frameAcc.getMeanVect()[j];

			if (reduce) {
				f[j] /= sqrt(frameAcc.getCovVect()[j]);
			}
		}
		features.seekFeature(i);
		features.writeFeature(f);
	}
}

void normBySeg(ClusterServer & clusters, FeatureServer & features) {
	ClusterMap& clusterMap = clusters.getContainer();				
	ClusterMap::iterator iterMap = clusterMap.begin();
	for(;iterMap != clusterMap.end(); iterMap++) {
		//we get the current cluster
		Cluster cluster = iterMap->second;
					
		//we get the segments
		for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
			const Segment & seg = cluster.getCurrentSeg();
			computeNorm(clusters,seg,features);
		}
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
		
		//--sOutMask=path.uem.seg (sph/myFile.seg)
		cc.addStringParam("saveSegServerFileExtension", MANDATORY, ARG_REQUIRED);
		cc.addStringParam("extSegClustOut", MANDATORY, ARG_REQUIRED);
		
		//--tInMask=./gmm_ubm/ubm.gmm
		cc.addStringParam("loadMixtureFileExtension", OPTIONAL, ARG_REQUIRED);
		cc.addStringParam("mixtureFilesPath", OPTIONAL, ARG_REQUIRED);
		cc.addStringParam("world", OPTIONAL, ARG_REQUIRED);
		
		//--sOutFormat=0 (seg)
		cc.addStringParam("saveSegServerFileFormat", MANDATORY, ARG_REQUIRED); //RAW or XML only ? --> TODO .seg

		//method
		cc.addStringParam("method", MANDATORY, ARG_REQUIRED);
		
		//threshold
		cc.addFloatParam("thr", MANDATORY, ARG_REQUIRED);

		//kind (full or diagonal)
		cc.addStringParam("distribType", MANDATORY, ARG_REQUIRED);

		//Nb max of merge (for hierarchical)
		cc.addIntegerParam("mergeMax", MANDATORY, ARG_REQUIRED);

		//min speaker
		cc.addIntegerParam("minSpk", OPTIONAL, ARG_REQUIRED);

		//top
		cc.addIntegerParam("topDistribsCount", OPTIONAL, ARG_REQUIRED);

		//nbTrainIt
		cc.addIntegerParam("nbTrainIt", OPTIONAL, ARG_REQUIRED);
		
		//get the line command
		CmdLine cmdLine(argc, argv);
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "clust 1.0" << endl
			<< "Usage: clust [option]..." << endl << endl
			<< "       BIC based clustering software used a hierarchical bottom-up algorithm." << endl << endl
			<< cc.getParamList();
		} else {
			if (cmdLine.displayVersionRequired()) {
				cout << "Version: 1.0" << endl;
			} else {
				//define the configuration
				Config config;
				cmdLine.copyIntoConfig(config);

				//check all the arguments : if one is missing, program stop
				cc.check(config);
				
				if (config.getParam_debug()) {
					cout << "mode debug" << endl;
				}
				
				// performs job
				cout << endl << "-------------- Start clust -----------" << endl << endl;
				//create segments and cluster
				ClusterServer clusters(config);
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				clusters.read(filename);

				// Features
				FeatureServer features(config, clusters.getSourceNameInXLine());			
						
				// methods
				ClusterServer clustersRes(config);
				if (config.getParam("method") == /*Param::CLUST_H_BIC*/"h") {
					cerr << "[hclust] \t begin clust" << endl;
					clustersRes = hclust(clusters, features);
					cerr << "[hclust] \t end clust" << endl;
				} else if (config.getParam("method") == /*Param::CLUST_L_BIC*/"l") {
					cerr << "[lclust] \t begin clust" << endl;
					clustersRes = lclust(clusters, features);
					cerr << "[lclust] \t end clust" << endl;
				} else if (config.getParam("method") == /*Param::CLUST_H_CLR*/"c") {
					cerr << "[clrhclust] \t begin clust" << endl;
					cerr << "[clrhclust] \t norm feature" << endl;
					normBySeg(clusters,features);
					cerr << "[clrhclust] \t clustering" << endl;
					clustersRes = cclust(clusters, features);
					cerr << "[clrhclust] \t end clust" << endl;
				} else if (config.getParam("method") == /*Param::CLUST_H_NCLR*/"nclr") {
					cerr << "[clrhclust] \t begin clust" << endl;
					cerr << "[clrhclust] \t norm feature" << endl;
					normBySeg(clusters,features);
					cerr << "[clrhclust] \t clustering" << endl;
					clustersRes = cclust(clusters, features, true);
					cerr << "[clrhclust] \t end clust" << endl;
				} else {
					cerr << "error \t unknown method" << endl;
					return -1;
				}
				//write the clusters in the file
				filename = config.getParam("show")+config.getParam("extSegClustOut");
				clustersRes.write(filename);
				cout << endl << "-------------- End clust -----------" << endl << endl;
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
