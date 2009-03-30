/*!
\file segInit.cpp
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
#include "segment.h"
#include "cluster.h"
#include <math.h>

using namespace alize;

/* segInit
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
		
		//--sOutFormat=0 (seg)
		cc.addStringParam("saveSegServerFileFormat", MANDATORY, ARG_REQUIRED); //RAW or XML only ? --> TODO .seg
		
		//get the line command
		CmdLine cmdLine(argc, argv);
		//two differents cases : help or version display
		if (cmdLine.displayHelpRequired()) {
			cout << "segInit 1.0" << endl
			<< "Usage: segInit [option]..." << endl << endl
			<< "       Perform two safety checks on a given feature file. "
			<< endl << endl
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
				cout << endl << "-------------- Start segInit -----------" << endl << endl;

				#if !defined NDEBUG
				cout << "debug mode" << endl << endl;
				#else
				cout << "release mode" << endl << endl;
				#endif

				//create segments and cluster
				ClusterServer inputClusterServer(config);
				String filename = config.getParam("show")+config.getParam("extSegClustIn");
				inputClusterServer.read(filename);
				
				//To do : get the file list from the command line
				FeatureServer fs(config, inputClusterServer.getSourceNameInXLine());

				Feature currentFeature;
				Feature previousFeature;
				
				ClusterServer outputClusterServer(config);

				//we get all the segments extracted from the .seg file : clusters + segments
				//iteration
				ClusterMap& clusterMap = inputClusterServer.getContainer();
				cerr << "trace[mSegInit] \t Begin compare features..." << endl;
				
				ClusterMap::iterator iterMap = clusterMap.begin();
				for(;iterMap != clusterMap.end();iterMap++) {
					//we get the current cluster
					Cluster clust = iterMap->second;
					//we copy him into the new server
					Cluster resultSegCluster(clust.getId());
					
					//we add the name of locutor
					long id = resultSegCluster.getId();
					String locutorName = inputClusterServer.getName(id);

					outputClusterServer.addClusterName(id,locutorName);
					
					for(clust.seekBeginSeg(); clust.isEof() == false; clust.seekNextSeg()){
						const Segment & seg = clust.getCurrentSeg();
						
						//we check the name of the source
						String sourceName = inputClusterServer.getSourceName(seg.getIdxSrcName());
					
						// we add the source name
						SrcNameMap & mapOfSrc = outputClusterServer.getSrcNameMap();
						mapOfSrc[sourceName] = seg.getIdxSrcName();

						//in order to manage multiple stream source, we add to the begin and the end
						//the position in feature of the source
						unsigned long position = fs.getFirstFeatureIndexOfASource(sourceName);
						
						//we get the minimum size between the feature and the seg
						unsigned long start = seg.getBegin()+position; // begin of the segment
						unsigned long end = start+seg.getLength();
						cerr << "trace[mSegInit] \t src name : " << sourceName << " size file : " << end-position << " size feature : "<< fs.getFeatureCountOfASource(sourceName) << endl;
						unsigned long mini = min(end, (fs.getFeatureCountOfASource(sourceName)+position));
						
						if (end > mini) {
							cerr << "WARNING[seginit] \t segment end upper to features end" << endl;
						}
						
						//Start with the second feature because we compare the current feature with the previous
						for (unsigned long k = start+1; k < mini; k++) {
							fs.seekFeature(k);
							fs.readFeature(currentFeature);
							fs.seekFeature(k-1);
							fs.readFeature(previousFeature);

							//Use the overloaded == operator to compare two features
							//Doc for == : Two Feature objects are equal if their dimensions are equal and their validities are equal and 
							//	       all their parameters (vector) are the same.
							if(currentFeature != previousFeature) {
								//The two features are not the same
								//add them to the SegServer
						
								//we create a segment of 1 long
								Segment newSeg;
								newSeg.setIdxSrcName(seg.getIdxSrcName());
									
								//the same segment, except that we don't have the same begin and length
								//case of the first element
								if (k == (start+1)) {
									newSeg.setBegin(k-1-position);
									newSeg.setLength(1);
									resultSegCluster.insertSeg(newSeg);
								}
								newSeg.setBegin(k-position);
								newSeg.setLength(1);
							
								// we add it to the list of the cluster
								resultSegCluster.insertSeg(newSeg);
							}
							else {
								//The two features are the same, there is certainly a problem in the file
								//Skip the current feature
								cerr << "WARNING[mSeginit] \t feature equal : " << k <<endl;
							}
						}
					}//end of iteration on segments
					
					//we add the cluster in the map
					outputClusterServer.addCluster(id,resultSegCluster);
				}//end of iteration on clusters
				
				cerr << "trace[mSegInit] \t ...End compare features" << endl;
				cerr << "trace[mSegInit] \t Collapse segments in progress..."<<endl;
				outputClusterServer.collapse();
				cerr << "trace[mSegInit] \t ... End collapse segments"<<endl;
				
				//save into a file in format .seg
				//we try to write
				filename = config.getParam("show")+config.getParam("extSegClustOut");
				outputClusterServer.write(filename);


				cout << endl << "-------------- End mSegInit -------------" << endl <<endl;
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
	
#ifdef WIN32
	int xx;
	cin >> xx;
#endif
	return 0;
}
