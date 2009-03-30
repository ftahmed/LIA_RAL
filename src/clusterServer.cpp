/*! 
\file clusterServer.cpp
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

#ifndef _CLUSTER_SERVER_CPP
#define _CLUSTER_SERVER_CPP

#include "clusterServer.h"

#include <cstdio>
#include <iostream>

using namespace lium_seg;

ClusterServer::ClusterServer(void)  {
	maxIdxSourceName = 0;
	maxIdxClusterName = 0;
}

ClusterServer::ClusterServer(Config & config, const ClusterServer & _clusters)
{
	maxIdxSourceName = 0;
	maxIdxClusterName = 0;
	clusterMap = _clusters.clusterMap;
	srcNameMap = _clusters.srcNameMap;
	maxIdxSourceName = _clusters.maxIdxSourceName;
	clusterNameMap = _clusters.clusterNameMap;
	maxIdxClusterName = _clusters.maxIdxClusterName;
	pconfig = &config;
}

ClusterServer::ClusterServer(Config & config,const SrcNameMap & sm, long si, const ClusterNameMap & nm, long ni)
{
	srcNameMap = sm;
	maxIdxSourceName = si;
	clusterNameMap = nm;
	maxIdxClusterName = ni;
	pconfig = &config;
}

ClusterServer::ClusterServer(Config & config)
{
	maxIdxSourceName = 0;
	maxIdxClusterName = 0;
	pconfig = &config;
}

ClusterServer & ClusterServer::operator=(const ClusterServer & _clusters)
{
	clusterMap = _clusters.clusterMap;
	srcNameMap = _clusters.srcNameMap;
	maxIdxSourceName = _clusters.maxIdxSourceName;
	clusterNameMap = _clusters.clusterNameMap;
	maxIdxClusterName = _clusters.maxIdxClusterName;
	pconfig = _clusters.pconfig;
	return *this;
}

void ClusterServer::erase(long j)
{	
	String locutorName = getName(j);
	cerr << "Locutor erase : " << locutorName << " ";
	cerr << "Size of list : " << clusterNameMap.size() << endl;	
	clusterNameMap.erase(locutorName);
	clusterMap.erase(j);
}

void ClusterServer::merge(long i, long j)
{
	SegLst & lsti = clusterMap[i].getContainer();
	SegLst & lstj = clusterMap[j].getContainer();
	lsti.insert(lstj.begin(), lstj.end());
	erase(j);
}

// container
ClusterMap & ClusterServer::getContainer(void)
{
	return clusterMap;
}

const ClusterMap & ClusterServer::getContainer(void) const
{
	return clusterMap;
}

void ClusterServer::collapse(void)
{
	ClusterMap::iterator it = clusterMap.begin();
	for(;it != clusterMap.end(); it++){
		it->second.collapse();
	}
}

String ClusterServer::getName(long id)
{
	String result = "";
	//iterator
	ClusterNameMap::iterator it;
	for(it = clusterNameMap.begin(); it != clusterNameMap.end(); it++){
		if (it->second == id){
			result = it->first;
			break;
		}
	}
	return result;
}

String ClusterServer::getSourceName(long id)
{
	String result = "";
	//iterator
	SrcNameMap::iterator it;
	for(it = srcNameMap.begin(); it != srcNameMap.end(); it++){
		if (it->second == id){
			result = it->first;
			break;
		}
	}
	return result;
}

//debug
/*void ClusterServer::debug(void)
{
	ClusterMap::iterator it;
	for (it = clusterMap.begin(); it != clusterMap.end(); it++) {
		it->second.debug();
	}
}*/

XLine & ClusterServer::getSourceNameInXLine()
{
	XLine& myList = XLine::create();
	SrcNameMap::iterator it = srcNameMap.begin();
	for(;it != srcNameMap.end(); it++){
		myList.addElement(it->first);
	}
	return myList;
}

// reader, writer
/// \todo remove show[2048]
/// \todo use ioFile class
void ClusterServer::read(FILE * f)
{
	cerr << "trace[readFile] \t Begin read..." << endl;
	//read a line
	char line[5001];

	size_t n = 5000;
	char sourceName[2048];
	char name[2048];
	line[0]='\0';
	fgets(line, n, f);
	while (!feof(f)) {
		if (line[0] == '\n') {
			continue;
		}						// empty line 
		if (line[0] == '#') {
			continue;
		}						// empty line 
		if ((line[0] == ';') && (line[1] == ';')) {
			continue;
		}						// rem line
		Segment seg;
		long segChannel;
		char segGender[2];
		char segBand[2];
		char segEnv[2];

 		long begin;
		long length;

		int res = sscanf(line,"%s %ld %ld %ld %s %s %s %s\n", sourceName, &(segChannel), &(begin), &(length), segGender, segBand, segEnv, name);
		
		// we insert the source name
		long showSrcName = setIdxSrcName(sourceName);
		// we put the num of source name in segment
		seg.setIdxSrcName(showSrcName);

		//we get the name of the locutor
		long clusterNameIdx = setIdxClusterName(name);
		
		// all the informations of a segment
		String index1="channel";
		String index2="gender";
		String index3="band";
		String index4="environment";
		String strGender = segGender;
		String strBand = segBand;
		String strEnv = segEnv;
		String strChannel = String::valueOf(segChannel);

		seg.setBegin(begin);
		seg.setLength(length);

		seg.getInfos().addLine(index1,strChannel);// add all the informations we have : first :  channel
		seg.getInfos().addLine(index2,strGender); // gender
		seg.getInfos().addLine(index3,strBand); // band
		seg.getInfos().addLine(index4,strEnv);// environnement

		// creation of the cluster
		// search if a cluster exist : if true, we add his id, else we create it and add his id
		ClusterMap::iterator iterMap = clusterMap.find(clusterNameIdx); 
		if (iterMap == clusterMap.end()) {		
			Cluster * clustTemp = new Cluster(clusterNameIdx);
			clusterMap[clusterNameIdx] = *clustTemp;
		}
		//we insert the segment
		Cluster & cluster = clusterMap[clusterNameIdx];
		cluster.getContainer().insert(seg);
		if (res != 8) {
			cerr << "clusters: read() error seg segmentation read error" << endl;
			break;
		}
		line[0]='\0';
		fgets(line, n, f);
	}
	cerr << "trace[readFile] \t ...End read" << endl;
	//collapse();
	//free(line);
}

void ClusterServer::read(const String & sourceName)
{
	String segInputMask = "";
	cerr << "************* type : " << pconfig->getParam_loadSegServerFileFormat() << endl;
	switch (pconfig->getParam_loadSegServerFileFormat())
	{
		case SegServerFileReaderFormat_LIUM:
		{
			segInputMask = ".seg";
			break;
		}
		default:
			throw Exception("Unable to read this seg server file format",
			__FILE__, __LINE__);
	}
	String path = pconfig->getParam_segServerFilesPath();

	String concat = path + sourceName + segInputMask;
	
	cerr << "trace[readFile] \t Name file : " << concat << endl;

	FILE *f = fopen(concat.c_str(), "r");
	if (f == NULL) {
		cerr << "clusters: read(lst) could not read file " +
					  concat << endl;
	}
	read(f);
	fclose(f);	
}

void ClusterServer::write(const String& sourceName)
{
	String segOutputMask = ".seg";
	
	String path = pconfig->getParam_segServerFilesPath();

	String concat = path + sourceName + segOutputMask;
	cerr << "trace[writeFile] \t Name file : " << concat << endl;

	ClusterMap::iterator it;
	FILE *f = fopen(concat.c_str(), "w");
	
	write(f);
	fclose(f);
}

void ClusterServer::write(FILE * f)
{
	#if !defined NDEBUG
	cerr << "trace[writeFile] \t Begin write..." << endl;
	#endif

	ClusterMap::iterator iterMap = clusterMap.begin();
	for(;iterMap != clusterMap.end();iterMap++) {
		Cluster myCluster = iterMap->second;
		String nameLocutor = getName(myCluster.getId());
		SegLst& mySegments = myCluster.getContainer();
		SegLst::iterator iterSeg = mySegments.begin();
		for(;iterSeg != mySegments.end();iterSeg++) {
			Segment mySeg = *iterSeg;
			String index1 = "channel";
			String index2 = "gender";
			String index3 = "band";
			String index4 = "environment";
			String sourceName = getSourceName(mySeg.getIdxSrcName());
			//to get all the XList infos : mySeg.getInfos().searchValue(index1)
			fprintf(f, "%s %s %ld %ld %s %s %s %s\n", sourceName.c_str(), mySeg.getInfos().searchValue(index1).c_str(), mySeg.getBegin(), mySeg.getLength(), mySeg.getInfos().searchValue(index2).c_str(), mySeg.getInfos().searchValue(index3).c_str(), mySeg.getInfos().searchValue(index4).c_str(), nameLocutor.c_str());
		}		
	}
	#if !defined NDEBUG
	cerr << "trace[writeFile] \t ...End write" << endl;
	#endif
}

SegLst ClusterServer::getSegLst(bool collapse)
{
	SegLst segLst;
	ClusterMap::iterator it = getContainer().begin();
	//Compute Model
	while(it != getContainer().end()) {
		Cluster & c = it->second;
		SegLst & sLst = c.getContainer();
		segLst.insert(sLst.begin(), sLst.end());
		it++;
	}
	if (collapse) Cluster::collapse(segLst);
	return segLst;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/*MapSeg clustersToFrames(ClusterServer & c)
{
	MapSeg res;
	ClusterMap::iterator itCluster;
	ClusterMap & clusterMap = c.getContainer();
	for (itCluster = clusterMap.begin(); itCluster != clusterMap.end(); itCluster++) {
		SegLst & segLst = itCluster->second.getContainer();
		SegLst::iterator itSeg;
		string eti = itCluster->second.getName();
		for (itSeg = segLst.begin(); itSeg != segLst.end(); itSeg++) {
			long start = itSeg->getStart();
			long len =	itSeg->getLen();
			for(long i = start; i < (start + len); i++){
				Seg se = (*itSeg);
				se.setStart(i);
				se.setLen(1);
				se.setTmpName(eti);
				se.setGender(itCluster->second.getGender());
				res[i] = se;				
			}
		}
	}
	return res;
}*/

/*VectorSeg clustersToVector(ClusterServer & c)
{
	SegLst segLst = c.getSegLstWithName();
	VectorSeg vSeg;
	SegLst::iterator it = segLst.begin();
	for (; it != segLst.end(); it++){
		vSeg.push_back(*it);
	}
	return vSeg;
}

ClusterServer vectorToClusters(VectorSeg & vSeg, ClusterServer & clusters, bool collapse = false)
{
	ClusterServer res = clusters;
	ClusterMap & clusterMap = res.getContainer();
	ClusterMap::iterator itCluster = clusterMap.begin();
	for (; itCluster != clusterMap.end(); itCluster++) {
		itCluster->second.getContainer().clear();
	}
	long size = vSeg.size();
	for(long i = 0; i < size; i++) {
		string name = vSeg[i].getTmpName();
		long nameIdx = res.setIdxName(name);
		Cluster & c = res.getContainer()[nameIdx];
		c.getContainer().insert(vSeg[i]);
	}
	if (collapse) {
		res.collapse();
	}
	return res;
}*/
#endif
