/*! 
\file bichclust.cpp
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

#ifndef _BICHCLUST_CPP
#define _BICHCLUST_CPP

#include "bichclust.h"
#include "hclust.h"

using namespace lium_seg;
using namespace alize;

BICHClust::BICHClust(ClusterServer & _clusters, FeatureServer & _features,
					 double _alpha, String & _kind):HClust(_clusters, _features, 
					 _alpha, _kind)
{
	double d = features.getVectSize();
	BICCst = Distance::BICCst(kind, d, alpha);
}

BICHClust::~BICHClust(void)
{
	for(unsigned long i = 0; i < vectorClusterAcc.size(); i++) {
		delete vectorClusterAcc[i];
	}
}
// TODO: make the free of FrameAcc !
void BICHClust::trainClusters()
{
	cerr << "[bichclust] \t Begin trainClusters " << endl; 
	long s = idVect.size();
	vectorClusterAcc.reserve(s);
	len = 0;
	//mxServer.reset();
	for (long i = 0; i < s; i++) {		
		// iteration in all the segments of the cluster
		Cluster & clust = clusters.getContainer()[idVect[i]];

		FrameAcc * accu; 
		if(clusters.getConfig().getParam("distribType") == "GF") {
			accu = new FrameAccGF();
		}
		else {
			accu = new FrameAccGD();
		}
		//we get the segments
		for(clust.seekBeginSeg(); clust.isEof() == false; clust.seekNextSeg()) {
			//we get the segment
			const Segment & seg = clust.getCurrentSeg();
			// we want the features of the source of the segment
			String sourceName = clusters.getSourceName(seg.getIdxSrcName());
			//in order to deal with multiple stream source, we add to the begin and the end
			//the position in feature of the source
			unsigned long position = features.getFirstFeatureIndexOfASource(sourceName);
			//we get the minimum size between the feature and the seg
			unsigned long start = seg.getBegin() + position; // begin of the segment
			unsigned long end = start + seg.getLength();
			Feature f;
			for (unsigned long j = start; j < end; j++) {
				features.seekFeature(j);
				features.readFeature(f);
				accu->accumulate(f);
			}
		}
		vectorClusterAcc.push_back(accu);
	}
	cerr << "[bichclust] \t End trainClusters " << endl; 
}

double BICHClust::computeDistance(const long i, const long j)
{
	double v = 0.0;
	FrameAcc * gi = vectorClusterAcc[i];
	FrameAcc * gj = vectorClusterAcc[j];
	FrameAcc * gij = NULL;
	if(clusters.getConfig().getParam("distribType") == "GF") {
		gij = new FrameAccGF();
		((FrameAccGF*)gij)->add((FrameAccGF&)(*gi));
		((FrameAccGF*)gij)->add((FrameAccGF&)(*gj));
	}
	else {
		gij = new FrameAccGD();
		((FrameAccGD*)gij)->add((FrameAccGD&)(*gi));
		((FrameAccGD*)gij)->add((FrameAccGD&)(*gj));
	}
	
	v = Distance::BICLocal(*gi, *gj, *gij, BICCst);
	delete gij;
	return v;
}

void BICHClust::init(void)
{
	cerr << "[bichclust] \t Begin init " << endl; 
	trainClusters();
	
	values.setAllValues(0.0);
	for (long i = 0; i < (long)idVect.size(); i++) {
		for (long j = i + 1; j < (long)idVect.size(); j++) {
			values(i, j) = computeDistance(i, j);
		}
	}
	cerr << "[bichclust] \t End init " << endl; 
}

void BICHClust::updateModels(long ci, long cj)
{

	FrameAcc * gi = vectorClusterAcc[ci];
	FrameAcc * gj = vectorClusterAcc[cj];
	if(clusters.getConfig().getParam("distribType") == "GF") {
		((FrameAccGF*)gi)->add((FrameAccGF&)(*gj));
	}
	else {
		((FrameAccGD*)gi)->add((FrameAccGD&)(*gj));
	}

	vectorClusterAcc.erase(vectorClusterAcc.begin() + cj);
}

#endif
