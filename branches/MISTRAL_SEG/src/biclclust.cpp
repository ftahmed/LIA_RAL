/*! 
\file biclclust.cpp
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

#ifndef _LCLUST_CPP
#define _LCLUST_CPP

#include "biclclust.h"
#include "float.h"

using namespace lium_seg;
using namespace alize;

BICLClust::BICLClust(ClusterServer & _clusters, FeatureServer & _features,
					 double _alpha, String & _kind):
					 BICHClust(_clusters,_features, _alpha, _kind)
{
}

void BICLClust::init(void)
{
	trainClusters();
}

long BICLClust::getClusterEnd(long &ci) 
{
	Cluster & cluster = clusters.getContainer()[idVect[ci]];
	long s = cluster.size();
	if (s > 1) {
		cerr << "WARNING : more than 1 segment in the cluster size=" << s << endl;
	}
	cluster.seekBeginSeg();
	const Segment & seg = cluster.getCurrentSeg();
	return seg.getBegin() + seg.getLength();
}

long BICLClust::getClusterStart(long &ci) 
{
	Cluster & cluster = clusters.getContainer()[idVect[ci]];
	cluster.seekBeginSeg();
	const Segment & seg = cluster.getCurrentSeg();
	return seg.getBegin();
}

void BICLClust::next(double &min, long &ci, long &cj)
{
	if (cj == (long)idVect.size()) {
		min = DBL_MAX;
	} else {
		long end = getClusterEnd(ci);
		long start = getClusterStart(cj);
		if ((start - end) > 0) {
			cerr << "WARNING : there is a hole between segments" << endl;
			min = 1.0;
		} else {			
			min = computeDistance(ci, cj);
		}
	}
}

void BICLClust::update(long ci, long cj)
{
	updateSwap(ci, cj);
	updateCluster(ci, cj);
	clusters.getContainer()[idVect[ci]].collapse();
	updateIdVect(cj);
	updateModels(ci, cj);
}

#endif
