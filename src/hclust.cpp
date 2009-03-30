/*! 
\file hclust.cpp
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

#ifndef _HCLUST_CPP
#define _HCLUST_CPP

#include "hclust.h"
#include "float.h"

using namespace alize;
using namespace lium_seg;

HClust::HClust(ClusterServer & aClusters, FeatureServer & aFeatures, double aAlpha, String & aKind):
	features(aFeatures),mixtureServer(aClusters.getConfig()),
	statServer(aClusters.getConfig(), mixtureServer),kind(aKind),
	config(aClusters.getConfig())
{
	alpha = aAlpha;
	clusters = aClusters;
	ClusterMap & map = clusters.getContainer();
	size = map.size();
	values = DoubleSquareMatrix(size);
	ClusterMap::iterator it;
	for (it = map.begin(); it != map.end(); it++) {
		idVect.push_back(it->first);
	}
	mixtureServer.reset();
}

HClust::~HClust(void)
{
}

void HClust::reset(void)
{
	mixtureServer.reset();
	values.setSize(0, true);
	statServer.reset();
	
}


void HClust::updateCluster(long ci, long cj)
{
	clusters.merge(idVect[ci], idVect[cj]);
}

void HClust::updateValuesReduce(long ci, long cj)
{
	// scores ci is cij, cj delete
	for (long i = cj + 1; i < (long) idVect.size(); i++) {
		for (long j = 0; j < (long) idVect.size(); j++) {
			values(i - 1, j) = values(i, j);
		}
	}
	for (long i = cj + 1; i < (long) idVect.size(); i++) {
		for (long j = 0; j < (long) idVect.size(); j++) {
			values(j, i - 1) = values(j, i);
		}
	}
	for (long i = 0; i < (long) idVect.size(); i++) {
		values(i, idVect.size() - 1) = DBL_MAX;
		values(idVect.size() - 1, i) = DBL_MAX;
	}
}

void HClust::updateIdVect(long cj)
{
	// update idVect, delete cj
	idVect.erase(idVect.begin() + cj);
}

void HClust::updateValues(long ci, long cj)
{
	for (long i = 0; i < ci; i++) {
		values(i, ci) = computeDistance(i, ci);
		cerr << "score (" << i << ", " << ci << ") = " << values(i, ci) << endl;
	}
	for (long i = ci + 1; i < (long) idVect.size(); i++) {
		values(ci, i) = computeDistance(ci, i);
		cerr << "score (" << ci << ", " << i << ") = " << values(ci, i) << endl;
	}
}

void HClust::updateSwap(long &ci, long &cj)
{
	if (ci > cj) {
		long tmp = ci;
		ci = cj;
		cj = tmp;
	}
}

void HClust::update(long ci, long cj)
{
	updateSwap(ci, cj);
	updateCluster(ci, cj);
	updateValuesReduce(ci, cj);
	updateIdVect(cj);
	updateModels(ci, cj);
	updateValues(ci, cj);
}

/*Clusters HClust::get()
{
	return clusters;
}*/

ClusterServer & HClust::getContainer()
{
	return clusters;
}

/*void HClust::debug(void)
{
	cerr << "debug[HClust] \t BIC score" << endl;
	for (long i = 0; i < values.getDim(); i++) {
		for (long j = 0; j < values.getDim(); j++) {
			cerr << " " << values(i, j);
		}
		cerr << endl;
	}
}
**/
void HClust::next(double &min, long &ci, long &cj)
{
	min = DBL_MAX;
	ci = -1;
	cj = -1;
	if (idVect.size() > 1) {
		ci = -1;
		cj = -1;
		for (long i = 0; i < (long) idVect.size(); i++) {
			for (long j = i + 1; j < (long) idVect.size(); j++) {
				if (values(i, j)  < min) {
					ci = i;
					cj = j;
					min = values(i, j);
				}
			}
		}
	}
}

#endif
