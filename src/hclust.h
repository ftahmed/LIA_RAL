/*! 
\file hclust.h

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

\brief abstract hierarchical clustering class

*/

#ifndef _HCLUST_H
#define _HCLUST_H

#include <string>
#include <vector>

#include "clusterServer.h"
#include "distance.h"
#include <alize.h>

using namespace std;
using namespace alize;
using namespace lium_seg;

/*#include "clusters.h"
#include "model.h"
#include "modelGMM.h"
#include "modelGauss.h"
#include "modelGaussFull.h"
#include "modelGaussDiag.h"
#include "modelTools.h"
#include "tools.h"*/

namespace lium_seg
{
/// \class HClust
/// \brief Abstract Hierarchical clustering class.
class HClust {
public:
	explicit HClust(ClusterServer & clusters, FeatureServer & features, double alpha, String & kind);

	virtual ~ HClust(void);

	/// \brief initialize the clustering.
	///
	/// Train the model of each cluster and compute the distances between the clusters.
	virtual void init(void) = 0;

	/// \brief Find the next candidate to a merge.
	virtual void next(double &min, long &ci, long &cj);
	
	/// \brief Merge two clusters and update model and distances.
	virtual void update(long ci, long cj);
	/*
	/// \brief Get the resulting segmentation
	Clusters get();*/

	/// \brief Access to the inner #clusters
	ClusterServer & getContainer();

	/// \brief Return the name of the clusters
	String getClusterName(long i) { return clusters.getName(i); }

	/// \brief uninitialize the clustering.
	virtual void reset(void);

	/*virtual void debug(void);
	*/
protected:
	
	virtual double computeDistance(const long i, const long j) = 0;
	
	/// \brief merge \e ci and \e cj clusters in #clusters. The new cluster is \e ci.
	virtual void updateCluster(long ci, long cj);

	/// \brief Resize the distance matrix.
	virtual void updateValuesReduce(long ci, long cj);

	/// \brief Remove \e cj index in #idVect.
	virtual void updateIdVect(long cj);
	
	/// \brief Compute the resulting merged model and remove \e cj model.  
	virtual void updateModels(long ci, long cj) = 0;

	/// \brief Compute the new distances.
	virtual void updateValues(long ci, long cj);

	/// \brief Swap \e ci and \e cj if \e ci > \e cj.
	virtual void updateSwap(long &ci, long &cj);

	FeatureServer & features;
	long size; //!< Current number of clusters.
	ClusterServer clusters; //!< Inner clusters.
	MixtureServer mixtureServer;
	StatServer statServer;
	DoubleSquareMatrix values; //!< Matrix of distances.
	vector<long> idVect; //!< List of index (name of a cluster).
	double alpha;
	String & kind;
	Config & config;
private:
	HClust(const HClust&); /*!Not implemented*/
	const HClust& operator=(const HClust&);/*! not implemented */
	bool operator==(const HClust&) const;/*! not implemented */
	bool operator!=(const HClust&) const;/*! not implemented */
};

}
#endif
