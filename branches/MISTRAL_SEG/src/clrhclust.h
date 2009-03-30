/*! 
\file clrhclust.h

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

\brief hierarchical clustering class based on clr distance and GMM

*/

#ifndef _CLRHCLUST_H
#define _CLRHCLUST_H

using namespace std;

#include "hclust.h"
#include "liaToolsWrapper.h"

namespace lium_seg
{
/// \class CLRHClust
/// \brief Hierarchical clustering class based on clr distance and GMM.
///
/// The clustering is based upon a bottom-up hierarchical clustering. 
///	Each segment is associated to a cluster providing the initial set of clusters. 
///	The two closest clusters \f$c_i\f$ and \f$c_j\f$ are merged at each iteration 
///	until a stop criterion is met. 
///
///	
class CLRHClust:public HClust {
public:
	explicit CLRHClust(ClusterServer & _clusters, FeatureServer & _features,
			 double _alpha, String & kind, bool _nclr);

	virtual ~CLRHClust(void);

	/// Train the model of each cluster and compute the distances between the clusters.
	virtual void init(void);

	/// \brief uninitialize the clustering.
	virtual void reset(void);

	/// \brief Merge two clusters and update model and distances.
	virtual void update(long ci, long cj);
protected:
	unsigned long getMAP(Mixture & client, Mixture & world, MixtureStat &clientAcc, Cluster & cluster);

	/// \brief Train a cluster. 
	virtual void trainCluster(const long i);

	/// /see #trainCluster
	void trainClusters();

	/// \brief compute the log likelihood mean of the world model for each cluster
	void computeScoreWld();

	/// \brief Compute a BIC distance between \e i and \e j.
	virtual double computeDistance(const long i, const long j);

	/// \brief update Models ci and cj by computing a new one ci using ci and cj data
	virtual void updateModels(long ci, long cj);

	/// \brief update the log likelihood mean of the world model of the cluster new ci (merge of ci and cj)
	virtual void updateScoreWld(long ci, long cj);

	double cst; //!< Constant in BIC.
	long len; //!< Number of features in segmentation.
	vector<double> scoreWorld;
	MixtureServer mixtureServerWorld;
	StatServer statServerWorld;
	bool useTop;
	bool nclr;
private:
	CLRHClust(const CLRHClust&); /*!Not implemented*/
	const CLRHClust& operator=(const CLRHClust&);/*! not implemented */
	bool operator==(const CLRHClust&) const;/*! not implemented */
	bool operator!=(const CLRHClust&) const;/*! not implemented */
};
}
#endif
