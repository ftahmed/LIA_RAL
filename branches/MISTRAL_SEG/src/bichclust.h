/*! 
\file bichclust.h

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

\brief hierarchical clustering class based on a BIC distance computed over mono-gaussian

*/

#ifndef _BICHCLUST_H
#define _BICHCLUST_H

using namespace std;

#include "hclust.h"

namespace lium_seg
{
/*! \class BICHClust
 \brief BICHierarchical clustering class based on a BIC distance computed 
	over mono Gaussian

 The clustering is based upon a bottom-up hierarchical clustering. 
	Each segment is associated to a cluster providing the initial set of clusters. 
	The two closest clusters \f$c_i\f$ and \f$c_j\f$ are merged at each iteration 
	until a stop criterion is met. Various metric and stop criterion have been 
	proposed in [Siu1991,Solomonoff1998,Chen1998a,Reynolds1998]. One of 
	those is \f$\Delta BIC\f$ metric employed to select the clusters to group as 
	well as to stop the merge process. The candidate clusters \f$c_i\f$ and \f$c_j\f$ 
	are selected according to:
 \f[\min_{i>j, i \in \{1, \ldots, K\}}\Delta BIC_{i,j}\f]
 with \f$K\f$ the number of cluster. Moreover the merge process is 
	stopped when \f$\Delta BIC_{i,j} > 0\f$. 
 Two kinds of BIC metrics differing by the penalty factor are proposed: 
 \li the local BIC penalty factor \f$P_l\f$ with \f$n_i\f$ and \f$n_j\f$ 
	respectively the sum of segment length in \f$c_i\f$ and \f$c_j\f$:
 \f[P_l=\frac{1}{2}\left(d + \frac{d(d+1)}{2}\right) log(n_1+n_2)\f]
 \li or the global BIC penalty factor \f$P_g\f$ with \f$N\f$ the length of 
	the whole signal:
 \f[P_g=\frac{1}{2}\left(d + \frac{d(d+1)}{2}\right) log(N)\f]
 The LIMSI proved that local BIC penalty factor gives better result in 
 [Barras2004] and this is confirmed in [Tranter2004].
 Provide method to perform a hierarchical clustering.
 2 methods are provide: 
 \li to select the next candidates: the call of method #next returns the 
	next couples of clusters to merge.
 \li the merge are perform by the call of the #update method ;
	the couples of clusters are merge and models, distances are updated.\n
 then the end of the clustering could be control out side of the class.
*/

typedef vector<FrameAcc*> VectorFrameAcc;


class BICHClust:public HClust {
public:
	explicit BICHClust(ClusterServer & _clusters, FeatureServer & _features,
			 double _alpha, String & kind);
	
	virtual ~BICHClust(void);
	
	/// \brief initialize the clustering.
	///
	/// Train the model of each cluster and compute the distances between the clusters.
	virtual void init(void);
protected:
	/// \brief Train clusters. 
	void trainClusters();
	
	/// \brief Compute a BIC distance between \e i and \e j.
	virtual double computeDistance(const long i, const long j);
	
	virtual void updateModels(long ci, long cj);
	//
	double BICCst; //!< Constant in BIC.
	long len; //!< Number of features in segmentation.
	VectorFrameAcc vectorClusterAcc;
private:
	BICHClust(const BICHClust&); /*!Not implemented*/
	const BICHClust& operator=(const BICHClust&);/*! not implemented */
	bool operator==(const BICHClust&) const;/*! not implemented */
	bool operator!=(const BICHClust&) const;/*! not implemented */
};
}
#endif
