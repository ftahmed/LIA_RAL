/*! 
\file distance.h

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

\brief tool classes for models.
*/

#ifndef _DISTANCE_H
#define _DISTANCE_H

#include <string>
#include <alize.h>
#include "segment.h"
#include "cluster.h"
#include "clusterServer.h"

using namespace std;
using namespace alize;
using namespace lium_seg;

/// \brief Distance class, ie GLR, BIC, KL2, GD, CLR for Gaussian or GMM.
class Distance{
	public:
		/**
	/// \brief Get a GLR score for GMMs, merged model is trained by EM. 
	static double GLR_GMM(GMM & gi, GMM & gj, SegLst & si,
		SegLst & sj, Features & features, Param::TrainInitMethod trainInitMethod, 
		long minIt, long maxIt, double minGain, double flooring, double ceilling, 
		bool trace);
	
	/// \brief Get a GLR score for GMMs, merged model is trained by MAP. 
	static double GLR_GMM_MAP(GMM & gi, GMM & gj, SegLst & si,
		SegLst & sj, Features & features, const GMM & wld, double alpha, bool m, 
		bool c, bool w, long minIt, long maxIt, double minGain, double flooring, 
		double ceilling, bool trace);
**/
	/// \brief Compute partial GLR score form an accumulator. 
	static real_t computePartialGLR(FrameAcc & acc);

	/// \brief Get a GLR score for Gaussians. 
	static real_t GLR(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij);

	/// \brief Get a BIC score for Gaussians given a constant and the length. 
	static real_t BIC(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij, real_t cst, real_t len);

	/// \brief Get a BIC score for Gaussians using a length of the clusters. 
	static real_t BICLocal(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij, real_t cst);
/**	
	/// \brief Get a BIC score for EM-GMMs using a length of the clusters. 
	static double BICLocal_GMM_EM(GMM & gi, GMM & gj, SegLst & si, SegLst & sj, 
		Features & features, Param::TrainInitMethod trainInitMethod, long minIt, 
		long maxIt, double minGain,	double cst, double flooring, 
		double ceilling, bool trace);
	
	/// \brief Get a BIC score for MAP-GMMs using a length of the clusters. 
	static double BICLocal_GMM_MAP(GMM & gi, GMM & gj, SegLst & si, SegLst & sj, 
		Features & features, const GMM & wld, double alpha, bool m, bool c, 
		bool w, long minIt, long maxIt, double minGain, double cst, 
		double flooring, double ceilling, bool trace);	
	**/
	/// \brief Get a BIC constant. 
	static real_t BICCst(String kind, real_t dim, real_t alpha);

	/// \brief Get a KL2 score for diagonal Gaussians. 
	static real_t KL2(FrameAccGD & accGi, FrameAccGD & accGj);

	/// \brief Get a GD score for diagonal gaussians. 
	static real_t GD(FrameAccGD & g1, FrameAccGD & g2);
	
	/// \brief Get a H2 score for diagonale Gaussians. 
	static real_t H2(FrameAccGD & gi, FrameAccGD & gj);
	
	/// \brief Get a CLR distance with UBM score precomputed
	static real_t CLR(MixtureStat & msi, MixtureStat & msj, 
		Cluster & ci, Cluster & cj, ClusterServer & clusters,
		double vw_si, double vw_sj, FeatureServer & fs, bool useTop = false);

	/// \brief Get a NCLR distance
	static real_t NCLR(MixtureStat & msi, MixtureStat & msj, 
		Cluster & ci, Cluster & cj, ClusterServer & clusters,
		FeatureServer & fs, bool useTop = false);

	static real_t getScore(FeatureServer & fs, MixtureStat &ms, 
		Cluster & c, ClusterServer & clusters, long & count, bool useTop = false);
	/// \brief Get the log-likelihood of a GMM over a list of segments.
/*	static double getScoreSetTop(Mixture & g, long nbTop, SegLst & s,
		FeatureServer & features, long & count);*/
	protected:

};


#endif
