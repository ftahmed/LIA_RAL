/*!
\file distance.cpp
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


#ifndef _DISTANCE_CPP
#define _DISTANCE_CPP

#include <vector>
using namespace std;

#include "distance.h"

real_t Distance::computePartialGLR(FrameAcc & acc) 
{
//	cerr << "det " << log(acc.getDet()) << " count " <<acc.getCount() << " mean" << acc.getMeanVect()[0] << endl;
	return (real_t)0.5 * acc.getCount() * log(acc.getDet());
}

real_t Distance::GLR(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij)
{
	real_t res = computePartialGLR(accGij) - computePartialGLR(accGi) - computePartialGLR(accGj);
//	cerr << "glr " << res << endl;
	return res;
}

real_t Distance::BIC(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij, real_t cst, real_t len)
{
	real_t res =  GLR(accGi, accGj, accGij) - cst * log(len);
	return res;	
}

real_t Distance::BICLocal(FrameAcc & accGi, FrameAcc & accGj, FrameAcc & accGij, real_t cst)
{
	real_t res = BIC(accGi, accGj, accGij, cst, accGij.getCount());
	return res;
}

real_t Distance::BICCst(String kind, real_t dim, real_t alpha)
{
	
	if (kind == "GF") {
		return 0.5 * alpha * (dim + 0.5 * ((dim + 1) * dim));
	}
	
	return 0.5 * alpha * (dim + dim);
}

real_t Distance::KL2(FrameAccGD & accGi, FrameAccGD & accGj)
{
	real_t s = 0.0;

	long dim = accGi.getVectSize();
	const DoubleVector & meanVectGi = accGi.getMeanVect();
	const DoubleVector & meanVectGj = accGj.getMeanVect();
	const DoubleVector & covVectGi = accGi.getCovVect();
	const DoubleVector & covVectGj = accGj.getCovVect();
	
	for (long j = 0; j < dim - 1; j++) {
		real_t m1 = meanVectGi[j];
		real_t m2 = meanVectGj[j];
		
		real_t dmean = m1 - m2;
		
		real_t v1 = covVectGi[j];
		real_t v2 = covVectGj[j];
		s += 0.25 * ((v1 / v2 + v2 / v1) +
					 dmean * dmean * (1.0 / v1 + 1.0 / v2) - 2.0);
	}
	s /= (real_t) dim;
	return s;
}

real_t Distance::GD(FrameAccGD & accGi, FrameAccGD & accGj)
{
	real_t s = 0.0;
	long dim = accGi.getVectSize();
	const DoubleVector & meanVectGi = accGi.getMeanVect();
	const DoubleVector & meanVectGj = accGj.getMeanVect();
	const DoubleVector & covVectGi = accGi.getCovVect();
	const DoubleVector & covVectGj = accGj.getCovVect();

	for (long j = 0; j < dim - 1; j++) {
		real_t dmean = meanVectGi[j] - meanVectGj[j];
		real_t v = sqrt(covVectGi[j]) * sqrt(covVectGj[j]);
		if (v < 0) {
			cerr << "waring[Distance] \t GD: variance problem" << endl;
			v = 1e-8;
		}
		s += (dmean * dmean) / v;
	}
	s /= (real_t) dim;
	return s;
}

real_t Distance::H2(FrameAccGD & accGi, FrameAccGD & accGj)
{
	real_t s = 0.0;
	long dim = accGj.getVectSize();
	FrameAccGD accGij = accGi;
	accGij.add(accGj);
	const DoubleVector & meanVectGi = accGi.getMeanVect();
	const DoubleVector & meanVectGj = accGj.getMeanVect();
	const DoubleVector & covVectGij = accGij.getCovVect();
	for (long j = 0; j < dim - 1; j++) {
		real_t dmean = meanVectGi[j] - meanVectGj[j];
		real_t v = covVectGij[j];
		if (v < 0) {
			cerr << "waring \t H2: variance problem" << endl;
			v = 1e-8;
		}
		s += (dmean * dmean) / v;
	}
	s /= (real_t) dim;
	return s;
}

real_t Distance::CLR(MixtureStat & msi, MixtureStat & msj, 
	Cluster & ci, Cluster & cj, ClusterServer & clusters,
	double vw_si, double vw_sj, FeatureServer & fs, bool useTop)
{
	long count_i = 0;
	long count_j = 0;
	real_t vgi_sj = getScore(fs, msi, cj, clusters, count_j, useTop);
	real_t vgj_si = getScore(fs, msj, ci, clusters, count_i, useTop);
	cerr << "CLR" << vgi_sj << " " << vw_sj << " " << vgj_si << " " << vw_si << endl;
	return (vgi_sj - vw_sj) + (vgj_si - vw_si);
}

real_t Distance::NCLR(MixtureStat & msi, MixtureStat & msj, 
	Cluster & ci, Cluster & cj, ClusterServer & clusters,
	FeatureServer & fs, bool useTop)
{
	long count_i = 0;
	long count_j = 0;
	real_t vgi_si = getScore(fs, msi, ci, clusters, count_i, useTop);
	real_t vgj_sj = getScore(fs, msj, cj, clusters, count_j, useTop);
	real_t vgi_sj = getScore(fs, msi, cj, clusters, count_j, useTop);
	real_t vgj_si = getScore(fs, msj, ci, clusters, count_i, useTop);
	cerr << "NCLR " << vgi_sj << " " << vgj_sj << " " << vgj_si << " " << vgi_si << endl;
	return (vgi_sj - vgj_sj) + (vgj_si - vgi_si);
}

real_t Distance::getScore(FeatureServer & fs, MixtureStat &ms, 
Cluster & c, ClusterServer & clusters, long & count, bool useTop)
{
	//we get the segments
	ms.resetLLK();
	count = 0;
	for(c.seekBeginSeg(); c.isEof() == false; c.seekNextSeg()){
		const Segment & seg = c.getCurrentSeg();
		const vector<LKVector> & top = seg.getTop();
		// Find the index of the first frame of the file in the buffer
		unsigned long begin = seg.getBegin() + 
			fs.getFirstFeatureIndexOfASource(clusters.getSourceName(seg.getIdxSrcName())); 
		
		// go to the frame in the buffer (and load it if needed)
		fs.seekFeature(begin);                                           
		for (unsigned long n=0; n < seg.getLength(); n++){
			Feature f;
			if(fs.readFeature(f)==false) cout<<"No more features"<<endl;
			if (useTop == true) {
				ms.computeAndAccumulateLLK(f, top[n]);
			} else {
				ms.computeAndAccumulateLLK(f);
			}
			count++;
		}
	}
	return ms.getMeanLLK();
}
#endif
