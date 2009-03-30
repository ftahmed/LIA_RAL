/*! 
\file clrhclust.cpp
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

#ifndef _CLRHCLUST_CPP
#define _CLRHCLUST_CPP

#include "clrhclust.h"
using namespace lium_seg;

CLRHClust::CLRHClust(ClusterServer & _clusters, FeatureServer & _features,
					double _alpha, String & _kind, bool _nclr):
					HClust(_clusters, _features, _alpha, _kind),
					mixtureServerWorld(_clusters.getConfig()),
					statServerWorld(_clusters.getConfig(), mixtureServerWorld)					 
{
	config.setParam("computeLLKWithTopDistribs", "PARTIAL");//COMPLETE or PARTIAL @see documentation
	double d = features.getVectSize();
	cst = Distance::BICCst(kind, d, alpha);
//	useTop = false;
	useTop = (config.getParam_topDistribsCount() > 0 ? true : false);
	nclr = _nclr;
	//get world model
	String inputFilename = config.getParam("world");
	mixtureServerWorld.load(inputFilename);
	if (mixtureServerWorld.getMixtureCount() != 1) {
		cerr << "error : CLRHClust, bad number of wld" << endl;
		exit(0);
	}
	Mixture & world = mixtureServerWorld.getMixture(0);
	statServerWorld.createAndStoreMixtureStat(world);
	world.computeAll();
	computeScoreWld();
}

CLRHClust::~CLRHClust(void)
{
}

void CLRHClust::reset(void)
{
	HClust::reset();
	scoreWorld.clear();
	mixtureServerWorld.reset();
	statServerWorld.reset();

	for (long i = 0; i < (long) idVect.size(); i++) {
		Cluster & cluster = clusters.getContainer()[idVect[i]];
		for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
			Segment seg = cluster.getCurrentSeg();
			vector<LKVector> & topVect = seg.getTop();
			topVect.clear();
			topVect.resize(0);
		}
	}
}

unsigned long CLRHClust::getMAP(Mixture & client, Mixture & world, MixtureStat &clientAcc, Cluster & cluster)
{
	long nbIt = config.getIntegerParam("nbTrainIt");
	//real_t minGainEM = config.getFloatParam("minGainEM");
	unsigned long frameCount = 0;
	
	real_t llkOld = 0.0;
	for (long trainIt = 0; trainIt < nbIt; trainIt++)
	{
		clientAcc.resetEM();
		
		real_t llk = 0.0;
		unsigned long sumLen = 0;
		cluster.seekBeginSeg();
		for(; cluster.isEof() == false; cluster.seekNextSeg()) {
			const Segment & seg = cluster.getCurrentSeg();
			
			String sourceName = clusters.getSourceName(seg.getIdxSrcName());
			llk += LiaToolsWrapper::accStatEM(statServer, features, clientAcc, seg, sourceName, config);
			sumLen += seg.getLength();
		}
		llk /= (double) sumLen;
		if (trainIt != 0) {
			real_t diff = llk - llkOld;
			cerr << "it = " << trainIt << " llk = " << llk << " delta = " << diff << endl;
//			if (diff < minGainEM) {
//				break;
//			}
		} else {
			cerr << "it = " << trainIt << " llk = " << llk << endl;
		}
		llkOld = llk;
		client = clientAcc.getEM();
		long frameCount = (unsigned long) clientAcc.getEMFeatureCount();
		//MAP
		LiaToolsWrapper::makeMAP(mixtureServer, (MixtureGD&)world, (MixtureGD&)client, frameCount, config);
		//TODO : variance flooring and model normalization
	}
	return frameCount;
}

void CLRHClust::trainCluster(const long i)
{
	Cluster & cluster = clusters.getContainer()[idVect[i]];
	
	//init
	Mixture & world = mixtureServerWorld.getMixture(0);
	Mixture & client = mixtureServer.duplicateMixture(world,DUPL_DISTRIB);
	client.computeAll();
	
	MixtureStat &clientAcc = statServer.createAndStoreMixtureStat(client);

	//calculate a map
	cerr << "trace[CLRHClust] \t mTrainMAP cluster=" << cluster.getId() << endl;
	getMAP(client, world, clientAcc, cluster);	
}


void CLRHClust::trainClusters()
{
	long s = idVect.size();
	len = 0;
	statServer.deleteAllMixtureStat();
	cerr << "trace ------------------------------------" << endl;
	for (long i = 0; i < s; i++) {
		trainCluster(i);
	}
}

void CLRHClust::computeScoreWld()
{
	unsigned long nTop = config.getParam_topDistribsCount();
	
	MixtureStat &worldMixtureStat = statServerWorld.getMixtureStat(0);
	scoreWorld.reserve(idVect.size());
	for (long i = 0; i < (long)idVect.size(); i++) {
		Cluster & cluster = clusters.getContainer()[idVect[i]];
		worldMixtureStat.resetLLK();
		real_t llk = 0.0;
		unsigned long lenGlobal = 0;
		for(cluster.seekBeginSeg(); cluster.isEof() == false; cluster.seekNextSeg()){
			Segment seg = cluster.getCurrentSeg();
			cluster.eraseCurrentSeg();
			vector<LKVector> & topVect = seg.getTop();
			topVect.clear();
			//A revoir
			unsigned long idxBeginFrame = seg.getBegin() + 
				features.getFirstFeatureIndexOfASource(clusters.getSourceName(seg.getIdxSrcName())); 
			Feature f;
			unsigned long len = seg.getLength();
			lenGlobal += len;
			for (unsigned long idxFrame = idxBeginFrame; idxFrame < (len + idxBeginFrame); idxFrame++){
				features.seekFeature(idxFrame);
	
				features.readFeature(f);
				llk += worldMixtureStat.computeAndAccumulateLLK(f, 1.0, DETERMINE_TOP_DISTRIBS);
				if (useTop) {
					LKVector top = statServerWorld.getTopDistribIndexVector();
					top.pack(nTop);
					topVect.push_back(top);
				}
			}
			cluster.insertSeg(seg);
		}
		llk = Distance::getScore(features, worldMixtureStat, cluster, clusters, (long&)lenGlobal, useTop);
		//llk /= (real_t)lenGlobal;
		cerr << "ubm " << i << "=" << llk << endl;
		scoreWorld.push_back(llk);
	}
}

double CLRHClust::computeDistance(const long i, const long j)
{
	double v = 0.0;
	Cluster & ci = clusters.getContainer()[idVect[i]];
	Cluster & cj = clusters.getContainer()[idVect[j]];
	MixtureStat & msi = statServer.getMixtureStat(i);
	MixtureStat & msj = statServer.getMixtureStat(j);
	
	if (nclr == false) {
		v = -1.0 * Distance::CLR(msi, msj, ci, cj, clusters, scoreWorld[i], scoreWorld[j], features, useTop);
	} else {
		v = -1.0 * Distance::NCLR(msi, msj, ci, cj, clusters,  features, useTop);
	}
	return v;		
}

void CLRHClust::init(void)
{
	trainClusters();
	values.setAllValues(0.0);
	for (long i = 0; i < (long) idVect.size(); i++) {
		for (long j = i + 1; j < (long) idVect.size(); j++) {
			values(i, j) = computeDistance(i, j);
			cerr << "score (" << i << ", " << j << ") = " << values(i, j) << endl;
		}
	}
}

void CLRHClust::updateModels(long ci, long cj)
{
	Cluster & cluster = clusters.getContainer()[idVect[ci]];

	//init
	Mixture & world = mixtureServerWorld.getMixture(0);
	MixtureStat & clientStat = statServer.getMixtureStat(ci);
	Mixture & client = mixtureServer.getMixture(ci);
	client = world;
	//MAP
	getMAP(client, world, clientStat, cluster);
	//delete old model
	mixtureServer.deleteMixtures(cj, cj);
	mixtureServer.deleteUnusedDistribs();
	statServer.deleteMixtureStat(cj, cj);
}

void CLRHClust::updateScoreWld(long ci, long cj)
{
	Cluster & c_ci = clusters.getContainer()[idVect[ci]];
	Cluster & c_cj = clusters.getContainer()[idVect[cj]];
	real_t nb_ci = (real_t)c_ci.getLength();
	real_t nb_cj = (real_t)c_cj.getLength();
	cerr << "update world ci=" << ci << " s=" << scoreWorld[ci] << " nb=" << nb_ci  << " cj=" << cj << " s=" << scoreWorld[cj] << "nb=" << nb_cj << endl;
	scoreWorld[ci] = (scoreWorld[ci] * nb_ci + scoreWorld[cj] * nb_cj) / (nb_ci + nb_cj);
	scoreWorld.erase(scoreWorld.begin() + cj);
}

void CLRHClust::update(long ci, long cj)
{
	updateSwap(ci, cj);
	updateScoreWld(ci, cj);
	updateCluster(ci, cj);
	updateValuesReduce(ci, cj);
	updateIdVect(cj);
	updateModels(ci, cj);
	updateValues(ci, cj);
}

#endif
