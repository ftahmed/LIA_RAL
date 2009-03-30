/*!
\file liaToolsWrapper.cpp
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


#ifndef _LIATOOLSWRAPPER_CPP
#define _LIATOOLSWRAPPER_CPP

#include "liaToolsWrapper.h"

real_t LiaToolsWrapper::accStatEM(StatServer & ss, FeatureServer &fs,MixtureStat &emAcc,
	const Segment& seg,String & sourceName,Config &config)
{
 	// Find the index of the first frame of the file in the buffer
	unsigned long begin = seg.getBegin()+fs.getFirstFeatureIndexOfASource(sourceName);              
	unsigned long len = seg.getLength();
	real_t v = accumulateStatEM(ss, fs, emAcc, begin, len, config);
	return v;
	
}

void LiaToolsWrapper::makeMAP(MixtureServer &ms,const MixtureGD& world,MixtureGD &client,
	unsigned long frameCount,Config &config)
{
	computeMAP(ms, world, client, frameCount, config);
}

void LiaToolsWrapper::checkModel(MixtureGD& world, DoubleVector &globalCov, Config &config)
{
	real_t varianceFlooring = config.getFloatParam("varFloor");
	real_t varianceCeiling  = config.getFloatParam("varCel");   

	varianceControl(world, varianceFlooring, varianceCeiling, globalCov);
    if (config.getParam("nomModel") == "True") normalizeMixture(world, config);	

}


#endif
