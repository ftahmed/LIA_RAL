/*! 
\file trainTools.h

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

#ifndef _LIATOOLSWRAPPER_H
#define _LIATOOLSWRAPPER_H

#include <string>
#include <alize.h>
#include <liatools.h>
#include "segment.h"
#include "cluster.h"
#include "clusterServer.h"


using namespace std;
using namespace alize;
using namespace lium_seg;

/// \brief LIA Train tools warpper
class LiaToolsWrapper{
	public:
	static real_t accStatEM(StatServer & ss, FeatureServer &fs, 
		MixtureStat &emAcc,const Segment& seg, String& sourceName,Config &config);

	static void makeMAP(MixtureServer &ms, const MixtureGD& world,
		MixtureGD &client, unsigned long frameCount, Config &config);

	static void checkModel(MixtureGD& world, DoubleVector &globalCov, Config &config);
};


#endif
