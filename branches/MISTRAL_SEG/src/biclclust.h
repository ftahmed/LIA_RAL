/*! 
\file biclclust.h

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

\date 30 Nov 2006
\version 0.0.2

\brief linear clustering, ie a clustering that merge only adjacent segment.

\history
-version 0.0.1, 12 Oct 2006 : initial version
-version 0.0.2, 30 Nov 2006 : the segment must be adjacent

*/

#ifndef _BICLCLUST_H
#define _BICLCLUST_H

#include "bichclust.h"

namespace lium_seg
{
/// \brief linear clustering class
///	
/// The computation cost of a hierarchical clustering depends of the number of 
///	input segments. The computation cost could be reduced by applying first a 
///	clustering that merged only adjacent segments \f$s_i\f$ and \f$s_{i+1}\f$ into 
///	\f$s_i\f$ from \f$i=1\f$ to \f$i=K\f$ that satisfy \f$\Delta BIC_{i,j} < 0\f$. 
///	if \f$\Delta BIC_{i,j} \geq 0\f$ the next candidates are \f$s_{i+1}\f$ and 
///	\f$s_{i+2}\f$ else the next candidates are the new \f$s_i\f$ and \f$s_{i+2}\f$. 
///	This method is denoted linear clustering.
class BICLClust:public BICHClust {
  public:
	BICLClust(ClusterServer & _clusters, FeatureServer & _features,
			 double _alpha, String & kind);

	virtual void init(void);

	virtual void next(double &min, long &ci, long &cj);

	virtual void update(long ci, long cj);
	
	long getClusterEnd(long &ci);
	
	long getClusterStart(long &ci);
	
	void debug(void) 
	{
	}
private:
	BICLClust(const BICLClust&); /*!Not implemented*/
	const BICLClust& operator=(const BICLClust&);/*! not implemented */
	bool operator==(const BICLClust&) const;/*! not implemented */
	bool operator!=(const BICLClust&) const;/*! not implemented */
};
}
#endif
