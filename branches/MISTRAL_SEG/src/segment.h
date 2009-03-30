/*! 
\file segment.h

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

#ifndef _SEGMENT_H
#define _SEGMENT_H

#include <map>
#include <vector>
#include <alize.h>

using namespace std;
using namespace alize;

namespace lium_seg{

/// \brief List of top Gaussians.
typedef vector<long> TopVect;	
//typedef vector<double> DblVect;	
typedef vector<TopVect> VectTopVect;	
	
/// \brief Segment
class Segment {
public:
	Segment(void);

	Segment(const Segment & _s);

	//operator
	Segment & operator=(const Segment & _s);

	bool operator() (const Segment & seg1, const Segment & seg2);

	//setter
	/// \brief Set index of the source name.
	void setIdxSrcName(const long &c);

	/// \brief Set start of the segment, ie a feature index.
	void setBegin(const unsigned long &c);

	/// \brief Set the length of a segment in number of features.
	void setLength(const unsigned long &c);

	/// \brief Set the length of a segment in number of features.
	void setInfos(XList list);

	//getter
	/// \brief Get the index of the show.
	long getIdxSrcName() const;

	/// \brief Get the start index of the segment.
	unsigned long getBegin() const;

	/// \brief Get the informations.
	XList getInfos() const;

	/// \brief Get the length of the segment.
	unsigned long getLength() const;

	//get the top
	const vector<LKVector> & getTop() const { return top; };

	vector<LKVector> & getTop() { return top; };
private:
	long idxSrcName; //!< Index of the SrcName.
	unsigned long begin; //!< Start time in feature.
	unsigned long length; //!< Length of the segment
	XList info;
	vector<LKVector> top;
};
}
#endif
