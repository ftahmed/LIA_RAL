/*!
\file segment.cpp
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

#ifndef _SEGMENT_CPP
#define _SEGMENT_CPP

#include "segment.h"
using namespace lium_seg;

#include <cstdio>
#include <iostream>
using namespace std;

#include <cluster.h>

//constructor
Segment::Segment(void):top(0,0)
{
	//index in cluster
	idxSrcName = -1;

	//some informations
	String index1="channel";
	String index2="gender";
	String index3="band";
	String index4="environment";
	String strGender = "U";
	String strBand = "U";
	String strEnv = "U";
	String strChannel = "1"	;
	info.addLine(index1,strChannel);// add all the informations we have : first :  channel
	info.addLine(index2,strGender); // gender
	info.addLine(index3,strBand); // band
	info.addLine(index4,strEnv);// environnement

	//begin and length
	begin = 0;
	length = 0;
}

Segment::Segment(const Segment & _s)
{
	idxSrcName = _s.idxSrcName;
	info = _s.info;
	begin = _s.begin;
	length = _s.length;
	top = _s.top;	
}

//operator
Segment & Segment::operator=(const Segment & _s)
{
	idxSrcName = _s.idxSrcName;
	info = _s.info;
	begin = _s.begin;
	length = _s.length;
	top = _s.top;	
	return *this;
}

bool Segment::operator() (const Segment & seg1, const Segment & seg2) {
	if (seg1.idxSrcName == seg2.idxSrcName) {
		return (seg1.begin < seg2.begin);
	}
//	cerr << "Warning Seg::operator()" << endl;
	return (seg1.idxSrcName < seg2.idxSrcName);
}


void Segment::setIdxSrcName(const long &c)
{
	idxSrcName = c;
}

void Segment::setBegin(const unsigned long &c)
{
	begin = c;
}

void Segment::setLength(const unsigned long &c)
{
	length = c;
}

void Segment::setInfos(XList list)
{
	info = list;
}

//
long Segment::getIdxSrcName()  const
{
	return idxSrcName;
}

unsigned long Segment::getBegin() const
{
	return begin;
}

unsigned long Segment::getLength() const
{
	return length;
}

XList Segment::getInfos() const
{
	return info;
}
#endif
