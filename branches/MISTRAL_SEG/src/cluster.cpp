/*! 
\file cluster.cpp
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

#ifndef _CLUSTER_CPP
#define _CLUSTER_CPP

#include "cluster.h"
#include "clusterServer.h"

using namespace std;
using namespace lium_seg;

//constructor
Cluster::Cluster(const long & _id)
{
	setId(_id);
}

Cluster::Cluster(const Cluster & _c)
{
	lst = _c.lst;
	info = _c.info;
	id = _c.id;
}

Cluster::Cluster()
{
}

//operator
Cluster & Cluster::operator=(const Cluster & _c)
{
	lst = _c.lst;
	info = _c.info;
	id = _c.id;
	return *this;
}

//setter getter
void Cluster::setId(const long &_id)
{
	id = _id;
}

long Cluster::getId(void)
{
	return id;
}

// container
SegLst & Cluster::getContainer(void)
{
//	cerr << "Cluster::getContainer(void)" << endl;
	return lst;
}

const SegLst & Cluster::getContainer(void) const
{
//	cerr << "const Cluster::getContainer(void)" << endl;
	return lst;
}

void Cluster::collapse(void)
{
	//cerr << "collapse" << endl;
	collapse(lst);
}

void Cluster::collapse(SegLst & segLst)
{
	SegLst::iterator its = segLst.begin();
	SegLst::iterator prev = its;
	if (its != segLst.end()) its++;
	while(its != segLst.end()){
		if (its->getIdxSrcName() == prev->getIdxSrcName()) {
			if ((prev->getBegin()+prev->getLength()) >= its->getBegin()){
				Segment s = *prev;
				s.setLength(its->getBegin() - s.getBegin() + its->getLength());
				segLst.erase(*prev);
				prev = segLst.insert(s).first;
				segLst.erase(*its);
				its = prev;
			}
		}
		prev = its;		
		its++;
	}
}

bool Cluster::find(unsigned long idx)
{
	SegLst & lst = getContainer();
	SegLst::iterator it = lst.begin();
	for(; it != lst.end(); it++){
		unsigned long end = it->getBegin()+it->getLength();
		if ((idx>=it->getBegin())&&(idx < end)){
			return true;
		}
	}
	return false;
}

//debug
/*void Cluster::debug(void)
{
	for (SegLst::iterator it = lst.begin(); it != lst.end(); it++) {
		cerr << "debug[Cluster] \t seg show= " << it->getShow() << " idxShow= " << it->getIdxShow();
		cerr << " id= " << id << " name= " << name;
		cerr << " channel= " << it->getChannel();
		cerr << " start= " << it->getStart() << " len= ";
		cerr << it->getLen() << " gender= " << gender;
		cerr << " band= " << it->getBand() << " env= " << it->getEnv() << endl;
	}
}*/

#endif
