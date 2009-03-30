/*!
\file cluster.h

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

\brief cluster class

 */

#ifndef _CLUSTER_H
#define _CLUSTER_H

#include <set>
#include "segment.h"
#include <alize.h>
#include <cstdio>
#include <iostream>

using namespace std;
using namespace alize;
using namespace lium_seg;

namespace lium_seg {
class ClusterServer;
	/// Sorted set of seg.
	typedef set	< Segment, Segment > SegLst;
	
	/// \brief Container for the storage of segments.
	///
	/// - A cluster generally corresponds to speaker.
	/// - This is a container of stl container (a set).
	///
	/// \warning the segment is added directly in the container,
	/// ie no add, del, etc function.
	///
	/// \todo use enum in setter and getter
	/// \e getContainer() in protected, do method to add/remove segments.
	
	
class Cluster {
		friend class ClusterServer; 
public:
		
		Cluster(const long & _id);

		Cluster(const Cluster & _c);

		Cluster();

		Cluster & operator=(const Cluster & _c);

		/// \brief Set the cluster index.
		/// \param _id the \# of the cluster.
		void setId(const long & _id);

		/// \brief Get the cluster index.
		/// \return the \# of the cluster.
		long getId(void);

		/// \brief Collapse segments form a cluster that are contiguous.
		/// Very useful when all segments have a length of 1 frame.
		void collapse(void);

		static void collapse(SegLst & segLst);

		bool find(unsigned long idx);

		/// \brief Information of the cluster.
		//void debug(void);
		
		unsigned long getLength() { 
			unsigned long len = 0;
			for(seekBeginSeg(); isEof() == false; seekNextSeg()){
				len += getCurrentSeg().getLength();
			}
			return len;
		}
		
		void seekBeginSeg() { 
			it = lst.begin();
		}

		void seekEndSeg() { 
			
			it = --lst.end();
		}

		bool isEof() { 
			return (it == lst.end());
		}

		void seekNextSeg() { 
			it++;
		}

		void seekPrevSeg() { 
			it--;
		}

		void eraseCurrentSeg() { 
			SegLst::iterator sv_it = it;
			sv_it--;
			lst.erase(it);
			it = sv_it;
		}
		
		Segment copyCurrentSeg() { 
			return *it; 
		}
		const Segment & getCurrentSeg() { 
			return *it; 
		}
		void insertSeg(const Segment & s) {
			it = lst.insert(s).first;
		}
		
		const long size(void) { return lst.size();}

protected:
		SegLst lst;		//!< sorted list (a set) of segments, sorted by start time
		XList info;		//!< list of informations
		long id;		//!< # of cluster
		SegLst::iterator it;
		/// \brief Get the container.
		///
		/// Useful to add, del segments.
		/// \return the container, ie the list of the segments.
		SegLst & getContainer(void);
		const SegLst & getContainer(void) const;
private:
		bool operator==(const Cluster&) const; /*! not implemented */
		bool operator!=(const Cluster&) const; /*! not implemented */
	};
}

#endif
