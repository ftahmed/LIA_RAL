/*!
\file clusters.h

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

\brief clusters class

*/

#ifndef _CLUSTER_SERVER_H
#define _CLUSTER_SERVER_H

#include <string>
#include <map>
#include <alize.h>
#include "cluster.h"

using namespace std;
using namespace alize;
using namespace lium_seg;

namespace lium_seg {

	/// Map container for a cluster \e long get the \e Cluster.
	typedef map < long, Cluster > ClusterMap;
	/// Map container to store the list of the shows.
	typedef map < String, long > SrcNameMap;
	/// Map container to store the list of the cluster name.
	typedef map < String, long > ClusterNameMap;

/// \brief Container for the storage of clusters.
///
/// \li A clusters corresponds to segmentation.
/// \li This is a container of stl container (a map).
/// \li it maintains also a map of shows. 
///
/// \warning the cluster is added directly in the container,
/// ie no add, del, etc functions are provided, only a merge
/// function between two clusters.
///
/// \todo 
///	\li move \e getContainer() in protected, do method to add/remove segments.
/// \li read/write: use \e IoFile instead of direct \e IOFILE *.


// show = srcName
// name = clusterName

class ClusterServer {
public:
		//Constructors
		ClusterServer(Config & config);
		ClusterServer(void);
		ClusterServer(Config & config, const ClusterServer & _clusters);
		ClusterServer(Config & config, const SrcNameMap & sm, long si, const ClusterNameMap & nm, long ni);

		ClusterServer & operator=(const ClusterServer & _clusters);
		ClusterMap & getContainer(void);
		const ClusterMap & getContainer(void) const;

		///  \brief Get max Idx Show
		long getMaxIdxSrcName(void) const { return maxIdxSourceName; }
		
		/// \brief Get the container show map.
		const SrcNameMap & getSrcNameMap(void) const {
			return srcNameMap;
		}

		/// \brief Get the container show map.
		SrcNameMap & getSrcNameMap(void) {
			return srcNameMap;
		}

		/// \brief Get the absolute index of the show.
		long getIdxSrcName(String& srcName) {
			return srcNameMap[srcName];
		}

		/// \brief Store a show and automatically assign a index number.
		long setIdxSrcName(String srcName) {
			SrcNameMap::iterator res = srcNameMap.find(srcName); 
			if (res == srcNameMap.end()) {
				srcNameMap[srcName] = maxIdxSourceName++;
				res = srcNameMap.find(srcName); 
			}
			return res->second;
		}

		///  \brief Get max Idx name
		long getMaxIdxClusterName(void) const { return maxIdxClusterName; }
		
		/// \brief Get the container name map.
		const ClusterNameMap & getClusterNameMap(void) const {
			return clusterNameMap;
		}
		
		/// \brief Get the container name map.
		ClusterNameMap & getClusterNameMap(void) {
			return clusterNameMap;
		}
		
		/// \brief Get the absolute index of the show.
		long getIdxClusterName(String name) {
			return clusterNameMap[name];
		}

		/// \brief Get the config.
		Config& getConfig() {
			return *pconfig;
		}

		/// \brief Store a name and automatically assign a index number.
		long setIdxClusterName(String name) {
			ClusterNameMap::iterator res = clusterNameMap.find(name); 
			if (res == clusterNameMap.end()) {
				clusterNameMap[name] = maxIdxClusterName++;
				res = clusterNameMap.find(name); 
			}
			return res->second;
		}

		String getName(long i);

		String getSourceName(long i);

		/// \brief Insert the name of the locutor in the map.
		/// \param id the id of the cluster.
		/// \param name the name of the locutor.
		void addClusterName(long id, String name) {
			clusterNameMap[name]=id;
		}

		/// \brief Add a cluster int the map of clusters.
		/// \param id the id of the cluster.
		/// \param name the name of the locutor.
		void addCluster(long id, Cluster cluster) {
			cluster.setId(id);
			clusterMap[id]=cluster;
		}

		/// \brief Store the config.
		void setConfig(Config& config) {
			pconfig = &config;
		}

		/// \brief Erase a cluster.
		/// \param j the cluster to erase.
		void erase(long j);

		/// \brief Merge two clusters.
		/// \param i the first cluster to merge.
		/// \param j the second cluster to merge.
		void merge(long i, long j);

		/// \brief Collapse segments form a cluster that are contiguous.
		/// \see Cluster::collapse(void).
		void collapse(void);
		
		/// \brief Put the name of a source in the List		
		XLine & getSourceNameInXLine();

		//void debug(void);

		/// \brief Read a list of segmentation file.
		void read(const String & sourceName);

		/// \brief Read a segmentation.
		void read(FILE * f);

		/// \brief Write a segmentation.
		void write(const String & sourceName);

		/// \brief Write a segmentation.
		void write(FILE * f);

		SegLst getSegLst(bool collapse=true);

protected:
		Config * pconfig;
		ClusterMap clusterMap; 	//!< container of the clusters.
		SrcNameMap srcNameMap;	//!< container of the show.
		ClusterNameMap clusterNameMap;	//!< container of the locutor name.	
		long maxIdxSourceName;	//!< maximum of the absolute show index.
		long maxIdxClusterName;	//!< maximum of the absolute name index.
private:
		bool operator==(const ClusterServer&) const; /*! not implemented */
		bool operator!=(const ClusterServer&) const; /*! not implemented */
	};
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// TODO if we need
//MapSeg clustersToFrames(ClusterServer & c);
//VectorSeg clustersToVector(ClusterServer & c);
//ClusterServer vectorToClusters(VectorSeg & vSeg, ClusterServer & clusters, bool collapse);
	
}
#endif
