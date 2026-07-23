#ifndef ClusterMerger_INCLUDED
#define ClusterMerger_INCLUDED

/** \file ClusterMerger.h
 *  \brief Incremental cluster merger over a PeriodicCellList (merges particles
 *         into clusters as a connection distance grows).  Header-only; split
 *         out of PeriodicCellList.h for readability. */

#include <list>
#include "PeriodicCellList.h"   // PeriodicCellList<Empty>, Empty

/** \brief Incrementally merges particles of a PeriodicCellList into clusters as a connection distance grows */
//This class initializes by assigining each particle in a PeriodicCellList to its own cluster
//Then, as the user call member function Proceed(double r), particles within distance r are merged to a single cluster
//The user can then query various member objects of this class for various aspects of clusters:
//member clusters is a list of all clusters. Each cluster is represented as a list of its members' index
//affiliation[i] is the cluster particle i currently belong to
class ClusterMerger
{
protected:
	typedef std::list<size_t> cluster;
private:
	/** \brief Forbid copy
	 *  \param src Source object (copy is forbidden; declaration only). */
	//forbid copy
	ClusterMerger(const ClusterMerger & src);
	struct Bond
	{
		double l2;
		size_t i, j;
	};
	std::vector<Bond> vBond;
	std::vector<Bond>::iterator next;
	size_t _N;
	std::list<cluster> _clusters;
	std::vector<std::list<cluster>::iterator> _affiliation;
	size_t _LargestClusterSize;
public:
	const size_t & N=_N;

//clusters is a list of all clusters. Each cluster is represented as a list of its members' index
	const std::list<cluster> & clusters=_clusters;

//affiliation[i] is the cluster particle i currently belong to
	const std::vector<std::list<cluster>::iterator> & affiliation=_affiliation;
	const size_t & LargestClusterSize=_LargestClusterSize;


public:
	/** \brief Build the bond list of a configuration, placing each particle in its own cluster.
	 *  \param c    Configuration whose particles are merged into clusters.
	 *  \param Rmax Maximum bond distance considered when building the bond list. */
	ClusterMerger(const PeriodicCellList<Empty> & c, double Rmax)
	{
		_N = c.NumParticle();
		for (size_t i = 0; i < _N; i++)
		{
			//at the beginning, each particle consist its own cluster
			cluster temp;
			temp.insert(temp.begin(), i);
			_clusters.insert(_clusters.begin(), temp);
			_affiliation.push_back(_clusters.begin());
			std::vector<Bond> & vb = vBond;
			//generate a vector of bonds for future computation
			c.IterateThroughNeighbors(c.GetRelativeCoordinates(i), Rmax, [&vb, &i](const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) ->void
			{
				if (i < Sourceparticle)
				{
					Bond temp;
					temp.i = i;
					temp.j = Sourceparticle;
					temp.l2 = shift.Modulus2();
					vb.push_back(temp);
				}
			});
		}
		std::sort(vBond.begin(), vBond.end(), [](const Bond & left, const Bond & right) ->bool
		{
			return left.l2 < right.l2;
		});
		_LargestClusterSize = (_N == 0) ? 0 : 1;
		next = vBond.begin();
	}

	/** \brief Merge _clusters until distances of any cluster is larger than r
	 *  \param r Connection distance up to which particles are merged into clusters. */
	//merge _clusters until distances of any cluster is larger than r
	void Proceed(double r)
	{
		double r2 = r*r;
		while (next != vBond.end() && next->l2 < r2)
		{
			std::list<cluster>::iterator pi = _affiliation[next->i], pj = _affiliation[next->j];
			if (pi != pj)
			{
				//merge the two _clusters
				cluster & cj = *_affiliation[next->j];
				cluster & ci = *_affiliation[next->i];
				//1. change the _affiliation of all particles in the same cluster as j
				for (auto iter = cj.begin(); iter != cj.end(); ++iter)
				{
					_affiliation[*iter] = _affiliation[next->i];
				}
				//2. merge the two _clusters
				ci.splice(ci.begin(), cj);
				_clusters.erase(pj);
				_LargestClusterSize = std::max(_LargestClusterSize, ci.size());
			}
			next++;
		}
	}
};

#endif
