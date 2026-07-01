#include "PeriodicCellList.h"
#include <exception>




/** \brief Tile a Configuration NumberInEachSide times along every box direction. 
 *  \param src             Source configuration to tile.
 *  \param NumberInEachSide Number of replicas along each box direction.
 *  \return                The tiled configuration. */
Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide)
{
	double dNumberInEachSide=static_cast<double>(NumberInEachSide);
	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i))*(1.0/dNumberInEachSide);
		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide)
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide;
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide)
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}

/** \brief Tile a Configuration with per-direction replication counts. 
 *  \param src             Source configuration to tile.
 *  \param NumberInEachSide Per-direction number of replicas.
 *  \return                The tiled configuration. */
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide)
{
	assert(NumberInEachSide.size()==src.GetDimension());
	std::vector<double> dNumberInEachSide;
	for(auto iter=NumberInEachSide.begin(); iter!=NumberInEachSide.end(); iter++)
		dNumberInEachSide.push_back((double)(*iter));

	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide[i]);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i));
		for(int i=0; i< src.GetDimension(); i++)
			OriginalRelative.x[i]/=dNumberInEachSide[i];

		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide[dim-1])
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide[j];
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide[j])
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}

/** \brief Return an empty unit cubic box configuration of dimension d. 
 *  \param d        Spatial dimension of the box.
 *  \param CellSize Target linear cell size.
 *  \return         An empty unit cubic box configuration. */
Configuration GetUnitCubicBox(DimensionType d, double CellSize)
{
	std::vector<GeometryVector> vbas;
	for(DimensionType i=0; i<d; i++)
	{
		GeometryVector t(d);
		t.x[i]=1.0;
		vbas.push_back(t);
	}
	Configuration result(d, &vbas[0], CellSize);

	return result;
}



/** \brief Open (or attach to) the pack file identified by prefix. 
 *  \param prefix File name prefix identifying the pack files. */
void ConfigurationPack::Open(const std::string & prefix)
{
	this->IndexName = prefix + ".ConfigPackIndex";
	this->PackName = prefix + ".ConfigPack";
	std::fstream ifile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	if (ifile.good())
	{
		ifile.read((char*)(&NConfig), sizeof(NConfig));
		if (ifile.eof())
			NConfig = 0;
	}
	else
		NConfig = 0;
}

/** \brief Return the number of configurations currently in the pack. */
long long ConfigurationPack::NumConfig(void)
{
	return NConfig;
}

/** \brief Append a configuration to the pack files. 
 *  \param c Configuration to append. */
void ConfigurationPack::AddConfig(const Configuration & c)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary);
	if (indexFile.good() == false)
	{
		indexFile.clear();
		indexFile.close();
		indexFile.clear();
		indexFile.open(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	}
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::ate);
	if (packFile.good() == false)
	{
		packFile.clear();
		packFile.close();
		packFile.clear();
		packFile.open(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	}
	this->NConfig++;
	indexFile.seekp(0, std::fstream::beg);
	indexFile.write((char*)(&NConfig), sizeof(NConfig));
	long long loc = packFile.tellp();
	indexFile.seekp(0, std::fstream::end);
	indexFile.write((char*)(&loc), sizeof(loc));
	c.WriteBinary(packFile);
	indexFile.close();
	if (indexFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Index file failing!\n";
	packFile.close();
	if (packFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Pack file failing!\n";
}

/** \brief Read and return the i-th configuration from the pack. 
 *  \param i Index of the configuration to read.
 *  \return  The i-th configuration in the pack. */
Configuration ConfigurationPack::GetConfig(long long i)
{
	if (i >= this->NConfig)
	{
		std::cerr << "Error in ConfigurationPack::GetConfig : Index is larger than number of configuration!\n";
		exit(1);
	}
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::binary);
	long long loc;
	indexFile.seekg(sizeof(long long)*(i + 1));
	indexFile.read((char*)(&loc), sizeof(loc));
	packFile.seekg(std::streampos(loc));
	return packFile.good() ? Configuration(packFile) : Configuration();
}

/** \brief Discard all configurations that are already in the file. */
void ConfigurationPack::Clear(void)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	this->NConfig = 0;
}

/** \brief Pixelizing the configuration to MeshSide pixes in each side. 
 *  \param[in]  c        Configuration to pixelize.
 *  \param[in]  radius   Sphere radius placed on each particle.
 *  \param[in]  MeshSide Number of pixels along each box direction.
 *  \param[out] occu     Filled with the occupancy grid (1 = occupied). */
//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & occu)
{
	DimensionType dim = c.GetDimension();
	size_t MeshTotal = std::pow(MeshSide, dim);
	struct VicinityMesh
	{
		signed long DeltaIndex[::MaxDimension] = {};
	};

	//get the list of mesh points occupied by a sphere
	std::vector<VicinityMesh> mesh;

	GeometryVector newbas[::MaxDimension] = {};
	for (int i = 0; i<dim; i++)
		newbas[i] = c.GetBasisVector(i)*(1.0 / MeshSide);
	Configuration * pt = new Configuration(dim, newbas, ::MaxDistance, false);
	pt->Insert("A", GeometryVector(dim));
	pt->IterateThroughNeighbors(GeometryVector(dim), radius, [&](const GeometryVector & s, const GeometryVector & l, const signed long * p, size_t source) ->void
	{
		VicinityMesh temp;
		for (int i = 0; i<dim; i++)
			temp.DeltaIndex[i] = p[i];
		mesh.push_back(temp);
	});
	delete pt;

	occu.assign(MeshTotal, 0);
	signed long m = MeshSide;
	int num = c.NumParticle();
	for (int i = 0; i<num; i++)
	{
		GeometryVector rel = c.GetRelativeCoordinates(i);
		for (auto iter = mesh.begin(); iter != mesh.end(); ++iter)
		{
			size_t index = 0;
			for (int j = 0; j<dim; j++)
			{
				index *= MeshSide;
				signed long temp = std::floor(rel.x[j] * MeshSide);
				temp += iter->DeltaIndex[j];
				while (temp >= m)
					temp -= m;
				while (temp<0)
					temp += m;
				index += temp;
			}
			occu[index] = 1;
		}
	}
}

/** \brief Calculate volume fraction if each particle is replaced with a hypersphere of certain radius; do this by pixelizing the configuration to MeshSide pixes in each side. 
 *  \param c        Configuration whose occupied volume is measured.
 *  \param radius   Sphere radius placed on each particle.
 *  \param MeshSide Number of pixels along each box direction.
 *  \param pbuffer  Optional scratch buffer for the occupancy grid; allocated internally if null.
 *  \return         The occupied (covered) volume. */
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer)
{
	size_t MeshTotal = std::pow(MeshSide, c.GetDimension());
	size_t InsideCount=0;
	if (pbuffer == nullptr)
	{
		std::vector<VoxelType> occu;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	else
	{
		std::vector<VoxelType> & occu = *pbuffer;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	return c.PeriodicVolume()*InsideCount/MeshTotal;
}


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <set>

/** \brief Return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected. 
 *  \param c Configuration of particles.
 *  \param a Index of the first particle.
 *  \param b Index of the second particle.
 *  \return  Minimum sphere diameter D that connects particles a and b. */
//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b)
{
	double clength = 0.0;
	double mlength = 2 * std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / c.GetDimension());
	std::set<size_t> cluster;
	struct toSearchStruct
	{
		double dist;
		size_t i;
		toSearchStruct(double dist, size_t i) : dist(dist), i(i)
		{}
	};
	typedef boost::multi_index_container<
		toSearchStruct,
		boost::multi_index::indexed_by<
		boost::multi_index::ordered_unique<boost::multi_index::member<toSearchStruct, size_t, &toSearchStruct::i> >,
		boost::multi_index::ordered_non_unique<boost::multi_index::member<toSearchStruct, double, &toSearchStruct::dist> >
		>
	> toSearch_Set;
	toSearch_Set toAdd;
	toSearch_Set::nth_index<0>::type & iview = toAdd.get<0>();
	toSearch_Set::nth_index<1>::type & dview = toAdd.get<1>();

	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		if (cluster.find(SourceAtom) == cluster.end())
		{
			double dist = std::sqrt(shift.Modulus2());
			auto iter = iview.find(SourceAtom);
			if (iter == iview.end())
				toAdd.insert(toSearchStruct(dist, SourceAtom));
			else if (iter->dist > dist)
			{
				toSearchStruct temp(dist, SourceAtom);
				toAdd.replace(iter, temp);
			}
		}
	};
	cluster.insert(a);
	GeometryVector rel = c.GetRelativeCoordinates(a);
	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
	for (;;)
	{
		if (toAdd.empty())
		{
			mlength *= 2;
			for (auto iter = cluster.begin(); iter != cluster.end(); iter++)
			{
				GeometryVector rel = c.GetRelativeCoordinates(*iter);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
			}
			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. toAdd size="<<toAdd.size()<<'\n';
		}
		else
		{
			auto ibegin = dview.begin();
			toSearchStruct temp = *ibegin;
			dview.erase(ibegin);
			clength = std::max(clength, temp.dist);
			size_t np = temp.i;
			if (np == b)
				return clength;
			if (cluster.find(np) == cluster.end())
			{
				cluster.insert(np);
				GeometryVector rel = c.GetRelativeCoordinates(np);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
				//std::cerr<<"Explore particle "<<np<<'\n';
			}
		}
	}

	//shouldn't reach this line
	std::cerr << "Error in MinConnectDiameter\n";
	return 0.0;
}
//double oldMinConnectDiameter(const Configuration & c, size_t a, size_t b)
//{
//	double clength=0.0;
//	double mlength=1.5*std::pow(c.PeriodicVolume()/c.NumParticle(), 1.0/c.GetDimension());
//	std::set<size_t> cluster;
//	std::multimap<double, size_t> toAdd;
//	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
//	{
//		if(cluster.find(SourceAtom)==cluster.end())
//		{
//			double dist=std::sqrt(shift.Modulus2());
//			toAdd.insert(std::make_pair(dist, SourceAtom));
//		}
//	};
//	cluster.insert(a);
//	GeometryVector rel=c.GetRelativeCoordinates(a);
//	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//	for(;;)
//	{
//		if(toAdd.empty())
//		{
//			mlength*=2;
//			for(auto iter=cluster.begin(); iter!=cluster.end(); iter++)
//			{
//				GeometryVector rel=c.GetRelativeCoordinates(*iter);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. \n";
//		}
//		else
//		{
//			auto ibegin=toAdd.begin();
//			std::pair<double, size_t> temp=*ibegin;
//			toAdd.erase(ibegin);
//			clength=std::max(clength, temp.first);
//			size_t np=temp.second;
//			if(np==b)
//				return clength;
//			if(cluster.find(np)==cluster.end())
//			{
//				cluster.insert(np);
//				GeometryVector rel=c.GetRelativeCoordinates(np);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//		}
//	}
//}

/** \brief Return the minimum D such that if each particle is replaced with a sphere with diameter D, particle n will connect to its periodic image in direction dir. 
 *  \param c   Configuration of particles.
 *  \param n   Index of the particle.
 *  \param dir Box direction along which percolation is tested.
 *  \return    Minimum sphere diameter D that connects particle n to its periodic image. */
//return the minimum D such that if each particle is replaced with a sphere with diameter D,
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir)
{
	std::vector<size_t> m(c.GetDimension(), 1);
	m[dir] = 2;
	Configuration t = MultiplicateStructure(c, m);
	return MinConnectDiameter(t, 2 * n, 2 * n + 1);
}

//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, 
//which means the cluster size is actually infinite rather than the returned value.
#include <list>
/** \brief Replace each particle in c with a hypersphere with radius, what's the cluster size InitParticle is in? FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, which means the cluster size is actually infinite rather than the returned value. 
 *  \param[in]  c                  Configuration of particles.
 *  \param[in]  InitParticle       Index of the seed particle whose cluster is measured.
 *  \param[in]  radius             Sphere radius placed on each particle.
 *  \param[out] FoundPeriodicImages Set true if more than one periodic image of a particle was found (cluster is effectively infinite).
 *  \return                        Number of particles in the cluster containing InitParticle. */
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages)
{
	FoundPeriodicImages = false;
	std::map<size_t, GeometryVector> Found;
	struct toSearchStruct
	{
		size_t i;
		GeometryVector LatticeShift;
		toSearchStruct(size_t i, const GeometryVector & LatticeShift) : i(i), LatticeShift(LatticeShift)
		{}
	};
	std::list<toSearchStruct> toSearch;
	GeometryVector alreadyLatticeShift = GeometryVector(c.GetDimension());
	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		auto iter = Found.find(SourceAtom);
		GeometryVector l = alreadyLatticeShift + LatticeShift;
		if (iter != Found.end())
		{
			if (!FoundPeriodicImages)
			{
				//this particle has been found before. Test if they belong to the same periodic cell
				GeometryVector deltaLatticeShift = l - iter->second;
				if (deltaLatticeShift.Modulus2() > LengthPrecision*LengthPrecision)
					FoundPeriodicImages = true;
			}
		}
		else
		{
			toSearch.push_back(toSearchStruct(SourceAtom, l));
			Found.insert(std::make_pair(SourceAtom, l));
		}
	};
	Found.insert(std::make_pair(InitParticle, GeometryVector(c.GetDimension())));
	c.IterateThroughNeighbors(c.GetRelativeCoordinates(InitParticle), radius, IterateFunction);
	while (!toSearch.empty())
	{
		toSearchStruct s = toSearch.front();
		toSearch.pop_front();
		alreadyLatticeShift = s.LatticeShift;
		c.IterateThroughNeighbors(c.GetRelativeCoordinates(s.i), radius, IterateFunction);
	}
	return Found.size();
}

/** \brief Find all reciprocal lattice vectors (wavevectors) of a given periodic simulation box. 
 *  \param tempList          A periodic simulation box.
 *  \param CircularKMax      Upper magnitude bound for the circularly-selected wavevectors.
 *  \param LinearKMax        Upper magnitude bound for the integer multiples of the selected wavevectors.
 *  \param SampleProbability When <1, randomly keep wavevectors within CircularKMax with this probability.
 *  \param use_iteratethrough True = use IterateThroughNeighbors() (may surge memory use).
 *  \param CircularKMin      Lower magnitude bound for the circularly-selected wavevectors (0 by default).
 *  \return                  The selected reciprocal lattice vectors (wavevectors). */
std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability, bool use_iteratethrough, double CircularKMin)
{
	//Exclude probablematic cases.
	if (CircularKMax <= 0.0) {
		std::cerr << "GetKs():: CircularKMax should be positive." << std::endl;
		return std::vector<GeometryVector>();
	}
	if (LinearKMax <= 0.0) {
		std::cerr << "GetKs():: Linear KMax should be positive." << std::endl;
		return std::vector<GeometryVector>();
	}

	RandomGenerator gen(98765);
	std::vector<GeometryVector> results;
	DimensionType dim = tempList.GetDimension();
	if (dim == 0) {
		std::cerr << "Dimension of PeriodicCellList must be positive!\n";
		exit(1);
	}
	else if (dim > 1)
	{
		std::vector<GeometryVector> bs;
		for (DimensionType i = 0; i<dim; i++)
			bs.push_back(tempList.GetReciprocalBasisVector(i));
		
		//a periodic list of reciprocal lattice
		PeriodicCellList<Empty> reciprocal(dim, &bs[0], std::sqrt(bs[0].Modulus2()));
		reciprocal.Insert(Empty(), GeometryVector(dim));
		std::vector<GeometryVector> ks;
		std::vector<GeometryVector> preKs;
		double KMin2 = CircularKMin*CircularKMin;

		if (use_iteratethrough)
		{
			reciprocal.IterateThroughNeighbors(GeometryVector(dim), CircularKMax, [&ks, &dim, &KMin2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{
			bool Trivial = true;
			for (DimensionType d = 0; d<dim; d++)
			{
				if (PeriodicShift[d]>0)
				{
					Trivial = false;
					break;
				}
				else if (PeriodicShift[d]<0)
				{
					Trivial = true;
					break;
				}
			}
			if (!Trivial && (shift.Modulus2() > KMin2))
				ks.push_back(shift);
			});
		}
		else{
			/* Get wavevectors without IterateThroughNeighbors */
			size_t num_ks = (size_t) round(0.5*HyperSphere_Volume(dim,CircularKMax) / reciprocal.PeriodicVolume());
			ks.reserve(num_ks);

			std::vector<int> limits;
			std::vector<int> PeriodicShift;
			size_t iter = 1;
			for (size_t i =0; i < dim; i++){
				int val = int(std::ceil( CircularKMax* std::sqrt( tempList.GetBasisVector(i).Modulus2()) /(2.*::pi) ));
				limits.push_back(val);
				PeriodicShift.push_back(-val);
				iter *= (2*val + 1);
			}

			/* iterate over all reciprocal lattice points */
			double K2 = CircularKMax*CircularKMax;
			for (size_t i=0; i < iter; i++){
				bool Trivial = false;
				
				for (DimensionType d = 0; d<dim; d++)
				{
					if (PeriodicShift[d]>0)
					{
						Trivial = false;
						break;
					}
					else if (PeriodicShift[d]<0)
					{
						Trivial = true;
						break;
					}
				}
				if (!Trivial){
					/* Ignore inverted vectors */
					GeometryVector k(dim);
					for (DimensionType d = 0; d<dim; d++)
					{
						k = k + PeriodicShift[d] * bs[d];
					} 
					double k2 = k.Modulus2();
					if ((k2 < K2) && (k2 > KMin2))
						ks.emplace_back(k);
				}

				/* Change PeriodicShift */
				PeriodicShift[0] ++;
				for (DimensionType d = 0; d<dim; d++)
				{	
					if (PeriodicShift[d] > limits[d]){
						PeriodicShift[d] = - limits[d];
						PeriodicShift[d+1] ++ ; 
						/* Note that PeriodicShift[dim-1] cannot exceed limits[dim-1] */ 
					}
				} 
			}	
		}
		std::sort(ks.begin(), ks.end(), [](const GeometryVector & left, const GeometryVector & right) ->bool {return (left.Modulus2()<right.Modulus2()) && (left.x[0]<right.x[0]); });
		//ks contains K points in that circle

		if (LinearKMax == CircularKMax)
			return ks;

		std::vector<GeometryVector> tempBase;
		for (DimensionType i = 0; i<dim; i++)
		{
			GeometryVector t(dim);
			t.x[i] = 2;
			tempBase.push_back(t);
		}
		PeriodicCellList<Empty> KDirection(dim, &tempBase[0], std::sqrt(tempBase[0].Modulus2())*std::pow(1.0 / ks.size(), 1.0 / (dim))*10.0);//use this list to discover k points of the same direction
		for (auto iter = ks.begin(); iter != ks.end(); iter++)
		{
			double Length = std::sqrt(iter->Modulus2());
			if (Length == 0)
				continue;
			GeometryVector dir(dim);
			for (DimensionType i = 0; i<dim; i++)
				dir.x[i] = iter->x[i] / Length / 2.0;

			bool NeighborFound = false;
			KDirection.IterateThroughNeighbors(dir, 1e-10, [&NeighborFound](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{NeighborFound = true; });
			if (NeighborFound == false)
			{
				preKs.push_back(*iter);
				KDirection.Insert(Empty(), dir);
			}
		}
		//preKs contains K points in that circle with different directions

		ks.clear();
		for (auto iter = preKs.begin(); iter != preKs.end(); iter++)
		{
			for (size_t i = 1;; i++)
			{
				GeometryVector temp = static_cast<double>(i)*(*iter);
				if (temp.Modulus2()>LinearKMax*LinearKMax)
					break;
				else
					if(SampleProbability>=1.0 || gen.RandomDouble()<SampleProbability)
						results.push_back(temp);
			}
		}
		preKs.clear();
	}
	else
	{
		/* 1D case */
		//std::cout<< "GetKs(...):: 1D case. \n";
		GeometryVector b = tempList.GetReciprocalBasisVector(0);
		if (CircularKMin <= 0.0){
			for (double i = 1.0;; i += 1.0)
			{
				GeometryVector k = i*b;
				if (k.x[0] > LinearKMax)
					break;
				if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
					results.push_back(k);
			}
		}
		else{
			/* search wavevectors between CircularKMin and CircularKMax */
			//std::cout<< "GetKs(...):: CircularKMin is positive. \n";
			size_t num_inside_circular = 0;
			for (double i = 1.0;; i += 1.0)
			{
				GeometryVector k = i*b;
				if (k.x[0] > CircularKMax)
					break;
				if (k.x[0] > CircularKMin){
					if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
						results.push_back(k);
					num_inside_circular ++;
				}
			}
			
			if (results.size() > 0){
				/* consider integer multiples in results. */
				bool OutsideRange = false;
				for (double i = 2.0;; i+= 1.0){
					for (size_t j = 0; j<num_inside_circular; j++){
						GeometryVector k = i * results[j];
						if (k.x[0] > LinearKMax){
							OutsideRange = true;
							break;
						}
						if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
							results.push_back(k);
					}
					if (OutsideRange)
						break;
				}
			}
			else{
				std::cerr << "GetKs(...):: There is no wavevector bewteen CircularKMin and CircularKMax " <<std::endl;
			}

		}
	}
	return results;
}


#ifdef VORONOI_INCLUDED
#include "Voronoi.h"
#include <gsl/gsl_sf_legendre.h>
/** \brief Evaluate the spherical-harmonic-normalized associated Legendre function P_l^m(x), including negative m. 
 *  \param l Degree of the associated Legendre function.
 *  \param m Order of the associated Legendre function (may be negative).
 *  \param x Argument in [-1, 1].
 *  \return  The spherical-harmonic-normalized value P_l^m(x). */
double LegenderSphericalPlm(int l, int m, double x)
{
	if (m >= 0)
		return gsl_sf_legendre_sphPlm(l, m, x);
	else
	{
		int mm = (-1)*m;
		double result = gsl_sf_legendre_sphPlm(l, mm, x);
		if (mm % 2 == 1)
			result *= (-1.0);
		for (int t = l - mm + 1; t <= l + mm; t++)
			result /= t;
		return result;
	}
}




#endif
