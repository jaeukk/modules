#ifndef VORONOI_INCLUDED
#define VORONOI_INCLUDED

#include <set>
#include <PeriodicCellList.h>
#ifdef USE_MATHGL
#include <mgl2/mgl.h>
#endif
//calculate all vertexes of the voronoi cell of particle i, in vertices
//if pLinkage is not nullptr, it will be filled with linkage information of the format:
// (*pLinkage)[i*(Dimension+1)+j] is the index of a vertex linked with vertex i
// (*pLinkage)[i*(Dimension+1)+j] is -1 if vertex i have less than j linked vertexes.
void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage = nullptr, std::set<size_t> * pNeighbors = nullptr);

//calculate all vertexes of the voronoi cell of all particles under the periodic boundary condition.
//v is a set of all vertices of Voronoi cells
//pNeighbor[i] is a set of vertices indices of a cell that contains ith particle.
void GetVoronoiCell(const PeriodicCellList<Empty> & c, std::vector<GeometryVector> &v, std::vector< std::vector<size_t>> &pNeighbors);

double PolygonVolume(const std::vector<GeometryVector> & Vertices, DimensionType d);
//color = edge color
void PlotVoronoi2D(std::string prefix, const Configuration & List, const char * edgeColor = "r-", const std::string & remark = "");
//color = edge color
void PlotVoronoi2D(std::string prefix, const SpherePacking & List, const char * sphereColor = "b", const char * edgeColor = "r-" , const std::string & remark = ""); //plot hardsphere at the same time

//Num=-1 : display voronoi of all particles Num>=0 : display voronoi of a specific particle
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius=0.0, long Num=-1);
void DisplayConfiguration_3D_WithVoronoi(const SpherePacking & list, double Radius=0.0, long Num = -1);


//Probability density function of the voronoi cell volume
//clear Result, then fill it with the PDF calculated from Config
//Result is filled with 4-dimensional GeometryVectors, the elements are:
//(r, p(r), \delta r, \delta p(r) )
//gets SampleSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize=500000, double ResolutionPreference=1.0);

void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result);

void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision);

/* Platform to use VORO++ routines
	Configurations are separated into several subdomains
	Perform Voronoi tessellations within each subdomain
*/
void VoronoiCellRoutine(const PeriodicCellList<Empty> &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids);
void VoronoiCellRoutine_radial(const SpherePacking &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids, bool radial = false);

bool transformBasis(const PeriodicCellList<Empty> & input, Configuration & output);
#if defined(USE_VORO)
#include <Discrepancy.h>
#if defined(DELLA)
#include <voro++.hh>
#else
#include <voro++/voro++.hh>
#endif
struct IDs {
	int blockId, order;
	IDs() : blockId(0), order(0)
	{}
};

//This is a class only for exploiting ::put_locate_block
class container_leak : public voro::container_periodic
{
public:
	Configuration c_;	//New fundamental cell
	Configuration c0;	//Fundamental cell of the original config;
	bool transform;//true :when basis vectors were changed
	std::vector<IDs> id_leak;


	container_leak(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int mem) : voro::container_periodic(bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, mem)
	{	}

	container_leak(const PeriodicCellList<Empty> & cInit, const std::vector<GeometryVector> & newBasis, bool trans, int nx, int ny, int nz, int mem=6);
	container_leak(const container_leak & source);

	//same as ordinary put, but now it return a block ID ijk;
	inline int put_ijk(int n, double x, double y, double z) {
		int ijk = 0;
		put_locate_block(ijk, x, y, z);
		id[ijk][co[ijk]] = n;
		double *pp = p[ijk] + 3 * co[ijk]++;
		*(pp++) = x; *(pp++) = y; *pp = z;
		return ijk;
	}
	
	void GetVORONOICELL(int i, std::vector<GeometryVector> &v, std::vector<signed long> * pLinkage = nullptr, std::set<size_t> * pNeighbors = nullptr);
};

void GetConversionParam(const PeriodicCellList<Empty> &c, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem=6);
void GetConversionParam(const PeriodicCellList<Empty> &c, int approxNumPrts, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem = 6);
container_leak Config2Container(const PeriodicCellList<Empty> &c);
#endif


#include <complex>
//return Psi_6 of particle j in Config
//in 2D, Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(6i \\theta) ).
//in 3D, an additional parameter, m, has to be provided. Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(m * i \\phi)*P_6^m(cos(\\theta)) ), where P_l^m is the normalized associated Legendre polynomial.
//definition of neighbors depend on NeighborType
//0: six nearest neighbors in 2D
//1: neighbor in Voronoi diagram
//where theta is the angle between the bond and some reference direction
std::complex<double> Psi6(const PeriodicCellList<Empty> & Config, size_t j, int NeighborType = 0, signed int m=0);




#endif