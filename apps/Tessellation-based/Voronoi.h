#ifndef VORONOI_INCLUDED
#define VORONOI_INCLUDED

#include <set>
#include <PeriodicCellList.h>
#ifdef USE_MATHGL
#include <mgl2/mgl.h>
#endif
/** \brief calculate all vertexes of the voronoi cell of particle i, in vertices
 *  if pLinkage is not nullptr, it will be filled with linkage information of the format:
 *  (*pLinkage)[i*(Dimension+1)+j] is the index of a vertex linked with vertex i;
 *  (*pLinkage)[i*(Dimension+1)+j] is -1 if vertex i have less than j linked vertexes.
 *  \param c          Configuration whose Voronoi cell is computed.
 *  \param i          Index of the particle whose cell is computed.
 *  \param[out] vertices  Vertices of the Voronoi cell of particle i.
 *  \param pLinkage   Optional output: vertex-linkage information (see above).
 *  \param pNeighbors Optional output: indices of the neighboring particles. */
void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage = nullptr, std::set<size_t> * pNeighbors = nullptr);

/** \brief calculate all vertexes of the voronoi cell of all particles under the periodic boundary condition.
 *  \param c          Configuration whose Voronoi cells are computed.
 *  \param[out] v     a set of all vertices of Voronoi cells.
 *  \param[out] pNeighbors  pNeighbor[i] is a set of vertices indices of a cell that contains ith particle. */
void GetVoronoiCell(const PeriodicCellList<Empty> & c, std::vector<GeometryVector> &v, std::vector< std::vector<size_t>> &pNeighbors);

/** \brief Volume (d-content) of the polytope spanned by the given vertices.
 *  \param Vertices  Vertices of the polytope.
 *  \param d         Spatial dimension.
 *  \return          Volume of the polytope. */
double PolygonVolume(const std::vector<GeometryVector> & Vertices, DimensionType d);

/** \brief Plot the 2D Voronoi diagram of a point configuration.
 *  \param prefix     Output filename prefix.
 *  \param List       Configuration to plot.
 *  \param edgeColor  edge color.
 *  \param remark     Annotation text drawn on the figure. */
void PlotVoronoi2D(std::string prefix, const Configuration & List, const char * edgeColor = "r-", const std::string & remark = "");

/** \brief Plot the 2D Voronoi diagram of a sphere packing (plot hardsphere at the same time).
 *  \param prefix       Output filename prefix.
 *  \param List         Sphere packing to plot.
 *  \param sphereColor  Color used for the spheres.
 *  \param edgeColor    edge color.
 *  \param remark       Annotation text drawn on the figure. */
void PlotVoronoi2D(std::string prefix, const SpherePacking & List, const char * sphereColor = "b", const char * edgeColor = "r-" , const std::string & remark = ""); //plot hardsphere at the same time

/** \brief Display a 3D configuration together with its Voronoi cells.
 *  \param list    Configuration to display.
 *  \param Radius  Sphere radius used for rendering particles (0 = automatic).
 *  \param Num     Num=-1 : display voronoi of all particles; Num>=0 : display voronoi of a specific particle. */
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius=0.0, long Num=-1);

/** \brief Display a 3D sphere packing together with its Voronoi cells.
 *  \param list    Sphere packing to display.
 *  \param Radius  Sphere radius used for rendering particles (0 = automatic).
 *  \param Num     Num=-1 : display voronoi of all particles; Num>=0 : display voronoi of a specific particle. */
void DisplayConfiguration_3D_WithVoronoi(const SpherePacking & list, double Radius=0.0, long Num = -1);


/** \brief Probability density function of the voronoi cell volume.
 *  clear Result, then fill it with the PDF calculated from Config.
 *  Result is filled with 4-dimensional GeometryVectors, the elements are
 *  (r, p(r), \delta r, \delta p(r) ); gets SampleSize number of sample pair
 *  distances, use it to generate bins, then count all pair distances.
 *  \param GetConfigsFunction  Callable returning configuration i.
 *  \param NumConfigs          Number of configurations in the ensemble.
 *  \param[out] Result         PDF samples (r, p(r), dr, dp(r)).
 *  \param SampleSize          Number of sample volumes used to build the bins.
 *  \param ResolutionPreference  >1 to get better r resolution, <1 to get better g_2 resolution. */
void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize=500000, double ResolutionPreference=1.0);

/** \brief Distribution of the number of sides of Voronoi cells over an ensemble.
 *  \param GetConfigsFunction  Callable returning configuration i.
 *  \param NumConfig           Number of configurations in the ensemble.
 *  \param[out] Result         (number of sides, probability) pairs. */
void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result);

/** \brief Two-point correlation of Voronoi-cell volumes, C(r)=<(V(x)-<V>)(V(x+r)-<V>)>.
 *  \param GetConfigsFunction  Callable returning configuration i.
 *  \param NumConfigs          Number of configurations in the ensemble.
 *  \param[out] Result         Correlation samples (r, C(r), dr, dC(r)).
 *  \param MaxDistance         Maximum pair distance considered.
 *  \param RPrecision          Bin width in r. */
void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision);

/** \brief Platform to use VORO++ routines.
 *  Configurations are separated into several subdomains; perform Voronoi
 *  tessellations within each subdomain.
 *  \param c           Configuration whose Voronoi cells are computed.
 *  \param numThreads  Number of threads used for the computation.
 *  \param cellVolume  Optional output: per-particle cell volumes.
 *  \param vertices    Optional output: per-particle vertex counts.
 *  \param faces       Optional output: per-particle face counts.
 *  \param centroids   Optional output: per-particle cell centroids. */
void VoronoiCellRoutine(const PeriodicCellList<Empty> &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids);

/** \brief Voronoi cell routine for a sphere packing, optionally radical (radial) tessellation.
 *  \param c           Sphere packing whose Voronoi cells are computed.
 *  \param numThreads  Number of threads used for the computation.
 *  \param cellVolume  Optional output: per-particle cell volumes.
 *  \param vertices    Optional output: per-particle vertex counts.
 *  \param faces       Optional output: per-particle face counts.
 *  \param centroids   Optional output: per-particle cell centroids.
 *  \param radial      If true, use the radical/radial tessellation (currently unimplemented). */
void VoronoiCellRoutine_radial(const SpherePacking &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids, bool radial = false);

/** \brief Construct a lower-triangular representation of the basis when rotation is required.
 * \param input Configuration providing the original basis.
 * \param[out] output Empty configuration receiving the transformed basis when the function returns true.
 * \return True if a transformed basis was written; false if no transformation was required or the dimension is unsupported. */
bool transformBasis(const PeriodicCellList<Empty> & input, Configuration & output);
#if defined(USE_VORO)
#if defined(DELLA)
#include <voro++.hh>
#else
#include <voro++/voro++.hh>
#endif
struct IDs {
	int blockId, order;
	/** \brief Construct with zero block id and order. */
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


	/** \brief Construct an empty periodic container with the given box and grid.
	 *  \param bx_   Box length along x.
	 *  \param bxy_  xy shear component of the box.
	 *  \param by_   Box length along y.
	 *  \param bxz_  xz shear component of the box.
	 *  \param byz_  yz shear component of the box.
	 *  \param bz_   Box length along z.
	 *  \param nx_   Number of grid blocks along x.
	 *  \param ny_   Number of grid blocks along y.
	 *  \param nz_   Number of grid blocks along z.
	 *  \param mem   Initial per-block memory allocation. */
	container_leak(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int mem) : voro::container_periodic(bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, mem)
	{	}

	/** \brief Build a container from a configuration, optionally with a rotated basis.
	 *  \param cInit     Source configuration providing the particles.
	 *  \param newBasis  Basis vectors of the container cell.
	 *  \param trans     true if the basis vectors were changed (rotated).
	 *  \param nx        Number of grid blocks along x.
	 *  \param ny        Number of grid blocks along y.
	 *  \param nz        Number of grid blocks along z.
	 *  \param mem       Initial per-block memory allocation. */
	container_leak(const PeriodicCellList<Empty> & cInit, const std::vector<GeometryVector> & newBasis, bool trans, int nx, int ny, int nz, int mem=6);

	/** \brief Copy constructor that re-inserts all particles of the source.
	 *  \param source  Container to copy. */
	container_leak(const container_leak & source);

	/** \brief same as ordinary put, but now it return a block ID ijk
	 *  \param n  Particle index to insert.
	 *  \param x  Cartesian x coordinate.
	 *  \param y  Cartesian y coordinate.
	 *  \param z  Cartesian z coordinate.
	 *  \return   Block id (ijk) into which the particle was placed. */
	inline int put_ijk(int n, double x, double y, double z) {
		int ijk = 0;
		put_locate_block(ijk, x, y, z);
		id[ijk][co[ijk]] = n;
		double *pp = p[ijk] + 3 * co[ijk]++;
		*(pp++) = x; *(pp++) = y; *pp = z;
		return ijk;
	}
	
	/** \brief Compute the Voronoi cell of stored particle i, undoing any basis rotation.
	 *  \param i          Index of the particle within this container.
	 *  \param[out] v     Vertices of the Voronoi cell of particle i.
	 *  \param pLinkage   Optional output: vertex-linkage information.
	 *  \param pNeighbors Optional output: indices of the neighboring particles. */
	void GetVORONOICELL(int i, std::vector<GeometryVector> &v, std::vector<signed long> * pLinkage = nullptr, std::set<size_t> * pNeighbors = nullptr);
};

/** \brief Derive the container basis and grid resolution from a configuration.
 *  \param c           Configuration to convert.
 *  \param[out] basis  3x3 basis vectors for the container cell.
 *  \param[out] trans  true if the basis vectors had to be rotated.
 *  \param[out] nx     Number of grid blocks along x.
 *  \param[out] ny     Number of grid blocks along y.
 *  \param[out] nz     Number of grid blocks along z.
 *  \param mem         Target number of particles per grid block. */
void GetConversionParam(const PeriodicCellList<Empty> &c, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem=6);

/** \brief Derive container basis and grid resolution for an approximate particle count.
 *  \param c            Configuration to convert.
 *  \param approxNumPrts  Approximate number of particles used to size the grid.
 *  \param[out] basis   3x3 basis vectors for the container cell.
 *  \param[out] trans   true if the basis vectors had to be rotated.
 *  \param[out] nx      Number of grid blocks along x.
 *  \param[out] ny      Number of grid blocks along y.
 *  \param[out] nz      Number of grid blocks along z.
 *  \param mem          Target number of particles per grid block. */
void GetConversionParam(const PeriodicCellList<Empty> &c, int approxNumPrts, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem = 6);

/** \brief Build a voro++ container_leak from a configuration.
 *  \param c  Configuration to load into the container.
 *  \return   Container holding the configuration's particles. */
container_leak Config2Container(const PeriodicCellList<Empty> &c);
#endif


#include <complex>
/** \brief return Psi_6 of particle j in Config.
 *  in 2D, Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(6i \\theta) ).
 *  in 3D, an additional parameter, m, has to be provided. Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(m * i \\phi)*P_6^m(cos(\\theta)) ), where P_l^m is the normalized associated Legendre polynomial.
 *  where theta is the angle between the bond and some reference direction.
 *  \param Config        Configuration containing the particle.
 *  \param j             Index of the particle whose bond-orientational order is computed.
 *  \param NeighborType  definition of neighbors: 0 = six nearest neighbors in 2D; 1 = neighbor in Voronoi diagram.
 *  \param m             Spherical-harmonic order m (3D only).
 *  \return              Complex Psi_6 order parameter of particle j. */
std::complex<double> Psi6(const PeriodicCellList<Empty> & Config, size_t j, int NeighborType = 0, signed int m=0);




#endif