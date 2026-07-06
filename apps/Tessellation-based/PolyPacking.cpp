#include "PolyPacking.h"
#include <RandomSequentialAddition.h>

//==============================================================================
// MakeMono2Poly / MakeMono2Poly_parallel sliced verbatim from the legacy
// etc1.cpp (lines 415-562).  They use GetVoronoiCell / PolygonVolume /
// VoronoiCellRoutine (Voronoi.h), all CGAL/voro++-only, so they stay guarded
// behind USE_CGAL exactly as in the legacy file.
//==============================================================================
#ifdef USE_CGAL
#include "Voronoi.h"

double MakeMono2Poly(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, std::vector<int> * numVertices) {
	//make resultConfig empty
	resultConfig = SpherePacking(config);
	resultConfig.RemoveParticles();

	if (numVertices != nullptr) {
		(*numVertices).resize(config.NumParticle(), 0);
	}

	DimensionType d = config.GetDimension();
	double v_1 = HyperSphere_Volume(d, config.GetMaxRadius());
	Volumes.erase(Volumes.begin(), Volumes.end());
	Volumes.resize(config.NumParticle());

	double max_cellVol = 0;
	double temp = 0;
	double r = 0;

	double factor = pow(targetPhi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	//------------- measure volumes of Voronoi cells
	std::vector<GeometryVector> vertices;
	//std::vector<std::vector<size_t>> pN;

	//GetVoronoiCell(config, vertices, pN);
	//for (int i = 0; i < config.NumParticle(); i++) {
	//	
	//}
	for (int i = 0; i < config.NumParticle(); i++) {
		GetVoronoiCell(config, i, vertices);
		if(numVertices != nullptr)
			numVertices->at(i) = (int)vertices.size();
		temp = PolygonVolume(vertices, d);
		max_cellVol = (max_cellVol > temp) ? max_cellVol : temp;

		if (d == 2)
			r = sqrt(temp);
		else if (d == 3)
			r = cbrt(temp);
		else { r = pow(temp, 1.0 / d); }

		r *= factor;
		Volumes.at(i) = temp;
		resultConfig.Insert(r, config.GetRelativeCoordinates(i));
	}

	//for (int i = 0; i < Volumes.size(); i++)
	//	Volumes.at(i) /= (double)config.NumParticle();

	return v_1 / max_cellVol;
}
//return the maximal packing fractions!
double MakeMono2Poly_parallel(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, int numThreads, std::vector<int> * numVertices) {

	resultConfig = SpherePacking(config);
	//resultConfig.RemoveParticles();

	if (numVertices != nullptr) {
		(*numVertices).resize(config.NumParticle(), 0);
	}

	DimensionType d = config.GetDimension();
	int numPrts = config.NumParticle();
	Volumes = std::vector<double>(numPrts, 0.0);
	double v_1 = HyperSphere_Volume(d, config.GetMaxRadius());
	double r = 0;
	double factor = pow(targetPhi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	
	std::vector<GeometryVector> centr;
	//------------- measure volumes of Voronoi cells

#if defined (USE_VORO)
	std::vector<double> cellvol;
	if(d>1)
		VoronoiCellRoutine(config, numThreads, &cellvol, nullptr, nullptr, nullptr);
	else if(d==1)
		VoronoiCellRoutine(config, numThreads, &cellvol, nullptr, nullptr, &centr);
	
	Volumes.swap(cellvol);

#else
#pragma omp parallel num_threads(numThreads)
	{
		SpherePacking c_copy(config);
		std::vector<GeometryVector> vertices;
		double temp = 0;
#pragma omp for
		for (int i = 0; i < numPrts; i++) {
			GetVoronoiCell(c_copy, i, vertices);
			if (numVertices != nullptr)
				numVertices->at(i) = (int)vertices.size();
			temp = PolygonVolume(vertices, d);
			//#pragma omp critical
			//			{
			//				max_cellVol = (max_cellVol > temp) ? max_cellVol : temp;
			//			}
			Volumes.at(i) = temp;
		}
	}
#endif
	double phi1, phi2;
	double temp, temp1, temp2;
	if (d == 1) {
		resultConfig.RemoveParticles();
		double dx = 0.0;
		GeometryVector centr_prt(static_cast<int> (1));
		RandomGenerator rng(std::time(NULL));
		for (size_t i = 0; i < numPrts; i++) {
			temp = Volumes[i];
			r = 0.5*targetPhi*temp;
			dx = (1.0 - targetPhi)*temp;
			centr_prt.x[0] = centr[i].x[0] + dx * (rng.RandomDouble() - 0.5);
			
			resultConfig.Insert(r, resultConfig.CartesianCoord2RelativeCoord(centr_prt));
		}
		phi1 = 1.0;
		phi2 = 1.0;
	}
	else if (d > 1) {
		for (size_t i = 0; i < numPrts; i++) {
			temp = Volumes[i];
			if (d == 2)
				r = sqrt(temp);
			else if (d == 3)
				r = cbrt(temp);
			else { r = pow(temp, 1.0 / d); }
			r *= factor;
			resultConfig.UpdateCharacteristics(r, i);
			//find the max. packing fraction
			//double r1 = config.GetCharacteristics(i);
			double r2 = 0.5*resultConfig.NearestParticleDistance(i);
			temp1 = HyperSphere_Volume(d, config.GetCharacteristics(i)) / temp;
			temp2 = HyperSphere_Volume(d, r2) / temp;
			if (i == 0) {
				phi1 = temp1;
				phi2 = temp2;
			}
			else {
				phi1 = std::min(phi1, temp1);
				phi2 = std::min(phi2, temp2);
			}
		}
	}

	std::cout << "Max Vol Fraction (1) = " << phi1 << ",\t (2) = " << phi2 << "\n";
	return phi2;


}
#endif // USE_CGAL

//-----------------------------------------------------------------------------
// RSA-substrate overload (etc1.cpp lines 627-636): always compiled, but it
// delegates to the CGAL-only MakeMono2Poly above.
//-----------------------------------------------------------------------------
double MakeMono2Poly(DimensionType d, size_t numParticles, double boxSize, double phi0, SpherePacking &resultConfig, double targetPhi, std::vector<double> &Volumes, int seed) {
#ifdef USE_CGAL
	RSA setting(d, numParticles);
	setting.gen.seed(std::time(NULL)+seed);
	double r = pow(phi0/((double) numParticles *HyperSphere_Volume(d, 1.0)), 1.0 / d);
	setting.radius = r;
	setting.RSA_I(numParticles);
	SpherePacking substrate(*setting.pSpheres);
	substrate.Rescale(boxSize);
	return MakeMono2Poly(substrate, resultConfig, targetPhi, Volumes);
#else
	(void)d; (void)numParticles; (void)boxSize; (void)phi0;
	(void)resultConfig; (void)targetPhi; (void)Volumes; (void)seed;
	std::cerr << "MakeMono2Poly (RSA-substrate overload) requires USE_CGAL.\n";
	return 0.0;
#endif
}
