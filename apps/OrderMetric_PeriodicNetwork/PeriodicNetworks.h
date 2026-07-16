/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file PeriodicNetworks.h
 * \brief Header file for computing the \overline{B}_V for periodic networks.
 * It takes advantage of the Fourier-space representation of local volume-fraction variance, and the spectral density of periodic networks.
 * At this moment, the dispersion class is incomplete to compute the spectral density
 */



#ifndef PERIODIC_NETWORKS_INCLUDED
#define PERIODIC_NETWORKS_INCLUDED

#if defined USE_SHAPES
#include <ParticleShapes.h>
#endif
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <GenerateConfigs.h>
#include <StructureFactor.h>
#include <Dispersions.h>

#include <vector>

/** \brief Shared throwaway target for optional output-vector parameters whose
 * default the caller does not care about. Replaces a `= std::vector<...>()`
 * default (binding a non-const reference to a temporary), which is a GNU
 * extension rejected as an error under -std=c++17. Runtime behavior is
 * unchanged: the default writes into a single discarded static. */
inline std::vector<GeometryVector> & _bv_discard_vec() {
	static std::vector<GeometryVector> v;
	return v;
}

#if defined USE_SHAPES

/** A function to generate certain periodic networks at packing fraction = 1. */
void GetPeriodicNetworks(std::istream & ifile, std::ostream & ofile, Dispersion & c, std::string & filename);

void GetPeriodicNetworks(const std::string & name, const std::string & scale, double val, Dispersion & c, std::string & filename);

/** A function to generate certain periodic spehre packings at the maximal packing fraction.
 * @param[in]	ifile	Input stream
 * @param[out]	ofile	Output stream
 * @param[out]	c	a periodic sphere packing in a fundamental cell with the maximal packing fraction. 
 * @param[out]	filename	a string for the prefix of file name. */
void GetPeriodicPackings(std::istream & ifile, std::ostream & ofile, SpherePacking & c, std::string & filename);

void GetPeriodicPackings(const std::string & NetworkName, const std::string & Lengthscale, double val, SpherePacking & sphC, std::string & filename);


/** \struct	A data structure for save */
struct PeriodicNetworkParameters {
	int NumThreads = 4;		/**< the number of threads in openMP */
	/** parameters for computations in the Fourier space */
	double kMax,	/**< the largest wavenumber up to which the summation is computed in the Fourier space. */
		kMin,		/**< the smallest wavenumber up to which B_V is computed. */
		dK;			/**< the change in wavenumber over which B_V is computed. */
	std::vector<GeometryVector> Ks;	/**< A list of all wavevectors whose magnitudes are smaller than kMax. */
	bool Init_Parameter_Fourier = false;	/**< true = parameters are initialized. */
	/** parameters for computations in the Direct space */
	double Rmax,	/**< the largest radius of observation windows. */
		dR,			/**< the difference between window radii. */
		Rmin,		/**< the smallest window radius up to which the running average is computed. */
		intv_R;		/**< the interval of window radius up to which the running average is computed. */
	int num_sampling_windows = 10000;	/**< the number of observation windows to compute the local volume-fraction variance. */
	bool Init_Parameter_Direct = false;	/**< true = parameters are initialized. */
	/** Initialize parameters required for computing B_V in the Fourier space. */
	void SetParam_Fourier(double kMax, double kMin, double dK) {
		this->kMax = kMax; this->kMin = kMin; this->dK = dK;
		Init_Parameter_Fourier = true;
	}
	/** Initialize parameters required for computing B_V in the direct space. */
	void SetParam_Direct(double R_max, double dR, double R_min, double intv, int num_windows = 10000) {
		this->Rmax = R_max; this->dR = dR; this->Rmin = R_min; this->intv_R = intv; this->num_sampling_windows = num_windows;
		Init_Parameter_Direct = true;
	}

	/** Compute the surface-coefficient term B_v for nonspherical particles from the summation in the Fourier space. 
	 * @param[in] network	A periodic fundamental cell of interest.
	 * @param[in] phi_void	Volume fraction of the void phase (1-packing fration)
	 * @param[out] result	result[i] = {wavenumber, B_V summed up to wavenumber}. */
	void ComputeBV_Fourier(const Dispersion & network, double phi_void, std::vector<GeometryVector> & result);
	
	void ComputeBV_Direct(const Dispersion & network, double phi_void, std::vector<GeometryVector> & result, std::vector<GeometryVector> & LV = _bv_discard_vec());

	/** Compute the surface-coefficient term B_v for spherical particles from the summation in the Fourier space.
 	 * @param[in] network	A periodic fundamental cell of interest.
	 * @param[in] packingfraction	Packing fraction
	 * @param[out] result	result[i] = {wavenumber, B_V summed up to wavenumber}. */
	void ComputeBV_Fourier(const SpherePacking & network, double packingfraction, std::vector<GeometryVector> & result);

	/** Compute the surface-coefficient term B_v for spherical particles from the local volume fraction variance.
	 * @param[in] network	A periodic fundamental cell of interest.
	 * @param[in] packingfraction	Packing fraction
	 * @param[out] result	result[i] = {R, B_V from the running average of the local volume-fraction variance up to R}. */
	void ComputeBV_Direct(const SpherePacking & network, double packingfraction, std::vector<GeometryVector> & result, std::vector<GeometryVector> & LV = _bv_discard_vec());

	/** Compute the local volume-fraction variance for spherical windows via the Fourier representations.
	 * @param[in]	network	A periodic fundamental cell of interest.
	 * @param[in]	phi_void	Volume fraction of the void phase (1-packing fration)
	 * @param[out]	result	result[i] = {R, sigma V}. */
	void ComputeVV_Fourier(const Dispersion & network, double phi_void,
		std::vector<GeometryVector> & result);

	/** Compute the local volume-fraction variance for spherical windows via the Fourier representations.
	 * @param[in]	packing	A periodic fundamental cell of interest.
	 * @param[in]	phi_void	Volume fraction of the void phase (1-packing fration)
	 * @param[out]	result	result[i] = {R, sigma V}. */
	void ComputeVV_Fourier(const SpherePacking & network, double phi_void,
		std::vector<GeometryVector> & result);

	/** Compute the local volume-fraction variance for spherical windows via the Fourier representations.
	 * @param[in]	packing	A periodic fundamental cell of interest.
	 * @param[in]	phi_void	Volume fraction of the void phase (1-packing fration)
	 * @param[out]	result	result[i] = {R, sigma V}. */
	void ComputeVV_Direct(const SpherePacking & network, double phi_void,
		std::vector<GeometryVector> & result);

	/** Compute the local volume-fraction variance for spherical windows via the Fourier representations.
	 * @param[in]	spectral_density	A spectral density...
	 * @param[in]	packing fraction
	 * @param[out]	result	result[i] = {R, sigma V}. */
	void ComputeVV_Fourier(const std::vector<GeometryVector> & spectral_density, double packingfraction,
		std::vector<GeometryVector> & result);
};










#endif


 
 
 

void PeriodicNetworkCLI(std::istream & ifile, std::ostream & ofile);




 /** \brief Header file to compute B_V coefficients for periodic two-phase systems (Sphere packings, networks). */
struct PeriodicNetwork{
	Configuration c;	/**< A fundamental cell */
	std::vector<Shape *> Particles;	/**< Particle shapes */
	double phi;			/**< volume fraction  */
	int num_threads;

	/** A default constructor */
	PeriodicNetwork() { c = Configuration(); this->phi = -1; this->num_threads = 0; };

	/** A constructor */
	PeriodicNetwork(const Configuration & FundamentalCell, std::vector<Shape *> & Prts) {
		this->c = Configuration(FundamentalCell);
		for (size_t i = 0; i < Prts.size(); i++)
			this->Particles.push_back(Prts[i]);
		this->phi = -1;
		this->CheckFundamentalCell();
		this->num_threads = 0;
	}

	/** Destructor */
	~PeriodicNetwork() {};	

	/** Check the consistentcy of a fundamental cell.
	 *	It checks the number of particles and centroids,
	 *  and volume fraction.
	 *  If the volume fraction is uninitialized, this function will initialize the volume fraction.
	 * @return True = all quantities are consistent.
	 */
	bool CheckFundamentalCell();

	/**	Change volume fraction of this network.
	 * @param volumefraction	Volume fraction of this periodic newtwork.
	 */
	void SetPhi(double volumefraction);

	/** Chnage set the number of threads 
	 * @param num	The number of threads
	 */
	void SetNumThreads(int num) { this->num_threads = num; };

	/**	Print fundamental cell, and all particles.
	 *	@param ofile	Output stream.
	 */
	void Print(std::ostream & ofile) const;

	double GetPhi() const { return phi; };			/** @return volume fraction. */
	int GetDimension() const { return c.GetDimension(); };	/** @return space dimension. */
	
	void ComputeBV_Fourier(double kmin, double kmax, double dk, std::vector<GeometryVector> & result, std::vector<GeometryVector> & ks = _bv_discard_vec());
	void ComputeBV_Direct(double Rmax, double dR);

private :
	/** Compute a squared form factor of particles in a fundamental cell
	 * @param k	A wavevector
	 * @return |\tilde{m}(k)|^2, where \tilde{m}(k) is the Fourier transform of particles in a fundamental cell 
	 */
	inline double GetSquaredFormFactor(const GeometryVector & k) const {
		std::complex<double> mtilde(0.0, 0.0);
		for (size_t i = 0; i < this->c.NumParticle(); i++) {
			mtilde += this->Particles[i]->GetFormFactor(k)*std::exp(std::complex<double>(0.0,k.Dot(this->c.GetCartesianCoordinates(i))) );
			/* The factor e^{i k*x} should be used because the form factor is computed from \int_particle e^{i k.x}\d{x}. */
		}
		
		//std::cout << k.x[0]<<"\t"<<k.x[1] << ":\t" << norm(mtilde) << "\n";
		return norm(mtilde);
	} 
};


void ComputeBV_CLI();
#endif