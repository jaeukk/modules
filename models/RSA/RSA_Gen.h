/** 
 *  Author	:	Jaeuk Kim
 *  Email	: phy000.kim@gmail.com
 *  Date	: Febrauray 2022
 */

#ifndef RSA_GEN_H
#define RSA_GEN_H

/* in $(core) */
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "RandomSequentialAddition.h"
#include <GenerateConfigs.h>

/** \brief Packing fractions of saturated RSA packings.
 *  Delegates to the single canonical table SaturatedPackingFraction() in
 *  RandomSequentialAddition.h (Reference: G. Zhang and S. Torquato,
 *  Phys. Rev. E 88, 053312 (2013)).
 *  \param	d	Space dimension.
 *  \return	the saturated packing fraction.
 */
inline double PhiC(size_t d) {
	return SaturatedPackingFraction((DimensionType)d);
}


class RSA_Gen : public ConfigGen {
private:
	double targetPhi = 0.0;	//target packing fraction
	bool saturated = true;	//indicate whether we generate a saturated RSA packing
	size_t numMax = 0;
public:

	/** \brief Default constructor: build an RSA generator with default parameters. */
	RSA_Gen() : ConfigGen::ConfigGen() {};

	/** \brief Construct an RSA generator for a target packing fraction.
	 *  set numberdensity = -1 if you don't want to rescale the system.
	 *  \param dimension        Spatial dimension.
	 *  \param numPrts          Target number of particles.
	 *  \param packingFraction  Target packing fraction.
	 *  \param numberdensity    Number density to rescale to (-1 to skip rescaling). */
	RSA_Gen(int dimension, int numPrts, double packingFraction, double numberdensity);

	/** \brief Construct an RSA generator that produces a saturated packing.
	 *  \param dimension     Spatial dimension.
	 *  \param numPrts       Target number of particles.
	 *  \param numberdensity Number density to rescale to (-1 to skip rescaling). */
	RSA_Gen(int dimension, int numPrts, double numberdensity);

	/** \brief Construct an RSA generator interactively (CLI form).
	 *  \param ifile  Input stream supplying the parameters.
	 *  \param ofile  Output stream for prompts and echoes. */
	RSA_Gen(std::istream & ifile, std::ostream & ofile);

	/** \brief Set the target particle count and recompute the working count from the saturation state.
	 *  Non-positive values are rejected without changing the generator.
	 *  \param n New target number of particles. */
	virtual void SetNumPrts(int n) {
		if (n <= 0) {
			std::cerr << "RSA_Gen::SetNumPrts: particle number must be positive.\n";
			return;
		}
		this->numPrts = (size_t)n;
		//Recompute the working count from the target/saturation state (as the constructors do).
		if (this->saturated)
			this->numMax = this->numPrts;
		else
			this->numMax = (size_t)round(PhiC(this->d) / this->targetPhi * this->numPrts);
	};
	/** \brief Set the number density used to rescale the generated packing.
	 *  \param numberdensity  New number density. */
	virtual void SetRho(double numberdensity) {
		this->rho = numberdensity;
	}

	/** \brief Set the target packing fraction and recompute the saturation state and working count.
	 *  Fractions at least 99 percent of the tabulated saturation fraction are treated as saturated.
	 *  \param phi New target packing fraction; must be positive. */
	virtual void SetPackingFraction(double phi) {
		this->targetPhi = phi;
		this->saturated = (this->targetPhi >= 0.99*PhiC(this->d)) ? true : false;
		if (this->saturated)
			this->numMax = numPrts;
		else
			this->numMax = (int)round(PhiC(this->d) / phi * this->numPrts);
	}

	/** \brief generate a point configuration; system will be rescaled to a given number density.
	 *  \param[out] c  Generated point configuration. */
	virtual void GenerateC(Configuration &c);

	/** \brief generate a sphere packing; system will be rescaled to a given number density.
	 *  \param[out] c  Generated sphere packing. */
	virtual void GenerateP(SpherePacking &c);
};


/****************************
 *	RSA_Gen: constructors
 ****************************/
//set numberdensity = -1 if you don't want to rescale the system.
inline RSA_Gen::RSA_Gen(int dimension, int numPrts, double packingFraction, double numberdensity) {
	this->d = dimension;	this->numPrts = numPrts;	this->rho = numberdensity;	this->targetPhi = packingFraction;
	this->numThreads = 1;
	this->saturated = (this->targetPhi >= 0.99*PhiC(this->d)) ? true : false;
	if (this->saturated)
		this->numMax = numPrts;
	else
		this->numMax = (int)round(PhiC(dimension) / packingFraction * numPrts);
}
//Saturated RSA
inline RSA_Gen::RSA_Gen(int dimension, int numPrts, double numberdensity) {
	this->d = dimension;	this->numPrts = numPrts;	this->rho = numberdensity;
	this->saturated = true;	this->numThreads = 1;
	this->numMax = numPrts;
}
//CLI form
inline RSA_Gen::RSA_Gen(std::istream & ifile, std::ostream & ofile) {
	ofile << "Dimensions = ";		ifile >> this->d;			ofile << this->d;
	ofile << "PackingFraction = ";	ifile >> this->targetPhi;	ofile << this->targetPhi;
	ofile << "Number density = ";	ifile >> this->rho;			ofile << this->rho;
	ofile << "Particle number = ";	ifile >> this->numPrts;		ofile << this->numPrts;

	this->numThreads = 1;
	this->saturated = (this->targetPhi >= 0.99*PhiC(this->d)) ? true : false;
	if (this->saturated)
		this->numMax = numPrts;
	else
		this->numMax = (int)round(PhiC(this->d) / this->targetPhi * numPrts);
}


 /****************************
 *	RSA_Gen: member functions
 ****************************/

inline void RSA_Gen::GenerateC(Configuration &c) {
	SpherePacking temp;
	this->GenerateP(temp);
	c = Configuration(temp, "a");
}

inline void RSA_Gen::GenerateP(SpherePacking &c) {
//#ifdef RandomSequentialAddition_Included
	bool parallel = (this->numThreads > 1) ? true : false;
	omp_set_num_threads((int)this->numThreads);
	int numCut = (this->saturated) ? (int)numPrts*1.2 : numPrts;
	
	c = GenerateRSAPacking(d, this->numMax, parallel ? (-3) : 3, 1000, 6000, 0.49, GetPreciseClock(), Verbosity > 5 ? std::cout : logfile, nullptr, false, numCut);
		
	//if (this->saturated) {
	//	c = GenerateRSAPacking(d, numPrts, parallel, (int)numPrts*1.2);
	//}
	//else {
	//	c = GenerateRSAPacking(d, this->numMax, parallel, numPrts);
	//}
	if(Verbosity >4)
		std::cout << "A RSA packing is generated: particle number = "<< c.NumParticle() <<"\n";

//#else
//	std::cout << "Please include \"RandomSequentialAdditions.h\"\n";
//#endif  
	//Rescale systems.
	//NOTE (semantics): the box is rescaled so that the REQUESTED count `numPrts`, not the
	//actual returned count c.NumParticle(), yields number density `rho`. For saturated
	//packings the generator may return a different count (up to ~1.2*numPrts is produced,
	//then cut to numCut), so `numPrts` is treated as the target count for this rescaling and
	//the realized number density can differ from `rho` when the counts differ.
	if (this->rho > 0.0) {
		double factor = pow(this->numPrts / this->rho, 1.0/this->d);
		c.Rescale(factor);
	}
}



#endif