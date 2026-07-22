/**
 *  Author	:	Jaeuk Kim
 *  Email	: phy000.kim@gmail.com
 *  Date	: Febrauray 2022
 */

/** \file HSF_Gen.h
 * \brief ConfigGen wrapper around the equilibrium hard-sphere-fluid sampler HSF.
 *
 * Ported from the legacy c++/Models/HSF/HSF_Gen.h and adapted to the code-reviewed
 * module HardSphereFluids API (unique_ptr pSpheres, per-thread deterministic
 * seeding, serial Evolve fallback). The constructor argument order was unified with
 * RSA_Gen -- (dimension, numPrts, packingFraction, numberdensity) -- rather than the
 * legacy (dimension, N, numberDensity, packingFraction), so RSA_Gen and HSF_Gen share
 * one calling convention.
 */

#ifndef HSF_Gen_H
#define HSF_Gen_H

/* in $(core) */
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "HardSphereFluids.h"
#include <GenerateConfigs.h>

/** \brief ConfigGen model that draws configurations from an equilibrium hard-sphere fluid. */
class HSF_Gen : public ConfigGen {
private:
	int samplingInterval = 1000;	//< number of MC cycles between saved configurations
public:
	HSF setting;					//< the underlying hard-sphere-fluid Monte-Carlo sampler

	/** \brief Default constructor: build an HSF generator with default parameters. */
	HSF_Gen() : ConfigGen::ConfigGen() {};

	/** \brief Construct an HSF generator for the given target.
	 *  \param dimension        Spatial dimension.
	 *  \param numPrts          Number of particles.
	 *  \param packingFraction  Target packing fraction of the fluid.
	 *  \param numberdensity    Number density to rescale the configuration to. */
	HSF_Gen(int dimension, int numPrts, double packingFraction, double numberdensity);

	/** \brief Construct an HSF generator interactively (CLI form).
	 *  \param ifile  Input stream supplying the parameters.
	 *  \param ofile  Output stream for prompts and echoes. */
	HSF_Gen(std::istream & ifile, std::ostream & ofile);

	/** \brief Unused override: HSF fixes the particle count at construction.
	 *  \param num  Ignored. */
	virtual void SetNumPrts(int num) {
		(void)num;
		std::cout << "SetNumPrts(int) is not supported in HSF\n";
	}
	/** \brief Rescale the configuration to a new number density.
	 *  \param rho  New number density. */
	virtual void SetRho(double rho) {
		this->setting.RescaleConfig(rho);
	}
	/** \brief Rescale particle radii to reach a new packing fraction.
	 *  \param phi  New packing fraction. */
	virtual void SetPackingFraction(double phi) {
		double factor = pow(phi / this->setting.PackingFraction, 1.0 / this->d);
		SpherePacking temp(*this->setting.pSpheres);
		temp.Rescale_prtOnly(factor);
		this->setting.ReplaceConfig(temp);
		this->setting.Prepare();
	}

	/** \brief Set the sampling interval between saved configurations.
	 *  \param sampleinterval  Number of MC cycles between samples. */
	void SetInterval(int sampleinterval) { this->samplingInterval = sampleinterval; }

	/** \brief Generate a point configuration by sampling the fluid.
	 *  \param[out] c  Generated point configuration. */
	virtual void GenerateC(Configuration &c);
	/** \brief Generate a sphere packing by sampling the fluid.
	 *  \param[out] c  Generated sphere packing. */
	virtual void GenerateP(SpherePacking &c);
};

inline HSF_Gen::HSF_Gen(int dimension, int numPrts, double packingFraction, double numberdensity) : ConfigGen(dimension, numPrts, numberdensity) {
	setting = HSF(dimension, numPrts, packingFraction);	setting.SetNumThreads(1);
	setting.RescaleConfig(numberdensity);
	//setting.Prepare();
	std::cout << "Default sampling interval = " << this->samplingInterval << " MCcycles\n";
}

inline HSF_Gen::HSF_Gen(std::istream & ifile, std::ostream & ofile) : ConfigGen(ifile, ofile) {
	double phi;
	ofile << "packing fraction = ";	Echo(ifile, ofile, phi);
	setting = HSF(this->d, (int)this->numPrts, phi);
	setting.RescaleConfig(this->rho);
	int seed, num_threads;
	ofile << "random seed = "; Echo(ifile, ofile, seed);
	setting.SetRandomSeed(seed);
	ofile << "number of threads = "; Echo(ifile, ofile, num_threads);
	setting.SetNumThreads(num_threads);
	ofile << "sampling interval (MC  cycles) = "; Echo(ifile, ofile, this->samplingInterval);
}

inline void HSF_Gen::GenerateP(SpherePacking &c) {
	if (!this->setting.IsEquilibrated()) {
		this->setting.Prepare();
		this->setting.Evolve(200);
	}
	this->setting.Evolve(this->samplingInterval);
	c = SpherePacking(*this->setting.pSpheres);
}

inline void HSF_Gen::GenerateC(Configuration &c) {
	SpherePacking temp;
	this->GenerateP(temp);
	c = Configuration(temp, "a");
}

#endif // HSF_Gen_H
