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

/** Packing fractions of saturated RSA packings.
 * @param	d	Space dimension.
 * @return	the saturated packing fraction. 
 */
inline double PhiC(size_t d) {
	/* Reference: 
	G. Zhang and S. Torquato, 
	Precise Algorithm to Generate Random Sequential Addition of Hard Hyperspheres at Saturation, Physical Review E, 88, 053312 (2013). */
	double volumeFraction = 0.0;
	switch (d) {
	case 1:
		volumeFraction = 0.7476;
		break;
	case 2:
		volumeFraction = 0.54707;
		break;
	case 3:
		volumeFraction = 0.38413;
		break;
	default:
		volumeFraction = 0.2;
		std::cout << "RandomSequentialAddition::Rc : out-of-range";
		break;
	}
	return volumeFraction;
}


class RSA_Gen : public ConfigGen {
private:
	double targetPhi;	//target packing fraction
	bool saturated;		//indicate whether we generate a saturated RSA packing
	size_t numMax;
public:

	RSA_Gen() : ConfigGen::ConfigGen() {};

	//set numberdensity = -1 if you don't want to rescale the system.
	RSA_Gen(int dimension, int numPrts, double packingFraction, double numberdensity);

	//Saturated RSA
	RSA_Gen(int dimension, int numPrts, double numberdensity);

	//CLI form
	RSA_Gen(std::istream & ifile, std::ostream & ofile);

	virtual void SetNumPrts(int n) {
		this->numMax = (int)round(n / this->numPrts);
		this->numPrts = n;
	};
	virtual void SetRho(double numberdensity) {
		this->rho = numberdensity;
	}
	virtual void SetPackingFraction(double phi) {
		this->targetPhi = phi;
		this->saturated = (this->targetPhi >= 0.99*PhiC(this->d)) ? true : false;
		if (this->saturated)
			this->numMax = numPrts;
		else
			this->numMax = (int)round(PhiC(this->d) / phi * this->numPrts);
	}
	//generate a point configuration/sphere packing
	//system will be rescaled to a given number density
	virtual void GenerateC(Configuration &c);
	virtual void GenerateP(SpherePacking &c);
};


/****************************
 *	RSA_Gen: constructors
 ****************************/
//set numberdensity = -1 if you don't want to rescale the system.
RSA_Gen::RSA_Gen(int dimension, int numPrts, double packingFraction, double numberdensity) {
	this->d = dimension;	this->numPrts = numPrts;	this->rho = numberdensity;	this->targetPhi = packingFraction;
	this->numThreads = 1;
	this->saturated = (this->targetPhi >= 0.99*PhiC(this->d)) ? true : false;
	if (this->saturated)
		this->numMax = numPrts;
	else
		this->numMax = (int)round(PhiC(dimension) / packingFraction * numPrts);
}
//Saturated RSA
RSA_Gen::RSA_Gen(int dimension, int numPrts, double numberdensity) {
	this->d = dimension;	this->numPrts = numPrts;	this->rho = numberdensity;
	this->saturated = true;	this->numThreads = 1;
	this->numMax = numPrts;
}
//CLI form
RSA_Gen::RSA_Gen(std::istream & ifile, std::ostream & ofile) {
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

void RSA_Gen::GenerateC(Configuration &c) {
	SpherePacking temp;
	this->GenerateP(temp);
	c = Configuration(temp, "a");
}

void RSA_Gen::GenerateP(SpherePacking &c) {
//#ifdef RandomSequentialAddition_Included
	bool parallel = (this->numThreads > 1) ? true : false;
	omp_set_num_threads(this->numThreads);
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
	//Rescale systems
	if (this->rho > 0.0) {
		double factor = pow(this->numPrts / this->rho, 1.0/this->d);
		c.Rescale(factor);
	}
}



#endif