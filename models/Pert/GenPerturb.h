/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file GenPerturb.h
 * \brief ConfigGen-derived transformations that turn one point process into
 *        another: GenPtsDefect (introduces vacancies / interstitials) and
 *        GenPerturb (stochastically displaces particles by a given distribution).
 *
 * Split out of GenerateConfigs.h during the per-model reorganization; the class
 * bodies are the code-reviewed (canonical) module versions. Every definition is
 * header-inline so this header is ODR-safe when pulled into several translation
 * units through the <GenerateConfigs.h> facade.
 */

#ifndef MODELS_GENPERTURB_H
#define MODELS_GENPERTURB_H

#include <functional>

#include <GenerateConfigs.h>

/********* Transformation **********
 * Converts a given point process to another by a specified method
 * Require an initial configuration.
 * Default initial system = d-dimensional cubic lattice
 **********************************/

/** /brief A struct to introduce point defects to a given system.
 *	It stores the mean of point defect concentration = (# of vacancies) / (total # of points ).
 *  Unnecessary information; dimension, number density, the number of particles */
class GenPtsDefect :ConfigGen {
protected:
	double MeanPointDefects;		//< mean fraction of vacancies
	int dist_type;					//< 0 = fixed number of point defects
									//< 1 = fluctuating number of point defects (follow Poisson distribution).
	int DefectType;					//< 0 = vacancy
									//< 1 = interstitials.
	RandomGenerator RNG;			//< a random number generator
	/** \brief A member function to get the total number of point defects.
	 *  \param TotalPrts	The total number of particles in the initial configuration.
	 *  \return The expected number of point defects.
	 */
	inline size_t GetNumDefects(size_t TotalPrts){
		double mean = TotalPrts * MeanPointDefects;
		if (dist_type == 0) {
			return (int)round(mean);
		}
		else if (dist_type == 1) {
			return RNG.PoissonRNG(mean);
		}
		else {
			std::cerr << "Undefined distributions\n";
			return 0;
		}
	}
public:
	/** \brief A default constructor.
	 *  \param mean				The mean of concentration of point defects (<1).
	 *  \param distribution_type	Distribution of the number of point defects (0 => fixed, 1 => Poisson distribution).
	 *  \param defect_type		Type of point defects (0 = vacancy, 1 = interstitials). 	 */
	GenPtsDefect(double mean, int distribution_type, int defect_type) : ConfigGen(), MeanPointDefects(mean), dist_type(distribution_type), DefectType(defect_type) {
		if (mean > 1.0) {
			std::cout << "GenPtsDefect:: mean is too high!" << std::endl;
			return;
		}
		if (distribution_type != 0 && distribution_type != 1) {
			std::cout << "GenPtsDefect:: undefined distributions" << std::endl;
			return;
		}
		if (defect_type != 0 && defect_type != 1) {
			std::cout << "GenPtsDefect:: undefined defect types" << std::endl;
			return;
		}
	};

	/** \brief A copy constructor.
	 *  \param source  The GenPtsDefect to copy. */
	GenPtsDefect(const GenPtsDefect & source) : GenPtsDefect(source.MeanPointDefects, source.dist_type, source.DefectType) {}

	/** \brief A constructor via streams.
	 *  \param ifile  Input stream supplying the parameters.
	 *  \param ofile  Output stream for prompts and echoes. */
	GenPtsDefect(std::istream & ifile, std::ostream & ofile);

	/** \brief A member function to reset.
	 *  \param numberdensity  Ignored (GenPtsDefect doesn't use number density). */
	virtual void SetRho(double numberdensity) { (void)numberdensity; std::cout << "GenPtsDefect doesn't use number density\n"; }
	/** \brief Unused override: GenPtsDefect doesn't use the number of particles.
	 *  \param numPrts  Ignored. */
	virtual void SetNumPrts(int numPrts) { (void)numPrts; std::cout << "GenPtsDefect doesn't use the number of particles\n"; }

	/** \brief A member function to introduce point defects in a Configuraiton object.
	 *  \param[in,out] c	A point configuration in which point defects are introduced.
	 */
	virtual void GenerateC(Configuration &c) {
		if (c.NumParticle() > 0) {
			size_t num_defects = this->GetNumDefects(c.NumParticle());
			if (this->DefectType == 0) {
				//Point vacancies;
				int idx = 0;
				for (size_t i = 0; i < num_defects; i++) {
					idx = (int)floor(c.NumParticle()*RNG.RandomDouble());
					c.DeleteParticle(idx);
				}
			}
			else if (this->DefectType == 1) {
				//Interstitial
				GeometryVector rel_pos(static_cast<int>(c.GetDimension()));
				for (size_t i = 0; i < num_defects; i++) {
					for (int j = 0; j < c.GetDimension(); j++)
						rel_pos.x[j] = RNG.RandomDouble();

					c.Insert("a", rel_pos);
				}
			}
			else {
				std::cerr << "Undefined defect_type" << std::endl;	return;
			}
		}
		else {
			std::cout << "c should have at least 1 particle." << std::endl;
			return;
		}
	}

	/** \brief A memebr function to introduce point defects in a SpherePacking object.
	 *  \param[out] c  Sphere packing (generation not implemented).
	 */
	virtual void GenerateP(SpherePacking & c) {
		(void)c;
		std::cout << "GenPtsDefects::GenerateP is not implemented.\n";
		return;
	}

	/** \brief Reset a seed.
	 *  \param seed  New seed for the random number generator. */
	void SetSeed(int seed) { this->RNG.seed(seed); }
};

/** \brief Construct a GenPtsDefect class via an input stream. */
inline GenPtsDefect::GenPtsDefect(std::istream & ifile, std::ostream & ofile) : ConfigGen() {
	//Initialize every member to a safe value before parsing so invalid input never leaves it indeterminate.
	this->MeanPointDefects = 0.0;
	this->dist_type = 0;
	this->DefectType = 0;
	int temp1;
	double temp2;
	ofile << "----------------------------------\n"
		  << "		      Point defects		    \n"
		  << "----------------------------------\n" << std::endl;
	ofile << "The expected concentration of point defects (# of point defects / total # of particles ) = ";
	ifile >> temp2;
	if (temp2 > 0 && temp2 < 1.0) {
		this->MeanPointDefects = temp2;
	}
	else {
		ofile << "This value should be between 0 and 1" << std::endl;
		return;
	}

	ofile << "Type of distribution of the number of defects (0 = fixed/ 1 = Poisson) = ";
	ifile >> temp1;
	if (temp1 == 0 || temp1 == 1)
		this->dist_type = temp1;
	else {
		std::cout << "Please choose either of 0 and 1\n";
	}

	ofile << "Type of point defects (0 =Vacancies/ 1 = interstitial)";
	ifile >> temp1;
	if (temp1 == 0 || temp1 == 1)
		this->DefectType = temp1;
	else {
		std::cout << "Please choose either of 0 and 1\n";
	}
}


/** \brief A class to stochastically and independently displace particle positions by a given distribution!
	This class doesn't require "the number of particles", "number density"
*/
class GenPerturb : public ConfigGen
{
protected:
	std::vector<double> parameters;			//< parameters for iCMD
	std::function<double (double) > iCMD;	//< the inverse of cummulative function of radial displacements.
											//< iCMD: [0,1]-> [0, infty].
	std::function<GeometryVector(RandomGenerator &)> GetRandomDispl;	//< A Lambda function to get a random displacement.
	RandomGenerator rng;					//< a number number generator
	int distType = 0;						//< 0 = radial via inverse-cumulative iCMD, 1 = uniform, 2 = gaussian
	bool ready = true;						//< false when no valid displacement strategy was configured

	/** \brief (Re)build GetRandomDispl so it captures THIS object's own state.
	 *  Centralizes the per-distribution callbacks; a copied object rebuilds via this
	 *  helper so its callback never captures the source object's `this`. */
	void BuildRandomDispl() {
		switch (this->distType) {
		case 1://uniform inside a ball of radius parameters[0]
			this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector {
				GeometryVector result = RandomUnitVector(this->d, g);
				double r;
				switch (this->d) {
				case 1: r = g.RandomDouble();	break;
				case 2: r = sqrt(g.RandomDouble());	break;
				case 3: r = cbrt(g.RandomDouble()); break;
				default: r = pow(g.RandomDouble(), 1.0 / this->d);
				}
				return this->parameters[0] * r * result;
			};
			break;
		case 2://isotropic Gaussian with standard deviation parameters[0]
			this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector {
				GeometryVector x(static_cast<int>(this->d));
				for (int i = 0; i < x.Dimension; i++)
					x.x[i] = g.RandomDouble_normal(this->parameters[0]);
				return x;
			};
			break;
		default://radial displacement whose magnitude follows the inverse cumulative iCMD
			this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector {
				return g.RandomDouble(this->iCMD) * RandomUnitVector(this->d, g);
			};
			break;
		}
	}
public:
	/** \brief A constructor.
	 *  \param dimension   Spatial dimension.
	 *  \param inverseCMD  Inverse cumulative distribution of the radial displacement.
	 *  \param seed        Seed for the random number generator. */
	GenPerturb(size_t dimension, const std::function <double(double)> & inverseCMD, int seed = 0) : ConfigGen() {
		this->d = dimension;
		this->iCMD = inverseCMD;
		this->rng.seed(seed);
		this->distType = 0;
		this->ready = true;
		this->BuildRandomDispl();
	}

	/** \brief A constructor via an input stream.
	 *  Predefined distributions = uniform, gaussian, powerlaw.
	 *  If Type of distribution = etc => it requries the inverse of cumulative distribution function.
	 *  \param ifile	Input stream supplying the parameters.
	 *  \param ofile	Output stream for prompts and echoes.
	 *  \param inverseCMD	The inverse function of cumulative distribution of radial displacement.  */
	GenPerturb(std::istream & ifile, std::ostream & ofile, const std::function <double(double)> & inverseCMD = nullptr) : ConfigGen() {
		ofile	<< " ----------------------------------\n "
				<< "             GenPerturb            \n"
				<< " ----------------------------------\n" << std::endl;
		ofile << "Dimension = ";
		ifile >> this->d;

		this->ready = false;
		std::string name;
		ofile << "Type of distribution = ";
		ifile >> name;
		std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch){ return (char)std::tolower(ch); });
		if (strcmp(name.c_str(), "uniform")==0) {
			double maxRadius;
			ofile << "Maximal radius (positive real) = ";
			ifile >> maxRadius;
			parameters.push_back(maxRadius);
			this->distType = 1;
			this->ready = true;
		}
		else if (strcmp(name.c_str(), "gaussian") == 0) {
			double standard_deviation;
			ofile << "Input the standard deviation (positive real)";
			ifile >> standard_deviation;
			parameters.push_back(standard_deviation);
			this->distType = 2;
			this->ready = true;
		}
		else if (strcmp(name.c_str(), "powerlaw") == 0|| strcmp(name.c_str(), "power-law") == 0) {
			double delta, alpha;
			ofile << "Input the length scale (positive real)";
			ifile >> delta;
			ofile << "Input the power exponent (positive real)";
			ifile >> alpha;
			this->parameters.push_back(delta);	this->parameters.push_back(alpha);

			//Capture the scalar parameters BY VALUE so a copied GenPerturb (whose copy
			//ctor delegates with source.iCMD) never dereferences the source object's
			//`this`/parameters and thus cannot dangle after the source is destroyed.
			this->iCMD = [delta, alpha](double y)->double {
				if (y  <  0.5 / (1.0 + alpha))
				{
					return  -delta / pow(2.0*(1.0 + alpha)*y, 1.0 / alpha);
				}
				else if (y  <  0.5*(1.0 + alpha / (1.0 + alpha)))
				{
					return  (2.0*y - 1)*(1 + alpha) / alpha * delta;
				}
				else
				{
					return  delta / pow(2.0*(1.0 - y)*(1.0 + alpha), 1.0 / alpha);
				}
			};
			this->distType = 0;
			this->ready = true;
		}
		else if (strcmp(name.c_str(), "etc") == 0) {
			if (inverseCMD != nullptr) {
				this->iCMD = inverseCMD;
				this->distType = 0;
				this->ready = true;
			}
			else {
				ofile << "Distribution is not defined!" << std::endl;
			}
		}
		else{
			ofile << "Distribution is not defined!" << std::endl;
		}

		if (this->ready)
			this->BuildRandomDispl();
		else
			this->GetRandomDispl = [](RandomGenerator &)->GeometryVector {
				std::cerr << "GenPerturb: no displacement distribution configured.\n";
				return GeometryVector();
			};
	}

	/** \brief A copy constructor.
	 *  \param source  The GenPerturb to copy. */
	GenPerturb(const GenPerturb & source) : GenPerturb(source.d, source.iCMD){
		this->parameters = std::vector<double>(source.parameters);
		this->distType = source.distType;
		this->ready = source.ready;
		//Rebuild the callback so it captures THIS object's state, never the source's `this`
		//(copying source.GetRandomDispl would alias the source and dangle if it dies first).
		this->BuildRandomDispl();
	}


	/** \brief Reset random seed.
	 *  \param seed  New seed for the random number generator. */
	void SetSeed(int seed) { this->rng.seed(seed); 	}

	/** \brief Stochastically displace particles in c.
	 *  \param[in,out] c	A point configuration that will be perturbed. It should have at least one particle, and positive volume.	 */
	virtual void GenerateC(Configuration & c) {
		if (!this->ready) {
			std::cerr << "GenPerturb::GenerateC: no displacement distribution configured.\n";
			return;
		}
		if (c.GetDimension() == this->d) {
			if (c.NumParticle() > 0) {
				if (c.PeriodicVolume() > 0) {
					GeometryVector Displacement_Rel(static_cast<int> (this->d));

					for (size_t i = 0; i < c.NumParticle(); i++) {
						//A random and isotropic displacement in Cartesian coordinates.
						Displacement_Rel = this->GetRandomDispl(this->rng);//this->rng.RandomDouble(this->iCMD) * RandomUnitVector(this->d, this->rng);
						//Convert it into relative coordinates.
						Displacement_Rel = c.CartesianCoord2RelativeCoord(Displacement_Rel);
						//Displace particle i by this dispalcement.
						c.MoveParticle(i, c.GetRelativeCoordinates(i) + Displacement_Rel);
					}
				}
				else {
					std::cout << "c should have positive volume" << std::endl;
					return;
				}
			}
			else {
				std::cout << "c should have at least one particle." << std::endl;
				return;
			}
		}
		else {
			std::cout << "The dimension of c should be identical to the dimension of GenPerturb" << std::endl;
			return;
		}
	}


	/** \brief Stochastically displace spherical particles in c.
	 *  If particle i is overlapping with another at a new position, trial movements are made until particle i is not overlapping with others.
	 *  \param[in,out] c	A SpherePacking that will be perutrbed. This member function ends up with a perturbed packing.*/
	virtual void GenerateP(SpherePacking & c) {
		if (!this->ready) {
			std::cerr << "GenPerturb::GenerateP: no displacement distribution configured.\n";
			return;
		}
		if (c.GetDimension() == this->d) {
			if (c.NumParticle() > 0) {
				if (c.PeriodicVolume() > 0) {
					GeometryVector	TrialPos_Rel(static_cast<int> (this->d)),	//Trial position in relative coordinates.
									OldPos_Rel(static_cast<int> (this->d));		//Initial position in relative coordinates.

					const size_t maxTrials = 1000000;	//per-particle trial cap to avoid an unbounded rejection loop
					double totalTrials = 0.0;			//total number of trial displacements attempted
					for (size_t i = 0; i < c.NumParticle(); i++) {
						OldPos_Rel = c.GetRelativeCoordinates(i);	//Initial particle position.
						double r = c.GetCharacteristics(i);			//Particle radius.

						if (c.NumParticle() > 1) {	//There are multiple particles!!
							//Temporarily set aside particle i to either particle 0 or particle i-1.
							int nextPrtIdx = (i == 0) ? 1 : i - 1;
							c.MoveParticle(i, c.GetRelativeCoordinates(nextPrtIdx));
							c.UpdateCharacteristics(0.0, i);
							//Try to move particle i until it is not overlapping with others (bounded).
							size_t trials = 0;
							bool placed = false;
							do {
								totalTrials++;	trials++;
								//A random and isotropic displacement in Cartesian coordinates.
								TrialPos_Rel = this->GetRandomDispl(this->rng);
								//Convert it into relative coordinates, and add initial position.
								TrialPos_Rel = c.CartesianCoord2RelativeCoord(TrialPos_Rel) + OldPos_Rel;
								if (!c.CheckOverlap(TrialPos_Rel, r)) { placed = true; break; }
							} while (trials < maxTrials);
							if (!placed) {
								std::cout << "GenPerturb::GenerateP: exceeded " << maxTrials
										  << " trials for particle " << i << "; keeping its original position.\n";
								TrialPos_Rel = OldPos_Rel;
							}
						}
						else {	//There is only one particle.
							totalTrials++;
							//A random and isotropic displacement in Cartesian coordinates.
							TrialPos_Rel = this->GetRandomDispl(this->rng);
							//Convert it into relative coordinates, and add initial position.
							TrialPos_Rel = c.CartesianCoord2RelativeCoord(TrialPos_Rel) + OldPos_Rel;
						}
						//Move particle i into a new position and recover its radius.
						c.MoveParticle(i, TrialPos_Rel);
						c.UpdateCharacteristics(r, i);
					}
					if (totalTrials > 0.0)
						std::cout << "Success rate of displacing particles = " << c.NumParticle() / totalTrials << std::endl;
				}
				else {
					std::cout << "c should have positive volume" << std::endl;
					return;
				}
			}
			else {
				std::cout << "c should have at least one particle." << std::endl;
				return;
			}
		}
		else {
			std::cout << "The dimension of c should be identical to the dimension of GenPerturb" << std::endl;
			return;
		}
	}


	/** \brief Undefined member functions !
	 *  \param rho  Ignored. */
	virtual void SetRho(double rho) { (void)rho; std::cout << "GenPerturb::SetRho is not defined\n"; };
	/** \brief Unused override: GenPerturb does not use the number of particles.
	 *  \param num  Ignored. */
	virtual void SetNumPrts(int num) { (void)num; std::cout << "GenPerturb::SetNumPrts is not defined\n"; };
	/** \brief Unused override: GenPerturb does not use the packing fraction.
	 *  \param phi  Ignored. */
	virtual void SetPackingFraction(double phi) { (void)phi; std::cout << "GenPerturb::SetPackingFraction is not defined\n"; };
};

#endif // MODELS_GENPERTURB_H
