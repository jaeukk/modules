#include "HardSphereFluids.h"
#include <GenerateConfigs.h>
void HSF::LatticePacking(int numPrts, double PackingFraction, std::vector<GeometryVector> & basis, SpherePacking & config) {
	int dimension = basis.size();
	int L = (int)ceil(pow(numPrts, 1.0 / dimension));
	double inverseL = 1.0 / (double)L;
	std::vector<GeometryVector> basis_cell;
	for (int i = 0; i < dimension; i++) {
		basis_cell.emplace_back((double)L*basis[i]);
	}
	
	config = SpherePacking(dimension, &basis_cell[0], 1.0);

	//Check the realizability!
	//Determine the particle radius = half the nearest lattice spacing (avoid overlaps)
	double r = 0.5*sqrt(basis[0].Modulus2());
	for (int i = 1; i < dimension; i++)
		r = std::min(r, 0.5*sqrt(basis[i].Modulus2()));

	double MaxPackingFraction = HyperSphere_Volume(dimension, r)*numPrts / config.PeriodicVolume();
	if (MaxPackingFraction >= PackingFraction) {
		GeometryVector x;
		//Generate a perfect lattice
		for (int i = 0; i < L; i++) {
			if (dimension > 1) {
				for (int j = 0; j < L; j++) {
					if (dimension > 2) {
						for (int k = 0; k < L; k++) {
							if (dimension > 3) {}
							else{
								config.Insert(r, inverseL*(i*basis[0] + j * basis[1] + k * basis[2]));
							}
						}
					}
					else {	config.Insert(r, inverseL*(i*basis[0] + j * basis[1]));}
				}
			}
			else { config.Insert(r, inverseL*(i*basis[0])); }
		}
		//Remove particles
		{
			size_t numRemoval = (int)floor((MaxPackingFraction - PackingFraction)*config.PeriodicVolume() / HyperSphere_Volume(dimension, r));
			if (numRemoval >= config.NumParticle() - numPrts) {
				for (int i = config.NumParticle() ; i > numPrts; i--)
					config.DeleteParticle((int)floor(config.NumParticle()*(this->rng[0].RandomDouble())));
			}

		}
		//Shrink particle radius
		{
			double factor = pow(PackingFraction / config.PackingFraction(), 1.0 / dimension);
			config.Rescale_prtOnly(factor);
			std::cout << "initial lattice packing is successfully generated: packing fraction = " << config.PackingFraction() <<std::endl;
			//config.UpdateMaxRadius();
		}
		//Verify the generated packing is overlap-free before returning it.
		if (config.CheckOverlap()) {
			std::cout << "LatticePacking::the generated packing has overlaps; returning an empty packing.\n";
			config = SpherePacking();
		}
	}
	else
		std::cout << "Target packing fraction is too high"<<std::endl;
		
	//return config;
}



void HSF::LatticePacking(int numPrts, double PackingFraction, const std::string &lattice, SpherePacking & config) {
	std::vector<GeometryVector> basis;
	////Determine basis vectors
	//if (strcmp(lattice.c_str(), "Integer") == 0) {
	//	basis.emplace_back(GeometryVector(static_cast<double>(1.0)));
	//}
	//else if (strcmp(lattice.c_str(), "Square") == 0) {
	//	basis.emplace_back(GeometryVector(1.0, 0.0));
	//	basis.emplace_back(GeometryVector(0.0, 1.0));
	//}
	//else if (strcmp(lattice.c_str(), "Triangular") == 0) {
	//	basis.emplace_back(GeometryVector(1.0, 0.0));
	//	basis.emplace_back(GeometryVector(cos(pi/3.0), sin(pi/3.0)));
	//}
	//else if (strcmp(lattice.c_str(), "Cubic") == 0) {
	//	basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
	//	basis.emplace_back(GeometryVector(0.0, 1.0, 0.0));
	//	basis.emplace_back(GeometryVector(0.0, 0.0, 1.0));
	//}
	////else if (strcmp(lattice.c_str(), "FCC") == 0) {
	////	
	////}
	//else {
	//	std::cout << "Please use different lattice" << std::endl;
	//}
	GetPredefinedBasis(lattice, basis);
	LatticePacking(numPrts, PackingFraction, basis, config);
}

HSF::HSF() {
	this->Dimension = 0;
	this->NumPrts = 0;
	this->NumThreads = 1;
	this->PackingFraction = 0;
	this->MaxDispl = 1.0;
	this->Diameter = 0;
	this->pSpheres.reset(new SpherePacking());
	this->rng.emplace_back(RandomGenerator(0));//guarantee one generator so serial Evolve() never reads rng[0] out of bounds
	adjusted = false;	equilibrated = false;
}

HSF::HSF(int Dimension, int NumPrts, double PackingFraction): HSF() {
	std::string initLattice;
	switch (Dimension) {
	case 1:	
		initLattice = std::string("Integer");	 break;
	case 2:
		initLattice = std::string("Square"); /* std::string("Triangular");*/ break;
	case 3:
		initLattice = std::string("Cubic"); break;//initLattice = std::string("fcc_cubic"); break;//std::string("FCC_cubic");	 break;
	default:
		std::cout << "I don't know\n";
		break;
	}

	SpherePacking Config;
	this->ensembleType = 0;
	{
		std::vector<GeometryVector> basis;	std::vector<size_t> nums;
		double L = ceil( pow(NumPrts, 1.0/Dimension));

		LatticeGen temp(initLattice);
		temp.SetRho(1.0);
		temp.GetBasisVectors(basis);

		for (int i = 0; i < Dimension; i++)
			nums.push_back((size_t)floor(L/fabs(basis[i].x[i])));

		/*for (int i = 0; i < Dimension; i++)
			basis[i] = (L/(double)nums[i])*basis[i];*/

		temp.SetNumCells(nums);
		SpherePacking config;
		temp.GenerateP(config);

		/*if (config.CheckOverlap())
			std::cout << "particles are overlapping\n";*/

		//remove excessive particles (deterministically drop the trailing particles)
		for (int i = config.NumParticle(); i > NumPrts; i--)
			config.DeleteParticle(config.NumParticle()-1);
		
		//Generate a cubic box
		Configuration temp2 = GetUnitCubicBox(Dimension);
		temp2.Rescale(L);
		Config = SpherePacking(temp2, 0.0);

		//Move particles into a cubic box
		for (int i = 0; i < NumPrts; i++)
			Config.Insert(config.GetCharacteristics(i), Config.CartesianCoord2RelativeCoord(config.GetCartesianCoordinates(i)));
		
		//scale particle size to "PackingFraction"
		double factor = pow(PackingFraction/Config.PackingFraction(), 1.0/Dimension);
		Config.Rescale_prtOnly(factor);
	}

	this->ReplaceConfig(Config);
}

void HSF::ReplaceConfig(const SpherePacking & newPacking) {
	this->pSpheres.reset(new SpherePacking(newPacking));
	this->Dimension = newPacking.GetDimension();
	this->NumPrts = newPacking.NumParticle();
	this->PackingFraction = this->pSpheres->PackingFraction();
	this->Diameter = 2.0*(this->pSpheres->GetMaxRadius());
	this->MaxDispl = 1.0/pow(this->NumPrts, 1.0 / this->Dimension);
	switch (this->Dimension) {
	case 2:
		this->MaxDispl *= (1.0 - this->PackingFraction / (pi*sqrt(3.0) / 6.0));	break;
	case 3:
		this->MaxDispl *= (1.0 - this->PackingFraction / (pi / (3.0*sqrt(2))));	break;
	default:
		this->MaxDispl *= (1.0 - this->PackingFraction);	break;
	}
	adjusted = false;
	equilibrated = false;
	//A replaced packing is treated as monodisperse until Switch2Polydisperse() is called again.
	this->isMonodisperse = true;
	this->pSphereRadii.clear();

	//if (this->NumThreads > 1)
	//	this->pSpheres->PrepareIterateThroughNeighbors(3.0*this->Diameter);
}

void HSF::RescaleConfig(double rho) {
	double rho_prev = (double)this->NumPrts / this->pSpheres->PeriodicVolume();
	double factor = pow(rho_prev / rho, 1.0/this->Dimension);
	this->pSpheres->Rescale(factor);
	this->Diameter *= factor;
}

double HSF::Evolve() {
	if (this->rng.empty()) {
		std::cerr << "HSF::Evolve: no random generator available.\n";
		return 0;
	}
	if (this->NumThreads > 1) {
		//Parallel evolution is not implemented; fall back to serial with a one-time warning.
		//atomic exchange makes the "warn once" race-free under concurrent calls.
		static std::atomic<bool> warned{false};
		if (!warned.exchange(true)) {
			std::cout << "HSF::Evolve: parallel evolution is not implemented; falling back to serial.\n";
		}
	}
	return Evolve_Serial();
}

double HSF::Evolve_Serial() {
	double acceptance = 0.0;
	double D = this->Diameter;
	double D2 = D * D;
	if (this->isMonodisperse){
		if (ensembleType == 0) {
			GeometryVector dis;
			for (size_t i = 0; i < this->NumPrts; i++) {
				bool overlap = false;
				TrialMove(dis, this->rng[0]);
				dis = dis + this->pSpheres->GetRelativeCoordinates(i);
				this->pSpheres->IterateThroughNeighbors(dis, D, 
					[D2, &overlap, i](const GeometryVector & shift, const GeometryVector & /*latticeShift*/, const signed long * /*PeriodicShift*/, size_t sourceParticles)->void{
					if (sourceParticles != i) { overlap = (shift.Modulus2() < D2); }
				}, &overlap);
				if (!overlap) {
					acceptance++;
					this->pSpheres->MoveParticle(i, dis);
				}
			}
			acceptance /= (double)this->NumPrts;
		}
		else if(ensembleType == 1){
			std::cout << "NPT is not ready yet\n";
		}
	}
	else{
		/* Polydisperse case!! */
		if (ensembleType == 0) {
			GeometryVector dis;
			for (size_t i = 0; i < this->NumPrts; i++) {
				bool overlap = false;
				TrialMove(dis, this->rng[0]);
				dis = dis + this->pSpheres->GetRelativeCoordinates(i);
				this->pSpheres->IterateThroughNeighbors(dis, D, 
					[&overlap, &i, this](const GeometryVector & shift, const GeometryVector & /*latticeShift*/, const signed long * /*PeriodicShift*/, size_t sourceParticles)->void{
					if (sourceParticles != i) { 
						double Dij = this->pSphereRadii[i] + this->pSphereRadii[sourceParticles];
						overlap = (shift.Modulus2() < Dij*Dij); }
				}, &overlap);
				if (!overlap) {
					acceptance++;
					this->pSpheres->MoveParticle(i, dis);
				}
			}
			acceptance /= (double)this->NumPrts;
		}
		else if(ensembleType == 1){
			std::cout << "NPT is not ready yet\n";
		}		

	}

	return acceptance;
}

/*
*	Measure g_2 within a range of D<r<D+dD
*	dD is proportional to .1 * 10^4/ N
*/
double HSF::g2_at_D() {
	double result = 0.0, D_;
	if (this->NumPrts > 10000)
		D_ = this->Diameter*(1.0 + 0.1* (double)10000 / this->NumPrts);
	else
		D_ = this->Diameter*1.2;
	double rho = this->NumPrts / this->pSpheres->PeriodicVolume();
//#pragma omp parallel for 
	for (size_t i = 0; i < this->NumPrts; i++) {
		this->pSpheres->IterateThroughNeighbors(this->pSpheres->GetRelativeCoordinates(i), D_,
			[&result,i](const GeometryVector & /*shift*/, const GeometryVector & /*latticeShift*/, const signed long * /*PeriodicShift*/, size_t Source)->void {
			if(i!=Source)
				result++;
		}
		);
	}

	result = result / (this->NumPrts*rho* (HyperSphere_Volume(this->Dimension, D_) - HyperSphere_Volume(this->Dimension, this->Diameter)));
	return result;
}

double HSF::Pressure() {
	double phi = this->PackingFraction;
	double factor = pow(2, this->Dimension - 1) * phi;
	if (Verbosity > 3) {
		double g2 = 0.0;
		//std::cout << "Scaled-particle\tPercus-Yevick\tnumerical\n";

		//scaled-particle theory
		if (this->Dimension == 1)
			g2 = 1.0 / (1.0 - phi);
		else if (this->Dimension == 2)
			g2 = (1.0 - phi / 2.0) / (1.0 - phi) / (1.0 - phi);
		else if (this->Dimension == 3) {
			g2 = (1.0 - phi / 2.0 + phi * phi / 4.0) / (1.0 - phi) / (1.0 - phi) / (1.0 - phi);
		}
		std::cout << 1 + factor *g2 << "\t";//" : scaled-particle theory\n";

		//PY thoery
		g2 = 0.0;
		if (this->Dimension == 1)
			g2 = 1.0 / (1.0 - phi);
		else if (this->Dimension == 3) {
			g2 = (1.0 + phi / 2.0 + phi * phi / 4.0) / (1.0 - phi) / (1.0 - phi) ;
		}
		else if (this->Dimension == 5) {
			g2 = (pow(1.0 + 18.0*phi + 6.0*phi*phi, 1.5) - 1.0 + 33.0*phi + 87.0*phi*phi + 6.0*phi*phi*phi) / (60.0*phi*(1.0 - phi)*(1.0 - phi)*(1.0 - phi));
		}
		std::cout << 1.0 + factor *g2 << "\t";
	}
	double result = 1.0 + factor*this->g2_at_D();
	return result;
}

double HSF::AdjustMaxDisp() {
	int numCycles = 10;
	int totalCycles = 0;
	if (Verbosity > 2)
		std::cout << "HSF::start optimizing the max. disp.\n";
	double acceptance = 0.0, devi = 0.0;
	if(! this->adjusted){
		do {
			acceptance = this->Evolve(numCycles);
			totalCycles++;
			devi = acceptance - targetAcceptance;
			if (Verbosity>2) {
				std::cout << "Max disp.: " << this->MaxDispl;
				std::cout << ",\tacceptance: " << acceptance << std::endl;
			}
			if (devi < -tolerance)
				this->MaxDispl *= 0.99;
			else if (devi > tolerance)
				this->MaxDispl *= 1.05;
		} while (std::fabs(devi) > tolerance && this->MaxDispl<0.1);
	}
	else {
		std::cout << "Already adjusted\n";
	}
	if (Verbosity > 2) {
		std::cout << "HSF::end optimizing the max. disp.\n";
		totalCycles *= numCycles;
		std::cout << "\t time spent : " << totalCycles << " MC cycles\n";
	}
	adjusted = true;
	return MaxDispl;
}
/*	Consider that the configuration is in equilibrium 
*	when g_2(r) becomes stationary.
*/
int HSF::Equilibrate() {
	int numTest = 20, numCycles = 10, totalCycles = 0;
	if (pSpheres->NumParticle() < 4000)
		numTest = 200;
	std::vector<double> g2D_(numTest, 0.0);
	
	//The estimated number of MC cylces to translate particles by 10*(interparticle distance)
	int numMoves = 500;//(int)ceil(1.0/this->targetAcceptance * pow(10.0 / this->AdjustMaxDisp() / pow(this->NumPrts,1.0/(double)this->pSpheres->GetDimension()), 2));
	if (Verbosity > 2)
		std::cout << "Estimated MC cycles = " << numMoves << "\n";
	if (this->ensembleType == 0) {
		if (Verbosity > 2)
			std::cout << "Equilibration started: NVT ensemble\n";
		auto stat = [&g2D_](double &av, double &std)->void {
			av = 0.0;	std = 0.0;
			for (size_t idx = 0; idx < g2D_.size(); idx++) {
				av += g2D_[idx];	std += g2D_[idx] * g2D_[idx];
			}
			av /= (double)g2D_.size();	std = fabs(std / (double)g2D_.size() - av * av) / ((double)g2D_.size() - 1.0);
			std = sqrt(std);
			for (size_t idx = 1; idx < g2D_.size(); idx++)
				g2D_[idx - 1] = g2D_[idx];
		};
		for (int i = 0; i < numTest-1; i++) {
			this->Evolve(numCycles);
			totalCycles++;
			g2D_[i] = this->g2_at_D();
			this->Evolve(numCycles);
			totalCycles++;
		}

		double relativeError, av; //av = average
		do {
			this->Evolve(numCycles);
			totalCycles += numCycles;
			g2D_[numTest - 1] = this->g2_at_D();
			stat(av, relativeError);
			relativeError /= av;
			if (pSpheres->NumParticle() < 1000)
				relativeError /= 10.0;
			if (Verbosity > 3)
				std::cout << "g2 = " << av << " + " << relativeError << "(%)\n";
		} while (relativeError > error || av <1.0 || totalCycles < numMoves);
		//Repeat Evolve until 
		this->equilibrated = true;//	totalCycles *= numCycles;
		if (Verbosity > 2) {
			std::cout << "HSF::Equilibration ended\n";
			std::cout << "\t time spent : " << totalCycles << " MC cycles\n";
		}
	}

	return totalCycles;
}

void HSF::Prepare() {
	this->pSpheres->SetCellSize(3.0);
	this->Evolve(100);		//Disturb an initial config
	this->AdjustMaxDisp();	//Optimize the max dispalcement

	{
		auto start = std::time(nullptr);
		int numCycles = this->Equilibrate();		//Equilibrate the system
		auto time = std::time(nullptr) - start;
		std::cout << "time spent: " << time << "seconds\n";
		std::cout << "N = " << this->NumPrts << "\t" << time / (double)numCycles << "seconds/MC cycle \n";
	}

}