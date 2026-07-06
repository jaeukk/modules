#include "SpacefillingPacking.h"
#include <HardSphereFluids.h>

#ifdef USE_MATHGL
#include <mgl2/mgl.h>
void PlotPacking2D(std::string prefix, const SpherePacking & List, const char * sphereColor, const std::string & remark);
void PlotCircle2D(mglGraph & graph, double x, double y, double r, const char * MGLOptions, size_t NumLineSegments);
void PlotDisk2D(mglGraph & graph, double x, double y, double r, const char * MGLOptions, size_t NumLineSegments);
#endif

//==============================================================================
// Body sliced verbatim from the legacy TilingProblem.cpp (lines 16-340):
// the SpacefillingPacking methods + getSpacefillingPacking / ...Packing2 drivers.
//==============================================================================
SpacefillingPacking::SpacefillingPacking(const SpherePacking& config) : RSA(config, 0.0)
{
	phi = 0.0;
	N = 0;
	gen.seed(std::time(NULL));
	spherevolume = HyperSphere_Volume(dimension, 1.0);
	voxelSize = 0.01;
}

SpacefillingPacking::SpacefillingPacking(const SpherePacking& config, const std::vector<double> & param, const char * fnName) : RSA(config, 0.0)
{
	phi = -1.0;
	gen.seed(std::time(NULL));
	spherevolume = HyperSphere_Volume(dimension, 1.0);
	std::vector<double> new_para;
	if(param[0]>0.0 && param[1]>1.0){//Determine N and maxR
		double V = this->pSpheres->PeriodicVolume();
		double sum = 0.0;
		if (strcmp(fnName, "Zeta") == 0) {
			sum = boost::math::zeta(param[1]);
			//std::cout << "GetSpacefillingSpherePackings::Zeta param[1] out of range\n";
		}
		else if (strcmp(fnName, "Power") == 0) {
			sum = param[1] / (param[1] - 1.0);
			//std::cout << "GetSpacefillingSpherePackings::Power param[1] out of range\n";
		}
		else {
			std::cout << "GetSpacefillingSpherePackings::unknown scaling function\n";
			return;
		}
		//The number of spheres of the same radius
		double V_prt_max = param[0];
		this-> N = (int)ceil((1.0-config.PackingFraction())*V / (V_prt_max * sum));
		V_prt_max = V / ((double)this->N * sum);
		{//Debug
			std::cout << "error:: " << (1.0 - config.PackingFraction())*V - (N*V_prt_max*sum) << std::endl;
		}
		//Maximal particle radius
		this->maxR = pow(V_prt_max / this->spherevolume, 1.0/(double)this->dimension );

		this->phi = config.PackingFraction() ;
		new_para.push_back(V_prt_max);
		new_para.push_back(param[1]);
	}
	this->setScalingFunction(new_para, fnName);
	this->voxelSize = 0.01;
	this->pSpheres->SetCellSize(2.0*this->maxR);
}
/*
SpacefillingPacking::~SpacefillingPacking(){//Identical to that of RSA
	if (this->pSpheres != nullptr) delete this->pSpheres;
	if (this->pVoxels != nullptr) delete this->pVoxels;
};*/
void SpacefillingPacking::setN(int n) { this->N = n; }
void SpacefillingPacking::setScalingFunction(const std::vector<double> & param, const char * fnName) {
	double r = param[1];
	if (r > 1.0) {
		if (strcmp(fnName, "Zeta") == 0)
			this->f = [r](int n)->double {	return 1.0 / pow(round((double)n), r); };
		else if (strcmp(fnName, "Power") == 0)
			this->f = [r](int n)->double { return 1.0 / pow(r, (double)n-1.0); };
		else {
			this->f = nullptr;
			std::cout << "unknown scaling function!\n";
		}
	}
	else {
		std::cout << "Inappropriate parameter \n";
	}
}

void SpacefillingPacking::setScalingFunction(const std::function <double(int n)> & function) {
	this->f = function;
}

void SpacefillingPacking::setRmax(double r) {
	this->maxR = r;
}


void SpacefillingPacking::RepeatAddition(int numSteps, int n0) {
	double r = 0;
	for (int n = n0; n < n0 + numSteps; n++) {
		//std::cout << n << ":\n";
		r = this->maxR*pow(this->f(n), 1.0/(double)this->dimension );
		this->AddSpheres(r, this->N);
	}
}


void SpacefillingPacking::AddSpheres(double r, int numAdd, int numThreads) {
	bool parallel = (numThreads == 1) ? false : true;
	int targetVoxelNumPerSide = 100;
	double VoxelRatio = 1.0;// 0.33;
	const long OriginalTrial = std::max(3 * (int)std::pow(10, this->dimension), numAdd * 10);
	long trial = OriginalTrial;

	int numAdded = 0, numRemain = numAdd - numAdded;
	int numInit = this->pSpheres->NumParticle();
	this->radius = r;	//Change the current radius
	//saturated packing fraction of RSA 
	double threshold = HyperSphere_Volume(this->dimension, Rc(this->dimension, 1.0 ));
	threshold *= 0.7;

	double expectedCoveringFraction =  this->phi + HyperSphere_Volume(this->dimension, r)*numAdd / this->pSpheres->PeriodicVolume();
	//When a significant portion of space is going to be covered.
	if (expectedCoveringFraction  >= threshold) {
		//Set [VoxelRatio] s.t. VoxelNum ~ targetVoxel Number
		VoxelRatio = (sqrt(this->pSpheres->GetBasisVector(0).Modulus2())/2.0/this->radius/ (double)targetVoxelNumPerSide);
		this->GetVoxelList(VoxelRatio);
		int SplitCount = 0, numVoxel = this->pVoxels->NumVoxel();
		//If the number of voxel is much larger than the number of trials
		if (Verbosity > 3) {
			std::cout << "\tthe number of voxels = " << numVoxel << " : remaining prts = "<<numRemain<<std::endl;
		}
		double acceptance = 0.0;
		bool largeVoxel, repeat, smallVoxel = false;
		int currNum, chang_in_num;
		while (numVoxel != 0 && numRemain > 0)
		{
			largeVoxel = this->pVoxels->GetVoxelSize() > 2.0*this->radius;
			//When voxel size is much larger than the system.
			if (largeVoxel)
				trial = 3 * pow(10, this->dimension);
			else if (acceptance < 0.015 && smallVoxel)
				trial = 0.1*numVoxel;
			else
				trial = std::min((unsigned long) 4*numRemain, this->pVoxels->NumVoxel());
			//we no longer split the voxel but make more trials!!
			

			if (Verbosity > 3) {
				std::cout << "\tVoxel size = " << this->pVoxels->GetVoxelSize() << "\n";
				std::cout << "\tthe number of trials = " << trial << "\n";
			}

			//Try to add (numAdd) particles by using a crude Voxel list
			currNum = this->pSpheres->NumParticle();
			//int numLimit = std::min(numRemain, 1000);
			int minInsertedNum2SplitVoxel = std::max((int)round(0.015*trial),1);
			this->RSA_II(minInsertedNum2SplitVoxel, trial, parallel);
			//Remove particles when there are too many particles are added
			numAdded = this->pSpheres->NumParticle() - numInit;
			while (numAdded > numAdd) {
				this->pSpheres->DeleteParticle(this->pSpheres->NumParticle() - 1);
				numAdded--;
			}
			chang_in_num = this->pSpheres->NumParticle() - currNum;
			acceptance = (double)chang_in_num / trial;
			numVoxel = this->pVoxels->NumVoxel();
			if (Verbosity > 3) {
				std::cout << "\tInserted particles = "<<chang_in_num<<"\n";
				std::cout << "\tAcceptance = " << acceptance << "\n";
			}
			numRemain = numAdd - numAdded;
			/*if (numRemain > 0 && acceptance < 0.015) {
				this->SplitVoxel(parallel, 1);
				if (smallParticles && numAdded < numAdd)
					this->voxelSize = this->pVoxels->GetVoxelSize();
				SplitCount++;
				numVoxel = this->pVoxels->NumVoxel();

				if (Verbosity>3)
					std::cout << "at time" << std::time(NULL) - ::ProgramStart << ", " << "current voxel:" << this->pVoxels->NumVoxel() << '\n';
			}
			else {
				if(Verbosity >3)
					std::cout << "\n";
			}*/
			/*
			Don't splt the voxel
			When the acceptance rate is high, or the voxel size is very small and there are many empty voxels
			*/
			smallVoxel = (this->pVoxels->GetVoxelSize() < r*0.01) && (numVoxel > MAX_NUM_VOXEL);
			repeat = smallVoxel;// (acceptance > 0.015) || smallVoxel;
			if (numVoxel != 0 && numRemain > 0) {
				if (repeat) {
				}
				else {//Split the voxel
					this->SplitVoxel(parallel, 1);
					SplitCount++;
					numVoxel = this->pVoxels->NumVoxel();
					if (Verbosity>3)
						std::cout << "Voxel splitted at time" << std::time(NULL) - ::ProgramStart << ", " << "current voxel:" << numVoxel << '\n';
				}
			}
		}
		//std::cout << "at time" << std::time(NULL) - ::ProgramStart << ", " << this->N - numAdded << " Spheres lacked\n";
		//If particles are large, we need to udpate the voxel list in every iterations.
		//if (largeParticles)
		//	this->pVoxels = nullptr;
	}
	//Slightly covered => check and insert particles in brutal force!
	else {
		//RSA_II_Serial(numAdd, trial);
		this->RSA_I(numAdd+numInit);
	}
	numAdded = this->pSpheres->NumParticle() - numInit;

	//Update covering fraction
	double particleVol_added = this->spherevolume * pow(r, (double)this->dimension) *(double)numAdded;//(double)numAdded;
	this->phi += particleVol_added / this->pSpheres->PeriodicVolume();
}

double SpacefillingPacking::GetTheoreticalPackingFraction(int stage) {
	double volume_fund_cell = this->pSpheres->PeriodicVolume(), volume_max_prt = HyperSphere_Volume(this->dimension, this->maxR);
	double sum = 0.0;

	for (int i = 1; i <= stage; i++)
		sum += this->f(i);

	return volume_max_prt * this->N / volume_fund_cell * sum;
}

double SpacefillingPacking::getPackingFraction() {
	if (phi < 0)
		this->phi = this->pSpheres->PackingFraction();

	return this->phi;
}


double getSpacefillingPacking(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, const std::vector<int> & outputIndex) {
	std::vector<double> p;	p.push_back(maxPrtVol);	p.push_back(parameter);
	SpacefillingPacking info(c, p, fnName);
	
	//What is the expected covering fraction?
	double phiMax = HyperSphere_Volume(c.GetDimension(), info.maxR);
	{
		double sum = 0.0;
		for (int i = 1; i <= numStages; i++) {
			sum += info.f(i);
		}
		phiMax *= sum;
	}

	double r = 0.0;
	size_t tick = 0;
	for (int j = 1; j < numStages+1; j++) {
		r = info.maxR*pow(info.f(j), 1.0 / (double)info.dimension);
		info.AddSpheres(r, info.N);

		if (outputIndex.size() > tick && outputIndex[tick] == j) {
			char name[300] = {};
			sprintf(name, "%s__%d", prefix,j);
			//ofile <<j<<":\t"<<info.getPackingFraction() << "\n";
			WriteConfiguration(*info.pSpheres, name);
			tick++;
		}

	}
	{
		char name[300] = {};
		sprintf(name, "%s__%d", prefix, numStages);
		//ofile <<j<<":\t"<<info.getPackingFraction() << "\n";
		std::cout << name << "\n";
		WriteConfiguration(*info.pSpheres, name);
	}

	double error = info.getPackingFraction() - phiMax;
	return error;
}


double getSpacefillingPacking2(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads, const std::vector<int> & outputIndex) {
	return getSpacefillingPacking2(c, 1, numStages, maxPrtVol, parameter, fnName, prefix, randomSeed, numThreads, outputIndex);
}


double getSpacefillingPacking2(SpherePacking & c, int initStages, int finStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads, const std::vector<int> & outputIndex) {
	std::vector<double> p;	p.push_back(maxPrtVol);	p.push_back(parameter);
	SpherePacking tempConfig(c);	tempConfig.RemoveParticles();
	SpacefillingPacking info(tempConfig, p, fnName);
	info.gen.seed(randomSeed);

	//EXPECTED VOLUME FRACTION
	//note that phi(0) = 0.0;
	//The packing fraction at the final stage in theory
	double phiMax = info.GetTheoreticalPackingFraction(finStages);
		
	info.ReplaceConfig(c);
	double r = 0.0;
	size_t tick = 0;

	double PHI_THEORY = info.GetTheoreticalPackingFraction(initStages-1);
	double PHI_RESIDUE = 0.0;
	int numMorePrts = 0;
	for (int j = initStages; j < finStages + 1; j++) {
		r = info.maxR*pow(info.f(j), 1.0 / (double)info.dimension);
		PHI_RESIDUE = PHI_THEORY - info.getPackingFraction();
		numMorePrts = (int)floor(PHI_RESIDUE / (info.spherevolume*pow(r, (double)info.dimension)));
		
		if(Verbosity>3){//Debug
		std::cout << j << " stage: r=" << r << "\tN=" <<info.N+numMorePrts<< "\tphi="<<info.getPackingFraction()<<std::endl;
		}
		info.AddSpheres(r, info.N + numMorePrts, numThreads);

		//OUTPUT configurations at certain indices
		if (outputIndex.size() > tick && outputIndex[tick] == j) {
			char name[300] = {};
			if (tick >0) {//delete the temporary configurations
				sprintf(name, "%s__%d.txt", prefix, outputIndex[tick-1]);
				std::remove(name);
			}
			sprintf(name, "%s__%d", prefix, j);
			//ofile <<j<<":\t"<<info.getPackingFraction() << "\n";
			WriteConfiguration(*info.pSpheres, name);
			tick++;
		}
		PHI_THEORY = info.GetTheoreticalPackingFraction(j);
	}
	/*OUPUT index is not specified => print the final config
	if (outputIndex.empty()) {
		char name[300] = {};
		sprintf(name, "%s__%d", prefix, finStages);
		//ofile <<j<<":\t"<<info.getPackingFraction() << "\n";
		//std::cout << name << "\n";
		WriteConfiguration(*info.pSpheres, name);
	}*/
	c = SpherePacking(*info.pSpheres);
	double error = info.getPackingFraction() - phiMax;
	
	return error;
}
