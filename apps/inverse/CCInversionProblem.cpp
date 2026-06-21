#include "CCInversionProblem.h"

/** \brief Generate a Poisson (uniformly random) point configuration in a unit cubic box scaled to N particles. (CLAUDE)
 *  \param dim  Spatial dimension.
 *  \param N    Number of points to insert (the box is resized so that the number density is unity).
 *  \param g    Random generator used to draw the relative coordinates.
 *  \return     The generated Poisson configuration. */
Configuration GeneratePoisson(int dim, int N, RandomGenerator & g){
	Configuration c = GetUnitCubicBox(dim);
	c.Resize((double)N);
	GeometryVector p(static_cast<int>(dim));
	for (int i = 0; i < N; i++){
		for (int j = 0; j < dim ; j++)
			p.x[j] = g.RandomDouble();

		c.Insert("a", p);
	}
	return c;
}


bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t numConstraints, std::vector<double> &results, std::fstream & outFile, double Tolerance, size_t NKPointsSkip, bool RealImagParts, char *method) {
	DimensionType d = target.GetDimension();
	size_t NumParticle = target.NumParticle();
	RandomGenerator gen((int)std::time(NULL));

	//if RealImagPart = true (consider them together) => numKpoints  = numConstraints
	//if RealImagPart = false (consider two parts separately) => numKpoints = 0.5* numConstraints
	size_t N_kPoints = (size_t)round(numConstraints*((RealImagParts) ? 1.0 : 0.5) + 0.1);
	//generate k point list
	double KLimit = 2 * pi*std::pow(1.2*N_kPoints / target.PeriodicVolume() / HyperSphere_Volume(d, 1.0), 1.0 / d);
	std::vector<GeometryVector> kPoints;
	while (kPoints.size() < N_kPoints + NKPointsSkip)
	{
		kPoints = GetKs(target, KLimit, KLimit);
		KLimit *= 1.5;
	}
	auto compareFunc = [](const GeometryVector & left, const GeometryVector & right) -> bool
	{ return left.Modulus2() < right.Modulus2(); };
	std::sort(kPoints.begin(), kPoints.end(), compareFunc);
	size_t num = numConstraints*((RealImagParts) ? 2 : 1);
	FrankInvertCCPotential_re pot(d, num);	pot.ParallelNumThread = 6;
	for (size_t i = 0; i < N_kPoints; i++)
	{
		GeometryVector k = kPoints[i + NKPointsSkip];
		std::complex<double> rho;
		for (size_t i = 0; i<NumParticle; i++)
			rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));
		pot.constraints.push_back(FrankInvertCCPotential_re::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
	}

	double maxDr = 0.0;
	double maxE = 0;
	int count = 0;
	//generate configurations and test
	for (size_t i = 0; i < Nc; i++)
	{
		double E = 100.0;
		Configuration temp = target;
		while (E > 1e-20)
		{
			temp.RemoveParticles();
			for (size_t j = 0; j < NumParticle; j++)
				temp.Insert("A", gen);
			if (strcmp(method, "NLOP+MINOP") == 0) {
				RelaxStructure_NLOPT(temp, pot, 0.0, 0, 0.0);
				RelaxStructure_MINOP(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=NLOP+MINOP\n";
			}
			else if (strcmp(method, "SteepestDescent") == 0) {
				RelaxStructure_SteepestDescent(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=Steepest Descent\n";
			}
			else if (strcmp(method, "ConjugateGradient") == 0) {
				RelaxStructure_ConjugateGradient(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=ConjugateGradient\n";
			}
			else if (strcmp(method, "LocalGradientDescent") == 0) {
				RelaxStructure_LocalGradientDescent(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=LocalGradientDescent\n";
			}
			else {
				RelaxStructure_NLOPT(temp, pot, 0.0, 0, 0.0);
				RelaxStructure_MINOP(temp, pot, 0.0, 0, 0.0);
				//std::cout <<"default method=NLOP+MINOP\n";
			}

			pot.SetConfiguration(temp);
			E = pot.Energy();
			count++;
		}
		maxE = (E > maxE) ? E : maxE;

		for (size_t j = 0; j < NumParticle; j++)
		{
			GeometryVector r = target.GetRelativeCoordinates(j);
			double dr = temp.NearestParticleDistance(r);
			maxDr = std::max(dr, maxDr);

			//debug temp
			//std::cout << "E=" << E << std::endl;
			//if (maxDr > Tolerance)
			//            std::cout << "nonuniqueness found\n";
		}
	}
	//outFile << "max dr=" << maxDr<<", max E="<<maxE<<", numTrials="<<count;
	results.resize(0);
	results.push_back(NumParticle);
	results.push_back(numConstraints);
	results.push_back(maxDr);
	results.push_back(maxE);
	results.push_back(count);

	if (maxDr > Tolerance)
	{
		return false;
	}
	return true;
}


bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::fstream & outFile, double Tolerance, size_t NKPointsSkip) {
	std::vector<double> dummy;
	return FrankInvCCUniqueTest(target, Nc, N_kPoints, dummy, outFile, Tolerance, NKPointsSkip, true);
}


bool FrankInvCCUniqueTest2(const Configuration & target, size_t Nc, size_t numConstraints, std::vector<double> &results, std::fstream & outFile, double Tolerance, size_t NKPointsSkip, bool RealImagParts, char *method) {
	DimensionType d = target.GetDimension();
	size_t NumParticle = target.NumParticle();
	RandomGenerator gen(std::time(NULL) + (int)floor(target.GetCartesianCoordinates(1).x[0]));

	//if RealImagPart = true (consider them together) => numKpoints  = numConstraints
	//if RealImagPart = false (consider two parts separately) => numKpoints = 0.5* numConstraints
	size_t N_kPoints = (size_t)round(numConstraints*((RealImagParts) ? 1.0 : 0.5) + 0.1);
	//generate k point list
	double KLimit = 2 * pi*std::pow(1.2*N_kPoints / target.PeriodicVolume() / HyperSphere_Volume(d, 1.0), 1.0 / d);
	std::vector<GeometryVector> kPoints;
	while (kPoints.size() < N_kPoints + NKPointsSkip)
	{
		kPoints = GetKs(target, KLimit, KLimit);
		KLimit *= 1.5;
	}
	auto compareFunc = [](const GeometryVector & left, const GeometryVector & right) -> bool
	{ return left.Modulus2() < right.Modulus2(); };
	std::sort(kPoints.begin(), kPoints.end(), compareFunc);
	size_t num = numConstraints*((RealImagParts) ? 2 : 1);
	FrankInvertCCPotential_re pot(d, num);	pot.ParallelNumThread = omp_get_num_threads();
	for (size_t i = 0; i < N_kPoints; i++)
	{
		GeometryVector k = kPoints[i + NKPointsSkip];
		std::complex<double> rho;
		for (size_t i = 0; i<NumParticle; i++)
			rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));
		pot.constraints.push_back(FrankInvertCCPotential_re::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
	}

	double maxDr = 0.0;
	double maxE = 0;
	int count_trials = 0;
	int count_suc_reconst = 0;
	//generate configurations and test
	for (size_t i = 0; i < Nc; i++)
	{
		double E = 100.0;
		Configuration temp = target;
		while (E > 1e-20)
		{
			temp.RemoveParticles();
			for (size_t j = 0; j < NumParticle; j++)
				temp.Insert("A", gen);
			if (strcmp(method, "NLOP+MINOP") == 0) {
				RelaxStructure_NLOPT(temp, pot, 0.0, 0, 0.0);
				RelaxStructure_MINOP(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=NLOP+MINOP\n";
			}
			else if (strcmp(method, "SteepestDescent") == 0) {
				RelaxStructure_SteepestDescent(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=Steepest Descent\n";
			}
			else if (strcmp(method, "ConjugateGradient") == 0) {
				RelaxStructure_ConjugateGradient(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=ConjugateGradient\n";
			}
			else if (strcmp(method, "LocalGradientDescent") == 0) {
				RelaxStructure_LocalGradientDescent(temp, pot, 0.0, 0, 0.0);
				//std::cout << "method=LocalGradientDescent\n";
			}
			else {
				RelaxStructure_NLOPT(temp, pot, 0.0, 0, 0.0);
				RelaxStructure_MINOP(temp, pot, 0.0, 0, 0.0);
				//std::cout <<"default method=NLOP+MINOP\n";
			}

			pot.SetConfiguration(temp);
			E = pot.Energy();
			count_trials++;
		}
		maxE = (E > maxE) ? E : maxE;
		double maxDr_temp = 0;
		for (size_t j = 0; j < NumParticle; j++)
		{
			GeometryVector r = target.GetRelativeCoordinates(j);
			double dr = temp.NearestParticleDistance(r);
			maxDr_temp = std::max(dr, maxDr_temp);
			//debug temp
			//std::cout << "E=" << E << std::endl;
			//if (maxDr > Tolerance)
			//            std::cout << "nonuniqueness found\n";
		}
		if (maxDr_temp < Tolerance) {
			count_suc_reconst++;
		}
		maxDr = std::max(maxDr_temp, maxDr);

	}
	//if () {
	//	outFile << "max dr=" << maxDr << ", max E=" << maxE << ", numTrials=" << count_trials;
	//}

	results.resize(0);
	results.push_back(NumParticle);
	results.push_back(numConstraints);
	results.push_back(maxDr);
	results.push_back(maxE);
	results.push_back(count_trials);
	results.push_back(count_suc_reconst);


	if (maxDr > Tolerance)
	{
		return false;
	}
	return true;
}


//constType =	"both" (default): numConstraints = num kPoints 
//				"real" : numConstraints = num real and imaginary parts in collective coordinates, consider real part if numConstraints is odd.
//			    "image": numConstrainst = num real and imaginary parts in collective coordinates, consider imaginary part if numConstraints is odd.

bool FrankInvCCUniqueTest3(const Configuration & target, size_t Nc, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance, double EnergyTolerance, size_t NKPointsSkip, char * constType, int method, std::vector<std::vector<double>> & energy, const std::vector<Configuration> & initConfigs, std::vector<std::vector<double>> & proximity) {
	DimensionType d = target.GetDimension();
	size_t NumParticle = target.NumParticle();
	RandomGenerator gen(std::time(NULL) + (int)floor(target.GetCartesianCoordinates(1).x[0]));
	bool doCalculateProximity = !proximity.empty();
	if (doCalculateProximity) {
		proximity.resize(0);
		//proximity[0] = list of proximity metrics of all solutions w.r.t. the target 
		//proximity[1] = list of proximity metrics for trivial solutions
		//proximity[2] = list of proximity metrics for nontrivial solutions
		proximity.resize(2);
	}
	

	sols.resize(0);
	{
		//first solutions....
		Configuration targetC(target);
		Merge_Sort(targetC, 0);
		sols.push_back(targetC);
		//-------
	}
	energy.resize(0);
	{
		std::vector<double> ttt;	//ttt.push_back(E);
		energy.push_back(ttt);
	}
	char input[100] = {};
	if (strcmp(constType, "real") == 0) {
		sprintf(input, "real");
	}
	else if (strcmp(constType, "image") == 0) {
		sprintf(input, "image");
	}
	else {
		sprintf(input, "both");
	}

	FrankInvertCCPotential2 pot_(d, numConstraints);
	pot_.SetConstraintType(constType);
	size_t N_kPoints = pot_.GetnumKPoints();
	pot_.ParallelNumThread = 1;
	//generate k point list
	double KLimit = 2 * pi*std::pow(1.2*N_kPoints / target.PeriodicVolume() / HyperSphere_Volume(d, 1.0), 1.0 / d);
	std::vector<GeometryVector> kPoints;
	while (kPoints.size() < N_kPoints + NKPointsSkip)
	{
		kPoints = GetKs(target, KLimit, KLimit);
		KLimit *= 1.5;
	}
	auto compareFunc = [](const GeometryVector & left, const GeometryVector & right) -> bool {
		return left.Modulus2() < right.Modulus2();
	};
	std::sort(kPoints.begin(), kPoints.end(), compareFunc);
	//pot_.ParallelNumThread = omp_get_num_threads();

	for (size_t i = 0; i < N_kPoints; i++)
	{
		GeometryVector k = kPoints[i + NKPointsSkip];
		std::complex<double> rho;
		for (size_t i = 0; i<NumParticle; i++)
			rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));
		pot_.constraints.push_back(FrankInvertCCPotential2::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
	}

	double maxDr = 0.0;
	double maxE = 0;
	double minE = EnergyTolerance;
	int count_ground = 0;
	int count_suc_reconst = 0;

	//For sols[i], E_sols[i].x[0] = min E, E_sols[i].x[1] = max E, E_sols[i].x[2] = max dx
	std::vector<GeometryVector> E_sols;
	E_sols.push_back(GeometryVector(EnergyTolerance, 0.0, 0.0));
	//generate configurations and test
	
	for (size_t i = 0; i < Nc; i++)
	{
		double E = 0;//100.0;
		Configuration temp = target;
		//while (E > 1e-20)
		//{
		temp.RemoveParticles();
		if (initConfigs.size()!=0) {
			temp = initConfigs[i];
			//ConfigurationPack initConfigsSet(initConfigs);
			//if (initConfigsSet.NumConfig() > Nc) {//there are enough number of input configs
			//	temp = initConfigsSet.GetConfig(i);
			//	didYouLoad = true;
			//}
		}
		else {//there is no input configs or num of input configs is not enough...
			for (size_t j = 0; j < NumParticle; j++)
				temp.Insert("A", gen);
		}


		if (method == 0){ //strcmp(method, "BFGS+MINOP") == 0) { //0
			RelaxStructure_NLOPT(temp, pot_, 0.0, 0, 0.0);
			RelaxStructure_MINOP(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=NLOP+MINOP\n";
		}
		else if (method == 1){//strcmp(method, "Steepest Descent") == 0) { //1
			RelaxStructure_SteepestDescent(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=Steepest Descent\n";
		}
		else if (method == 2) {//strcmp(method, "ConjugateGradient") == 0) { //2
			RelaxStructure_ConjugateGradient(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=ConjugateGradient\n";
		}
		else if (method == 3) {//strcmp(method, "LocalGradientDescent") == 0) { //3
			RelaxStructure_LocalGradientDescent(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=LocalGradientDescent\n";
		}
		else {
			RelaxStructure_NLOPT(temp, pot_, 0.0, 0, 0.0);
			RelaxStructure_MINOP(temp, pot_, 0.0, 0, 0.0);
			//std::cout <<"default method=NLOP+MINOP\n";
		}
		pot_.SetConfiguration(temp);
		E = pot_.Energy();
		//count_trials++;
		//}


		if (E < EnergyTolerance) {//ground state!!
			count_ground++;
			maxE = std::max(E, maxE);
			minE = std::min(E, minE);
			bool isDegenerated = false;
			size_t identindex = 0;
			double maxDr_temp = 0;
			double proximityMetric = 0;

			for (int n = 0; n < sols.size() && !isDegenerated; n++) {
				maxDr_temp = 0;
				//proximityMetric = 0;
				for (size_t j = 0; j < NumParticle; j++)
				{
					GeometryVector r = sols[n].GetRelativeCoordinates(j);
					double dr = temp.NearestParticleDistance(r);
					maxDr_temp = std::max(dr, maxDr_temp);
				
					//if (doCalculateProximity)
					if (doCalculateProximity && n == 0)
						proximityMetric += dr*dr; //Proximity measure with respect to the target!!
					//debug temp
					//std::cout << "E=" << E << std::endl;
					//if (maxDr > Tolerance)
					//            std::cout << "nonuniqueness found\n";
				}

				if (doCalculateProximity)
					proximityMetric /= (double)NumParticle;

				bool areTheySame = (maxDr_temp < PositionTolerance);
				if (areTheySame) {
					if (n == 0) {//same as the target configuration
						count_suc_reconst++;
						maxDr = std::max(maxDr_temp, maxDr);
					}
					//energy distribution of n-th solution config
					energy[n].push_back(E); 
					
					E_sols[n].x[0] = std::min(E_sols[n].x[0], E);
					E_sols[n].x[1] = std::max(E_sols[n].x[1], E);
					E_sols[n].x[2] = std::max(E_sols[n].x[2], maxDr_temp);

				}

				if (doCalculateProximity) {
					if (n == 0) {//proximity metric of all solutions with respect to the target 
						proximity[0].push_back(proximityMetric);
					}

					if (areTheySame) {//proximity metric of n-th ground states with respect to the target.
						proximity[n + 1].push_back(proximityMetric);
					}
				}

				isDegenerated = isDegenerated || areTheySame;

			}
			if (!isDegenerated) {//when the obtained solution is different from all of known solutions.
				std::vector<double> t_E;	t_E.push_back(E);
				energy.push_back(t_E);
				Configuration tempConfig(temp);
				Merge_Sort(tempConfig, 0);
				sols.push_back(tempConfig);
				E_sols.push_back(GeometryVector(E, E, 0.0));
				
				if (doCalculateProximity) {
					std::vector<double> dummy;
					proximity.push_back(dummy);
				}

					
			}

		}

	}
	//if () {
	//	outFile << "max dr=" << maxDr << ", max E=" << maxE << ", numTrials=" << count_trials;
	//}

	results.resize(0);
	results.push_back(NumParticle); //0
	results.push_back(numConstraints);	//1
	results.push_back(maxDr);	//2
	results.push_back(maxE);	//3
	results.push_back(count_ground);	//4	
	results.push_back(count_suc_reconst);//5
	for (size_t i = 0; i < E_sols.size(); i++) {
		results.push_back(E_sols[i].x[0]); //Sol[i]:min energy
		results.push_back(E_sols[i].x[1]); //Sol[i]:max energy
		results.push_back(E_sols[i].x[2]); //Sol[i]:max dx
	}

	//----SORT ENERGY
	//energy[0]: for trivial solution
	//energy[1]: for secondary
	//energy[i>1]: for others 
	{
		int index_secSol = 1;
		int numSols = 0;
		for (int i = 0; i < energy.size(); i++) {//order energy[1]: secondary...
			MergeSort(energy[i]);
			if (i != 0 && energy[i].size() >= numSols) {
				index_secSol = i;	numSols = energy[i].size();
			}
		}
		if (energy.size() >= 4) {//there are more than 2 other solutions; combine unimportant ones into a single data
			for (int i = 3; i < energy.size(); i++) {
				for (int j = 0; j < energy[i].size(); j++) {
					energy[2].push_back(energy[i].at(j));
				}
			}
			while (energy.size()>3) {
				energy.pop_back();
			}
			MergeSort(energy[2]);
		}
	}


	if (maxDr > PositionTolerance)
	{
		return false;
	}
	return true;
}

//linear constraints: numConstraints = lenght of each lines
//spherical constraints: numConstraints = the total number of constraints.
bool FrankInvCCUniqueTest4(const Configuration & target, size_t Nc, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance, double EnergyTolerance, size_t NKPointsSkip, char * constType, char *method, char * initConfigs) {
	DimensionType d = target.GetDimension();
	size_t NumParticle = target.NumParticle();
	size_t M = 0; //the number of available constraints
	RandomGenerator gen(std::time(NULL) + (int)floor(target.GetCartesianCoordinates(1).x[0]));

	sols.resize(0);
	{//first solutions....
		Configuration targetC(target);
		Merge_Sort(targetC, 0);
		sols.push_back(targetC);
	}


	char input[100] = {};
	FrankInvertCCPotential pot_(d);
	if (strcmp(constType, "linear") == 0) {
		sprintf(input, "linear");

		M = (d + d*(d - 1) / 2)*numConstraints; //d + Binomial(d,2)
		pot_.constraints.reserve(M);
		std::vector<GeometryVector> basis; //basis vectors in the Fourier space
		for (int i = 0; i < d; i++) {
			basis.push_back(target.GetReciprocalBasisVector(i));
		}
		//set constraints
		for (int i = 0; i < d; i++) {
			for (int j = 0; j <= i; j++) {
				GeometryVector k;
				k.Dimension = d;
				for (int L = 1 + (int)NKPointsSkip; L <= (int)NKPointsSkip + (int)numConstraints; L++) {
					if (i == j) //k_i axes
						k = (double)L*basis[i];
					else //k_i + k_j axes
						k = (double)L*(basis[i] + basis[j]);
				}
				std::complex<double> rho;
				for (size_t i = 0; i<NumParticle; i++)
					rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));

				pot_.constraints.push_back(FrankInvertCCPotential2::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
			}
		}

	}
	else {
		sprintf(input, "spherical");
		M = numConstraints;
		double KLimit = 2 * pi*std::pow(1.2*M / target.PeriodicVolume() / HyperSphere_Volume(d, 1.0), 1.0 / d);
		std::vector<GeometryVector> kPoints;
		while (kPoints.size() < M + NKPointsSkip)
		{
			kPoints = GetKs(target, KLimit, KLimit);
			KLimit *= 1.5;
		}
		auto compareFunc = [](const GeometryVector & left, const GeometryVector & right) -> bool {
			return left.Modulus2() < right.Modulus2();
		};
		std::sort(kPoints.begin(), kPoints.end(), compareFunc);
		//pot_.ParallelNumThread = omp_get_num_threads();

		for (size_t i = 0; i < M; i++)
		{
			GeometryVector k = kPoints[i + NKPointsSkip];
			std::complex<double> rho;
			for (size_t i = 0; i<NumParticle; i++)
				rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));
			pot_.constraints.push_back(FrankInvertCCPotential2::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
		}
		M = pot_.constraints.size();
	}
	pot_.ParallelNumThread = 1;
	//generate k point list


	double maxDr = 0.0;
	double maxE = 0;
	double minE = EnergyTolerance;
	int count_ground = 0;
	int count_suc_reconst = 0;

	//For sols[i], E_sols[i].x[0] = min E, E_sols[i].x[1] = max E, E_sols[i].x[2] = max dx
	std::vector<GeometryVector> E_sols;
	E_sols.push_back(GeometryVector(EnergyTolerance, 0.0, 0.0));
	//generate configurations and test
	for (size_t i = 0; i < Nc; i++)
	{
		double E = 0;//100.0;
		Configuration temp = target;
		//while (E > 1e-20)
		//{
		bool didYouLoad = false;
		temp.RemoveParticles();
		if (initConfigs != nullptr) {
			ConfigurationPack initConfigsSet(initConfigs);
			if (initConfigsSet.NumConfig() > Nc) {//there are enough number of input configs
				temp = initConfigsSet.GetConfig(i);
				didYouLoad = true;
			}
		}
		if (!didYouLoad) {//there is no input configs or num of input configs is not enough...
			for (size_t j = 0; j < NumParticle; j++)
				temp.Insert("A", gen);
		}

		if (strcmp(method, "NLOP+MINOP") == 0) {
			RelaxStructure_NLOPT(temp, pot_, 0.0, 0, 0.0);
			RelaxStructure_MINOP(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=NLOP+MINOP\n";
		}
		else if (strcmp(method, "SteepestDescent") == 0) {
			RelaxStructure_SteepestDescent(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=Steepest Descent\n";
		}
		else if (strcmp(method, "ConjugateGradient") == 0) {
			RelaxStructure_ConjugateGradient(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=ConjugateGradient\n";
		}
		else if (strcmp(method, "LocalGradientDescent") == 0) {
			RelaxStructure_LocalGradientDescent(temp, pot_, 0.0, 0, 0.0);
			//std::cout << "method=LocalGradientDescent\n";
		}
		else {
			RelaxStructure_NLOPT(temp, pot_, 0.0, 0, 0.0);
			RelaxStructure_MINOP(temp, pot_, 0.0, 0, 0.0);
			//std::cout <<"default method=NLOP+MINOP\n";
		}
		pot_.SetConfiguration(temp);
		E = pot_.Energy();
		//count_trials++;
		//}


		if (E < EnergyTolerance) {//ground state!!
			count_ground++;
			maxE = std::max(E, maxE);
			minE = std::min(E, minE);
			bool isDegenerated = false;
			size_t n = 0;
			double maxDr_temp = 0;
			for (n = 0; n < sols.size() && !isDegenerated; n++) {
				maxDr_temp = 0;

				for (size_t j = 0; j < NumParticle; j++)
				{
					GeometryVector r = sols[n].GetRelativeCoordinates(j);
					double dr = temp.NearestParticleDistance(r);
					maxDr_temp = std::max(dr, maxDr_temp);
					//debug temp
					//std::cout << "E=" << E << std::endl;
					//if (maxDr > Tolerance)
					//            std::cout << "nonuniqueness found\n";
				}

				bool areTheySame = (maxDr_temp < PositionTolerance);
				if (areTheySame) {//same as the first solution (target configuration)
					if (n == 0) {
						count_suc_reconst++;
						maxDr = std::max(maxDr_temp, maxDr);
					}
					E_sols[n].x[0] = std::min(E_sols[n].x[0], E);
					E_sols[n].x[1] = std::max(E_sols[n].x[1], E);
					E_sols[n].x[2] = std::max(E_sols[n].x[2], maxDr_temp);
				}

				isDegenerated = isDegenerated || areTheySame;

			}
			if (!isDegenerated) {//when the obtained solution is different from all of known solutions.
				Configuration tempConfig(temp);
				Merge_Sort(tempConfig, 0);
				sols.push_back(tempConfig);
				E_sols.push_back(GeometryVector(E, E, 0.0));
			}

		}

	}
	//if () {
	//	outFile << "max dr=" << maxDr << ", max E=" << maxE << ", numTrials=" << count_trials;
	//}

	results.resize(0);
	results.push_back(NumParticle); //0
	results.push_back((double)M / (d*(NumParticle - 1)));	//1
	results.push_back(maxDr);	//2
	results.push_back(maxE);	//3
	results.push_back(count_ground);	//4	
	results.push_back(count_suc_reconst);//5
	for (size_t i = 0; i < E_sols.size(); i++) {
		results.push_back(E_sols[i].x[0]);
		results.push_back(E_sols[i].x[1]);
		results.push_back(E_sols[i].x[2]);
	}

	if (maxDr > PositionTolerance)
	{
		return false;
	}
	return true;
}

// dim=1, i=0 : perturbed lattices
// dim=1, i=1 : Poisson 
// dim=1, i=2 : non-Bravais lattice
// dim=2, i=0 :
/** \brief Build a target configuration of a requested type (perturbed lattice, Poisson, or non-Bravais lattice) for the inversion tests. (CLAUDE)
 *  \param i          Index of the target (used to seed the random generator so successive targets differ).
 *  \param type       Target family: 0=perturbed lattice, 1=Poisson, 2=non-Bravais lattice (d=1 only).
 *  \param N          Number of particles in the configuration.
 *  \param basis      Lattice basis vectors; their count sets the spatial dimension.
 *  \param parameter  Perturbation magnitude for type 0 (maximum random displacement); unused by other types.
 *  \param[out] systN  If non-null, receives a short name tag describing the generated system (e.g. "_PL_0.100", "_Poisson", "_nonBravais").
 *  \return           The generated target configuration. */
Configuration GetTarget(size_t i, int type, int N, std::vector<GeometryVector> & basis, double parameter, char * systN) {
	Configuration c;
	int dim = (int)basis.size();
	GeometryVector temp;    temp.SetDimension(dim);
	RandomGenerator g(std::time(NULL) + i);
	if (dim == 1) {
		if (type == 0) {
			c = GetLattice_config(basis, N);
			for (int j = 0; j < N; j++) {
				//int j = (int)floor(N*g.RandomDouble());
				RandomVector(dim, g, temp);
				temp = parameter*temp;
				c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
			}
			if(systN != nullptr)
				sprintf(systN, "_PL_%0.3f", parameter); //"%dD_PL_%0.3f"
		}
		else if (type == 1) {
			c = GeneratePoisson(dim, N, g);
			c.Resize(N);
			if (systN != nullptr)
				sprintf(systN, "_Poisson");
		}
		else if (type == 2) {
			GeometryVector *A = new GeometryVector[1];
			A[0].Dimension = 1;
			A[0].x[0] = N;
			c = Configuration(dim, A, 1.0);
			int numBasisAtom = Find_LeastCommonDivider(N);
			std::vector<GeometryVector> atomBasis;
			GeometryVector X(static_cast<int> (1));
			RandomGenerator g(std::time(NULL) + i);
			for (int i = 0; i < numBasisAtom; i++) {
				X.x[0] = (double)numBasisAtom*g.RandomDouble();
				atomBasis.push_back(X);
			}
			X.x[0] = 0;
			for (int i = 0; i < N; i++) {
				int j = i%numBasisAtom;
				int k = (int)floor(i / numBasisAtom);

				X.x[0] = k*numBasisAtom;
				c.Insert("a", c.CartesianCoord2RelativeCoord(X + atomBasis[j]));
			}
			if (systN != nullptr)
				sprintf(systN, "_nonBravais");
		}
	}
		
	else {//dim >1
		if (type == 0) {
			c = GetLattice_config(basis, (int)pow(N, 1.0 / (double)dim));
			for (int j = 0; j < N; j++) {
				RandomVector(dim, g, temp);
				temp = parameter*temp;
				c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
			}
			if (systN != nullptr)
				sprintf(systN, "_PL_%0.3f", parameter); //dim
 		}
		else if (type == 1) {
			c = Configuration(dim, &basis[0],1.0);
			for(int j=0; j<N; j++)
				c.Insert("a", g);

			if (systN != nullptr)
				sprintf(systN, "_Poisson");
		}
		else {
			std::cout << "dim>1, type>1 are not defined yet";
		}
	}
	return c;
}

#include "Interfaces.h"
#include "Plots.h"
void CCInversionProblem(bool outConfig) {
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char dir0[100] = "results/Reconst/";
	char dir[100] = {};
	char temp[200] = {};
	char name[200] = {};
	
	//---- general parameters
	DimensionType dim = 1;
	int Nc = 100; //the number of trial initial configs
	int numConfigs = 1000; //the number of target configs
	int NKSkip = 0;
	int numThreads = 0;
	char constType[100] = {};
	
	//---- swap parameters
	int N1 = 10;//min system size
	int N2 = 40;//max system size
	int systType1;
	int systType2;
	int dN = 1;
	int methodNum1 = 0;
	int methodNum2 = 0;
	int M1 = -1;
	int M2 = +1;
	double EnergyTolerance = 1e-20;
	double PositionTolerance = 1e-6;
	double delta1;
	double delta2;
	double d_delta;
	bool doProxMetric;
	double proximityUB = PositionTolerance*PositionTolerance;
	{
		ofile << "NKSkip=";
		ifile >> NKSkip;
		ofile << "dim=";
		ifile >> dim;
		ofile << "min N=";
		ifile >> N1;
		ofile << "max N=";
		ifile >> N2;
		ofile << "dN=";
		ifile >> dN;
		ofile << "numConfigs=";
		ifile >> numConfigs;
		ofile << "Nc=";
		ifile >> Nc;
		ofile << "system type1:";
		ifile >> systType1;
		ofile << "system type2:";
		ifile >> systType2;
		ofile << "initial method num: ";
		ifile >> methodNum1;
		ofile << "final method num: ";
		ifile >> methodNum2;
		ofile << "Energ Tolerance = ";
		ifile >> EnergyTolerance;
		ofile << "delta1 = ";
		ifile >> delta1;
		ofile << "delta2 = ";
		ifile >> delta2;
		ofile << "d_delta = ";
		ifile >> d_delta;
		ofile << "M1 = ";
		ifile >> M1;
		ofile << "M2 = ";
		ifile >> M2;
		ofile << "num threads = ";
		ifile >> numThreads;
		ofile << "Do we compute the proximity metric?  ";
		ifile >> doProxMetric;
	}


	{//SHOW PARAMETERS
		ofile << "NKSkip=";
		ofile << NKSkip;
		ofile << "\ndim=";
		ofile << dim;
		ofile << "\nmin N=";
		ofile << N1;
		ofile << "\nmax N=";
		ofile << N2;
		ofile << "\ndN=";
		ofile << dN;
		ofile << "\nnumConfigs=";
		ofile << numConfigs;
		ofile << "\nNc=";
		ofile << Nc;
		ofile << "\nsystem type1:";
		ofile << systType1;
		ofile << "\nsystem type2:";
		ofile << systType2;
		ofile << "\ninitial method num: ";
		ofile << methodNum1;
		ofile << "\nfinal method num: ";
		ofile << methodNum2;
		ofile << "\nEnerg Tolerance = ";
		ofile << EnergyTolerance;
		ofile << "\ndelta1 = ";
		ofile << delta1;
		ofile << "\ndelta2 = ";
		ofile << delta2;
		ofile << "\nd_delta = ";
		ofile << d_delta ;
		ofile << "\nM1 = ";
		ofile << M1 ;
		ofile << "\nM2 = ";
		ofile << M2 ;
		ofile << "\nnum threads = ";
		ofile << numThreads;
		ofile << "\nDo we compute the proximity metric?  ";
		ofile << doProxMetric << std::endl;
	
	}

	if (NKSkip > 0)
		std::cout << "Skip the central " << NKSkip << " k points.\n";

	sprintf(dir, "%sdE_%.0e/Nk_skip_%d/", dir0, EnergyTolerance, NKSkip);
	ofile << dir<<"\n";
	std::vector<GeometryVector> basis;
	for (int j = 0; j < dim;j++) {
		basis.push_back(GeometryVector());
		basis[j].SetDimension(dim);
		basis[j].x[j] = 1.0;
	}
	//----energy distribution----
	std::vector<double> bin;
	std::vector<double> dE;
	std::vector<double> av_X2; //for proximity metrics
	double Energy0 = 1e-34;
	double xdE = 2;
	double Energy = Energy0; //start energy
 	while (Energy < EnergyTolerance) {
		bin.push_back(Energy);
		Energy *= xdE;
		if (bin.size() >= 2) {
			dE.push_back(bin[bin.size()-1]-bin[bin.size()-2]);
		}
	}

	double x0 = 1e-38;
	double d_x2 = 2;
	if (doProxMetric) {
		double x2 = x0;
		while (x2 < 10.0) {
			x2 *= d_x2;
			av_X2.push_back(x2);
		}
	}


	for (double delta = delta1; delta <= delta2; delta += d_delta) {
		for (int method_num = methodNum1; method_num <= methodNum2; method_num++) {
			char method[100] = {};
			switch (method_num)
			{
			case 0:
				sprintf(method, "BFGS+MINOP"); break;
			case 1:
				sprintf(method, "Steepest Descent"); break;
			case 2:
				sprintf(method, "ConjugateGradient"); break;
			case 3:
				sprintf(method, "LocalGradientDescent"); break;
			default:
				break;
			}
			ofile << "method = " << method << std::endl;

			for (int type = systType1; type <= systType2; type++) {
				char systName[200] = {};
				GetTarget(0, type, 1, basis, delta, systName); //just get system name;;;
				ofile << "system type = " << systName << std::endl;
				ofile << "N\tM\tType\tmaxSols\tminSols\tProb_DiffSols\tav.NumSols\tstdSols\tMaxE_ground\tMaxdx\n";

				for (int N = N1; N <= N2; N += dN) {
					char filename[400] = {};
					sprintf(filename, "%s%dD%s__%s__N_%d.txt", dir, dim, systName, method, N);
					std::fstream outfile(filename, std::fstream::out);//
					std::fstream outfile1;
					auto getTargetC = [type,N, &basis, delta](size_t i)->Configuration {
						//RandomGenerator g(std::time(NULL)+i);
						return GetTarget(i, type, N, basis, delta, nullptr);
					};
					
					std::vector<GeometryVector> summaryReal1;
					std::vector<GeometryVector> summaryReal2;
					std::vector<GeometryVector> summaryImag1;
					std::vector<GeometryVector> summaryImag2;
					
					if (dim == 1) {
						for (size_t M = N + M1; M <= N + M2; M++) {
							outfile << "-----------------------------------------------------------------" << std::endl;
							outfile << "Testing N=" << N << "\tM=" << M << std::endl << std::endl;
							outfile << "num\tconstr\tconfigs";
							for (int z = 0; z < N; z++)
								outfile << "\t\t";
							outfile << "\tminE\t\tmaxE\t\tmaxDx" << std::endl;


							//-----------------------final results!!------------------------------------------------
							//results[i].x[0] <= real type / .x[1] <= imag type
							//[0].				maxSols
							//[1].				minSols
							//[2].				probability to gain nontrivial solutions
							//[3].				average number of solutions
							//[4].				std of sols
							//[5].				max energy of ground states
							//[6].				max position errors
							//[7].				success rate: ground states
							//[8].				error_success rate: ground states
							//[9].				success rate: correct solutions
							//[10].				error_success rate: correct solutions
							std::vector<GeometryVector> results; 
							{
								GeometryVector dump(static_cast<int> (2));
								results.assign(11, dump); 
								results[1].x[0] = Nc+1;
								results[1].x[1] = Nc+1;
							}
							int count = 0;
							//energy distributions...
							std::vector<std::vector<GeometryVector>> EnergyDist;
							std::vector<std::vector<GeometryVector>> ProximityDist; //measures of proximity metrics
							
							{
								//n = 0: real, n=1: imag
								//EnergyDist[n+0]: trivial sols
								//EnergyDist[n+1]: nontrivial sols
								//EnergyDist[n+2]: others

								//ProximityMetric[n+0]: trivial sols
								//ProximityMetric[n+1]: nontrivial sols
								//ProximityMetric[n+2]: others
								//ProximityMetric[n+3]: for all solutions w.r.t. the targets

								std::vector<GeometryVector> table;
								GeometryVector element(0, 0, 0, 0);	element.Dimension = 4;
								for (int i = 1; i < bin.size(); i++) {
									element.x[0] = 0.5*(bin[i] + bin[i - 1]);
									table.push_back(element);
								}
								for (int i = 0; i < 6; i++) {
									EnergyDist.push_back(table);
								}

								{
									table.resize(0);
									for (int i = 0; i < av_X2.size(); i++) {
										element.x[0] = (1.0 + 1.0 / d_x2)*0.5*av_X2[i];
										table.push_back(element);
									}
									for (int i = 0; i < 8; i++) {
										ProximityDist.push_back(table);
									}
								}

							}

							
							//temporary storage! version 1.1  buffers!!
							std::vector<std::vector<double>> tab_results; //store all results here!!
							std::vector<std::vector<std::vector<double>>> tab_energy;
							std::vector<std::vector<std::vector<double>>> tab_proxy;
							std::vector<int> numSols1;
							std::vector<int> numSols2;
							int numLegd = 4;
							{//fill storage so that the size of containers doesn't change during the for-loops
								//------------for i-th target configuration

								//tab_results[2*i] =	results with "real" constraints
								//tab_results[2*i+1] =	results with "imag" constraints

								//tab_energy[2*i] =		energy distribution of "real" constraints
								//tab_energy[2*i+1] =	energy distribution of "imag" constraints

								//tab_proximity[2*i] =	proximity distribution of "real" constraints
								//tab_proximity[2*i+1] = proximity distribution of "imag" constraints
								tab_results.resize(0);
								tab_energy.resize(0);
								tab_proxy.resize(0);
								numSols1.resize(0);
								numSols2.resize(0);

								std::vector<double> temp_result;
								std::vector<std::vector<double>> temp_box;
								
								tab_results.resize(2 * numConfigs, temp_result);
								tab_energy.resize(2 * numConfigs, temp_box);
								tab_proxy.resize(2 * numConfigs, temp_box);
								numSols1.resize(numConfigs, 0);
								numSols2.resize(numConfigs, 0);
							}

							
														

#pragma omp parallel for schedule (guided) num_threads(numThreads) shared( tab_results, tab_energy, tab_proxy, numSols1, numSols2)//shared(results, EnergyDist, ProximityDist) //shared(numSols, numDiffNumSols, maxSols, minSols, count, maxE_ground1, maxE_ground2, maxDx1, maxDx2)
							for (int i = 0; i < numConfigs; i++) {
								
								Configuration c = getTargetC(i); //TARGET
								std::vector<Configuration> initConfig;
								//initial configurations
								{
									RandomGenerator g(i + std::time(NULL));
									Configuration templateC(c);	templateC.RemoveParticles();
									for (int j = 0; j < Nc; j++) {
										Configuration t(templateC);
										for (int k = 0; k < N; k++) {
											t.Insert("a", g);
										}
										initConfig.push_back(t);
									}
								}

								//results_temp1[0] = N
								//results_temp1[1] = M
								//results_temp1[2] = maxDr
								//results_temp1[3] = maxE
								//results_temp1[4] = # ground states
								//results_temp1[5] = # correct solutions
								//results_temp1[6+3*i+0] = min energy: sol[i]
								//results_temp1[6+3*i+1] = max energy: sol[i]
								//results_temp1[6+3*i+2] = max dx: sol[i]
								std::vector<double> results_temp1;			//MISC INFO REAL CONSTRAINT
								std::vector<double> results_temp2;			//MISC INFO IMAG CONSTRAINT
								std::vector<Configuration> sols1;			//SOLUTION CONFIG REAL CONSTRAINT
								std::vector<Configuration> sols2;			//SOLUTION CONFIG IMAG CONSTRAINT
								std::vector<std::vector<double>> energy1;	//energy1[i] = list of energy of this solution
								std::vector<std::vector<double>> energy2;	//energy2[i] = list of energy of this solution

								std::vector<std::vector<double>> proximity1;	//proximity1[i] = list of proximity metrics of solution i
								std::vector<std::vector<double>> proximity2;	//proximity2[i] = list of proximity metrics of solution i
								if (doProxMetric) {//just give signal to initiate calculations
									std::vector<double> dummy;
									proximity1.push_back(dummy);
									proximity2.push_back(dummy);
								}
								//-------GENERATE INFO

								FrankInvCCUniqueTest3(c, Nc, M, results_temp1, sols1, outfile1, PositionTolerance, EnergyTolerance, NKSkip, "real", method_num, energy1, initConfig, proximity1);
								for (int j = 0; j < energy1.size(); j++)
									MergeSort(energy1[j]);
								

								if (doProxMetric && proximity1.size()>numLegd) {
									//if there are more than 2 nontrivial solutions
									//merge ternary solutions
									for (int l = proximity1.size() - 1; l >= numLegd; l--) {
										for (int k = 0; k < proximity1[l].size(); k++) {
											proximity1[numLegd - 1].push_back(proximity1[l].at(k));
										}
									}
									proximity1.resize(numLegd);
									for (int j = 0; j < proximity1.size(); j++) {
										MergeSort(proximity1[j]);
									}
								}

								if (M % 2 == 0) {//number of constraints = even => real and imaginary constraints are identical
									results_temp2 = results_temp1;
									sols2 = sols1;
									energy2 = energy1;
									proximity2 = proximity1;
								}
								else {//number of constraints = odd
									FrankInvCCUniqueTest3(c, Nc, M, results_temp2, sols2, outfile1, PositionTolerance, EnergyTolerance, NKSkip, "image", method_num, energy2, initConfig, proximity2);
									for (int j = 0; j < energy2.size(); j++)
										MergeSort(energy2[j]);

									if (doProxMetric && proximity2.size()>numLegd) {
											//if there are more than 2 nontrivial solutions
											//merge ternary solutions
											for (int l = proximity2.size() - 1; l >= numLegd; l--) {
												for (int k = 0; k < proximity2[l].size(); k++) {
													proximity2[numLegd - 1].push_back(proximity2[l].at(k));
												}
											}
											proximity2.resize(numLegd);
											for (int j = 0; j < proximity2.size(); j++) {
												MergeSort(proximity2[j]);
											}
									}
								}
								//--------GATHER INFO
								//SUCCESS RATE: GROUND STATES
								//SUCCESS RATE: UNIQUE SOLUTION
								//MAX SOLS/TARGET
								//MIN SOLS/TARGET
								//MAX ENERGY
								//MIN ENERGY
								//ENERGY DISTRIBUTION: TRIVIAL/SECONDARY/OTHERS....

								


#pragma omp critical
								{
									numSols1[i] = sols1.size();
									numSols2[i] = sols2.size();
								}


//#pragma omp atomic
//								numSols1[i] = sols1.size();
//#pragma omp atomic
//								numSols2[i] = sols2.size();

								//organize output!!
								if (M >= N && outConfig) {
									//TARGET
									Merge_Sort(sols1[0], 0);
									//SOLUTIONS: REAL CONDITION
									//if (sols1.size() > 1) {//multiple solutions in real constraints
									//	for (int t = 1; t < sols1.size(); t++) {
									//		Merge_Sort(sols1[t], 0);
									//	}
									//}

									//SOLUTIONS: IMAG CONDITION
									sols2[0] = sols1[0];
									if (sols2.size() > 1) {//multiple solutions in real constraints
										if (M % 2 == 0) {//when REAL and IMAG are identical
											sols2 = sols1;
										}
										//else {//when REAL and IMAG are distinct
										//	for (int t = 1; t < sols2.size(); t++) {
										//		Merge_Sort(sols2[t], 0);
										//	}
										//}

									}

#pragma omp critical
									{//OUTPUT:SOLS
										outfile << count << "\tTarget\t";
										{//TARGET
											
											for (int z = 0; z < sols1[0].NumParticle(); z++)
												outfile << std::setprecision(7) << sols1[0].GetCartesianCoordinates(z).x[0] << ",\t";

											outfile << "\n";
										}
									
										//REAL CONDITION
										if (sols1.size() == 1) {//UNIQUE SOLUTION
											outfile << count << "\tReal\t";
											if (results_temp1[5] > 0) { outfile << "unique"; }
											else { outfile << "fail"; }
											for (int m = 0; m < N; m++)
												outfile << "\t\t";

											outfile << results_temp1[6] << "\t" << results_temp1[7] << "\t" << results_temp1[8] << std::endl;
										}
										else {//MULTIPLE SOLUTIONS
											for (int t = 0; t < sols1.size(); t++) {
												outfile << count << "\tReal\t";
												for (int y = 0; y < sols1[t].NumParticle(); y++)
													outfile << sols1[t].GetCartesianCoordinates(y).x[0] << ",\t";

												outfile << results_temp1[3 * t + 6] << "\t" << results_temp1[3 * t + 7] << "\t" << results_temp1[3 * t + 8] << std::endl;
											}
										}

										//IMAG
										if (sols2.size() == 1) {//UNIQUE SOL
											outfile << count << "\tImag\t";
											if (results_temp2[5] > 0) { outfile << "unique"; }
											else { outfile << "fail"; }
											for (int m = 0; m < N; m++)
												outfile << "\t\t";

											outfile << results_temp2[6] << "\t" << results_temp2[7] << "\t" << results_temp2[8] << std::endl;
										}
										else {//multiple solutions in real constraints
											for (int t = 0; t < sols2.size(); t++) {
												outfile << count << "\tImag\t";
												for (int y = 0; y < sols2[t].NumParticle(); y++)
													outfile << sols2[t].GetCartesianCoordinates(y).x[0] << ",\t";

												outfile << results_temp2[3 * t + 6] << "\t" << results_temp2[3 * t + 7] << "\t" << results_temp2[3 * t + 8] << std::endl;
											}
										}
										count++;
									}
									

									
									
								}
																		
#pragma omp critical
								{
									tab_results[2 * i].swap(results_temp1);
									tab_results[2 * i + 1].swap(results_temp2);

									tab_energy[2 * i].swap(energy1);
									tab_energy[2 * i + 1].swap(energy2);

									tab_proxy[2 * i].swap(proximity1);
									tab_proxy[2 * i + 1].swap(proximity2);
								}
																	
							}

							
								//SUMMARY!!
//#pragma omp parallel
//							{
//							int priv_max_sols = 0;
//							int priv_min_sols = numConfigs + 1;
//							double rate_diff_sols = 0;
//							double av_numSols = 0;
//							double sq_numSols = 0;
//							double rate
//							}
//#pragma omp parallel for shared(tab_results, results, EnergyDist, numSols1, numSols2)
								for (int i = 0; i < numConfigs; i++) {
									//REAL CONDITION
									double suc_rate_ground1 = tab_results[2 * i][4] / (double)Nc;
									double suc_rate_uniq1 = (tab_results[2 * i][4] > 0.5) ? tab_results[2 * i][5] / tab_results[2 * i][4] : 0;

//#pragma omp atomic
									results[2].x[0] += (1.0 - suc_rate_uniq1);
//#pragma omp atomic
									results[3].x[0] += (double)numSols1[i];
//#pragma omp atomic
									results[4].x[0] += (double)numSols1[i] * numSols1[i];
//#pragma omp atomic
									results[7].x[0] += suc_rate_ground1;
//#pragma omp atomic
									results[8].x[0] += suc_rate_ground1*suc_rate_ground1;
//#pragma omp atomic
									results[9].x[0] += suc_rate_uniq1;
//#pragma omp atomic
									results[10].x[0] += suc_rate_uniq1*suc_rate_uniq1;
						//Max Sols _ C(i) condition			
//#pragma omp flush (results)
									if (results[0].x[0] < numSols1[i]) {
//#pragma omp critical
										{
											if (results[0].x[0] < numSols1[i])	
												results[0].x[0] = numSols1[i];
										}
									}
						//MIN sols _ C(i) condition
//#pragma omp flush (results)
									if (results[1].x[0] > numSols1[i]) {
//#pragma omp critical
										{
											if (results[1].x[0] > numSols1[i])	
												results[1].x[0] = numSols1[i];
										}
									}
						//Max E _ C(i) condition
//#pragma omp flush (results)
									if (results[5].x[0] < tab_results[2 * i][3]) {
//#pragma omp critical
										{
											if (results[5].x[0] <  tab_results[2 * i][3])	results[5].x[0] = tab_results[2 * i][3];
										}
									}
									//Max dx _ C(i) condition

//#pragma omp flush (results)
									if (results[6].x[0] < tab_results[2 * i][2]) {
//#pragma omp critical
										{
											if (results[6].x[0] < tab_results[2 * i][2])	results[6].x[0] = tab_results[2 * i][2];
										}
									}
	
									//IMAG CONDITION
									double suc_rate_ground2 = tab_results[2 * i + 1][4] / (double)Nc;
									double suc_rate_uniq2 = (tab_results[2 * i + 1][4] > 0.5) ? tab_results[2 * i + 1][5] / tab_results[2 * i + 1][4] : 0;

//#pragma omp atomic
									results[2].x[1] += (1.0 - suc_rate_uniq2);
//#pragma omp atomic
									results[3].x[1] += (double)numSols2[i];
//#pragma omp atomic
									results[4].x[1] += (double)numSols2[i] * numSols2[i];
//#pragma omp atomic
									results[7].x[1] += suc_rate_ground2;
//#pragma omp atomic
									results[8].x[1] += suc_rate_ground2*suc_rate_ground2;
//#pragma omp atomic
									results[9].x[1] += suc_rate_uniq2;
//#pragma omp atomic
									results[10].x[1] += suc_rate_uniq2*suc_rate_uniq2;
									//Max Sols _ C(ii) condition			
//#pragma omp flush (results)
									if (results[0].x[1] < numSols2[i]) {
//#pragma omp critical
										{
											if (results[0].x[1] < numSols2[i])	results[0].x[1] = numSols2[i];
										}
									}
									//MIN sols _ C(ii) condition
//#pragma omp flush (results)
									if (results[1].x[1] > numSols2[i]) {
//#pragma omp critical
										{
											if (results[1].x[1] > numSols2[i])	results[1].x[1] = numSols2[i];
										}
									}
									//Max E _ C(ii) condition
//#pragma omp flush (results)
									if (results[5].x[1] < tab_results[2 * i+1][3]) {
//#pragma omp critical
										{
											if (results[5].x[1] <  tab_results[2 * i+1][3])	results[5].x[1] = tab_results[2 * i+1][3];
										}
									}
									//Max dx _ C(ii) condition

//#pragma omp flush (results)
									if (results[6].x[1] < tab_results[2 * i+1][2]) {
//#pragma omp critical
										{
											if (results[6].x[1] < tab_results[2 * i+1][2])	results[6].x[1] = tab_results[2 * i+1][2];
										}
									}

								}

								//OUTPUT:SUMMARY
								for (int i = 0; i < 2; i++) {
									results[2].x[i] /= (double)numConfigs;
									results[3].x[i] /= (double)numConfigs;
									results[4].x[i] = sqrt(fabs(results[4].x[i] / ((double)numConfigs) - results[3].x[i] * results[3].x[i]) / (double)numConfigs);
									results[7].x[i] /= (double)numConfigs;
									results[8].x[i] = sqrt(fabs(results[8].x[i] / ((double)numConfigs) - results[7].x[i] * results[7].x[i]) / (double)numConfigs);
									results[9].x[i] /= (double)numConfigs;
									results[10].x[i] = sqrt(fabs(results[10].x[i] / ((double)numConfigs) - results[9].x[i] * results[9].x[i]) / (double)numConfigs);
								}
								ofile << N << "\t" << M << "\t" << "real" << "\t" << results[0].x[0] << "\t" << results[1].x[0] << "\t" << results[2].x[0] << "\t" << results[3].x[0] << "\t" << results[4].x[0] << "\t" << results[5].x[0] << "\t" << results[6].x[0] << "\n";
								ofile << N << "\t" << M << "\t" << "imag" << "\t" << results[0].x[1] << "\t" << results[1].x[1] << "\t" << results[2].x[1] << "\t" << results[3].x[1] << "\t" << results[4].x[1] << "\t" << results[5].x[1] << "\t" << results[6].x[1] << std::endl;

								{
									GeometryVector temp;	temp.Dimension = 4;
									temp = GeometryVector(M, results[7].x[0], results[8].x[0]);
									summaryReal1.push_back(temp);
									temp = GeometryVector(M, results[9].x[0], results[10].x[0]);
									summaryReal2.push_back(temp);
									temp = GeometryVector(M, results[7].x[1], results[8].x[1]);
									summaryImag1.push_back(temp);
									temp = GeometryVector(M, results[9].x[1], results[10].x[1]);
									summaryImag2.push_back(temp);
								}

								//Energy Distributions: preparation
#pragma omp parallel for  shared(tab_energy,Energy0, dE, EnergyDist)
								for (int i = 0; i < numConfigs; i++) {//For each TARGET
									
								  //C(i) condition: REAL
									for (int l = 0; l < tab_energy[2 * i].size(); l++) {//For each TRIAL initial configs
										int n = 0;
										std::vector<double> list(dE.size());
										for (int m = 0; m < tab_energy[2 * i][l].size(); m++) {
											n = (int)floor(log(tab_energy[2 * i][l][m] / Energy0) / log(xdE));
											if (n >= 0 && n < list.size())
												list.at(n)++;
										}

										double num = (double)tab_energy[2 * i][l].size();
										if (num > 0.5) {
											for (int m = 0; m < list.size(); m++) { //probability distributions
												list.at(m) /= num;
											}
										}
										
#pragma omp atomic
										EnergyDist[l][0].x[3]++; // the number of target configurations yielding l-th solution.

										
										for (int m = 0; m < list.size(); m++) { //Average of distributions!
#pragma omp atomic
												EnergyDist[l][m].x[1] += list[m];
#pragma omp atomic
												EnergyDist[l][m].x[2] += list[m] * list[m];
										}
										
									}


									//C(ii) condition: IMAG
									if (M % 2 != 0) {//real \neq imag
										for (int l = 0; l < tab_energy[2 * i + 1].size(); l++) {//For each TRIAL initial configs
											int n = 0;
											std::vector<double> list(dE.size());
											for (int m = 0; m < tab_energy[2 * i+1][l].size(); m++) {
												n = (int)floor(log(tab_energy[2 * i+1][l][m] / Energy0) / log(xdE));
												if (n >= 0 && n < list.size())
													list.at(n)++;
											}
											double num = (double)tab_energy[2 * i + 1][l].size();
											if (num > 0.5) {
												for (int m = 0; m < list.size(); m++) { //probability distributions
													list.at(m) /= num;
												}
											}

#pragma omp atomic
											EnergyDist[3+l][0].x[3]++; // the number of target configurations yielding l-th solution.

											
											for (int m = 0; m < list.size(); m++) { //Average of distributions!
#pragma omp atomic
													EnergyDist[3+l][m].x[1] += list[m];
#pragma omp atomic
													EnergyDist[3+l][m].x[2] += list[m] * list[m];
											}
											
										}
									}
									
								}

								//OUPUT: Energy Dist (not multi-threaded)
								for (int i = 0; i < 6; i++) {
									//The number of i-th solutions for all target configurations
									//i=0	trivial solution
									//i=1	nontrivial
									//i=2	secondary nontrivial
									double num = EnergyDist[i][0].x[3]; 
									

									if (num > 0.5) {//Only when there are more than a single datum
										for (int j = 0; j < EnergyDist[i].size(); j++) {
											EnergyDist[i][j].x[2] = sqrt(fabs(EnergyDist[i][j].x[2] / num - EnergyDist[i][j].x[1] / num * EnergyDist[i][j].x[1] / num) / num);
											EnergyDist[i][j].x[1] /= num; //*dE[j];
																		  //EnergyDist[i][j].x[2] /= num; //*dE[j];
										}

										if (i < 3) {
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Edist_real_%d", dir, dim, systName, method, N, (int)M, i);
										}
										else{
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Edist_imag_%d", dir, dim, systName, method, N, (int)M, i - 3);
										}
										EnergyDist[i][0].x[3] = 0;

										if (i >= 3 && M % 2 == 0) {//When C(i) and C(ii) are the same 
											WriteFunction(EnergyDist[i-3], filename);
										}
										else {
											WriteFunction(EnergyDist[i], filename);
										}
									}
								}

								tab_energy.clear();


								//Proximity Distributions
								if (doProxMetric) {
#pragma omp parallel for shared(ProximityDist, av_X2, tab_proxy)
									for (int i = 0; i < numConfigs; i++) {
										//original data are processed s.t. the size of proximity is not larger than 4 (total, trivial, and leading nontrivial, and nontrivials).
										
										//REAL; C(I) condition
										for (int l = 0; l < tab_proxy[2*i].size(); l++) {// l th solution!
											int n = 0;
											std::vector<double> list(av_X2.size()); //histogram
											for (int m = 0; m < tab_proxy[2*i][l].size(); m++) {
												n = (int)floor(log(tab_proxy[2*i][l][m] / x0) / log(d_x2));
												if (n >= 0 && n < list.size())
													list.at(n)++;
											}
											double deno = (double)tab_proxy[2 * i][l].size();
											if (deno > 0.5) {
												for (int m = 0; m < list.size(); m++) 											
													list.at(m) /= deno;//*dE[m];
											}

#pragma omp atomic
											ProximityDist[l][0].x[3]++;	//The occurence of l-th solution
											//double num_inv = 1.0 / (double)ProximityDist[l][0].x[3];
											//	std::cout << "num = " << num << ",\t";
											for (int m = 0; m < list.size(); m++) {
#pragma omp atomic
												ProximityDist[l][m].x[1] += list[m];
#pragma omp atomic
												ProximityDist[l][m].x[2] += list[m] * list[m];
												//ProximityDist[l][m].x[1] = (1.0 - num_inv)*ProximityDist[l][m].x[1] + num_inv * list[m];
												//ProximityDist[l][m].x[2] = (1.0 - num_inv)*ProximityDist[l][m].x[2] + num_inv * list[m] * list[m];
											}
										}
									
										//IMAG; C(ii) condition
										if (M % 2 != 0) {//a.k.a. C(i) and C(ii) are different
											for (int l = 0; l < tab_proxy[2 * i+1].size(); l++) {
												int n = 0;
												std::vector<double> list(av_X2.size()); //histogram
												for (int m = 0; m < tab_proxy[2 * i+1][l].size(); m++) {
													n = (int)floor(log(tab_proxy[2 * i+1][l][m] / x0) / log(d_x2));
													if (n >= 0 && n < list.size())
														list.at(n)++;
												}
												double deno = (double)tab_proxy[2 * i + 1][l].size();
												if (deno > 0.5) {
													for (int m = 0; m < list.size(); m++)
														list.at(m) /= deno;//*dE[m];  : list = probability
													
												}

#pragma omp atomic
												ProximityDist[l+ numLegd][0].x[3]++; 	//The occurence of l-th solution

												for (int m = 0; m < list.size(); m++) {
#pragma omp atomic
													ProximityDist[l + numLegd][m].x[1] += list[m];
#pragma omp atomic
													ProximityDist[l + numLegd][m].x[2] += list[m] * list[m];
												}
											}
										}

									}
								
									//OUTPUT: Proximity Measure

									for (int i = 0; i < 8; i++) {
										double num = ProximityDist[i][0].x[3];	//The occurence of ith solution
										if (num > 0) {
											for (int j = 0; j < ProximityDist[i].size(); j++) {
												ProximityDist[i][j].x[1] /= num;
												ProximityDist[i][j].x[2] = sqrt(fabs(ProximityDist[i][j].x[2]/num - ProximityDist[i][j].x[1] * ProximityDist[i][j].x[1]) / num);
												//ProximityDist[i][j].x[1] /= (1.0 - 1.0 / d_x2)*av_X2[j];
												//ProximityDist[i][j].x[2] /= (1.0 - 1.0 / d_x2)*av_X2[j];
											}


											if (i == 0) {
												sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_real_total", dir, dim, systName, method, N, (int)M, i);
											}
											else if (i < 4) {
												sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_real_%d", dir, dim, systName, method, N, (int)M, i - 1);
											}
											else if (i == 4) {
												sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_imag_total", dir, dim, systName, method, N, (int)M, i);
											}
											else {
												sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_imag_%d", dir, dim, systName, method, N, (int)M, i - 5);
											}
											ProximityDist[i][0].x[3] = 0;
											WriteFunction(ProximityDist[i], filename);
										}
									}
								}

								
								
								tab_proxy.clear();

		

						}
								



//#pragma omp critical
//								{
								//REAL CONSTRAINTS
//									double suc_rate_ground_real = (double)results_temp1[4] / Nc;
//									double suc_rate_uniq_real = (results_temp1[4]>0) ? (double)results_temp1[5] / results_temp1[4] : 0;
//									
//									results[0].x[0] = std::max(results[0].x[0], (double)sols1.size());
//									results[1].x[0] = std::min(results[1].x[0], (double)sols1.size());
//									results[5].x[0] = std::max(results[5].x[0], results_temp1[3]);
//									results[6].x[0] = std::max(results[6].x[0], results_temp1[2]);
//									
//
//									results[2].x[0] += (1.0 - suc_rate_uniq);
//									results[3].x[0] += (double)sols1.size();
//									results[4].x[0] += (double)sols1.size()*sols1.size();
//									results[7].x[0] += suc_rate_ground;
//									results[8].x[0] += suc_rate_ground*suc_rate_ground;
//									results[9].x[0] += suc_rate_uniq;
//									results[10].x[0] += suc_rate_uniq*suc_rate_uniq;
//
//									//IMAGE CONSTRAINTS
//									suc_rate_ground = (double)results_temp2[4] / Nc;
//									suc_rate_uniq = (results_temp2[4]>0) ? (double)results_temp2[5] / results_temp2[4] : 0;
//									results[0].x[1] = std::max(results[0].x[1], (double)sols2.size());
//									results[1].x[1] = std::min(results[1].x[1], (double)sols2.size());
//									results[5].x[1] = std::max(results[5].x[1], results_temp2[3]);
//									results[6].x[1] = std::max(results[6].x[1], results_temp2[2]);
//
//									
//									results[2].x[1] += (1.0 - suc_rate_uniq);
//									results[3].x[1] += (double)sols2.size();
//									results[4].x[1] += (double)sols2.size()*sols2.size();
//									results[7].x[1] += suc_rate_ground;
//									results[8].x[1] += suc_rate_ground*suc_rate_ground;
//									results[9].x[1] += suc_rate_uniq;
//									results[10].x[1] += suc_rate_uniq*suc_rate_uniq;
//
//									//if (M % 2 == 0) {//.x[0] = .x[1]
//									//	for (int l = 0; l < results.size(); l++)
//									//		results[l].x[1] = results[l].x[0];
//									//}
//									//else { //odd number of constraints
//									//	
//									//}
//
//									//output solutions---------------------------
//									if (M >= N && outConfig) {
//										outfile << count << "\tTarget\t";
//										{//target configs
//											Configuration tempC(sols1[0]);
//											Merge_Sort(tempC, 0);
//											for (int z = 0; z < tempC.NumParticle(); z++)
//												outfile << std::setprecision(7) << tempC.GetCartesianCoordinates(z).x[0] << ",\t";
//
//											outfile << "\n";
//										}
//
//										if (sols1.size() == 1) {//unique solution in real constraints
//											outfile << count << "\tReal\tunique";
//											for (int m = 0; m < N; m++)
//												outfile << "\t\t";
//
//											outfile << results_temp1[6] << "\t" << results_temp1[7] << "\t" << results_temp1[8] << std::endl;
//										}
//										else {//multiple solutions in real constraints
//											for (int t = 0; t < sols1.size(); t++) {
//												outfile << count << "\tReal\t";
//												Configuration tempC(sols1[t]);
//												Merge_Sort(tempC, 0);
//												for (int y = 0; y < tempC.NumParticle(); y++)
//													outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";
//
//												outfile << results_temp1[3 * t + 6] << "\t" << results_temp1[3 * t + 7] << "\t" << results_temp1[3 * t + 8] << std::endl;
//												//if (t == 1) {
//												//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
//												//}
//											}
//										}
//										if (M % 2 == 0) {//real = image
//											if (sols1.size() == 1) {//unique solution in real constraints
//												outfile << count << "\tImag\tunique";
//												for (int m = 0; m < N; m++)
//													outfile << "\t\t";
//
//												outfile << results_temp1[6] << "\t" << results_temp1[7] << "\t" << results_temp1[8] << std::endl;
//											}
//											else {//multiple solutions in real constraints
//												for (int t = 0; t < sols1.size(); t++) {
//													outfile << count << "\tImag\t";
//													Configuration tempC(sols1[t]);
//													Merge_Sort(tempC, 0);
//													for (int y = 0; y < tempC.NumParticle(); y++)
//														outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";
//
//													outfile << results_temp1[3 * t + 6] << "\t" << results_temp1[3 * t + 7] << "\t" << results_temp1[3 * t + 8] << std::endl;
//													//if (t == 1) {
//													//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
//													//}
//												}
//											}
//										}
//										else {
//											if (sols2.size() == 1) {//unique solution in real constraints
//												outfile << count << "\tImag\tunique";
//												for (int m = 0; m < N; m++)
//													outfile << "\t\t";
//
//												outfile << results_temp2[6] << "\t" << results_temp2[7] << "\t" << results_temp2[8] << std::endl;
//											}
//											else {//multiple solutions in real constraints
//												for (int t = 0; t < sols2.size(); t++) {
//													outfile << count << "\tImag\t";
//													Configuration tempC(sols2[t]);
//													Merge_Sort(tempC, 0);
//													for (int y = 0; y < tempC.NumParticle(); y++)
//														outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";
//
//													outfile << results_temp2[3 * t + 6] << "\t" << results_temp2[3 * t + 7] << "\t" << results_temp2[3 * t + 8] << std::endl;
//													//if (t == 1) {
//													//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
//													//}
//												}
//											}
//										}
//									}
//									count++;
//
//									{//------------Energy distribution------------
//										for (int l = 0; l < energy1.size(); l++) {
//											int n = 0;
//											std::vector<double> list(dE.size());
//											for (int m = 0; m < energy1[l].size(); m++) {
//												n = (int)floor(log(energy1[l][m] / Energy0) / log(xdE));
//												if (n >= 0 && n < list.size())
//													list.at(n)++;
//											}
//
//											for (int m = 0; m < list.size(); m++) {
//												list.at(m) /= (double)Nc;//*dE[m];
//											}
//
//											double num = ++EnergyDist[l][0].x[3];
//											//	std::cout << "num = " << num << ",\t";
//											for (int m = 0; m < list.size(); m++) {
//												EnergyDist[l][m].x[1] = (1.0 - 1.0 / num)*EnergyDist[l][m].x[1] + (1.0 / num)*list[m];
//												EnergyDist[l][m].x[2] = (1.0 - 1.0 / num)*EnergyDist[l][m].x[2] + (1.0 / num) *list[m] * list[m];
//											}
//
//										}
//										if (M % 2 != 0) {//real \neq imag
//											for (int l = 0; l < energy2.size(); l++) {
//												int n = 0;
//												std::vector<double> list(dE.size());
//												for (int m = 0; m < energy2[l].size(); m++) {
//													n = (int)floor(log(energy2[l][m] / Energy0) / log(xdE));
//													if (n >= 0 && n < list.size())
//														list.at(n)++;
//												}
//
//												for (int m = 0; m < list.size(); m++) {
//													list.at(m) /= (double)Nc;//*dE[m];
//												}
//
//												double num = ++EnergyDist[3 + l][0].x[3];
//												//	std::cout << "num = " << num << ",\t";
//												for (int m = 0; m < list.size(); m++) {
//													EnergyDist[3 + l][m].x[1] = (1.0 - 1.0 / num)*EnergyDist[3 + l][m].x[1] + (1.0 / num) * list[m];
//													EnergyDist[3 + l][m].x[2] = (1.0 - 1.0 / num)*EnergyDist[3 + l][m].x[2] + (1.0 / num) * list[m] * list[m];
//												}
//
//											}
//										}
//									}
//
//									if (doProxMetric) {//calculate proximity metrics
//										{
//											int numSols = 4;
//											if (proximity1.size()>numSols) {
//												//if there are more than 2 nontrivial solutions
//												//merge ternary solutions
//												for (int l = proximity1.size() - 1; l >= numSols; l--) {
//													for (int k = 0; k < proximity1[l].size(); k++) {
//														proximity1[numSols-1].push_back(proximity1[l].at(k));
//													}
//												}
//												proximity1.resize(numSols);
//											}
//
//											//therefore, the size of proximity1 is now not larger than 4 (total, trivial, and leading nontrivial, and nontrivials).
//											for (int l = 0; l < proximity1.size(); l++) {
//												int n = 0;
//												std::vector<double> list(av_X2.size()); //histogram
//												for (int m = 0; m < proximity1[l].size(); m++) {
//													n = (int)floor(log(proximity1[l][m] / x0) / log(d_x2));
//													if (n >= 0 && n < list.size())
//														list.at(n)++;
//												}
//
//												for (int m = 0; m < list.size(); m++) {
//													list.at(m) /= (double)proximity1[l].size();//*dE[m];
//												}
//
//												ProximityDist[l][0].x[3]++;	double num_inv = 1.0 / (double)ProximityDist[l][0].x[3];
//												//	std::cout << "num = " << num << ",\t";
//												for (int m = 0; m < list.size(); m++) {
//													ProximityDist[l][m].x[1] = (1.0 - num_inv)*ProximityDist[l][m].x[1] + num_inv * list[m];
//													ProximityDist[l][m].x[2] = (1.0 - num_inv)*ProximityDist[l][m].x[2] + num_inv * list[m] * list[m];
//												}
//											}
//										}
//
//										if (M % 2 != 0) {//when proximity1 is different from another.
//											int numSols = 4;
//											if (proximity2.size()>numSols) {
//												//if there are more than 2 nontrivial solutions
//												//merge ternary solutions
//												for (int l = proximity2.size() - 1; l >= numSols; l--) {
//													for (int k = 0; k < proximity2[l].size(); k++) {
//														proximity2[numSols - 1].push_back(proximity2[l].at(k));
//													}
//												}
//												proximity2.resize(numSols);
//											}
//										}
//										else {
//											proximity2 = proximity1;
//
//										}
//	
//										//therefore, the size of proximity2 is now not larger than 4 (total, trivial, and leading nontrivial, and nontrivials).
//										for (int l = 0; l < proximity2.size(); l++) {
//											int numSols = 4;
//											int n = 0;
//											std::vector<double> list(av_X2.size());
//											for (int m = 0; m < proximity2[l].size(); m++) {
//												n = (int)floor(log(proximity2[l][m] / x0) / log(d_x2));
//												if (n >= 0 && n < list.size())
//													list.at(n)++;
//											}
//
//											for (int m = 0; m < list.size(); m++) {
//												list.at(m) /= (double)proximity2[l].size();//*dE[m];
//											}
//
//											ProximityDist[l+numSols][0].x[3]++;	double num_inv = 1.0 / (double)ProximityDist[l+numSols][0].x[3];
//											//	std::cout << "num = " << num << ",\t";
//											for (int m = 0; m < list.size(); m++) {
//												ProximityDist[l + numSols][m].x[1] = (1.0 - num_inv)*ProximityDist[l + numSols][m].x[1] + num_inv * list[m];
//												ProximityDist[l + numSols][m].x[2] = (1.0 - num_inv)*ProximityDist[l + numSols][m].x[2] + num_inv * list[m] * list[m];
//											}
//										}
//
//									}
//								}
//
//							}

											
							

						}

						else if (dim == 2) {
							//on construction....
							ofile << "not yet...\n";
						}


						sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "real");
						WriteFunction(summaryReal1, filename);
						sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "real");
						WriteFunction(summaryReal2, filename);

						sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "imag");
						WriteFunction(summaryImag1, filename);
						sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "imag");
						WriteFunction(summaryImag2, filename);
					}
					

					

				}


			}

		}

}





void CCInversionProblem1(bool outConfig) {
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char dir0[100];// = "results/Reconst/";
	sprintf(dir0, "%s/results/Reconst/", masterDir.c_str());
	char dir[100] = {};
	char temp[200] = {};
	char name[200] = {};

	//---- general parameters
	DimensionType dim = 1;
	int Nc = 100; //the number of trial initial configs
	int numConfigs = 1000; //the number of target configs
	int NKSkip = 0;
	int numThreads = 0;
	char constType[100] = {};

	//---- swap parameters
	int N1 = 10;//min system size
	int N2 = 40;//max system size
	int systType1;
	int systType2;
	int dN = 1;
	int methodNum1 = 0;
	int methodNum2 = 0;
	int M1 = -1;
	int M2 = +1;
	double EnergyTolerance = 1e-20;
	double PositionTolerance = 1e-6;
	double delta1;
	double delta2;
	double d_delta;
	bool doProxMetric;
	double proximityUB = PositionTolerance * PositionTolerance;
	{
		ofile << "NKSkip=";
		ifile >> NKSkip;
		ofile << "dim=";
		ifile >> dim;
		ofile << "min N=";
		ifile >> N1;
		ofile << "max N=";
		ifile >> N2;
		ofile << "dN=";
		ifile >> dN;
		ofile << "numConfigs=";
		ifile >> numConfigs;
		ofile << "Nc=";
		ifile >> Nc;
		ofile << "system type1:";
		ifile >> systType1;
		ofile << "system type2:";
		ifile >> systType2;
		ofile << "initial method num: ";
		ifile >> methodNum1;
		ofile << "final method num: ";
		ifile >> methodNum2;
		ofile << "Energ Tolerance = ";
		ifile >> EnergyTolerance;
		ofile << "delta1 = ";
		ifile >> delta1;
		ofile << "delta2 = ";
		ifile >> delta2;
		ofile << "d_delta = ";
		ifile >> d_delta;
		ofile << "M1 = ";
		ifile >> M1;
		ofile << "M2 = ";
		ifile >> M2;
		ofile << "num threads = ";
		ifile >> numThreads;
		ofile << "Do we compute the proximity metric?  ";
		ifile >> doProxMetric;
	}
	{//SHOW PARAMETERS
		ofile << "NKSkip=";
		ofile << NKSkip;
		ofile << "\ndim=";
		ofile << dim;
		ofile << "\nmin N=";
		ofile << N1;
		ofile << "\nmax N=";
		ofile << N2;
		ofile << "\ndN=";
		ofile << dN;
		ofile << "\nnumConfigs=";
		ofile << numConfigs;
		ofile << "\nNc=";
		ofile << Nc;
		ofile << "\nsystem type1:";
		ofile << systType1;
		ofile << "\nsystem type2:";
		ofile << systType2;
		ofile << "\ninitial method num: ";
		ofile << methodNum1;
		ofile << "\nfinal method num: ";
		ofile << methodNum2;
		ofile << "\nEnerg Tolerance = ";
		ofile << EnergyTolerance;
		ofile << "\ndelta1 = ";
		ofile << delta1;
		ofile << "\ndelta2 = ";
		ofile << delta2;
		ofile << "\nd_delta = ";
		ofile << d_delta;
		ofile << "\nM1 = ";
		ofile << M1;
		ofile << "\nM2 = ";
		ofile << M2;
		ofile << "\nnum threads = ";
		ofile << numThreads;
		ofile << "\nDo we compute the proximity metric?  ";
		ofile << doProxMetric << std::endl;

	}

	if (NKSkip > 0)
		std::cout << "Skip the central " << NKSkip << " k points.\n";

	sprintf(dir, "%sdE_%.0e/Nk_skip_%d/", dir0, EnergyTolerance, NKSkip);
	ofile << dir << "\n";
	std::vector<GeometryVector> basis;
	for (int j = 0; j < dim; j++) {
		basis.push_back(GeometryVector());
		basis[j].SetDimension(dim);
		basis[j].x[j] = 1.0;
	}
	//----energy distribution----
	std::vector<double> bin;
	std::vector<double> dE;
	std::vector<double> av_X2; //for proximity metrics
	double Energy0 = 1e-34;
	double xdE = 2;
	double Energy = Energy0; //start energy
	while (Energy < EnergyTolerance) {
		bin.push_back(Energy);
		Energy *= xdE;
		if (bin.size() >= 2) {
			dE.push_back(bin[bin.size() - 1] - bin[bin.size() - 2]);
		}
	}

	double x0 = 1e-38;
	double d_x2 = 2;
	if (doProxMetric) {
		double x2 = x0;
		while (x2 < 10.0) {
			x2 *= d_x2;
			av_X2.push_back(x2);
		}
	}


	for (double delta = delta1; delta <= delta2; delta += d_delta) {
		for (int type = systType1; type <= systType2; type++) {
			char systName[200] = {};
			GetTarget(0, type, 1, basis, delta, systName); //just get system name;;;
			ofile << "system type = " << systName << std::endl;
			ofile << "N\tM\tType\tmaxSols\tminSols\tProb_DiffSols\tav.NumSols\tstdSols\tMaxE_ground\tMaxdx\tMethods\n";

			for (int N = N1; N <= N2; N += dN) {
				
				std::vector<Configuration> targetList;	targetList.reserve(numConfigs);
				for (int i = 0; i < numConfigs; i++)
					targetList.emplace_back( GetTarget(i, type, N, basis, delta, nullptr));


				auto getTargetC = [&targetList](size_t i)->Configuration {
					//RandomGenerator g(std::time(NULL)+i);
					return targetList[i];
				};

				for (int method_num = methodNum1; method_num <= methodNum2; method_num++) {
					char method[100] = {};
					switch (method_num)
					{
					case 0:
						sprintf(method, "BFGS+MINOP"); break;
					case 1:
						sprintf(method, "Steepest Descent"); break;
					case 2:
						sprintf(method, "ConjugateGradient"); break;
					case 3:
						sprintf(method, "LocalGradientDescent"); break;
					default:
						break;
					}

					char filename[400] = {};
					sprintf(filename, "%s%dD%s__%s__N_%d.txt", dir, dim, systName, method, N);
					std::fstream outfile(filename, std::fstream::out);//
					std::fstream outfile1;
			//ofile << "method = " << method << std::endl;

					std::vector<GeometryVector> summaryReal1;
					std::vector<GeometryVector> summaryReal2;
					std::vector<GeometryVector> summaryImag1;
					std::vector<GeometryVector> summaryImag2;

					if (dim == 1) {
						for (size_t M = N + M1; M <= N + M2; M++) {
							outfile << "-----------------------------------------------------------------\n";
							outfile << "Testing N=" << N << "\tM=" << M << "\n\n";
							outfile << "num\tconstr\tconfigs";
							for (int z = 0; z < N; z++)
								outfile << "\t";
							outfile << "\tminE\t\tmaxE\t\tmaxDx\tProbability" << std::endl;


							//-----------------------final results!!------------------------------------------------
							//results[i].x[0] <= real type / .x[1] <= imag type
							//[0].				maxSols
							//[1].				minSols
							//[2].				probability to gain nontrivial solutions
							//[3].				average number of solutions
							//[4].				std of sols
							//[5].				max energy of ground states
							//[6].				max position errors
							//[7].				success rate: ground states
							//[8].				error_success rate: ground states
							//[9].				success rate: correct solutions
							//[10].				error_success rate: correct solutions
							std::vector<GeometryVector> results;
							{
								GeometryVector dump(static_cast<int> (2));
								results.assign(11, dump);
								results[1].x[0] = Nc + 1;
								results[1].x[1] = Nc + 1;
							}
							int count = 0;
							//energy distributions...
							std::vector<std::vector<GeometryVector>> EnergyDist;
							std::vector<std::vector<GeometryVector>> ProximityDist; //measures of proximity metrics

							{
								//n = 0: real, n=1: imag
								//EnergyDist[n+0]: trivial sols
								//EnergyDist[n+1]: nontrivial sols
								//EnergyDist[n+2]: others

								//ProximityMetric[n+0]: trivial sols
								//ProximityMetric[n+1]: nontrivial sols
								//ProximityMetric[n+2]: others
								//ProximityMetric[n+3]: for all solutions w.r.t. the targets

								std::vector<GeometryVector> table;
								GeometryVector element(0, 0, 0, 0);	element.Dimension = 4;
								for (int i = 1; i < bin.size(); i++) {
									element.x[0] = 0.5*(bin[i] + bin[i - 1]);
									table.push_back(element);
								}
								for (int i = 0; i < 6; i++) {
									EnergyDist.push_back(table);
								}

								{
									table.resize(0);
									for (int i = 0; i < av_X2.size(); i++) {
										element.x[0] = (1.0 + 1.0 / d_x2)*0.5*av_X2[i];
										table.push_back(element);
									}
									for (int i = 0; i < 8; i++) {
										ProximityDist.push_back(table);
									}
								}

							}


							//temporary storage! version 1.1  buffers!!
							std::vector<std::vector<double>> tab_results; //store all results here!!
							std::vector<std::vector<std::vector<double>>> tab_energy;
							std::vector<std::vector<std::vector<double>>> tab_proxy;
							std::vector<int> numSols1;
							std::vector<int> numSols2;
							int numLegd = 4;
							{//fill storage so that the size of containers doesn't change during the for-loops
							 //------------for i-th target configuration

							 //tab_results[2*i] =	results with "real" constraints
							 //tab_results[2*i+1] =	results with "imag" constraints

							 //tab_energy[2*i] =		energy distribution of "real" constraints
							 //tab_energy[2*i+1] =	energy distribution of "imag" constraints

							 //tab_proximity[2*i] =	proximity distribution of "real" constraints
							 //tab_proximity[2*i+1] = proximity distribution of "imag" constraints
								tab_results.resize(0);
								tab_energy.resize(0);
								tab_proxy.resize(0);
								numSols1.resize(0);
								numSols2.resize(0);

								std::vector<double> temp_result;
								std::vector<std::vector<double>> temp_box;

								tab_results.resize(2 * numConfigs, temp_result);
								tab_energy.resize(2 * numConfigs, temp_box);
								tab_proxy.resize(2 * numConfigs, temp_box);
								numSols1.resize(numConfigs, 0);
								numSols2.resize(numConfigs, 0);
							}




#pragma omp parallel for schedule (guided) num_threads(numThreads) shared( tab_results, tab_energy, tab_proxy, numSols1, numSols2, targetList)//shared(results, EnergyDist, ProximityDist) //shared(numSols, numDiffNumSols, maxSols, minSols, count, maxE_ground1, maxE_ground2, maxDx1, maxDx2)
							for (int i = 0; i < numConfigs; i++) {

								Configuration c = getTargetC(i); //TARGET
								std::vector<Configuration> initConfig;
								//initial configurations
								{
									RandomGenerator g(i + std::time(NULL));
									Configuration templateC(c);	templateC.RemoveParticles();
									for (int j = 0; j < Nc; j++) {
										Configuration t(templateC);
										for (int k = 0; k < N; k++) {
											t.Insert("a", g);
										}
										initConfig.push_back(t);
									}
								}

								//results_temp1[0] = N
								//results_temp1[1] = M
								//results_temp1[2] = maxDr
								//results_temp1[3] = maxE
								//results_temp1[4] = # ground states
								//results_temp1[5] = # correct solutions
								//results_temp1[6+3*i+0] = min energy: sol[i]
								//results_temp1[6+3*i+1] = max energy: sol[i]
								//results_temp1[6+3*i+2] = max dx: sol[i]
								std::vector<double> results_temp1;			//MISC INFO REAL CONSTRAINT
								std::vector<double> results_temp2;			//MISC INFO IMAG CONSTRAINT
								std::vector<Configuration> sols1;			//SOLUTION CONFIG REAL CONSTRAINT
								std::vector<Configuration> sols2;			//SOLUTION CONFIG IMAG CONSTRAINT
								std::vector<std::vector<double>> energy1;	//energy1[i] = list of energy of this solution
								std::vector<std::vector<double>> energy2;	//energy2[i] = list of energy of this solution

								std::vector<std::vector<double>> proximity1;	//proximity1[i] = list of proximity metrics of solution i
								std::vector<std::vector<double>> proximity2;	//proximity2[i] = list of proximity metrics of solution i
								if (doProxMetric) {//just give signal to initiate calculations
									std::vector<double> dummy;
									proximity1.push_back(dummy);
									proximity2.push_back(dummy);
								}
								//-------GENERATE INFO

								FrankInvCCUniqueTest3(c, Nc, M, results_temp1, sols1, outfile1, PositionTolerance, EnergyTolerance, NKSkip, "real", method_num, energy1, initConfig, proximity1);
								for (int j = 0; j < energy1.size(); j++)
									MergeSort(energy1[j]);


								if (doProxMetric && proximity1.size()>numLegd) {
									//if there are more than 2 nontrivial solutions
									//merge ternary solutions
									for (int l = proximity1.size() - 1; l >= numLegd; l--) {
										for (int k = 0; k < proximity1[l].size(); k++) {
											proximity1[numLegd - 1].push_back(proximity1[l].at(k));
										}
									}
									proximity1.resize(numLegd);
									for (int j = 0; j < proximity1.size(); j++) {
										MergeSort(proximity1[j]);
									}
								}

								if (M % 2 == 0) {//number of constraints = even => real and imaginary constraints are identical
									results_temp2 = results_temp1;
									sols2 = sols1;
									energy2 = energy1;
									proximity2 = proximity1;
								}
								else {//number of constraints = odd
									FrankInvCCUniqueTest3(c, Nc, M, results_temp2, sols2, outfile1, PositionTolerance, EnergyTolerance, NKSkip, "image", method_num, energy2, initConfig, proximity2);
									for (int j = 0; j < energy2.size(); j++)
										MergeSort(energy2[j]);

									if (doProxMetric && proximity2.size()>numLegd) {
										//if there are more than 2 nontrivial solutions
										//merge ternary solutions
										for (int l = proximity2.size() - 1; l >= numLegd; l--) {
											for (int k = 0; k < proximity2[l].size(); k++) {
												proximity2[numLegd - 1].push_back(proximity2[l].at(k));
											}
										}
										proximity2.resize(numLegd);
										for (int j = 0; j < proximity2.size(); j++) {
											MergeSort(proximity2[j]);
										}
									}
								}
								//--------GATHER INFO
								//SUCCESS RATE: GROUND STATES
								//SUCCESS RATE: UNIQUE SOLUTION
								//MAX SOLS/TARGET
								//MIN SOLS/TARGET
								//MAX ENERGY
								//MIN ENERGY
								//ENERGY DISTRIBUTION: TRIVIAL/SECONDARY/OTHERS....




#pragma omp critical
								{
									numSols1[i] = sols1.size();
									numSols2[i] = sols2.size();
								}


								//#pragma omp atomic
								//								numSols1[i] = sols1.size();
								//#pragma omp atomic
								//								numSols2[i] = sols2.size();

								//organize output!!
								if (M >= N && outConfig) {
									//TARGET
									Merge_Sort(sols1[0], 0);
									//SOLUTIONS: REAL CONDITION
									//if (sols1.size() > 1) {//multiple solutions in real constraints
									//	for (int t = 1; t < sols1.size(); t++) {
									//		Merge_Sort(sols1[t], 0);
									//	}
									//}

									//SOLUTIONS: IMAG CONDITION
									sols2[0] = sols1[0];
									if (sols2.size() > 1) {//multiple solutions in real constraints
										if (M % 2 == 0) {//when REAL and IMAG are identical
											sols2 = sols1;
										}
										//else {//when REAL and IMAG are distinct
										//	for (int t = 1; t < sols2.size(); t++) {
										//		Merge_Sort(sols2[t], 0);
										//	}
										//}

									}

#pragma omp critical
									{//OUTPUT:SOLS
										count = i;
										outfile << count << "\tTarget\t";
										{//TARGET

											for (int z = 0; z < sols1[0].NumParticle(); z++)
												outfile << std::setprecision(7) << sols1[0].GetCartesianCoordinates(z).x[0] << ",\t";

											outfile << "\n";
										}

										//REAL CONDITION
										if (sols1.size() == 1) {//UNIQUE SOLUTION
											outfile << count << "\tReal\t";
											if (results_temp1[5] > 0) { outfile << "unique"; }
											else { outfile << "fail"; }
											for (int m = 0; m < N; m++)
												outfile << "\t\t";

											outfile << results_temp1[6] << "\t" << results_temp1[7] << "\t" << results_temp1[8] << std::endl;
										}
										else {//MULTIPLE SOLUTIONS
											for (int t = 0; t < sols1.size(); t++) {
												outfile << count << "\tReal\t";
												for (int y = 0; y < sols1[t].NumParticle(); y++)
													outfile << sols1[t].GetCartesianCoordinates(y).x[0] << ",\t";

												outfile << results_temp1[3 * t + 6] << "\t" << results_temp1[3 * t + 7] << "\t" << results_temp1[3 * t + 8] << "\t";
												if (t < 2) {
													double prob = (results_temp1[4] > 0.5) ? (double)energy1[t].size() / results_temp1[4] : 0;
													outfile  << prob;
												}

													outfile << std::endl;

											}
										}

										//IMAG
										if (sols2.size() == 1) {//UNIQUE SOL
											outfile << count << "\tImag\t";
											if (results_temp2[5] > 0) { outfile << "unique"; }
											else { outfile << "fail"; }
											for (int m = 0; m < N; m++)
												outfile << "\t\t";

											outfile << results_temp2[6] << "\t" << results_temp2[7] << "\t" << results_temp2[8] << std::endl;
										}
										else {//multiple solutions in real constraints
											for (int t = 0; t < sols2.size(); t++) {
												outfile << count << "\tImag\t";
												for (int y = 0; y < sols2[t].NumParticle(); y++)
													outfile << sols2[t].GetCartesianCoordinates(y).x[0] << ",\t";

												outfile << results_temp2[3 * t + 6] << "\t" << results_temp2[3 * t + 7] << "\t" << results_temp2[3 * t + 8] << "\t"; 
												if (t < 2) {
													double prob = (results_temp1[4] > 0.5) ? (double)energy1[t].size() / results_temp1[4] : 0;
													outfile << prob;
												}

												outfile << std::endl;
											}
										}
										//count++;
									}




								}

#pragma omp critical
								{
									tab_results[2 * i].swap(results_temp1);
									tab_results[2 * i + 1].swap(results_temp2);

									tab_energy[2 * i].swap(energy1);
									tab_energy[2 * i + 1].swap(energy2);

									tab_proxy[2 * i].swap(proximity1);
									tab_proxy[2 * i + 1].swap(proximity2);
								}

							}


							//SUMMARY!!
							//#pragma omp parallel
							//							{
							//							int priv_max_sols = 0;
							//							int priv_min_sols = numConfigs + 1;
							//							double rate_diff_sols = 0;
							//							double av_numSols = 0;
							//							double sq_numSols = 0;
							//							double rate
							//							}
							//#pragma omp parallel for shared(tab_results, results, EnergyDist, numSols1, numSols2)
							for (int i = 0; i < numConfigs; i++) {
								//REAL CONDITION
								double suc_rate_ground1 = tab_results[2 * i][4] / (double)Nc;
								double suc_rate_uniq1 = (tab_results[2 * i][4] > 0.5) ? tab_results[2 * i][5] / tab_results[2 * i][4] : 0.0;

								//#pragma omp atomic
								results[2].x[0] += (1.0 - suc_rate_uniq1);
								//#pragma omp atomic
								results[3].x[0] += (double)numSols1[i];
								//#pragma omp atomic
								results[4].x[0] += (double)numSols1[i] * numSols1[i];
								//#pragma omp atomic
								results[7].x[0] += suc_rate_ground1;
								//#pragma omp atomic
								results[8].x[0] += suc_rate_ground1 * suc_rate_ground1;
								//#pragma omp atomic
								results[9].x[0] += suc_rate_uniq1;
								//#pragma omp atomic
								results[10].x[0] += suc_rate_uniq1 * suc_rate_uniq1;
								//Max Sols _ C(i) condition			
								//#pragma omp flush (results)
								if (results[0].x[0] < numSols1[i]) {
									//#pragma omp critical
									{
										if (results[0].x[0] < numSols1[i])
											results[0].x[0] = numSols1[i];
									}
								}
								//MIN sols _ C(i) condition
								//#pragma omp flush (results)
								if (results[1].x[0] > numSols1[i]) {
									//#pragma omp critical
									{
										if (results[1].x[0] > numSols1[i])
											results[1].x[0] = numSols1[i];
									}
								}
								//Max E _ C(i) condition
								//#pragma omp flush (results)
								if (results[5].x[0] < tab_results[2 * i][3]) {
									//#pragma omp critical
									{
										if (results[5].x[0] <  tab_results[2 * i][3])	results[5].x[0] = tab_results[2 * i][3];
									}
								}
								//Max dx _ C(i) condition

								//#pragma omp flush (results)
								if (results[6].x[0] < tab_results[2 * i][2]) {
									//#pragma omp critical
									{
										if (results[6].x[0] < tab_results[2 * i][2])	results[6].x[0] = tab_results[2 * i][2];
									}
								}

								//IMAG CONDITION
								double suc_rate_ground2 = tab_results[2 * i + 1][4] / (double)Nc;
								double suc_rate_uniq2 = (tab_results[2 * i + 1][4] > 0.5) ? tab_results[2 * i + 1][5] / tab_results[2 * i + 1][4] : 0.0;

								//#pragma omp atomic
								results[2].x[1] += (1.0 - suc_rate_uniq2);
								//#pragma omp atomic
								results[3].x[1] += (double)numSols2[i];
								//#pragma omp atomic
								results[4].x[1] += (double)numSols2[i] * numSols2[i];
								//#pragma omp atomic
								results[7].x[1] += suc_rate_ground2;
								//#pragma omp atomic
								results[8].x[1] += suc_rate_ground2 * suc_rate_ground2;
								//#pragma omp atomic
								results[9].x[1] += suc_rate_uniq2;
								//#pragma omp atomic
								results[10].x[1] += suc_rate_uniq2 * suc_rate_uniq2;
								//Max Sols _ C(ii) condition			
								//#pragma omp flush (results)
								if (results[0].x[1] < numSols2[i]) {
									//#pragma omp critical
									{
										if (results[0].x[1] < numSols2[i])	results[0].x[1] = numSols2[i];
									}
								}
								//MIN sols _ C(ii) condition
								//#pragma omp flush (results)
								if (results[1].x[1] > numSols2[i]) {
									//#pragma omp critical
									{
										if (results[1].x[1] > numSols2[i])	results[1].x[1] = numSols2[i];
									}
								}
								//Max E _ C(ii) condition
								//#pragma omp flush (results)
								if (results[5].x[1] < tab_results[2 * i + 1][3]) {
									//#pragma omp critical
									{
										if (results[5].x[1] <  tab_results[2 * i + 1][3])	results[5].x[1] = tab_results[2 * i + 1][3];
									}
								}
								//Max dx _ C(ii) condition

								//#pragma omp flush (results)
								if (results[6].x[1] < tab_results[2 * i + 1][2]) {
									//#pragma omp critical
									{
										if (results[6].x[1] < tab_results[2 * i + 1][2])	results[6].x[1] = tab_results[2 * i + 1][2];
									}
								}

							}

							//OUTPUT:SUMMARY
							for (int i = 0; i < 2; i++) {
								results[2].x[i] /= (double)numConfigs;
								results[3].x[i] /= (double)numConfigs;
								results[4].x[i] = sqrt(fabs(results[4].x[i] / ((double)numConfigs) - results[3].x[i] * results[3].x[i]) / (double)numConfigs);
								results[7].x[i] /= (double)numConfigs;
								results[8].x[i] = sqrt(fabs(results[8].x[i] / ((double)numConfigs) - results[7].x[i] * results[7].x[i]) / (double)numConfigs);
								results[9].x[i] /= (double)numConfigs;
								results[10].x[i] = sqrt(fabs(results[10].x[i] / ((double)numConfigs) - results[9].x[i] * results[9].x[i]) / (double)numConfigs);
							}
							ofile << N << "\t" << M << "\t" << "real" << "\t" << results[0].x[0] << "\t" << results[1].x[0] << "\t" << results[2].x[0] << "\t" << results[3].x[0] << "\t" << results[4].x[0] << "\t" << results[5].x[0] << "\t" << results[6].x[0] <<"\t"<< method << "\n";
							ofile << N << "\t" << M << "\t" << "imag" << "\t" << results[0].x[1] << "\t" << results[1].x[1] << "\t" << results[2].x[1] << "\t" << results[3].x[1] << "\t" << results[4].x[1] << "\t" << results[5].x[1] << "\t" << results[6].x[1] << "\t" << method << std::endl;

							{
								GeometryVector temp;	temp.Dimension = 4;
								temp = GeometryVector(M, results[7].x[0], results[8].x[0]);
								summaryReal1.push_back(temp);
								temp = GeometryVector(M, results[9].x[0], results[10].x[0]);
								summaryReal2.push_back(temp);
								temp = GeometryVector(M, results[7].x[1], results[8].x[1]);
								summaryImag1.push_back(temp);
								temp = GeometryVector(M, results[9].x[1], results[10].x[1]);
								summaryImag2.push_back(temp);
							}

							//Energy Distributions: preparation
#pragma omp parallel for  shared(tab_energy,Energy0, dE, EnergyDist)
							for (int i = 0; i < numConfigs; i++) {//For each TARGET

																  //C(i) condition: REAL
								for (int l = 0; l < tab_energy[2 * i].size(); l++) {//For each TRIAL initial configs
									int n = 0;
									std::vector<double> list(dE.size());
									for (int m = 0; m < tab_energy[2 * i][l].size(); m++) {
										n = (int)floor(log(tab_energy[2 * i][l][m] / Energy0) / log(xdE));
										if (n >= 0 && n < list.size())
											list.at(n)++;
									}

									double num = (double)tab_energy[2 * i][l].size();
									if (num > 0.5) {
										for (int m = 0; m < list.size(); m++) { //probability distributions
											list.at(m) /= num;
										}
									}

#pragma omp atomic
									EnergyDist[l][0].x[3]++; // the number of target configurations yielding l-th solution.


									for (int m = 0; m < list.size(); m++) { //Average of distributions!
#pragma omp atomic
										EnergyDist[l][m].x[1] += list[m];
#pragma omp atomic
										EnergyDist[l][m].x[2] += list[m] * list[m];
									}

								}


								//C(ii) condition: IMAG
								if (M % 2 != 0) {//real \neq imag
									for (int l = 0; l < tab_energy[2 * i + 1].size(); l++) {//For each TRIAL initial configs
										int n = 0;
										std::vector<double> list(dE.size());
										for (int m = 0; m < tab_energy[2 * i + 1][l].size(); m++) {
											n = (int)floor(log(tab_energy[2 * i + 1][l][m] / Energy0) / log(xdE));
											if (n >= 0 && n < list.size())
												list.at(n)++;
										}
										double num = (double)tab_energy[2 * i + 1][l].size();
										if (num > 0.5) {
											for (int m = 0; m < list.size(); m++) { //probability distributions
												list.at(m) /= num;
											}
										}

#pragma omp atomic
										EnergyDist[3 + l][0].x[3]++; // the number of target configurations yielding l-th solution.


										for (int m = 0; m < list.size(); m++) { //Average of distributions!
#pragma omp atomic
											EnergyDist[3 + l][m].x[1] += list[m];
#pragma omp atomic
											EnergyDist[3 + l][m].x[2] += list[m] * list[m];
										}

									}
								}

							}

							//OUPUT: Energy Dist (not multi-threaded)
							for (int i = 0; i < 6; i++) {
								//The number of i-th solutions for all target configurations
								//i=0	trivial solution
								//i=1	nontrivial
								//i=2	secondary nontrivial
								double num = EnergyDist[i][0].x[3];


								if (num > 0.5) {//Only when there are more than a single datum
									for (int j = 0; j < EnergyDist[i].size(); j++) {
										EnergyDist[i][j].x[2] = sqrt(fabs(EnergyDist[i][j].x[2] / num - EnergyDist[i][j].x[1] / num * EnergyDist[i][j].x[1] / num) / num);
										EnergyDist[i][j].x[1] /= num; //*dE[j];
																	  //EnergyDist[i][j].x[2] /= num; //*dE[j];
									}

									if (i < 3) {
										sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Edist_real_%d", dir, dim, systName, method, N, (int)M, i);
									}
									else {
										sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Edist_imag_%d", dir, dim, systName, method, N, (int)M, i - 3);
									}
									EnergyDist[i][0].x[3] = 0;

									if (i >= 3 && M % 2 == 0) {//When C(i) and C(ii) are the same 
										WriteFunction(EnergyDist[i - 3], filename);
									}
									else {
										WriteFunction(EnergyDist[i], filename);
									}
								}
							}

							tab_energy.clear();


							//Proximity Distributions
							if (doProxMetric) {
#pragma omp parallel for shared(ProximityDist, av_X2, tab_proxy)
								for (int i = 0; i < numConfigs; i++) {
									//original data are processed s.t. the size of proximity is not larger than 4 (total, trivial, and leading nontrivial, and nontrivials).

									//REAL; C(I) condition
									for (int l = 0; l < tab_proxy[2 * i].size(); l++) {// l th solution!
										int n = 0;
										std::vector<double> list(av_X2.size()); //histogram
										for (int m = 0; m < tab_proxy[2 * i][l].size(); m++) {
											n = (int)floor(log(tab_proxy[2 * i][l][m] / x0) / log(d_x2));
											if (n >= 0 && n < list.size())
												list.at(n)++;
										}
										double deno = (double)tab_proxy[2 * i][l].size();
										if (deno > 0.5) {
											for (int m = 0; m < list.size(); m++)
												list.at(m) /= deno;//*dE[m];
										}

#pragma omp atomic
										ProximityDist[l][0].x[3]++;	//The occurence of l-th solution
																	//double num_inv = 1.0 / (double)ProximityDist[l][0].x[3];
																	//	std::cout << "num = " << num << ",\t";
										for (int m = 0; m < list.size(); m++) {
#pragma omp atomic
											ProximityDist[l][m].x[1] += list[m];
#pragma omp atomic
											ProximityDist[l][m].x[2] += list[m] * list[m];
											//ProximityDist[l][m].x[1] = (1.0 - num_inv)*ProximityDist[l][m].x[1] + num_inv * list[m];
											//ProximityDist[l][m].x[2] = (1.0 - num_inv)*ProximityDist[l][m].x[2] + num_inv * list[m] * list[m];
										}
									}

									//IMAG; C(ii) condition
									if (M % 2 != 0) {//a.k.a. C(i) and C(ii) are different
										for (int l = 0; l < tab_proxy[2 * i + 1].size(); l++) {
											int n = 0;
											std::vector<double> list(av_X2.size()); //histogram
											for (int m = 0; m < tab_proxy[2 * i + 1][l].size(); m++) {
												n = (int)floor(log(tab_proxy[2 * i + 1][l][m] / x0) / log(d_x2));
												if (n >= 0 && n < list.size())
													list.at(n)++;
											}
											double deno = (double)tab_proxy[2 * i + 1][l].size();
											if (deno > 0.5) {
												for (int m = 0; m < list.size(); m++)
													list.at(m) /= deno;//*dE[m];  : list = probability

											}

#pragma omp atomic
											ProximityDist[l + numLegd][0].x[3]++; 	//The occurence of l-th solution

											for (int m = 0; m < list.size(); m++) {
#pragma omp atomic
												ProximityDist[l + numLegd][m].x[1] += list[m];
#pragma omp atomic
												ProximityDist[l + numLegd][m].x[2] += list[m] * list[m];
											}
										}
									}

								}

								//OUTPUT: Proximity Measure

								for (int i = 0; i < 8; i++) {
									double num = ProximityDist[i][0].x[3];	//The occurence of ith solution
									if (num > 0) {
										for (int j = 0; j < ProximityDist[i].size(); j++) {
											ProximityDist[i][j].x[1] /= num;
											ProximityDist[i][j].x[2] = sqrt(fabs(ProximityDist[i][j].x[2] / num - ProximityDist[i][j].x[1] * ProximityDist[i][j].x[1]) / num);
											//ProximityDist[i][j].x[1] /= (1.0 - 1.0 / d_x2)*av_X2[j];
											//ProximityDist[i][j].x[2] /= (1.0 - 1.0 / d_x2)*av_X2[j];
										}


										if (i == 0) {
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_real_total", dir, dim, systName, method, N, (int)M, i);
										}
										else if (i < 4) {
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_real_%d", dir, dim, systName, method, N, (int)M, i - 1);
										}
										else if (i == 4) {
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_imag_total", dir, dim, systName, method, N, (int)M, i);
										}
										else {
											sprintf(filename, "%s%dD%s__%s__N_%d_M_%d_Proxydist_imag_%d", dir, dim, systName, method, N, (int)M, i - 5);
										}
										ProximityDist[i][0].x[3] = 0;
										WriteFunction(ProximityDist[i], filename);
									}
								}
							}



							tab_proxy.clear();



						}

					}

					else if (dim == 2) {
						//on construction....
						ofile << "not yet...\n";
					}


					sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "real");
					WriteFunction(summaryReal1, filename);
					sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "real");
					WriteFunction(summaryReal2, filename);

					sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "imag");
					WriteFunction(summaryImag1, filename);
					sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "imag");
					WriteFunction(summaryImag2, filename);
				}




			}


		}

	}

}