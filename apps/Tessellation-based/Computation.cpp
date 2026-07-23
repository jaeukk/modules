/** \file Computation.cpp
 *  \brief Translation unit for the Computation task hierarchy.
 *
 *  The substantial out-of-line-able method bodies (the large Compute / Write /
 *  AdditionalOptions / GetIndexList implementations) live here; the header
 *  carries their declarations together with the Doxygen documentation.  Method
 *  bodies that are trivial, templated, or entirely #ifdef-guarded
 *  (USE_VORO / USE_CGAL / USE_MATHGL) remain inline in Computation.h.
 */

#include "Computation.h"

//==============================================================================
// Out-of-line method bodies moved here from Computation.h.  The Doxygen
// documentation lives on the declarations in the header; these are plain
// definitions.  Bodies whose implementation is entirely #ifdef-guarded
// (USE_VORO / USE_CGAL) are intentionally left inline in the header.
//==============================================================================

void chivComputation::Compute(std::function<const SpherePacking(size_t i)> GetConfig, size_t numConfigs) {
	result.clear();
	omp_set_num_threads(this->numThreads);
	if (! this->lowStorage)
		::IsotropicSpectralDensity(GetConfig, numConfigs, CircularKMax, LinearKMax, result, KPrecision);
	else
		::IsotropicSpectralDensity_Lstorage(GetConfig, numConfigs, CircularKMax, LinearKMax, result, KPrecision);
}
void chivComputation::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_chiv");
	WriteFunction(result, name);
}
void chivComputation::AdditionalOptions(const std::string & option, std::istream & input, std::ostream & output) {
	if (option == "LowStorage") {
		//default = false; true => don't use map <,> in the calculation
		//input >> lowStorage;
		this->lowStorage = true;
	}
	else
		output << "Unrecognized command!\n";
}
void VoronoiCellDistribution::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_cellvolDist");
	WriteFunction(result, name);
}
void VoronoiCellVertexDistribution::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_VertexDist");
	WriteFunction(result, name);
}
void MaxPackingFraction::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_MaxPackingFraction");
	WriteFunction(maxPhi, name);
	std::ofstream add;	add.open(name + std::string(".txt"), std::fstream::app);
	add << "---------------------Summary-----------------------\n";
	add << result[0];
	add << result[1];
	add << result[2];
	add.close();
}
void LocalNumberVariance::Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	auto GetC = [&GetConfigs](size_t i)->Configuration {return Configuration(GetConfigs(i),"a"); };
	omp_set_num_threads(this->numThreads);
	MC_NV_Ensemble(GetC, numConfigs, result, Rmax, dR, (size_t)numSamplingPts);
}
void LocalNumberVariance::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_localNumVar");
	WriteFunction(result, name);
}
void LocalVolumeFractionVariance::Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	omp_set_num_threads(this->numThreads);
	MC_vfv(GetConfigs, numConfigs, result, Rmax, dR, numSamplingPts);
}
void LocalVolumeFractionVariance::Write(const std::string OutputPrefix) {
	std::string name = OutputPrefix + std::string("_localVolVar");
	WriteFunction(result, name);
}
void CentroidComputation::Write(const std::string OutputPrefix) {
	std::string fname = OutputPrefix + std::string("_SumDisplacements");
	WriteFunction(result, fname.c_str());
	if (this->output_each_config) {
		for (size_t i = 0; i < displacement_i.size(); i++) {
			fname = OutputPrefix + std::string("_SumDisplacements_") + std::to_string(i);
			WriteFunction(displacement_i[i], fname);
		}
	}
}
void MultiscaleComputation::GetIndexList(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	this->indices.resize(numConfigs, std::vector<int>());
	/*	R[i][j] = jth largest radii for configuration i
		it needs when the number of stages are not the same	*/
	std::vector<std::vector<double>> R(numConfigs, std::vector<double>());
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
	for (size_t i = 0; i < numConfigs; i++) {
		//packing fraction at each stage in THIS configuration
		SpherePacking c = GetConfigs(i);
		int numprts = c.NumParticle();
		int m = 1;	//stage
		double r_prev = c.GetCharacteristics(0), r_curr = 0.0;
		for (int j = 0; j < numprts; j++) {
			r_curr = c.GetCharacteristics(j);
			if (r_prev - r_curr > this->tolerance) {//entering to the next stage
				m++;
				R[i].push_back(r_prev);
				this->indices[i].push_back(j - 1);
				r_prev = r_curr;
			}
		}
		this->indices[i].push_back(numprts - 1);
	}
	/*-------
	for (int i = 0; i < 200; i++) {
		std::cout << i;
		for (int j = 0; j < 5; j++)
			if (i < R[j].size())
				std::cout << "\t" << R[j][i];
			else
				std::cout << "\t\t";
		std::cout << "\n";
	}*/


	bool isThereMissedStage = false;
	//check whether the number of stages is identical for all configurations
	for (size_t i = 1; i < numConfigs; i++)
		isThereMissedStage |= (R[i - 1].size() != R[i].size());
	if (isThereMissedStage) {
		/*	Particles of a certain radius are missed
			We make a complete list of particle radii.
			While doing so, we also remake the list of particle indices
			R[i][idx_stage[i]] will represent particle radius of configuration i that we are looking at
			*/
		std::vector<double> R_new;	
		std::vector<std::vector<int>> idx_list_new (numConfigs, std::vector<int> ());
		std::vector<size_t> idx_stage(numConfigs, 0);	
		
		bool conti;
		do {//repeat until idx_stage[i] == indices[i].size() - 1 for all i
			conti = false;
			bool radiiAreSame = true;
			std::vector<double> dev(numConfigs);
			double r0 = R[0][idx_stage[0]];	dev[0] = 0.0;
			for (size_t i = 1; i < numConfigs; i++) {
				dev[i] = R[i][idx_stage[i]];	//particle radius of config i
				radiiAreSame &= (fabs(dev[i] - r0) < this->tolerance);
			}
			
			if (!radiiAreSame) {
				double r_max = r0;	//find the largest radii at this stage
				for (size_t i = 1; i < numConfigs; i++)
					r_max = std::max(r_max, R[i][idx_stage[i]]);
													
				R_new.push_back(r_max);
				//update index only when particle radius is the same as the smallest one.
				for (size_t i = 0; i < numConfigs; i++) {
					idx_list_new[i].push_back(this->indices[i][idx_stage[i]]);
					if (r_max - R[i][idx_stage[i]] < this->tolerance && (idx_stage[i] < R[i].size() - 1)) {
						idx_stage[i]++;
					}
				}
			}
			else {
				R_new.push_back(r0);
				for (size_t i = 0; i < numConfigs; i++) {
					idx_list_new[i].push_back(this->indices[i][idx_stage[i]]);
					if (idx_stage[i] < R[i].size() - 1)
						idx_stage[i]++;
				}
			}
			
			
			//conti = true if not all of indices become its maximum size
			for (size_t i = 0; i < numConfigs; i++) 
				conti |= (idx_stage[i] < R[i].size() - 1);
		
		} while (conti);
		this->SetIndexList(idx_list_new);
		//this->indices=idx_list_new;
	}
}
void SpectralDensity_Multiscale::Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	std::cout << "Computing the spectral densities for multiscale coated-sphere model at every " << this->interval << ":\n";
	bool fillTheList = (this->indices.size() == 0);
	//find the index list
	if (fillTheList)
		this->GetIndexList(GetConfigs, numConfigs);

	for (size_t i = 1; i <= this->indices[0].size(); i++) {
		if (i==1 || i % this->interval == 0 || i == (size_t)this->indices[0].back())
			outputIndex.push_back(i-1);
	}

	struct SkBin
	{
		double Sum1, SumK1, SumK2, SumS, SumS2;
		SkBin() : Sum1(0), SumK1(0), SumK2(0), SumS(0), SumS2(0)
		{}
	};
	size_t NumBin = (size_t)std::floor(LinearKmax / KPrecision) + 1;

	//k-points at which we compute the spectral densities
	std::vector<GeometryVector> ks = GetKs(GetConfigs(0), CircularKmax, LinearKmax);
	//ensemble average of the spectral densities at each stage
	std::vector<std::vector<GeometryVector>> chiV_k(outputIndex.size(), std::vector<GeometryVector>(NumBin, GeometryVector(0.0, 0.0, 0.0, 0.0)));
	std::vector<double> sum1(NumBin, 0.0);
	
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
	for (size_t i = 0; i < numConfigs; i++) {
		//CONFIGURATIONS
		//Itilde[l] = chi_v(ks[l]): F.T. indicator fn of config i.
		std::vector<std::complex<double>> Itilde(ks.size(), std::complex<double>(0.0, 0.0));
		//chiV_k_m[m][idx_k] = spectral density at (idx_k) th bin in m-th stage
		std::vector<std::vector<SkBin>> chiV_k_m(outputIndex.size(), std::vector<SkBin>(NumBin, SkBin()));
		SpherePacking C = GetConfigs(i);
		double V = C.PeriodicVolume();

		int idx_init, idx_fin = -1, d=C.GetDimension();
		//STEPS: from stage outputIndex[m] to stage outputIndex[m+1]
		for (size_t m = 0; m < outputIndex.size(); m++) {
			//Set the particle indices from the stage "outputIndex[m-1]" to "outputIndex[m]"
			idx_init = idx_fin + 1;
			idx_fin = indices[i][outputIndex[m]];
			//Compute Itilde[l]'s
			
			for (size_t l = 0; l < ks.size(); l++) {//At a given wavevector
				std::map < double, double> mtildes;
				double k2 = ks[l].Modulus2();
				double k = std::sqrt(k2);
				double r = 0.0;
				double mt = 0.0;
				for (int n = idx_init; n <= idx_fin; n++) {//for prts
					r = C.GetCharacteristics(n);
					auto iter = mtildes.find(r);
					if (iter != mtildes.end())
						mt = iter->second;
					else
					{
						mt = Mtilde(k, r, d);
						mtildes.insert(std::make_pair(r, mt));
					}
					Itilde[l] += mt * std::exp(std::complex<double>(0.0, ks[l].Dot(C.GetCartesianCoordinates(n))));
				}

				size_t idx_k = (size_t)std::floor(k / KPrecision);
				double s = (Itilde[l].real() * Itilde[l].real() + Itilde[l].imag() * Itilde[l].imag()) / V;
				chiV_k_m[m][idx_k].Sum1 += 1.0;
				chiV_k_m[m][idx_k].SumK1 += k;
				chiV_k_m[m][idx_k].SumK2 += k2;
				chiV_k_m[m][idx_k].SumS += s;
				chiV_k_m[m][idx_k].SumS2 += s * s;
			}
		}

		/*To chiV_k*/
		for (size_t m = 0; m < outputIndex.size(); m++) {
			for (size_t j = 0; j < NumBin; j++) {
				if (chiV_k_m[m][j].Sum1 > 0.0) {
					if (i == 0) {
						//Wavevectors
						double k = chiV_k_m[m][j].SumK1 / chiV_k_m[m][j].Sum1;
						//var in Wavevectors
						double Dk = sqrt(fabs(chiV_k_m[m][j].SumK2 / chiV_k_m[m][j].Sum1 - k * k) / chiV_k_m[m][j].Sum1);
#pragma omp atomic
						chiV_k[m][j].x[0] += k;
#pragma omp atomic
						chiV_k[m][j].x[2] += Dk;
					}

					if (i == 0 && m == 0) {
#pragma omp atomic
						sum1[j] += chiV_k_m[m][j].Sum1;
					}

					double s = chiV_k_m[m][j].SumS / chiV_k_m[m][j].Sum1;
					double s2 = chiV_k_m[m][j].SumS2 / chiV_k_m[m][j].Sum1;
#pragma omp atomic
					chiV_k[m][j].x[1] += s;
#pragma omp atomic
					chiV_k[m][j].x[3] += s2;
				}
			}
		}
	}


	/*	Summarize and Plot

	*/
	for (size_t m = 0; m < outputIndex.size(); m++) {
		for (size_t k = 0; k < NumBin; k++) {
			if (sum1[k] > 0) {
				double s = chiV_k[m][k].x[1] / (double)numConfigs;
				double var_s = std::sqrt(std::fabs(chiV_k[m][k].x[3] / (double)numConfigs - s * s) / sum1[k]/(double)numConfigs);

				chiV_k[m][k].x[1] = s;
				chiV_k[m][k].x[3] = var_s;
			}
			else {
				chiV_k[m][k] = GeometryVector(0, 0, 0, 0);
			}
		}
	}
	this->result = chiV_k;
}
void SpectralDensity_Multiscale::Write(const std::string OutputPrefix) {
	std::string filename;
	if (aroundOriginOnly) {//Write chiv(0) vs m
		std::vector<GeometryVector> chiv0;
		for (size_t i = 0; i < outputIndex.size(); i++) {
			chiv0.emplace_back(GeometryVector(outputIndex[i]+1, result[i][0].x[1], result[i][0].x[3]));
		}
		filename = OutputPrefix + std::string("_chiv0");
		WriteFunction(chiv0, filename);
	}
	else {
		for (size_t i = 0; i < outputIndex.size(); i++) {
			filename = OutputPrefix + std::string("_chiv_") + std::to_string(outputIndex[i]+1);
			WriteFunction(result[i], filename);
		}
	}

}
void VolumeFraction_Multiscale::Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	std::cout << "Computing the volume fraction at every " << this->interval << " stages:\n";
	
	bool fillTheList = (this->indices.size() == 0);
	//find the index list
	if (fillTheList) 
		this->GetIndexList(GetConfigs, numConfigs);

	int binsize = (int)round((this->indices[0].size()) / interval);
	result.clear();
	{
		GeometryVector temp(0, 0, 0, 0);	result.resize(binsize, temp);
	}
	//Compute based on the index list
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
	for (size_t i = 0; i < numConfigs; i++) {
		//packing fraction at each stage in THIS configuration
		SpherePacking c = GetConfigs(i);
		int d = c.GetDimension();
		size_t m = 1; int idx_prev = -1, tick=0;	//stage
		double V = c.PeriodicVolume();
		double phi0 = 0.0;
		for (size_t j = 0; j < this->indices[i].size(); j++) {
			phi0 += HyperSphere_Volume(d, c.GetCharacteristics(this->indices[i][j]))*(this->indices[i][j] - idx_prev) / V;
			idx_prev = this->indices[i][j];
			if (m%interval == 0) {
#pragma omp atomic
				result[tick].x[1] += phi0;
#pragma omp atomic
				result[tick].x[2] += phi0 * phi0;

				tick++;
			}
			else if (m == this->indices[i].size()) {
#pragma omp atomic
				result[tick].x[1] += phi0;
#pragma omp atomic
				result[tick].x[2] += phi0 * phi0;
			}
			m++;
		}
	}
	//Summary
	for (size_t i = 0; i < result.size(); i++) {
		if (i < result.size() - 1)
			result[i].x[0] = (i + 1) * this->interval;
		else
			result[i].x[0] = (this->indices[0].size());
		result[i].x[1] /= (double)numConfigs;
		result[i].x[2] = fabs(result[i].x[2] / (double)numConfigs - result[i].x[1] * result[i].x[1]);
		result[i].x[2] = sqrt(result[i].x[2] / ((double)numConfigs - 1));
	}
}
void VolumeFraction_Multiscale::Write(const std::string OutputPrefix) {
	std::string filename = OutputPrefix + std::string("_PackingFractions");
	WriteFunction(result, filename);
}
void VolumeFractionVariance_Multiscale::Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
	std::cout << "Computing the local volume fraction variance at every " << this->interval << " stages:\n";
	bool fillTheList = (this->indices.size() == 0);
	//find the index list
	if (fillTheList)
		this->GetIndexList(GetConfigs, numConfigs);

	for (size_t i = 1; i <= this->indices[0].size(); i++) {
		if (i % this->interval == 0 || i == (size_t)this->indices[0].back())
			outputIndex.push_back(i);
	}

	for (auto iter = outputIndex.begin(); iter != outputIndex.end(); iter++) {
		std::vector<int> idx_prt; //idx_prt[i] = particle idx in the configuration i at the *iter stage
		for (size_t i = 0; i < numConfigs; i++)
			idx_prt.push_back(this->indices[i][*iter]);
		//Remove particles at stages higher than m
		auto getC = [&GetConfigs, &idx_prt](size_t j)->SpherePacking {
			SpherePacking c = GetConfigs(j);
			int num = c.NumParticle() - 1;
			for (int k = num; k > idx_prt[j]; k--) {
				c.DeleteParticle(k);
			}
			return c;
		};

		std::vector<GeometryVector> RESULT;
		MC_vfv(getC, numConfigs, RESULT, this->Rmax, this->dR, this->numSamplePts);
		this->result.push_back(RESULT);
	}
}
void VolumeFractionVariance_Multiscale::Write(const std::string OutputPrefix) {
	std::string filename;
	for (size_t i = 0; i < outputIndex.size(); i++) {
		filename = OutputPrefix + std::string("_LocalVolVar_") + std::to_string(outputIndex[i]);
		WriteFunction(result[i], filename);
	}
}
