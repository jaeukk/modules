/** \file TessellationAnalysisCLI.cpp
 *  \brief CLI for the tiling application, decoupled from the config generators.
 *
 *  A variant of main.cpp that is decoupled from the `models` generators
 *  (RSA_Gen, HardSphereFluids/HSF_Gen): the internal "GenerateConfig" branches
 *  are removed, so initial configurations are supplied only by loading files
 *  produced by other apps ("InitialConfiguration LoadPrefix ..."), in the spirit
 *  of the PairStatCLI driver. Everything else — the Voronoi-based and multiscale
 *  Computation tasks, and the polydisperse construction via MakeMono2Poly — is
 *  unchanged.
 *
 *  The heavy Voronoi analyses need CGAL / voro++ / MathGL (USE_CGAL / USE_VORO /
 *  USE_MATHGL); without them the corresponding branches are #ifdef-guarded out.
 *
 *  Usage:
 *      ./TessellationCLI.out [masterDir]
 *  then on stdin:
 *      Tiling           # enter the tiling CLI (see TilingProblemCLI)
 *      Exit
 */

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>

#include "Computation.h"        // Computation hierarchy (+ SpacefillingPacking, Voronoi)
#include "PolyPacking.h"        // MakeMono2Poly*

// Library globals owned by the executable (modular Verbosity convention).
size_t Verbosity = 2;
std::string masterDir = "./";

//==============================================================================
// TilingProblemCLI sliced verbatim from the legacy TilingProblem.cpp
// (lines 1364-2100).
//==============================================================================
/** \brief Interactive CLI for the tiling problems (Voronoi-based and multiscale packings). (CLAUDE)
 *  \return  0 on normal exit. */
int TilingProblemCLI() {
	std::string dir = "";
	std::string savePrefix = "";
	std::string loadPrefix = "";
	std::string loadFileType = "";
	std::string saveFileType(".txt");
	int idx_beg = 0, idx_end = 1;
	

	char tempstring[1000] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	std::function <SpherePacking(size_t i)> getInitConfig = nullptr;
	ConfigGen *initC = nullptr;
	SpherePacking * pSpheres = nullptr;
	int dim = 0;
	size_t numParticles = 0;
	size_t numConfigs = 0;
	int numThreads = 1;
	
	ofile << "Curr directory = " << masterDir << std::endl;
	int count_unrecognized_commands = 0;
	auto AutoExit = [&count_unrecognized_commands, &ofile](char * txt)->void {
		if (count_unrecognized_commands > 5) {
			ofile << "Too many unrecognized commands! AutoEuit the loop\n";
			sprintf(txt,"Exit");
		}
	};
	for (;;)
	{
		ifile >> tempstring;

		AutoExit(tempstring);
		if (strcmp(tempstring, "Exit") == 0) {
			if (pSpheres != nullptr)
				delete pSpheres;
			if (initC != nullptr)
				delete initC;

			return 0;
		}
		else if (strcmp(tempstring, "NumConfigs") == 0) {
			ifile >> numConfigs;
		}
		else if (strcmp(tempstring, "SavePrefix") == 0) {//If this is 
			ofile << "(default/ what you want)" << std::endl;
			ifile >> tempstring;	ofile << "\t" << tempstring << "\n";
			if (strcmp(tempstring, "default") == 0)
				savePrefix = "";
			else {
				savePrefix = masterDir + std::string(tempstring);
				ofile << "full OutputPrefix=" << savePrefix << std::endl;
			}
		}
		else if (strcmp(tempstring, "NumThreads") == 0) {
			ifile >> numThreads;
			omp_set_num_threads(numThreads);
		}
		else if (strcmp(tempstring, "dim") == 0) {
			ifile >> dim;	ofile << "\t" << dim << std::endl;
		}
		//Externally generate initial configurations:
		else if (strcmp(tempstring, "InitialConfiguration") == 0) {
			ifile >> tempstring;
			if (strcmp(tempstring, "LoadPrefix") == 0) {
				std::fstream readFile;

				bool files_are_good = true;	char tempName[300];
				do {//Repeat the input until acquring existent files 
					if (!files_are_good)
						ofile << "files do not exist!" << std::endl;

					files_are_good = true;
					ofile << "file name (includes \"%d\" for multiple files) = \t";
					ifile >> loadPrefix;	ofile << loadPrefix << std::endl;
					ofile << "file type = \t";
					ifile >> loadFileType; ofile << loadFileType << std::endl;
					
					ofile << "index::begin= \t";
					ifile >> idx_beg;	ofile << idx_beg << std::endl;
					ofile << "index::end = \t";
					ifile >> idx_end;	ofile << idx_end << std::endl;
					
					for (int i = idx_beg; i < idx_end; i++) {
						sprintf(tempName, (loadPrefix+ loadFileType).c_str(), i);
						ofile << tempName;
						readFile.open(tempName, std::fstream::in);
						files_are_good = files_are_good && readFile.good();
						readFile.close();
						if (readFile.good())
							ofile << ": good" << std::endl;
					}
						 
				} while (!files_are_good);
				//filePrefix = 
				numConfigs = idx_end - idx_beg;
				ofile << "The number of configuration is set to " << numConfigs << "\n";
				getInitConfig = [idx_beg, &loadFileType, &loadPrefix, numConfigs](size_t i)->SpherePacking {
					SpherePacking config;
					char tempName[300];
					if (strcmp(loadFileType.c_str(), ".txt") == 0) {
						sprintf(tempName, loadPrefix.c_str(), idx_beg+(i-idx_beg)%numConfigs);
						ReadConfiguration(config, tempName);
					}
					else if (strcmp(loadFileType.c_str(), ".SpherePack") == 0) {
					
					}
					return config;
				};
			}
			else if (strcmp(tempstring, "EmptyUnitBox") == 0) {
				if (dim <= 0) {
					ofile << "wrong dimension \n";
					break;
				}
				else {
					getInitConfig = [dim](size_t i)->SpherePacking {
						SpherePacking temp(GetUnitCubicBox(dim), 0.0);
						return temp;
					};
					ofile << "InitialConfiguration = EmptyUnitBox" << std::endl;
				}
			}


		}
		else if (strcmp(tempstring, "VoronoiBased") == 0) {
			const double phi_stored = 0.01;//packing fraction of the constructed dispersions
			
			ofile << "CLI for the tiling-based procedure with Voronoi tessellations\n";
			ofile << "Input \"Start\" after inputting tasks\n";
			std::vector<Computation *> tasks;
				
			bool defaultSavePrefix = true;
			bool construct = false;
			bool analyze = false;	bool overwrite = false;
			bool discardConfig = false;	int largestIndex = 1;
			bool mono = false;
			char systemType[30];
			double rho = 0.0, phi_init, phi_target = phi_stored;
			std::string fileName;
			for (;;) {
					
				ifile >> tempstring;
				AutoExit(tempstring);
				if (strcmp(tempstring, "Exit") == 0) {
					for (auto iter = tasks.begin(); iter != tasks.end(); iter++)
						delete *iter;
					ofile << "Exit from Voronoi-based tiling procedure\n";
					return 0;
				}
				if (strcmp(tempstring, "Construct")==0) {
					construct = true;	ofile << "Construction is on\n";
				}
				else if (strcmp(tempstring, "Analyze") == 0) {
					analyze = true; ofile << "Analysis is on\n";
				}
				else if (strcmp(tempstring, "OriginalVolumeFraction") == 0) {
					ifile >> phi_init;	ofile << "\tvolume fraction of the initial packing = " << phi_init<<"\n";
				}
				else if (strcmp(tempstring, "Start") == 0) {
					ofile << "Start!!"<<std::endl;
					if (getInitConfig != nullptr)
						break;
					else {
						ofile << "Initial packings are not defined!\n";
						ofile << "Input \"InitialConfiguration\"\n";
					}
				}
				else if (strcmp(tempstring, "NumConfigs") == 0) {ifile >> numConfigs;	ofile << "The number of configurations = \t" << numConfigs << "\n";	}
				else if (strcmp(tempstring, "InitialConfiguration") == 0) {
					ifile >> tempstring;
					if (strcmp(tempstring, "LoadPrefix") == 0) {
						std::fstream readFile;
						ofile << "Load Configurations\n";
						bool files_are_good = true;	char tempName[300];
						do {//Repeat the input until acquring existent files 
							if (!files_are_good)
								ofile << "files do not exist!" << std::endl;

							files_are_good = true;
							ofile << "file name (includes \"%d\" for multiple files) = \t";
							ifile >> loadPrefix;	ofile << loadPrefix << std::endl;
							ofile << "file type = \t";
							ifile >> loadFileType; ofile << loadFileType << std::endl;

							ofile << "index::begin= \t";
							ifile >> idx_beg;	ofile << idx_beg << std::endl;
							ofile << "index::end = \t";
							ifile >> idx_end;	ofile << idx_end << std::endl;

							for (int i = idx_beg; i < idx_end; i++) {
								sprintf(tempName, (loadPrefix + loadFileType).c_str(), i);
								if(Verbosity>3)
									ofile << tempName << std::endl;
								readFile.open(tempName, std::fstream::in);
								files_are_good = files_are_good && readFile.good();
								readFile.close();
							}

						} while (!files_are_good);
						//filePrefix = 

						getInitConfig = [idx_beg, &loadFileType, &loadPrefix](size_t i)->SpherePacking {
							SpherePacking config;
							char tempName[300];
							if (strcmp(loadFileType.c_str(), ".txt") == 0) {
								sprintf(tempName, loadPrefix.c_str(), i);
								ReadConfiguration(config, tempName);
							}
							else if (strcmp(loadFileType.c_str(), ".SpherePack") == 0) {

							}
							return config;
						};
					}
				}
				else if (strcmp(tempstring, "mono") == 0) {
					mono = true;	ofile << "\t monodisperse packings are considered\n";
				}
				else if (strcmp(tempstring, "SavePrefix") == 0) {
					ofile << "SavePrefix = ";
					ifile >> savePrefix; ofile << "\t" << savePrefix << "\n";
					defaultSavePrefix = false;
				}
				else if (strcmp(tempstring, "OverWrite") == 0) { overwrite = true; ofile << "Overwrite Configurations\n"; }
				else if (strcmp(tempstring, "NumThreads") == 0) {
					ofile << "The number of threads = ";
					ifile >> numThreads;	ofile << "\t" << numThreads << "\n";
					omp_set_num_threads(numThreads);
				}
				else if (strcmp(tempstring, "TargetVolumeFraction") == 0) {
					ofile << "Rescale to the Volume Fraction = ";
					ifile >> phi_target; ofile << "\t" << phi_target << "\n";
				}
				else if (strcmp(tempstring, "chivComputation") == 0) {
					tasks.emplace_back(new chivComputation(ifile, ofile));
					//ofile << "Additional Option? \t";
					//ifile >> tempstring;
					//(*tasks.end())->AdditionalOptions(tempstring, ifile, ofile);
				}
				else if (strcmp(tempstring, "MaxPackingFractionComputation") == 0) {
					tasks.emplace_back(new MaxPackingFraction(ifile, ofile));	
				}
				else if (strcmp(tempstring, "VoronoiCellVolumeDistribution") == 0) {
					tasks.emplace_back(new VoronoiCellDistribution(ifile, ofile));
				}
				else if (strcmp(tempstring, "VoronoiCellVertexDistribution") == 0) {
					tasks.emplace_back(new VoronoiCellVertexDistribution(ifile, ofile));
				}
				else if (strcmp(tempstring, "LocalNumberVarianceComputation") == 0) {
					tasks.emplace_back(new LocalNumberVariance(ifile, ofile));
				}
				else if (strcmp(tempstring, "LocalVolumeFractionVariance") == 0) {
					tasks.emplace_back(new LocalVolumeFractionVariance(ifile, ofile));
				}
				else if (strcmp(tempstring, "CentroidComputation") == 0) {
					tasks.emplace_back(new CentroidComputation(ifile, ofile));
				}
				else if (strcmp(tempstring, "dim") == 0) { ofile << "Dimension = "; ifile >> dim; ofile << "\t" << dim << "\n"; }
				else if (strcmp(tempstring, "DiscardConfigurations") == 0) {
					discardConfig = true;
					ofile << "Discard Configurations with indices lower than ...\t";
					ifile >> largestIndex;	ofile << "\t" << largestIndex << "\n";
				}
				else {
					ofile << "Unrecognized commands !!!\n";
					count_unrecognized_commands++;
				}
			}
			
			if (construct) {
				ofile << "-------------------------\n";
				ofile << "Construction starts\n";
				ofile << "-------------------------\n";
				if (defaultSavePrefix) {
					ofile << "default save prefix is used\n";
					savePrefix = masterDir;
					if (strcmp(systemType, "RSA") == 0)
						savePrefix += std::string("/results/RSA_configs/");
					else if (strcmp(systemType, "HSF") == 0)
						savePrefix += std::string("/results/HS_Fluids/");
					else if (strcmp(systemType, "ImLat") == 0)
						savePrefix += std::string("/results/ImperfectLattices/");
					else if (strcmp(systemType, "Poisson") == 0)
						savePrefix += std::string("/results/Poisson/");
					//e.g.:  2D_RSA_0.500__10000000
					char file[200];
					if (strcmp(systemType, "Poisson") == 0)
						sprintf(file, "%dD_%s_%0.4f__%d", dim, systemType, rho, (int)numParticles);
					else
						sprintf(file, "%dD_%s_%0.4f__%d", dim, systemType, phi_init, (int)numParticles);
					savePrefix += file;
				}
				/*
					Check existing configurations
				*/
				std::vector<int> config_idx;
				for (size_t i = 0; i < numConfigs; i++) {
					if (!overwrite) {
						char fname[300];	sprintf(fname, "%s_config_%04d___", savePrefix.c_str(), (int)i);
						std::fstream check(fname+std::string(".txt"), std::fstream::in);
						if (!check.good())
							config_idx.push_back(i);
					}
					else
						config_idx.push_back(i);
				}

				ofile << "We construct " << config_idx.size() << " configurations\n";
				ofile << "\t";
				for (size_t i = 0; i < config_idx.size(); i++)	ofile << config_idx[i] << ", ";
				ofile << "\n";

				/* Generate New Dispersions */
				auto t_init = std::time(NULL);
				for (size_t i = 0; i < config_idx.size(); i++) {
					int idx_c = config_idx[i];
					SpherePacking C = getInitConfig(idx_c);
					SpherePacking resultC;
					std::vector<double> vols;
					if (Verbosity > 2) 
						ofile << "Config " << idx_c << " is successfully constructed:\n\t";
					
#ifdef USE_CGAL
					double maxVolFraction = MakeMono2Poly_parallel(C, resultC, phi_target, vols, numThreads);
					(void)maxVolFraction;
#else
					(void)C; (void)resultC; (void)vols;
					ofile << "Construct requires USE_CGAL (polydisperse cell construction); skipping.\n";
					continue;
#endif
					//Save
					char fname[300];	sprintf(fname, "%s_config_%04d___", savePrefix.c_str(), idx_c);
					if (Verbosity > 2) {
						ofile << "\tsaved as "<<fname<<"\n";
					}
					WriteConfiguration(resultC, fname);
				}
				//if (Verbosity > 3)
				auto t_fin = std::time(NULL);
				ofile << "time = " << t_fin - t_init << " seconds ->\t" << (t_fin - t_init) / numConfigs << "seconds/config\n";
			}
			if (savePrefix.empty()) {
				ofile << "save prefix is empty. Use a default name = ";
				savePrefix = "./temp";
				ofile << savePrefix << "\n";
			}
			if (analyze) {
				ofile << "\n\n-------------------------\n";
				ofile << "Analysis starts\n";
				ofile << "-------------------------\n";
				std::function<SpherePacking(size_t i)> GetC;
				std::function<SpherePacking(size_t i)> GetC_mono;
				/*If new systems were constructed -> load them*/
				if (construct) {
					ofile << "Start analyzing the constructed configurations\n";
					GetC = [&savePrefix, phi_target](size_t i)->SpherePacking {
						SpherePacking c;	char name[300];	sprintf(name, "%s_config_%04d___", savePrefix.c_str(), (int)i);
						ReadConfiguration(c, name);
						double phi = c.PackingFraction();
						if (std::fabs(phi - phi_target) > 1e-5) {
							double factor = pow(phi_target / phi, 1.0 / c.GetDimension());
							c.Rescale_prtOnly(factor);
						}
						return c;
					};
				}

				/*exported configs*/
				else
					GetC = [&getInitConfig, phi_target](size_t i)->SpherePacking {
					SpherePacking c = getInitConfig(i);
					double phi = c.PackingFraction();
					if (std::fabs(phi - phi_target) > 1e-5) {
						double factor = pow(phi_target / phi, 1.0 / c.GetDimension());
						c.Rescale_prtOnly(factor);
					}
					return c;
				};

				if (mono) {
					GetC_mono = [&GetC, phi_target](size_t i)->SpherePacking {
						SpherePacking c = GetC(i);
						double r = pow(phi_target/(HyperSphere_Volume(c.GetDimension(), 1.0)*c.NumParticle()/c.PeriodicVolume()), 1.0 / c.GetDimension());
						for (size_t j = 0; j < c.NumParticle(); j++)
							c.UpdateCharacteristics(r, j);

						return c;
					};
				}


				for (auto iter = tasks.begin(); iter != tasks.end(); iter++){
					(*iter)->origPackingFraction = phi_init;
					(*iter)->SetNumThreads(numThreads);
					(*iter)->Compute(GetC, numConfigs);
					(*iter)->Write(savePrefix);
				}
				if (mono) {
					ofile << "Monodisperse packings are computed\n";
					for (auto iter = tasks.begin(); iter != tasks.end(); iter++) {
						if((*iter)->initialPackingFriendly){
							(*iter)->origPackingFraction = phi_init;
							(*iter)->Compute(GetC_mono, numConfigs);
							(*iter)->Write(savePrefix+std::string("_mono"));
						}
					}
				}

			}
			/* discard configurations whose indices are higher than "largestIndex"*/
			if (discardConfig) {
				char fname[300];
				int idx_begin = largestIndex, idx = idx_begin;
				int count = 0;
				
				bool check = true;
				while (check) {
					sprintf(fname, "%s_config_%04d___.txt", savePrefix.c_str(), (int)idx);
					std::fstream temp(fname, std::fstream::in);	check = temp.good();
					if (check) {
						std::remove(fname);
						count++; idx++;
					}
				}
				if (count > 0)
					ofile << "Config " << idx_begin << " - " << idx << ": successfully deleted\n";
				else
					ofile << "No such cnofigurations\n";
			}
			//delete pointers
			for (auto iter = tasks.begin(); iter != tasks.end(); iter++)
				delete *iter;
		}
			
		else if (strcmp(tempstring, "MultiscalePacking") == 0) {
			std::vector<double> p(2, 0.0);	//Parameters for the scaling function
			int ni = 0, nf = 0;				//initial & final stages
			char fnName[10];				//name of scaling function
			std::vector<int> saveIndex;// = std::vector<int>();
			bool multiThread_in_config = false; // true = each config is multi-threaded
												// false = loop is multi-threaded
			bool generate = false, analyze = false, overwrite =false;
			std::vector<MultiscaleComputation *> tasks;
			double targetPackingFraction = 1.0;
			ofile << "CLI for generating and analyzing the Multiscale sphere packings:\n";
			if (getInitConfig != nullptr) {
				for (;;) {//Input
					ifile >> tempstring;
					AutoExit(tempstring);
					if (strcmp(tempstring, "Start") == 0)
						break;
					if (strcmp(tempstring, "Construct") == 0)
						generate = true;
					else if (strcmp(tempstring, "Analyze") == 0) {
						analyze = true;
						ofile << "Analysis is on: the following options are available \n";
						ofile << "\tVolumeFraction\n\tSpectralDensities\n\tLocalVolumeFractionVariance\n";
						for (;;) {
							ifile >> tempstring;
							AutoExit(tempstring);

							if (strcmp(tempstring, "Exit") == 0) {
								ofile << "Exit from MultiscalePacking::Analysis options\n";
								break;
							}
							else if (strcmp(tempstring, "VolumeFraction") == 0) {
								tasks.emplace_back(new VolumeFraction_Multiscale(ifile, ofile));
							}
							else if (strcmp(tempstring, "SpectralDensities") == 0) {
								tasks.emplace_back(new SpectralDensity_Multiscale(ifile, ofile));
							}
							else if (strcmp(tempstring, "LocalVolumeFractionVariances") == 0) {
								//tasks.emplace_back(new LocalVolumeFractionVariance(ifile, ofile));
							}
							else {
								ofile << "Unrecognized Commands\n";
								count_unrecognized_commands++;
							}

						}
					}
					else if (strcmp(tempstring, "Exit") == 0) {
						ofile << "Exit from Multiscale-disk tilings\n";
						generate = false;	analyze = false;
						return 0;
					}
					else if (strcmp(tempstring, "InitStage") == 0) {
						ifile >> ni;	ofile << "\t" << ni << std::endl;
					}
					else if (strcmp(tempstring, "FinStage") == 0) {
						ifile >> nf;	ofile << "\t" << nf << std::endl;
					}
					else if (strcmp(tempstring, "ScalingFunction") == 0) {
						ofile << " (Power/Zeta) "; ifile >> fnName; ofile << "\t" << fnName << std::endl;
					}
					else if (strcmp(tempstring, "MaxParticleVolume") == 0) {
						ifile >> p[0]; ofile << "\t" << p[0] << std::endl;
					}
					else if (strcmp(tempstring, "ScalingParameter") == 0) {
						ifile >> p[1]; ofile << "\t" << p[1] << std::endl;
					}
					else if (strcmp(tempstring, "NumConfigs") == 0) {
						ifile >> numConfigs; ofile << "\t" << numConfigs << std::endl;
					}
					else if (strcmp(tempstring, "NumThreads") == 0) {
						ifile >> numThreads; ofile << "\t" << numThreads << std::endl;
					}
					else if (strcmp(tempstring, "AutosaveInterval") == 0) {
						if ((nf - ni) <= 0) {
							ofile << "Please input 'InitStage' and 'FinStage' first" << std::endl;
						}
						else {
							ofile << "Input a positive integer: " << std::endl;
							int interval = 0;
							ifile >> interval;	ofile << "\t" << interval << std::endl;
							if (interval > 0) {
								ofile << "Autosave indices are ";
								for (int idx = ni + interval - 1; idx < nf; idx += interval) {
									saveIndex.push_back(idx);
									ofile << idx << "\t";
								}
								ofile << std::endl;
							}

						}
					}
					else if (strcmp(tempstring, "MultiThread") == 0) {
						ofile << "Choose (0 = parallelize each config /1 = parallelize the loop)? \t";
						ifile >> multiThread_in_config; ofile << "\t" << multiThread_in_config << std::endl;
					}
					else if (strcmp(tempstring, "TargetVolumeFraction") == 0) {
						ofile << "Set the target volume fraction = \t";
						ifile >> targetPackingFraction; ofile << "\t" << targetPackingFraction << std::endl;
					}
					else if (strcmp(tempstring, "Overwrite") == 0) {
						ofile << "We overwrite configurations\n";	overwrite = true;
					}
					else {
						ofile << "Unrecognized Commands\n";	count_unrecognized_commands++;
					}
				}

				if (generate) {
					ofile << "--------------------------------------\n";
					ofile << "\tConstruction starts\n";
					ofile << "--------------------------------------\n";
					if (savePrefix.empty()) {//Use a default filename
						char filename[100];
						sprintf(filename, "%dD_%s_%0.2f_Vmax_%0.5f_%d", dim, fnName, p[1], p[0], nf);
						savePrefix = masterDir + std::string("/results/CoatedSpheres/") + std::string(filename);
						ofile << "Save file name is not specified\n the default format is used: " << savePrefix << std::endl;
					}

					/*
					Check existing configurations
					*/
					std::vector<int> config_idx;
					for (size_t i = 0; i < numConfigs; i++) {
						if (!overwrite) {
							char fname[300];	sprintf(fname, "%s_config_%04d", savePrefix.c_str(), (int)i);
							std::fstream check(fname + std::string(".txt"), std::fstream::in);
							if (!check.good())
								config_idx.push_back(i);
						}
						else
							config_idx.push_back(i);
					}

					ofile << "We will construct " << config_idx.size() << " configurations\n";
					ofile << "\t";
					for (size_t i = 0; i < config_idx.size(); i++)	ofile << config_idx[i] << ", ";
					ofile << "\n";

					if (multiThread_in_config) {//Generate multiple configurations at once.
#pragma omp parallel for schedule (guided) num_threads(numThreads)
						for (size_t i = 0; i < config_idx.size(); i++) {
							int idx = config_idx[i];
							SpherePacking c = getInitConfig(idx);
							double error = getSpacefillingPacking2(c, ni, nf, p[0], p[1], fnName, (savePrefix + (std::string("__") + std::to_string(idx))).c_str(), std::time(NULL)+i, 1, saveIndex);
							ofile << "config " << idx << "::completed error" << error << std::endl;
							
							char savename[300];
							sprintf(savename, "%s_config_%04d", savePrefix.c_str(), idx);
							if (strcmp(saveFileType.c_str(), ".txt") == 0) {
								WriteConfiguration(c, savename);
							}
							else {
							}
						}
					}
					else {//Generate a single configuration at once.
						for (size_t i = 0; i < config_idx.size(); i++) {
							int idx = config_idx[i];
							omp_set_num_threads(numThreads);
							SpherePacking c = getInitConfig(idx);
							double error = getSpacefillingPacking2(c, ni, nf, p[0], p[1], fnName, savePrefix.c_str(), std::time(NULL)+i, numThreads, saveIndex);
							ofile << "config " << idx << "::completed error" << error << std::endl;
							
							char savename[300];
							sprintf(savename, "%s_config_%04d", savePrefix.c_str(), idx);
							if (strcmp(saveFileType.c_str(), ".txt") == 0) {
								WriteConfiguration(c, savename);
							}
							else {
							}
						}
					}
				}

				if (analyze) {
					ofile << "--------------------------------------\n";
					ofile << "\tAnalysis starts\n";
					ofile << "--------------------------------------\n";
					std::function<SpherePacking(size_t i)> GetC;
					if (generate) {
						ofile << "Start analyzing the constructed configurations\n";
						GetC = [&savePrefix, targetPackingFraction](size_t i)->SpherePacking {
							SpherePacking c;	char name[300];	sprintf(name, "%s_config_%04d", savePrefix.c_str(), (int)i);
							ReadConfiguration(c, name);
							double factor = pow(targetPackingFraction, 1.0 / c.GetDimension());
							c.Rescale_prtOnly(factor);
							return c;
						};
					}
					/*exported configs*/
					else
						GetC = [&getInitConfig, targetPackingFraction](size_t i)->SpherePacking {
						SpherePacking c = getInitConfig(i);
						double factor = pow(targetPackingFraction, 1.0 / c.GetDimension());
						c.Rescale_prtOnly(factor);
						return c;
					};
					std::vector<std::vector<int>> indexList;
					for (auto iter = tasks.begin(); iter != tasks.end(); iter++) {
						(*iter)->SetNumThreads(numThreads);
						if (iter != tasks.begin())
							(*iter)->SetIndexList(indexList);
						(*iter)->Compute(GetC, numConfigs);
						if (iter == tasks.begin())
							(*iter)->OutputIndexList(indexList);
						(*iter)->Write(savePrefix);
					}
				}

			}
			else {
				ofile << "Need Initial Configurations\n";
				break;
			}

			for (auto iter = tasks.begin(); iter != tasks.end(); iter++)
				delete *iter;
		}
			
		else {
			ofile << "Invalid input!!" << std::endl;
			count_unrecognized_commands++;
		}
	}


	return 0;
}

//==============================================================================
// Entry point: the "Tiling" command extracted from the legacy Interfaces.h cli()
//==============================================================================
/** \brief Program entry point; reads top-level commands and dispatches "Tiling". (CLAUDE)
 *  \param argc  Argument count; argv[1] (if present) sets masterDir.
 *  \param argv  Argument vector.
 *  \return      0 on normal exit. */
int main(int argc, char ** argv)
{
	std::istream & ifile = std::cin;

	if (argc > 1)
		masterDir = argv[1];

	char tempstring[1000] = {};
	for (;;)
	{
		ifile >> tempstring;
		if (!ifile)
			return 0;                       // EOF -> done
		if (strcmp(tempstring, "Exit") == 0)
		{
			return 0;
		}
		else if (strcmp(tempstring, "Verbosity") == 0)
		{
			ifile >> Verbosity;
		}
		else if (strcmp(tempstring, "Tiling") == 0)
		{
			TilingProblemCLI();
		}
		else
			std::cerr << "Unrecognized command!\n";
	}
	return 0;
}
