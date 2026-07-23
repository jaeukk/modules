#include <list>
#include <set>
#include <Plots.h>
#include <MC_System.h>
#include <StructureOptimization.h>
#include <Interfaces.h>
#include <CollectiveCoordinatePotential.h>
#include <CCPotential_Kernels.h>
#include <PairCorrelation.h>
#include <StructureFactor.h>


int ReadProblemFile(char * Filename, size_t & NumParticles, DimensionType & dim, double & boxSize, Configuration * &pConfig, ShiftedCCPotential * &pPotential, int NumThread)
{
	/////////////
	// in : Filename, NumThread
	// out : NumKPoints, NumParticles, dim, boxSize, pConfig, pPotential
	if(pConfig!=nullptr)
	{
		delete pConfig;
		pConfig=nullptr;
	}
	if(pPotential!=nullptr)
	{
		delete pPotential;
		pPotential=nullptr;
	}

	signed HasWeight=-1;
	std::string s;
	//////////////////
	bool particlesInit = false;
	size_t NumKPoints = 0;
	bool dimInit = false;
	bool sizeInit = false;
	std::list< std::vector<int> > rawVec;
	std::list< double > rawWeight;
	std::fstream inputFile(Filename, std::fstream::in);
	while (!inputFile.fail()) 
	{
		getline(inputFile, s);

		if(!s.empty())
		{
			std::stringstream ss(s);

			if(s[0] == '#')
			{
				std::string token;

				ss >> token;

				if(token == "#function")
				{
					ss >> token;
					if(token=="stealth")
						HasWeight=0;
					else if(token=="weight")
						HasWeight=1;
				} 
				else if(token == "#particles")
				{
					ss >> NumParticles;
					particlesInit = true;
				} 
				else if(token == "#dimensions")
				{
					ss >> dim;
					dimInit = true;
				} 
				else if(token == "#size")
				{
					if(!dimInit)
					{
						std::cerr << "#dimensions parameter need to be before "
							<< "the #size parameter.\n";
						exit(1);
					}
					for(int d = 0; d < dim; d++)
					{
						double temp;
						ss >> temp;
						if(d==0)
							boxSize=temp;
						else if(boxSize!=temp)
						{
							std::cerr<<"Error in CollectiveCoordinateCLI : non-cubix box is not supported";
							return 1;
						}
					}
					sizeInit = true;
				}
			} 
			else 
			{

				std::vector<int> vec(dim);
				for(int d = 0; d < dim; d++)
				{
					ss >> vec[d];
				}

				rawVec.push_back(vec);
				if(HasWeight)
				{
					double w;
					ss >> w;
					rawWeight.push_back(w);
				}
				else
				{
					rawWeight.push_back(1.0);
				}

				NumKPoints++;
			}
		}
	}
	if(HasWeight==(-1))
	{
		std::cerr <<"Error in CollectiveCoordinateCLI : Problem file does not contain the #function parameter.\n";
		return 1;
	}
	if(!particlesInit)
	{
		std::cerr <<"Error in CollectiveCoordinateCLI : Problem file does not contain the #particles parameter.\n";
		return 1;
	}
	if(!dimInit)
	{
		std::cerr<<"Error in CollectiveCoordinateCLI : Problem file does not contain the #dimensions parameter.\n";
		return 1;
	}
	if(!sizeInit)
	{
		std::cerr << "Error in CollectiveCoordinateCLI : Problem file does not contain the #size parameter.\n";
		return 1;
	}
	double CellSize=1.0/std::pow(NumParticles, 1.0/(double)(dim) );
	pConfig=new Configuration(GetUnitCubicBox(dim, CellSize));
	pConfig->Rescale(boxSize);
	RandomGenerator gen(2);
	for(size_t i=0; i<NumParticles; i++)
		pConfig->Insert("A", gen);
	ShiftedCCPotential * MypPotential = new ShiftedCCPotential(dim);
	auto iter2=rawWeight.begin();
	for(auto iter=rawVec.begin(); iter!=rawVec.end()&&iter2!=rawWeight.end(); iter++, iter2++)
	{
		GeometryVector k(dim);
		for(DimensionType d=0; d<dim; d++)
			k.AddFrom(iter->at(d)*pConfig->GetReciprocalBasisVector(d));
		MypPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(k, *iter2));
	}
	pPotential=MypPotential;
	MypPotential->ParallelNumThread=NumThread;

	return 0;
}

//Type :
//0: stealth
//1: overlap
//2: (K-k)^6
//3: (K-k)^2
//4: ignore chi, set K=1, stealth
//5: ignore chi, set K=1, (K-k)^2
//6: ignore chi, set K=1, (K-k)^6
//7: ignore chi, set K=1, V(k)=|1-k|
//9: ignore chi, V(k<1)=1, V of the shell right after 1 is -1e-4
//10: ignore chi, (S-S0)^2 type potential, V(any k)=1, S0(k<1)=0, S0 of the shell right after 1 is S0Out
//11: ignore chi, (S-S0)^2 type potential, V(k<1)=1, S0(k<1)=1, "superideal gas"
//12: ignore chi, (S-S0)^2 type potential, V(k<1)=1, input some constant, then S0(k<1)=constant, "equiluminous"
//13: ignore chi, (S-S0)^2 type potential, V(k<1)=1, S0(k<1)=k/2
//14: ignore chi, (S-S0)^2 type potential, input S, V(k<1)=1, S0(0)=0, S0(1)=S, S0 is linear
//16: ignore chi, (S-S0)^2 type potential, input S, V(k<1)=1, S0(0)=0, S0(1)=S, S0 is quadratic
//17: ignore chi, (S-S0)^2 type potential, input S, input n, V(k<1)=1, S0(0)=0, S0(1)=S, S0 \prop k^n
//15: ignore chi, (<S>-S0)^2 type potential, input S, input Nc, V(k<1)=1, S0(0)=0, S0(1)=S, S0 is linear, treat the configuration as Nc sub-configurations
//100: K2K2 potential
//101: input n, KM3Kn potential
ShiftedCCPotential * GeneratePotential(Configuration * pConfig, double & chi, size_t NumThreads, int Type, std::istream & ifile)
{
	size_t NumP=pConfig->NumParticle();
#ifdef HAVETWOFIXEDPARTICLES
	NumP+=2;
#endif
	DimensionType Dimension=pConfig->GetDimension();
	ShiftedCCPotential * pCCPotential=new ShiftedCCPotential(Dimension);

	size_t min_k_points = std::ceil(chi * NumP * Dimension);
	// We will iterate over
	// a hypersphere, which is circumscribing to an hypercube whose
	// volume is equal to the number of k-points.
	double max_k;
	if(Type<4)
	{
		max_k= (std::sqrt((double)(Dimension)) / 2) * std::pow( (2 * (double)(Dimension) * min_k_points) , (1.0/(double)(Dimension)) ) *(2*pi/std::pow(pConfig->PeriodicVolume(), 1.0/(double)(Dimension) ));

		std::vector<GeometryVector> Reciprocal;
		for(DimensionType i=0; i<Dimension; i++)
			Reciprocal.push_back(pConfig->GetReciprocalBasisVector(i));
		PeriodicCellList<Empty> lReciprocal(Dimension, &Reciprocal[0], std::sqrt(Reciprocal[0].Modulus2())*1.1, false);
		lReciprocal.Insert(Empty(), GeometryVector(Dimension));

		//generate list of possible ks
		std::vector<GeometryVector> kVecs;
		lReciprocal.IterateThroughNeighbors(GeometryVector(Dimension), max_k, [&Dimension, &kVecs] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
		{
			bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
			for(DimensionType i=0; i<Dimension; i++)
			{
				if(PeriodicShift[i]>0)
				{
					NonTrivial=true;
					break;
				}
				else if(PeriodicShift[i]<0)
				{
					NonTrivial=false;
					break;
				}
			}
			if(NonTrivial)
				kVecs.push_back(shift);
		});

		std::sort(kVecs.begin(), kVecs.end(), [] (const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2();});
		size_t CountKPoints=0;
		double LastNorm2=0;
		while(CountKPoints<kVecs.size())
		{
			if( ! (CountKPoints<min_k_points || std::abs((LastNorm2-kVecs[CountKPoints].Modulus2())/LastNorm2)<1e-10) )
				break;
			LastNorm2=kVecs[CountKPoints].Modulus2();
			pCCPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(kVecs[CountKPoints], 1.0));
			CountKPoints++;
		}
		if(Type==1)
		{
			double K2=kVecs[CountKPoints].Modulus2();
			for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
			{
				double k2=iter->k.Modulus2();
				double k_K=std::sqrt(k2/K2);
				iter->V=(std::acos(k_K)-k_K*std::sqrt(1-k2/K2));
			}
		}
		else if(Type==2)
		{
			double K2=kVecs[CountKPoints].Modulus2();
			for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
			{
				double k2=iter->k.Modulus2();
				double k_K=std::sqrt(k2/K2);
				iter->V=std::pow( (1.0-k_K), 6.0);
			}
		}
		else if(Type==3)
		{
			double K2=kVecs[CountKPoints].Modulus2();
			for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
			{
				double k2=iter->k.Modulus2();
				double k_K=std::sqrt(k2/K2);
				iter->V=std::pow( (1.0-k_K), 2.0);
			}
		}
	}
	else if (4<=Type && Type<9)
	{
		max_k=1.0;
		std::vector<GeometryVector> Reciprocal;
		for(DimensionType i=0; i<Dimension; i++)
			Reciprocal.push_back(pConfig->GetReciprocalBasisVector(i));
		PeriodicCellList<Empty> lReciprocal(Dimension, &Reciprocal[0], std::sqrt(Reciprocal[0].Modulus2())*1.1, false);
		lReciprocal.Insert(Empty(), GeometryVector(Dimension));

		//generate list of possible ks
		std::vector<GeometryVector> kVecs;
		lReciprocal.IterateThroughNeighbors(GeometryVector(Dimension), max_k, [&Dimension, &kVecs] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
		{
			bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
			for(DimensionType i=0; i<Dimension; i++)
			{
				if(PeriodicShift[i]>0)
				{
					NonTrivial=true;
					break;
				}
				else if(PeriodicShift[i]<0)
				{
					NonTrivial=false;
					break;
				}
			}
			if(NonTrivial)
				kVecs.push_back(shift);
		});

		std::sort(kVecs.begin(), kVecs.end(), [] (const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2();});
		for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			pCCPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));

		if(Type!=4)
		{
			//find |k| of the smallest UNCONSTRAINED ring, let it be K
			for(double maxk=1.2; ;maxk*=1.2)
			{
				kVecs.clear();
				lReciprocal.IterateThroughNeighbors(GeometryVector(Dimension), maxk, [&Dimension, &kVecs] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
				{
					bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
					for(DimensionType i=0; i<Dimension; i++)
					{
						if(PeriodicShift[i]>0)
						{
							NonTrivial=true;
							break;
						}
						else if(PeriodicShift[i]<0)
						{
							NonTrivial=false;
							break;
						}
					}
					if(NonTrivial)
						kVecs.push_back(shift);
				});
				if(kVecs.size()>pCCPotential->constraints.size())
					break;
			}
			std::nth_element(kVecs.begin(), kVecs.begin()+pCCPotential->constraints.size(), kVecs.end(), [] (const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2();});
			double K2=kVecs[pCCPotential->constraints.size()].Modulus2();

			if(Type==5)
			{
				for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
				{
					double k2=iter->k.Modulus2();
					double k_K=std::sqrt(k2/K2);
					iter->V=std::pow( (1.0-k_K), 2.0);
				}
			}
			if(Type==6)
			{
				for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
				{
					double k2=iter->k.Modulus2();
					double k_K=std::sqrt(k2/K2);
					iter->V=std::pow( (1.0-k_K), 6.0);
				}
			}
			if(Type==7)
			{
				for(auto iter=pCCPotential->constraints.begin(); iter!=pCCPotential->constraints.end(); iter++)
				{
					double km=std::sqrt(iter->k.Modulus2());
					iter->V=1.0-km;
				}
			}
		}

	}
	else if(Type>=9 && Type <=17)
	{
		std::vector<GeometryVector> Reciprocal;
		for(DimensionType i=0; i<Dimension; i++)
			Reciprocal.push_back(pConfig->GetReciprocalBasisVector(i));
		PeriodicCellList<Empty> lReciprocal(Dimension, &Reciprocal[0], std::sqrt(Reciprocal[0].Modulus2())*1.1, false);
		lReciprocal.Insert(Empty(), GeometryVector(Dimension));

		//generate list of possible ks
		std::vector<GeometryVector> kVecs;

		max_k=1.1;
		while( (kVecs.empty())? true : kVecs.back().Modulus2()<1.0 )
		{
			lReciprocal.IterateThroughNeighbors(GeometryVector(Dimension), max_k, [&Dimension, &kVecs] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
			{
				bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
				for(DimensionType i=0; i<Dimension; i++)
				{
					if(PeriodicShift[i]>0)
					{
						NonTrivial=true;
						break;
					}
					else if(PeriodicShift[i]<0)
					{
						NonTrivial=false;
						break;
					}
				}
				if(NonTrivial)
					kVecs.push_back(shift);
			});

			std::sort(kVecs.begin(), kVecs.end(), [] (const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2();});
			max_k*=1.1;
		}
		if(Type==9)
		{
			double prevK2=0.0;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( prevK2>1.0 && (iter->Modulus2() - prevK2)/prevK2 > 1e-9 )
					break;
				double V= (iter->Modulus2()>1.0) ? -1e-4 : 1;
				pCCPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, V));
				prevK2=iter->Modulus2();
			}
		}
		else if (Type==10)
		{
			double S0Out;
			ifile>>S0Out;
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			double prevK2=0.0;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( prevK2>1.0 && (iter->Modulus2() - prevK2)/prevK2 > 1e-9 )
					break;
				double s0= (iter->Modulus2()>1.0) ? S0Out : 0;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
				prevK2=iter->Modulus2();
			}
		}
		else if (Type==11)
		{
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( iter->Modulus2()>1 )
					break;
				double s0= 1;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type==12)
		{
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			double s0= 1;
			std::cout<<"S_0=";
			std::cin>>s0;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( iter->Modulus2()>1 )
					break;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type==13)
		{
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( iter->Modulus2()>1 )
					break;
				double s0= std::sqrt(iter->Modulus2())/2;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type == 14)
		{
			double S0;
			ifile >> S0;
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			for (auto iter = kVecs.begin(); iter != kVecs.end(); iter++)
			{
				if (iter->Modulus2()>1)
					break;
				double s0 = std::sqrt(iter->Modulus2())*S0;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type == 16)
		{
			double S0;
			ifile >> S0;
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			for (auto iter = kVecs.begin(); iter != kVecs.end(); iter++)
			{
				if (iter->Modulus2()>1)
					break;
				double s0 = iter->Modulus2()*S0;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type == 17)
		{
			double S0;
			ifile >> S0;
			double n;
			ifile >> n;
			delete pCCPotential;
			ShiftedCCPotential_S2 * pS2Potential = new ShiftedCCPotential_S2(Dimension);
			pCCPotential = pS2Potential;
			for (auto iter = kVecs.begin(); iter != kVecs.end(); iter++)
			{
				if (iter->Modulus2()>1)
					break;
				double s0 = std::pow(iter->Modulus2(), n*0.5)*S0;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
		else if (Type == 15)
		{
			double S0;
			ifile>>S0;
			size_t Nc;
			ifile>>Nc;
			delete pCCPotential;
			ShiftedCCPotential_S2_MultiConfig * pS2Potential = new ShiftedCCPotential_S2_MultiConfig(Dimension, Nc);
			pCCPotential = pS2Potential;
			for(auto iter=kVecs.begin(); iter!=kVecs.end(); iter++)
			{
				if( iter->Modulus2()>1 )
					break;
				double s0= std::sqrt(iter->Modulus2())*S0;
				pS2Potential->constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, 1.0));
				pS2Potential->S0.push_back(s0);
			}
		}
	}
	else if(Type==100)
	{
		delete pCCPotential;
		pCCPotential = new K2K2Potential(Dimension);
		pCCPotential->SetConfiguration(*pConfig);
	}
	else if (Type == 101)
	{
		double n;
		ifile >> n;
		delete pCCPotential;
		pCCPotential = new KM3KnPotential(Dimension, n);
		pCCPotential->SetConfiguration(*pConfig);
	}
	else
		std::cerr<<"Unknown Potential Type!\n";

	pCCPotential->ParallelNumThread=NumThreads;
	chi=(double)(pCCPotential->constraints.size())/((NumP-1)*Dimension);


	return pCCPotential;
}

//DensestReciprocalBox = true : use triangle box in 2D, and BCC box in 3D
//DensestReciprocalBox = false : use cubic box in any dimension
Configuration * GenerateConfig(size_t NumParticles, double boxSize, DimensionType dim, RandomGenerator & gen, bool DensestReciprocalBox)
{
	Configuration * result=nullptr;
	if(NumParticles<=0)
		std::cerr<<"Specify NumParticles before Generating Configuration!\n";
	else if(boxSize<=0)
		std::cerr<<"Specify BoxSize before Generating Configuration!\n";
	else if(dim<=0)
		std::cerr<<"Specify Dimension before Generating Configuration!\n";
	else
	{
		if(DensestReciprocalBox)
		{
			if(dim==2)
			{
				std::vector<GeometryVector> bases;
				bases.push_back(GeometryVector(boxSize, 0.0));
				bases.push_back(GeometryVector((-0.5)*boxSize, std::sqrt(3)*0.5*boxSize));
				result=new Configuration(2, &bases[0], boxSize, false);
			}
			else if(dim==3)
			{
				std::vector<GeometryVector> bases;
				bases.push_back(GeometryVector(boxSize, 0.0, 0.0));
				bases.push_back(GeometryVector(0.0, boxSize, 0.0));
				bases.push_back(GeometryVector(0.5*boxSize, 0.5*boxSize, 0.5*boxSize));
				result=new Configuration(3, &bases[0], boxSize, false);
			}
			else
			{
				std::cerr<<"Error in GenerateConfig : unsupported dimension!\n";
				exit(1);
			}
		}
		else
		{
			double CellSize=1.0/std::pow(NumParticles, 1.0/(double)(dim) );
			result=new Configuration(GetUnitCubicBox(dim, CellSize));
			result->Rescale(boxSize);
		}
#ifdef HAVETWOFIXEDPARTICLES
		for(size_t i=0; i<NumParticles-2; i++)
			result->Insert("A", gen);
#else
		for(size_t i=0; i<NumParticles; i++)
			result->Insert("A", gen);
#endif
	}
	return result;
}

/** \brief Perform the energy optimization multiple times for a given potential energy.
	Obtain ground-state configurations from multi initial configurations.
	@param[in]	pConfig	A pointer of a point configuration.
	@param[in]	pPotential	A pointer of a Potential object.
	@param[in]	gen	A random number generator.
	@param[in]	chi	The fraction of restricted degrees of freedom (when the potential is in the form of collective-coordiante potential).
	@param[in]	Prefix	A prefix of output files.
						It generates three distinct ConfigurationPacks;
						*_InitConfig	Initial configurations
						*				Configurations after minimizing energy
						*_Success		Configurations whose energy is lower than 1e-11 (a subset of *)
	@param[in]	SampleNumber	The number of initial configurations
	@param[in]	TimeLimit		When running time becomes larger than this, this function automatically ends.
	@param[in]	Algorithm	Name of alogrithm for optimizations. The followings are allowed; LBFGS, LocalGradientDescent, ConjugateGradient, SteepestDescent, MINOP
	@return 0 = exit without any problem. */
int CollectiveCoordinateMultiRun(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double chi, std::string Prefix, size_t SampleNumber, time_t TimeLimit, const std::string & AlgorithmName)
{
	(void)chi;
	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	ConfigurationPack AfterRelaxPack(Prefix);
	ConfigurationPack SuccessPack(Prefix+"_Success");
	ConfigurationPack InitConfigPack(Prefix+"_InitConfig");
	if(AfterRelaxPack.NumConfig()>0)
		std::cout<<"Found "<<AfterRelaxPack.NumConfig()<<" configurations, continue MultiRun.\n";

	std::cout<<"MultiRun start."<<'\n';
	for(size_t i=0; i<SampleNumber; i++)
	{
		std::cout<<"at time"<<std::time(nullptr)-ProgramStart<<"generating config"<<i;
		logfile<<"at time"<<std::time(nullptr)-ProgramStart<<"generating config"<<i;
		Configuration result(*pConfig);
		for(size_t i=0; i<Num; i++)
		{
			GeometryVector temp(dim);
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			result.MoveParticle(i, temp);
		}
		InitConfigPack.AddConfig(result);
		if(AlgorithmName=="LBFGS")
			RelaxStructure_NLOPT(result, *pPotential, 0.0, 0, 0.0, 10000);
		else if(AlgorithmName=="LocalGradientDescent")
			RelaxStructure_LocalGradientDescent(result, *pPotential, 0.0, 0, 0.0, 10000);
		else if(AlgorithmName=="ConjugateGradient")
			RelaxStructure_ConjugateGradient(result, *pPotential, 0.0, 0, 0.0, 10000);
		else if(AlgorithmName=="SteepestDescent")
			RelaxStructure_SteepestDescent(result, *pPotential, 0.0, 0, 0.0, 10000);
		else if(AlgorithmName=="MINOP")
			RelaxStructure_MINOP_withoutPertubation(result, *pPotential, 0.0, 0, 0.0, 10000);
		AfterRelaxPack.AddConfig(result);
		pPotential->SetConfiguration(result);
		double E=pPotential->Energy();
		std::cout<<" \tE_relax="<<E<<" \n";
		logfile<<" \tE_relax="<<E<<" \n";
		if(E<1e-11)
		{
			std::cout<<"Add to success pack\n";
			logfile<<"Add to success pack\n";
			SuccessPack.AddConfig(result);
		}
		if(std::time(nullptr) > TimeLimit)
		{
			//std::fstream ofile2("continue.txt", std::fstream::out);
			//ofile2<<"Not Completed. Please run again\n";
			return 0;
		}
	}
	return 0;
}

#include <MD_System.h>
int CollectiveCoordinateMD(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double TimeStep, double Temperature, std::string Prefix, size_t SampleNumber, size_t StepPerSample, bool QuenchAfterMD, bool AllowRestore, time_t TimeLimit, size_t EquilibrateSamples, bool MDAutoTimeStep)
{
	if(TimeStep==0)
	{
		std::cout<<"Input TimeStep!\n";
		std::cin>>TimeStep;
	}
	if(Temperature==0)
	{
		std::cout<<"Input Temperature!\n";
		std::cin>>Temperature;
	}
	std::cin.clear();
	if(Prefix=="")
	{
		std::cout<<"Input Prefix!\n";
		std::cin>>Prefix;
	}
	if(SampleNumber==0)
	{
		std::cout<<"Input SampleNumber!\n";
		std::cin>>SampleNumber;
	}
	if(StepPerSample==0)
	{
		std::cout<<"Step per sample not specified. Use default value 3000!\n";
		StepPerSample=3000;
	}
	ConfigurationPack BeforeRelaxPack(Prefix+std::string("_BeforeRelax"));
	ConfigurationPack AfterRelaxPack(Prefix);

	size_t OneTenthStepPerSample = StepPerSample/10;
	if(OneTenthStepPerSample==0)
		OneTenthStepPerSample=1;

	RelaxStructure_NLOPT(*pConfig, *pPotential, 0.0, 0, 0.0, 1000);
	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	size_t dimTensor=dim*Num; (void)dimTensor;
	Potential * pPot=pPotential;
	//unsigned long tempNumThread=1;//use this when doing MD because there is parallelization in MD code, so no need to parallelize in the potential
	pPot->SetConfiguration(* pConfig);
	double E=pPot->Energy(); (void)E;
	ParticleMolecularDynamics * psystem = NULL;

	bool Restart = false;
	signed char stage=0;
	long long step=0;
	if(AllowRestore)
	{
		std::fstream ifile("CCMD_Dump.MDDump", std::fstream::in | std::fstream::binary);
		if(ifile.good())
		{
			ifile.read( (char*)(&stage), sizeof(stage) );
			ifile.read( (char*)(&step), sizeof(step) );
			psystem = new ParticleMolecularDynamics(ifile);
			Restart=true;
		}
	}
	if(Restart == false)
	{
		BeforeRelaxPack.Clear();
		AfterRelaxPack.Clear();
		psystem = new ParticleMolecularDynamics(*pConfig, TimeStep, 1.0); 
	}

	size_t NumExistConfig=0; (void)NumExistConfig;
	std::cout<<"CCMD start. Temperature="<<Temperature<<'\n';
	logfile<<"CCMD start. Temperature="<<Temperature<<'\n';
	
	if(stage==0)
	{
		stage++;
		step=0;
	}
	if(stage==1)
	{
		//stage 1: equilibration after adjusting time step
		for(long long i=step; (size_t)i<EquilibrateSamples; i++)
		{
			if (std::time(nullptr) > TimeLimit || std::time(nullptr) > ::TimeLimit)
			{
				std::fstream ofile("CCMD_Dump.MDDump", std::fstream::out | std::fstream::binary);
				ofile.write( (char*)(&stage), sizeof(stage) );
				ofile.write( (char*)(&i), sizeof(i) );
				psystem->WriteBinary(ofile);
				delete psystem;
				std::fstream ofile2("continue.txt", std::fstream::out);
				ofile2<<"Not Completed. Please run again\n";
				return 0;
			}
			double c0=psystem->Position.GetCartesianCoordinates(0).x[0]; (void)c0;

			for(size_t ii=0; ii<2; ii++)
			{
				psystem->SetRandomSpeed(Temperature, gen);
				//std::swap(tempNumThread, pPot->ParallelNumThread);
				if(MDAutoTimeStep)
					psystem->Evolve_AutoTimeStep(StepPerSample/2, *pPot, 0.0001/SampleNumber);
				else
					psystem->Evolve(StepPerSample/2, *pPot);
				//std::swap(tempNumThread, pPot->ParallelNumThread);
			}

			pPot->SetConfiguration(psystem->Position);
			std::cout<<"1:"<<i<<"/"<<(EquilibrateSamples)<<", x0="<<psystem->Position.GetCartesianCoordinates(0).x[0]<<", Ep="<<pPot->Energy()<<", Ek="<<psystem->GetKineticEnergy()<<", dt="<<psystem->TimeStep<<'\n';;
			logfile<<"1:"<<i<<"/"<<(EquilibrateSamples)<<", x0="<<psystem->Position.GetCartesianCoordinates(0).x[0]<<", Ep="<<pPot->Energy()<<", Ek="<<psystem->GetKineticEnergy()<<", dt="<<psystem->TimeStep<<'\n';;
			std::cout.flush();
		}
		stage++;
		step=0;
	}

	//stage 2: sample
	for(long long i=step; (size_t)i<SampleNumber; i++)
	{
		if (std::time(nullptr) > TimeLimit || std::time(nullptr) > ::TimeLimit)
		{
			std::fstream ofile("CCMD_Dump.MDDump", std::fstream::out | std::fstream::binary);
			ofile.write( (char*)(&stage), sizeof(stage) );
			ofile.write( (char*)(&i), sizeof(i) );
			psystem->WriteBinary(ofile);
			delete psystem;
			std::fstream ofile2("continue.txt", std::fstream::out);
			ofile2<<"Not Completed. Please run again\n";
			return 0;
		}
		std::cout<<"at time"<<std::time(nullptr)-ProgramStart;
		logfile<<"at time"<<std::time(nullptr)-ProgramStart;
		//std::swap(tempNumThread, pPot->ParallelNumThread);
		psystem->AndersonEvolve(StepPerSample/2, *pPot, Temperature, 0.01, gen);
		psystem->Evolve(StepPerSample/2, *pPot);
		//std::swap(tempNumThread, pPot->ParallelNumThread);
		Configuration result(psystem->Position);
		if(QuenchAfterMD)
			RelaxStructure_NLOPT(result, *pPot, 0.0, 0, 0.0, 100000);
		pPot->SetConfiguration(result);
		std::cout<<", 2:"<<i<<"/"<<(SampleNumber)<<" \tE_relax="<<pPot->Energy()<<" \t";
		logfile<<", 2:"<<i<<"/"<<(SampleNumber)<<" \tE_relax="<<pPot->Energy()<<" \t";
		pPot->SetConfiguration(psystem->Position);
		std::cout<<"E_p="<<pPot->Energy()<<" \t";
		logfile<<"E_p="<<pPot->Energy()<<" \t";
		std::cout<<"E_k="<<psystem->GetKineticEnergy()<<", dt="<<psystem->TimeStep<<'\n';
		logfile<<"E_k="<<psystem->GetKineticEnergy()<<", dt="<<psystem->TimeStep<<'\n';

		AfterRelaxPack.AddConfig(result);
		if(QuenchAfterMD)
			BeforeRelaxPack.AddConfig(psystem->Position);

		std::cout.flush();
	}
	::Output("FinalConfiguration", psystem->Position);
	delete psystem;

	return 0;
}
int CCMDAnneal(Configuration * pConfig, Potential * pPotential, RandomGenerator & /*gen*/, double AnnealingSchedule)
{
	//LogCCPotential pot(*dynamic_cast<ShiftedCCPotential *>(pPotential));
	Potential & pot = (*pPotential);
	ParticleMDAnnealing an(*pConfig, 12345, 0.0, pot, false, &std::cout, AnnealingSchedule); 
	std::cout<<"T_init=";
	std::cin>>an.Tinit;
	std::cout<<"T_cut=";
	std::cin>>an.Tcut;
	double TimeStep=nan("");
	std::cout<<"Time step(or auto)=";
	std::cin>>TimeStep;
	if(!isnan(TimeStep) && TimeStep>0 )
	{
		an.AutoTimeStep=false;
		an.pSys->TimeStep=TimeStep;
	}
	else
		an.pSys->TimeStep = 0.001;
	an.SameLocalMinimaDetection=false;
	an.Anneal(pot);

	::Verbosity=4;
	Configuration c(*an.pList);
	::RelaxStructure_NLOPT(c, pot, 0.0, 0, 0.0, 100000);

	::AlsoWriteEPS=true;
	::Plot("CCMDAnneal_Result", c);
	if(pConfig->GetDimension()==2)
		::PlotDirectionalStructureFactor2D(c, 40, "CCMDAnneal_Result");
	::Output("CCMDAnneal_Result", c);

	(*pConfig)=c;

	return 0;
}
//Serial code to calculate the eigen of 1 chi value
int CalculateEigen(double chi, size_t SideLength, const std::vector<std::string> & Structures)
{
	void GetHessianEigens(Configuration * pConfig, ShiftedCCPotential * pPotential, std::vector<double> & results, std::vector< std::vector<double> > * pEigenVectors = nullptr );
	
	for(auto iter=Structures.begin(); iter!=Structures.end(); iter++)
	{
		Configuration c=MultiplicateStructure(ReadPos(iter->c_str()), SideLength);
		double ch=chi;
		ShiftedCCPotential * ppot=GeneratePotential(&c, ch, 4, 0, std::cin);
		std::vector<double> eig;
		GetHessianEigens(&c, ppot, eig);

		size_t ZeroCount=0;
		for(auto iter=eig.begin(); iter!=eig.end(); iter++)
		{
			if(std::abs(*iter)<1e-5)
				ZeroCount++;
		}
		std::cout<<*iter<<" \t"<<ch<<" \t"<<ZeroCount<<" \n";

	}
	return 0;

}

#ifdef USE_EIGEN
#include <Eigen/Dense>
//only works for stealth potentials
//if pEigenVectors is not nullptr, then fill it with eigenvectors
void GetHessianEigens(Configuration * pConfig, ShiftedCCPotential * pPot, std::vector<double> & results, std::vector< std::vector<double> > * pEigenVectors = nullptr )
{
	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	size_t dimTensor=dim*Num;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian(dimTensor, dimTensor);

	pPot->SetConfiguration(*pConfig);
	pPot->GetRho();

//#pragma omp parallel for schedule(dynamic) num_threads(3)
	for(signed long i=0; i<Num; i++)
	{
		for(size_t j=0; j<=i; j++)
		{
			auto iter2=pPot->RhoReal.begin();
			auto iter3=pPot->RhoImag.begin();

			for(auto iter=pPot->constraints.begin(); iter!=pPot->constraints.end(); iter++, iter2++, iter3++)
			{
				GeometryVector ri=pConfig->GetCartesianCoordinates(i);
				GeometryVector rij=ri-pConfig->GetCartesianCoordinates(j);
				double coeff=std::cos(iter->k.Dot(rij));
				if(i==j)
					coeff-= ( (*iter2)*std::cos(iter->k.Dot(ri))+(*iter3)*std::sin(iter->k.Dot(ri)) );
				coeff=2*iter->V*coeff;

				for(DimensionType di=0; di<dim; di++)
					for(DimensionType dj=0; dj<dim; dj++)
						(Hessian(i*dim+di, j*dim+dj))+=coeff*iter->k.x[di]*iter->k.x[dj]; 

				if(i!=j)
					for(DimensionType di=0; di<dim; di++)
						for(DimensionType dj=0; dj<dim; dj++)
							(Hessian(j*dim+dj, i*dim+di))+=coeff*iter->k.x[di]*iter->k.x[dj]; 
			}
		}
	}

	Eigen::SelfAdjointEigenSolver<decltype(Hessian)> Solver(Hessian, (pEigenVectors==nullptr)?Eigen::EigenvaluesOnly:Eigen::ComputeEigenvectors);
	if (Solver.info() != Eigen::Success)
	{
		std::cout<<"Error in GetHessianEigens : Eigen failed\n";
		return ;
	}

	results.clear();
	auto eigval = Solver.eigenvalues();
	if(pEigenVectors!=nullptr)
	{
		pEigenVectors->clear();
		Hessian=Solver.eigenvectors();
	}

	for(size_t i=0; i<dimTensor; i++)
	{
		double dtemp=eigval(i);
		results.push_back(dtemp);
		if(pEigenVectors!=nullptr)
		{
			pEigenVectors->push_back(std::vector<double>());
			for(size_t j=0; j<dimTensor; j++)
				pEigenVectors->back().push_back( Hessian(j, i) );
		}
	}
}
void CCNewtonStep(Configuration * pConfig, ShiftedCCPotential * pPot, double * pStep=nullptr, double Gamma=1.0)
{
	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	size_t dimTensor=dim*Num;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian(dimTensor, dimTensor);

	pPot->SetConfiguration(*pConfig);
	pPot->GetRho();

//#pragma omp parallel for schedule(dynamic) num_threads(3)
	for(signed long i=0; i<Num; i++)
	{
		for(size_t j=0; j<=i; j++)
		{
			auto iter2=pPot->RhoReal.begin();
			auto iter3=pPot->RhoImag.begin();

			for(auto iter=pPot->constraints.begin(); iter!=pPot->constraints.end(); iter++, iter2++, iter3++)
			{
				GeometryVector ri=pConfig->GetCartesianCoordinates(i);
				GeometryVector rij=ri-pConfig->GetCartesianCoordinates(j);
				double coeff=std::cos(iter->k.Dot(rij));
				if(i==j)
					coeff-= ( (*iter2)*std::cos(iter->k.Dot(ri))+(*iter3)*std::sin(iter->k.Dot(ri)) );
				coeff=2*iter->V*coeff/pConfig->PeriodicVolume();

				for(DimensionType di=0; di<dim; di++)
					for(DimensionType dj=0; dj<dim; dj++)
						(Hessian(i*dim+di, j*dim+dj))+=coeff*iter->k.x[di]*iter->k.x[dj]; 

				if(i!=j)
					for(DimensionType di=0; di<dim; di++)
						for(DimensionType dj=0; dj<dim; dj++)
							(Hessian(j*dim+dj, i*dim+di))+=coeff*iter->k.x[di]*iter->k.x[dj]; 
			}
		}
	}

	Eigen::Matrix<double, Eigen::Dynamic, 1> NegativeDerivative(dimTensor);
	for(int i=0; i<Num; i++)
	{
		GeometryVector f;
		pPot->Force(f, i);
		for(int j=0; j<dim; j++)
			(NegativeDerivative(i*dim+j))=f.x[j];
	}
	//debug temp
	//{
	//	Configuration c2(*pConfig);
	//	GeometryVector r=c2.GetCartesianCoordinates(0);
	//	r.x[0]+=1e-6;
	//	r=c2.CartesianCoord2RelativeCoord(r);
	//	c2.MoveParticle(0, r);
	//	pPot->SetConfiguration(c2);
	//	for(int i=0; i<Num; i++)
	//	{
	//		GeometryVector f;
	//		pPot->Force(f, i);
	//		for(int j=0; j<dim; j++)
	//			std::cout<<Hessian(0, i*dim+j)<<" \t"<<(f.x[j]-NegativeDerivative(i*dim+j))/1e-6<<" \n";
	//	}
	//	exit(0);
	//}

	Eigen::Matrix<double, Eigen::Dynamic, 1> Step = Hessian.fullPivLu().solve(NegativeDerivative);
	if(pStep!=nullptr)
		(*pStep) = Step.norm();

	for(int i=0; i<Num; i++)
	{
		GeometryVector r=pConfig->GetCartesianCoordinates(i);
		for(int j=0; j<dim; j++)
			r.x[j]+=(Step(i*dim+j))*Gamma;
		r=pConfig->CartesianCoord2RelativeCoord(r);
		pConfig->MoveParticle(i, r);
	}

	if(Verbosity>=5)
	{
		pPot->SetConfiguration(*pConfig);
		std::cout<<"f="<<pPot->Energy()<<"|g|="<<NegativeDerivative.norm()<<"|s|="<<Step.norm()<<"\n";
	}

}
void CCNewtonOptimization(Configuration * pConfig, ShiftedCCPotential * pPot)
{
	const double MinStepSize=1e-15;
	double StepSize=1.0;
	double f;
	pPot->SetConfiguration(*pConfig);
	f=pPot->Energy();
	double gamma=1.0;
	while(StepSize>MinStepSize)
	{
		CCNewtonStep(pConfig, pPot, &StepSize, gamma);
		pPot->SetConfiguration(*pConfig);
		double e2=pPot->Energy();
		if(e2>f) gamma*=0.5;
		f=e2;
	}
}
#else
void GetHessianEigens(Configuration * /*pConfig*/, ShiftedCCPotential * /*pPot*/, std::vector<double> & /*results*/, std::vector< std::vector<double> > * /*pEigenVectors*/ = nullptr )
{
	std::cerr<<"Error in GetHessianEigens : Eigen not enabled!\n";
}
void CCNewtonOptimization(Configuration * /*pConfig*/, ShiftedCCPotential * /*pPot*/)
{
	std::cerr<<"Error in CCNewtonOptimization : Eigen not enabled!\n";
}
#endif




int CC_WLMC();
int CollectiveCoordinateCLI()
{

	Configuration * pConfig=nullptr;
	ShiftedCCPotential * pPotential=nullptr;
	size_t NumParticles = 0;
	DimensionType dim = 0;
	double boxSize = 1.0;
	int NumThreads=1;
	bool QuenchAfterMD = true;
	RandomGenerator gen;
	double chi=0.0;
	//////////////////
	double MDTimeStep=0, MDTemperature=0;
	size_t MDEquilibrateSamples=20;
	std::string MDPrefix="";
	size_t MDSampleNumber=0, MDStepPerSample=0;
	bool MDAllowRestore = false;
	bool MDAutoTimeStep = true;
	time_t MDTimeLimit = std::time(nullptr) + 8640000; //default time limit of 100 days
	//////////////////
	char tempstring[1000] = {};
	double AnnealingSchedule=1.0;
	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;


	for(;;)
	{
		ifile>>tempstring;
		if(strcmp(tempstring, "Exit")==0)
		{
			if(pPotential!=nullptr)
				delete pPotential;
			if(pConfig!=nullptr)
				delete pConfig;
			return 0;
		}
		else if(strcmp(tempstring, "MDTimeLimit")==0)
		{
			size_t deltaTime;
			ifile>>deltaTime;
			MDTimeLimit=std::time(nullptr)+deltaTime;
		}
		else if(strcmp(tempstring, "MDAllowRestore")==0)
		{
			ifile>>MDAllowRestore;
		}
		else if(strcmp(tempstring, "MDAutoTimeStep")==0)
		{
			ifile>>MDAutoTimeStep;
		}
		else if(strcmp(tempstring, "MDEquilibrateSamples")==0)
		{
			ifile>>MDEquilibrateSamples;
		}
		else if(strcmp(tempstring, "NumParticles")==0 || strcmp(tempstring, "NbrParticles")==0)
		{
			ifile>>NumParticles;
		}
		else if(strcmp(tempstring, "CoolingSchedule")==0 || strcmp(tempstring, "AnnealingSchedule")==0)
		{
			ifile>>AnnealingSchedule;
		}
		else if(strcmp(tempstring, "Anneal")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
				std::cerr<<"Initialize Configuration and Potential before annealing!\n";
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy before="<<pPotential->Energy()<<'\n';
				CCMDAnneal(pConfig, pPotential, gen, AnnealingSchedule);
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy after="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "NumThreads")==0 || strcmp(tempstring, "NbrThreads")==0)
		{
			ifile>>NumThreads;
			if(NumThreads<=0)
			{
				std::cerr<<"NumThreads<=0 does not make sense. Change it to NumThreads=1!\n";
				NumThreads=1;
			}
		}
		else if(strcmp(tempstring, "Dimension")==0)
		{
			ifile>>dim;
		}
		//else if(strcmp(tempstring, "CalculateEigen")==0)
		//{
		//	size_t SideLength;
		//	std::cout<<"SideLength=";
		//	std::cin>>SideLength;
		//	CalculateEigen(chi, SideLength);
		//}
		else if(strcmp(tempstring, "BoxSize")==0)
		{
			ifile>>boxSize;
		}
		else if(strcmp(tempstring, "Chi")==0)
		{
			ifile>>chi;
		}
		else if(strcmp(tempstring, "AutoRescale")==0)
		{
			//Rescale the configuration so that K=1
			//Do this before GeneratePotential

			//defined in PairCorrelation.cpp
			double HyperSphere_Volume(DimensionType n, double R);

			double rescale=2*pi/boxSize*std::pow(2*dim*chi*NumParticles/HyperSphere_Volume(dim, 1.0), 1.0/dim);

			if(pConfig==nullptr)
				std::cerr<<"Initialize Configuration before AutoRescale!\n";
			else
			{
				pConfig->Rescale(rescale);
				boxSize*=rescale;
			}
		}
		else if(strcmp(tempstring, "RandomSeed")==0)
		{
			int RandomSeed;
			ifile>>RandomSeed;
			gen.seed(RandomSeed);
		}
		else if(strcmp(tempstring, "ReadProblemFile")==0)
		{
			ifile>>tempstring;
			ReadProblemFile(tempstring, NumParticles, dim, boxSize, pConfig, pPotential, NumThreads);
		}
		else if(strcmp(tempstring, "WangLandauMonteCarlo")==0)
		{
			CC_WLMC();
		}
		else if(strcmp(tempstring, "ReadConfiguration")==0)
		{
			ifile>>tempstring;
			if(pConfig!=nullptr)
				delete pConfig;
			pConfig=new Configuration(ReadPos(tempstring));
#ifdef HAVETWOFIXEDPARTICLES
			pConfig->DeleteParticle(pConfig->NumParticle()-1);
			pConfig->DeleteParticle(pConfig->NumParticle()-1);
#endif
		}
		else if(strcmp(tempstring, "GenerateConfiguration")==0)
		{
			Configuration * pNewConfig=GenerateConfig(NumParticles, boxSize, dim, gen, false);
			if(pNewConfig!=nullptr)
			{
				if(pConfig!=nullptr)
					delete pConfig;
				pConfig=pNewConfig;
			}
		}
		else if(strcmp(tempstring, "LoadConfigPack")==0)
		{
			if(pConfig!=nullptr)
				delete pConfig;

			ofile << "File name = ";
			ifile >> tempstring;
			ConfigurationPack c(tempstring);
			ofile << "It contains "<<c.NumConfig()<<" Configurations \n";
			ofile << "Index = ";
			ifile >> tempstring;

			Configuration result = c.GetConfig(atoi(tempstring));
			WriteConfiguration(result,"test");
			pConfig= new Configuration(result);

		}
		else if(strcmp(tempstring, "GenerateLatticeConfiguration")==0)
		{
			//pConfig=GenerateConfig(NumParticles, boxSize, dim, gen, false);
			std::cout<<"Input lattice filename:";
			std::string latticename;
			std::cin>>latticename;
			Configuration l=ReadPos(latticename);
			size_t NumPerSide = std::floor( std::pow ( NumParticles/l.NumParticle(), 1.0/l.GetDimension() ) +0.5 );
			if( l.NumParticle()*std::pow(NumPerSide, l.GetDimension()) == NumParticles )
			{
				if(pConfig!=nullptr)
					delete pConfig;
				pConfig = new Configuration(MultiplicateStructure(l, NumPerSide));
				double factor = boxSize/std::sqrt(pConfig->GetBasisVector(0).Modulus2());
				pConfig->Rescale(factor);
			}
			else
			{
				std::cerr<<"Error in GenerateLatticeConfiguration : cannot find appropriate side length.\n";
			}
		}
		else if(strcmp(tempstring, "GenerateDensestReciprocalConfiguration")==0)
		{
			Configuration * pNewConfig=GenerateConfig(NumParticles, boxSize, dim, gen, true);
			if(pNewConfig!=nullptr)
			{
				if(pConfig!=nullptr)
					delete pConfig;
				pConfig=pNewConfig;
			}
		}
		else if(strcmp(tempstring, "GenerateFCCConfiguration")==0)
		{
			if(dim!=3)
				std::cout<<"GenerateFCCConfiguration only works for 3 dimensions!\n";
			else
			{
				if(pConfig!=nullptr)
					delete pConfig;

				std::vector<GeometryVector> bases;
				bases.push_back(GeometryVector(boxSize, 0.0, boxSize));
				bases.push_back(GeometryVector(boxSize, boxSize, 0.0));
				bases.push_back(GeometryVector(0.0, boxSize, boxSize));
				pConfig=new Configuration(3, &bases[0], boxSize, false);
				for(size_t i=0; i<NumParticles; i++)
					pConfig->Insert("A", gen);
			}
		}
		else if(strcmp(tempstring, "GeneratePotential")==0)
		{
			if(pConfig==nullptr)
			{
				std::cerr<<"Initialize Configuration before Generating Potential!\n";
			}
			else
			{
				if(chi==0.0)
				{
					ofile<<"Enter Expected chi:";
					ifile>>chi;
				}

				if(pPotential!=nullptr)
					delete pPotential;
				pPotential=GeneratePotential(pConfig, chi, NumThreads, 0, ifile);
				std::cout<<"Num of Constraints:"<<pPotential->constraints.size()<<'\n';
				std::cout<<"Final chi="<<chi<<'\n';
				std::cout<<"Constrained k points are:\n------------------------------------------\nweight \t k\n";
				for(auto iter=pPotential->constraints.begin(); iter!=pPotential->constraints.end(); iter++)
				{
					std::cout<<iter->V<<" \t";
					iter->k.OutputCoordinate(std::cout, pPotential->Dimension);
					std::cout<<'\n';
				}
				std::cout<<"------------------------------------------\n";
			}
		}
		else if(strcmp(tempstring, "GenerateOtherPotential")==0)
		{
			if(pConfig==nullptr)
			{
				std::cerr<<"Initialize Configuration before Generating Potential!\n";
			}
			else
			{
				ofile<<"Enter Expected chi:";
				ifile>>chi;
				ofile<<"Enter type:";
				int type;
				ifile>>type;

				if(pPotential!=nullptr)
					delete pPotential;
				pPotential=GeneratePotential(pConfig, chi, NumThreads, type, ifile);
				std::cout<<"Num of Constraints:"<<pPotential->constraints.size()<<'\n';
				std::cout<<"Final chi="<<chi<<'\n';
				std::cout<<"Constrained k points are:\n------------------------------------------\nweight \t k\n";
				if(Verbosity>3)
				for(auto iter=pPotential->constraints.begin(); iter!=pPotential->constraints.end(); iter++)
				{
					std::cout<<iter->V<<" \t";
					iter->k.OutputCoordinate(std::cout, pPotential->Dimension);
					std::cout<<'\n';
				}
				std::cout<<"------------------------------------------\n";
			}
		}
		else if(strcmp(tempstring, "CalculateEnergy")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "MultiRun")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				std::string algorithm;
				std::cout<<"Algorithm=";
				std::cin>>algorithm;
				CollectiveCoordinateMultiRun(pConfig, pPotential, gen, chi, MDPrefix, MDSampleNumber, MDTimeLimit, algorithm);
			}
		}
		else if(strcmp(tempstring, "Run")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy before="<<pPotential->Energy()<<'\n';
				RelaxStructure_NLOPT(* pConfig, * pPotential, 0.0, 0, 0.0, 5000);
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy after="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "RunNewton")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy before="<<pPotential->Energy()<<'\n';
				CCNewtonOptimization(pConfig, pPotential);
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy after="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "RunGSL")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy before="<<pPotential->Energy()<<'\n';
				RelaxStructure_ConjugateGradient(* pConfig, * pPotential, 0.0, 0, 0.0);
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy after="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "RunMINOP")==0)
		{
			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy before="<<pPotential->Energy()<<'\n';
				RelaxStructure_MINOP_withoutPertubation(* pConfig, * pPotential, 0.0, 0, 0.0);
				pPotential->SetConfiguration(*pConfig);
				std::cout<<"Energy after="<<pPotential->Energy()<<'\n';
			}
		}
		else if(strcmp(tempstring, "MDTimeStep")==0)
		{
			std::cin>>MDTimeStep;
		}
		else if(strcmp(tempstring, "MDTemperature")==0)
		{
			std::cin>>MDTemperature;
		}
		else if(strcmp(tempstring, "MDPrefix")==0)
		{
			std::cin>>MDPrefix;
		}
		else if(strcmp(tempstring, "Prefix")==0)
		{
			std::cin>>MDPrefix;
		}
		else if(strcmp(tempstring, "MDSampleNumber")==0)
		{
			std::cin>>MDSampleNumber;
		}
		else if(strcmp(tempstring, "SampleNumber")==0)
		{
			std::cin>>MDSampleNumber;
		}
		else if(strcmp(tempstring, "MDStepPerSample")==0)
		{
			std::cin>>MDStepPerSample;
		}
		else if(strcmp(tempstring, "QuenchAfterMD")==0)
		{
			std::cin>>QuenchAfterMD;
		}
		else if(strcmp(tempstring, "MD")==0)
		{
			int temp=1;

			if(pConfig==nullptr || pPotential==nullptr)
			{
				std::cerr<<"Initialize Configuration and Potential before running!\n";
			}
			else
			{
				std::swap(temp, NumThreads);
				CollectiveCoordinateMD(pConfig, pPotential, gen, MDTimeStep, MDTemperature, MDPrefix, MDSampleNumber, MDStepPerSample, QuenchAfterMD, MDAllowRestore, MDTimeLimit, MDEquilibrateSamples, MDAutoTimeStep);
				std::swap(temp, NumThreads);
			}
		}
		else if(strcmp(tempstring, "Output")==0)
		{
			if(pConfig==nullptr)
			{
				std::cerr<<"Initialize Configuration before outputting!\n";
			}
			else
			{
				ifile>>tempstring;
	#ifdef HAVETWOFIXEDPARTICLES
				pConfig->Insert("B", GeometryVector(dim));
				pConfig->Insert("B", GeometryVector(dim));
	#endif
				::Plot(tempstring, * pConfig);
				::Output(tempstring, * pConfig);
			}
		}
		else
		{
			ofile<<"Unrecognized command!\n";
			ifile.clear();
		}
		ifile.ignore(1000,'\n');
		if(ifile.eof()) return 2;
	}
	return 1;

}


//deprecated. Need to use Eigen instead of GSL
//int CollectiveCoordinateEigenMove(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen)
//{
//	DimensionType dim=pConfig->GetDimension();
//	size_t Num=pConfig->NumParticle();
//	size_t dimTensor=dim*Num;
//	ShiftedCCPotential * pPot= dynamic_cast<ShiftedCCPotential *>(pPotential);
//	gsl_matrix * S= gsl_matrix_calloc(dimTensor, dimTensor);
//	gsl_matrix * evals= gsl_matrix_calloc(dimTensor, dimTensor);
//	gsl_vector * deform= gsl_vector_calloc(dimTensor);
//	gsl_vector * PrevDeform= gsl_vector_calloc(dimTensor);
//	gsl_eigen_symmv_workspace * wor = gsl_eigen_symmv_alloc(dimTensor);
//	gsl_vector * eig = gsl_vector_calloc(dimTensor);
//	gsl_vector * vtemp = gsl_vector_calloc(dimTensor);
//
//	ConfigurationPack cPack("EigenMove");
//
//	for(size_t OCount=0; OCount<200; OCount++)
//	{
//		RelaxStructure_NLOPT(*pConfig, *pPotential, 0.0, 0, 0.0);
//		std::stringstream name;
//		std::string name2;
//		name<<"EigenMove/"<<OCount<<"_test";
//		name>>name2;
//		::Output(name2, *pConfig);
//		::Plot(name2, *pConfig);
//		cPack.AddConfig(*pConfig);
//		pPotential->SetConfiguration(*pConfig);
//		std::cout<<"Configuration "<<OCount<<" Completed!, E="<<pPotential->Energy()<<", Write as:"<<name2<<'\n';
//		for(size_t ICount=0; ICount<100; ICount++)
//		{
//			RelaxStructure_NLOPT_Emin(*pConfig, *pPotential, 0.0, 0, 0.0, 10000, 1e-6);
//			pPot->SetConfiguration(*pConfig);
//			auto iter2=pPot->RhoReal.begin();
//			auto iter3=pPot->RhoImag.begin();
//			
//			gsl_matrix_set_zero(S);
//			for(auto iter=pPot->constraints.begin(); iter!=pPot->constraints.end(); iter++, iter2++, iter3++)
//			{
//				for(size_t i=0; i<Num; i++)
//				{
//					for(size_t j=0; j<=i; j++)
//					{
//						GeometryVector ri=pConfig->GetCartesianCoordinates(i);
//						GeometryVector rij=ri-pConfig->GetCartesianCoordinates(j);
//						double coeff=std::cos(iter->k.Dot(rij));
//						if(i==j)
//							coeff-= ( (*iter2)*std::cos(iter->k.Dot(ri))+(*iter3)*std::sin(iter->k.Dot(ri)) );
//						coeff=2*iter->V*coeff;
//
//						for(DimensionType di=0; di<dim; di++)
//							for(DimensionType dj=0; dj<dim; dj++)
//								(*gsl_matrix_ptr(S, i*dim+di, j*dim+dj))+=coeff*iter->k.x[di]*iter->k.x[dj]; 
//
//						if(i!=j)
//							for(DimensionType di=0; di<dim; di++)
//								for(DimensionType dj=0; dj<dim; dj++)
//									(*gsl_matrix_ptr(S, j*dim+dj, i*dim+di))+=coeff*iter->k.x[di]*iter->k.x[dj]; 
//					}
//				}
//			}
//
//			gsl_eigen_symmv(S, eig, evals, wor);
//
//			gsl_vector_set_zero(deform);
//			//if(ICount==0 && OCount==0)
//			if(ICount==0)
//			//if(true)
//			{
//				//pick up a random direction
//				for(size_t i=0; i<dimTensor; i++)
//				{
//					double dtemp=gsl_vector_get(eig, i);
//					//std::cout<<dtemp<<" \t";
//					if(std::abs(dtemp)<1e-6)
//					{
//						gsl_matrix_get_col(vtemp, evals, i);
//						gsl_vector_scale(vtemp, gen.RandomDouble());
//						gsl_vector_add(deform, vtemp);
//						//break;
//					}
//				}
//			}
//			else
//			{
//				//pick up a direction that is as close to PrevDeform as possible
//				double MaxDot=0.0;
//				for(size_t i=0; i<dimTensor; i++)
//				{
//					double dtemp=gsl_vector_get(eig, i);
//					//std::cout<<dtemp<<" \t";
//					if(dtemp<1e-6)
//					{
//						gsl_matrix_get_col(vtemp, evals, i);
//						double dtemp2;
//						gsl_blas_ddot(vtemp, PrevDeform, &dtemp2);
//						dtemp2=dtemp2/gsl_blas_dnrm2(vtemp);
//						if(std::abs(dtemp2)>MaxDot)
//						{
//							MaxDot=std::abs(dtemp2);
//							gsl_vector_memcpy(deform, vtemp);
//							if(dtemp2<0.0)
//								gsl_vector_scale(deform, (-1.0));
//						}
//					}
//				}
//			}
//			if(gsl_vector_max(deform)==0.0)
//			{
//				std::cerr<<"Error in EigenMove : Unable to find move direction!\n";
//				std::cerr<<"i="<<ICount<<'\n';
//				for(size_t i=0; i<dimTensor; i++)
//				{
//					double dtemp=gsl_vector_get(eig, i);
//					std::cout<<dtemp<<" \t";
//				}
//				std::cout<<'\n';
//				assert(false);
//			}
//			
//			//for(size_t i=0; i<deform->size; i++)
//			//	(*gsl_vector_ptr(deform, i))=gen.RandomDouble();
//
//
//			//for(size_t i=0; i<dimTensor; i++)
//			//	std::cout<<gsl_vector_get(deform, i)<<" \t";
//			//std::cout<<'\n';
//			for(double i=1e-9; i<0.5; i*=2)
//			{
//				Configuration conf(*pConfig);
//				for(size_t j=0; j<Num; j++)
//				{
//					GeometryVector rj=conf.GetRelativeCoordinates(j);
//					for(DimensionType k=0; k<dim; k++)
//						rj.x[k]+=i*gsl_vector_get(deform, j*dim+k);
//
//					conf.MoveParticle(j, rj);
//				}
//
//				pPot->SetConfiguration(conf);
//				//std::cout<<i<<" \t"<<i*gsl_vector_max(deform)<<" \t"<<pPot->Energy()<<'\n';
//				if(pPot->Energy()>1e-4)
//				{
//					(*pConfig)=conf;
//					break;
//				}
//			}
//			gsl_vector_memcpy(PrevDeform, deform);
//		}
//	}
//
//	gsl_vector_free(eig);
//	gsl_vector_free(vtemp);
//	gsl_vector_free(deform);
//	gsl_eigen_symmv_free(wor);
//	gsl_matrix_free(S);
//	gsl_matrix_free(evals);
//	return 0;
//}

int RescaleEnergy()
{
	Configuration original = ReadPos("temp3");
	double N=original.NumParticle();

	for(double V=300; V<330; V+=0.1)
	{
		Configuration temp(original);
		temp.Rescale(V/temp.PeriodicVolume());

		double dtemp=0.0;
		ShiftedCCPotential * pPot=GeneratePotential(&temp, dtemp, 4, 7, std::cin);

		pPot->SetConfiguration(temp);
		double Eprime=pPot->Energy();
		double SumVk=0.0;
		for(auto iter=pPot->constraints.begin(); iter!=pPot->constraints.end(); iter++)
			SumVk+=iter->V;
		double E=Eprime/2.0/V-N/2.0/V*SumVk+N*(N-1)/2.0/V;
		//std::cout<<V<<" \t"<<Eprime<<" \t"<<SumVk<<" \t"<<E<<'\n';
		std::cout<<V<<" \t"<<E<<'\n';

		delete pPot;
	}

	return 0;
}
//{
//	//gradually increase chi, once it is above 0.5, it is trapped in local minimum!
//	for(double chi=0.45; chi<0.65; chi+=0.01)
//	{
//		double chi2=chi;
//		Potential * pPot = GeneratePotential(pConfig, chi2, 4, 0, std::cin);
//		RelaxStructure_NLOPT(*pConfig, *pPot, 0.0, 0, 0.0, 10000);
//		pPot->SetConfiguration(*pConfig);
//		std::cout<<chi<<" \t"<<pPot->Energy()<<'\n';
//		delete pPot;
//	}
//	return 0;
//}
//{
//	//add larger k first, not working either
//	double chi=0.65;
//	ShiftedCCPotential * pPot = GeneratePotential(pConfig, chi, 4, 0, std::cin);
//	std::vector<ShiftedCCPotential::KPointConstraint> KPs;
//	std::swap(KPs, pPot->constraints);
//	std::sort(KPs.begin(), KPs.end(), [] (const ShiftedCCPotential::KPointConstraint & left, const ShiftedCCPotential::KPointConstraint & right) ->bool
//	{
//		return left.k.Modulus2()>right.k.Modulus2();
//	});
//
//	for(auto iter=KPs.begin(); iter!=KPs.end(); iter++)
//	{
//		pPot->constraints.push_back(*iter);
//		double chi2=(double)(pPot->constraints.size())/pConfig->NumParticle()/pConfig->GetDimension();
//		if(chi2<0.45)
//			continue;
//		if( (iter-KPs.begin())%10 !=0)
//			continue;
//		RelaxStructure_NLOPT(*pConfig, *pPot, 0.0, 0, 0.0, 10000);
//		pPot->SetConfiguration(*pConfig);
//		std::cout<<chi2<<" \t"<<pPot->Energy()<<'\n';
//	}
//	delete pPot;
//	return 0;
//}
//{
//	//add larger k first, add k one by one. not working either
//	double chi=0.65;
//	ShiftedCCPotential * pPot = GeneratePotential(pConfig, chi, 4, 0, std::cin);
//	std::vector<ShiftedCCPotential::KPointConstraint> KPs;
//	std::swap(KPs, pPot->constraints);
//	std::sort(KPs.begin(), KPs.end(), [] (const ShiftedCCPotential::KPointConstraint & left, const ShiftedCCPotential::KPointConstraint & right) ->bool
//	{
//		return left.k.Modulus2()>right.k.Modulus2();
//	});
//
//	for(auto iter=KPs.begin(); iter!=KPs.end(); iter++)
//	{
//		pPot->constraints.push_back(*iter);
//		double chi2=(double)(pPot->constraints.size())/pConfig->NumParticle()/pConfig->GetDimension();
//		if(chi2<0.450)
//			continue;
//		RelaxStructure_NLOPT(*pConfig, *pPot, 0.0, 0, 0.0, 10000);
//		pPot->SetConfiguration(*pConfig);
//		std::cout<<chi2<<" \t"<<pPot->Energy()<<'\n';
//	}
//	delete pPot;
//	return 0;
//}



int VirialSum(Configuration * pConfig, Potential * pPotential, RandomGenerator & /*gen*/)
{
	ShiftedCCPotential * pPot = dynamic_cast<ShiftedCCPotential *>(pPotential);
	struct term
	{
		double cutoff;
		double rf;
	};
	std::vector<term> Terms;
	const size_t NumConfig=10;
	double Volume;
	DimensionType dim;
	for(int i=0; i<10000; i+=1000)
	{
		std::stringstream ss;
		ss<<"C:/cygwin64/home/gezhang/MD_fixK_Again_2D/540/1/_0.02_0.05_540_10000_"<<i;
		Configuration ci=ReadPos(ss.str());
		//debug temp
		ci=*pConfig;
		if(i==0)
		{
			Volume=ci.PeriodicVolume();
			dim=ci.GetDimension();
		}

		for(size_t j=0; j<ci.NumParticle(); j++)
		{
			//std::cout<<"doing particle"<<j<<'\n';
			GeometryVector rj=ci.GetRelativeCoordinates(j);
			ci.IterateThroughNeighbors(rj, 20, [&] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * /*PeriodicShift*/, const size_t Sourceparticle) ->void
			{
				size_t k=Sourceparticle;
				if(!(k<j))
					return;
				GeometryVector rk=ci.GetRelativeCoordinates(k); (void)rk;
				//calculate f(rjk)
				GeometryVector vf=SameDimensionZeroVector(rj);
				for(auto iter=pPot->constraints.begin(); iter!=pPot->constraints.end(); iter++)
					vf.AddFrom( (-2.0/Volume*iter->V*std::sin(iter->k.Dot(shift))) * iter->k );

				GeometryVector rjk=shift;

				term t;
				t.cutoff=std::sqrt(rjk.Modulus2());
				t.rf=rjk.Dot(vf);
				Terms.push_back(t);

				//debut temp
				//std::cout<<t.cutoff<<" \t"<<std::sqrt(vf.Modulus2())<<" \t"<<t.rf<<'\n';
				//std::cout<<vf;
			});
		}
	}

	std::sort(Terms.begin(), Terms.end(), [](const term & left, const term & right) ->bool
	{
		return left.cutoff<right.cutoff;
	});

	double Sum=0.0;
	size_t count=0;
	void PlotFunction_Grace(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::string & Title);
	std::vector<GeometryVector> result;
	for(auto iter=Terms.begin(); iter!=Terms.end(); iter++, count++)
	{
		Sum+=iter->rf;
		if(count%10000==999)
		{
			std::cout<<iter->cutoff<<" \t"<<(-1)*Sum/NumConfig/dim/Volume<<'\n';
			logfile<<iter->cutoff<<" \t"<<(-1)*Sum/NumConfig/dim/Volume<<'\n';
			result.push_back(GeometryVector(iter->cutoff, (-1)*Sum/NumConfig/dim/Volume));
		}
	}
	PlotFunction_Grace(result, "2D_0.05_VirialSum_K2", "cutoff", "virial pressure", "d=2, T\\sE\\N=0.02, \\xc\\f{}=0.05, N\\sp\\N=540, N\\sc\\N=10, V(k)=|K-k|\\S2");

	return 0;
}


//const size_t GridSize=1000;
//const double Radius =1.46;
//int CollectiveCoordinateDebug(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen)
//{
//	Configuration c=ReadPos("temp2");
//
//	//std::vector<GeometryVector> vresult;
//	//auto GetConfigFunction = [&c](size_t i) -> const Configuration
//	//{
//	//	return c;
//	//};
//	//IsotropicTwoPairCorrelation(GetConfigFunction, 1, 5, vresult);
//
//	std::fstream ofile("temp3.txt", std::fstream::out);
//	for(int i=0; i<GridSize; i++)
//	{
//		for(int j=0; j<GridSize; j++)
//		{
//			GeometryVector center( (double)(i)/GridSize, (double)(j)/GridSize );
//			bool covered = false;
//			c.IterateThroughNeighbors(center, Radius, [&covered](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
//			{
//				covered=true;
//			}, &covered);
//
//			ofile<<covered<<' ';
//		}
//		ofile<<'\n';
//	}
//	return 0;
//}

int QCMDAnneal(Configuration * /*pConfig*/, Potential * /*pPotential*/, RandomGenerator & gen)
{
	GeometryVector basis[2]={GeometryVector(0, 1.0), GeometryVector(1.0, 0)};
	Configuration c(2 , &basis[0], 0.01, true);
	for(int i=0; i<144; i++)
		c.Insert("A", gen);
	ShiftedCCPotential * p1 = new ShiftedCCPotential(2);
	GeometryVector b0=c.GetReciprocalBasisVector(0);
	GeometryVector b1=c.GetReciprocalBasisVector(1);
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(0.0*b0+14.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(8.0*b0+11.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(13.0*b0+4.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(8.0*b0-11.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(13.0*b0-14.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(0.0*b0+36.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(21.0*b0+29.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(34.0*b0+11.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(34.0*b0-11.0*b1, -1.0) );
	p1->constraints.push_back( ShiftedCCPotential::KPointConstraint(21.0*b0-29.0*b1, -1.0) );



	RnPotential * p2 = new RnPotential(2, 1e-13, -12, 0.25);

	CombinedPotential pc(p1, p2);

	ParticleMDAnnealing an(c, 12345, 0.0, pc, false, &std::cout, 0.1);
	an.Anneal(pc);


	::AlsoWriteEPS=true;
	::Plot("temp", *an.pList);
	::Output("temp", *an.pList);

	return 0;
}


//return the maximum K so that config is still stealthy (S(k)<SkLimit for all k<K)
double StealthyMaximumK(const Configuration & config, double SkLimit=1e-8)
{
	if(config.NumParticle()==0)
	{
		std::cout<<"Warning in StealthyMaximumK : no particle in configuration!\n";
		return 1e10;
	}
	DimensionType d=config.GetDimension();
	GeometryVector b[ ::MaxDimension];
	for(DimensionType i=0; i<d; i++)
		b[i]=config.GetReciprocalBasisVector(i);

	double KScale=std::sqrt(b[0].Modulus2());
	PeriodicCellList<Empty> RecipLattice(d, b, KScale, true);
	RecipLattice.Insert( Empty(), SameDimensionZeroVector(b[0]) );
	for(double rK=3*KScale; ; rK*=1.5)
	{
		std::vector<GeometryVector> KPoints;
		auto IterateFunction = [& KPoints] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * /*PeriodicShift*/, const size_t /*Sourceparticle*/) ->void
		{
			KPoints.push_back(shift);
		};
		RecipLattice.IterateThroughNeighbors(SameDimensionZeroVector(b[0]), rK, IterateFunction);
		if(KPoints.size()<2)
			continue;
		auto CompareFunction = [] (const GeometryVector & left, const GeometryVector & right) ->bool
		{
			return left.Modulus2()<right.Modulus2();
		};
		std::sort(KPoints.begin(), KPoints.end(), CompareFunction);

		for(auto iter=KPoints.begin()+1; iter!=KPoints.end(); iter++)
		{
			if( StructureFactor(config, *iter) > SkLimit )
			{
				return (1.0-1e-13)*std::sqrt( iter->Modulus2() );
			}
		}
	}
}
double GetChi( const Configuration & c )
{
	DimensionType d=c.GetDimension();
	double maxk = StealthyMaximumK(c, 1e-8);
	double ncons = HyperSphere_Volume(d, maxk)/std::pow(2*pi, (double)(d))*c.PeriodicVolume()/2.0;
	return ncons/d/c.NumParticle();
}

//when assigned a configuration, automatically generate a ShiftedCCPotential with correct K points list
//V(k)=(1-k/rK)^2
class StealthyPotential : public Potential
{
public:
	double rK;
	ShiftedCCPotential * pCCPot;
	GeometryVector CachedBasisVector[ ::MaxDimension];
	double TargetDensity;
	StealthyPotential(DimensionType d, double rK, double TargetDensity) : Potential(d), rK(rK), TargetDensity(TargetDensity)
	{
		this->pCCPot=nullptr;
	}
	~StealthyPotential()
	{
		if(this->pCCPot!=nullptr)
			delete this->pCCPot;
	}
	virtual double Energy()//return the total potential energy of the given configuration
	{
		return this->pCCPot->Energy();
	}
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		this->pCCPot->Force(result, i);
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	//{
	//	std::cerr<<"Error : StealthyPotential::EnergyDerivativeToBasisVectors is not implemented!\n";
	//}
	virtual void SetConfiguration(const Configuration & config)
	{
		DimensionType d=this->Potential::Dimension;
		if(this->pCCPot!=nullptr)
		{
			//test if basis vector has never changed
			//if so, we don't need to re-generate the constrained k points
			bool SameBasisVector = true;
			for(DimensionType i=0; i<d; i++)
			{
				if(this->CachedBasisVector[i]!=config.GetBasisVector(i))
				{
					SameBasisVector=false;
					break;
				}
			}
			if(SameBasisVector)
			{
				this->pCCPot->SetConfiguration(config);
				return;
			}
		}
		GeometryVector b[ ::MaxDimension];
		for(DimensionType i=0; i<d; i++)
			b[i]=config.GetReciprocalBasisVector(i);
		PeriodicCellList<Empty> RecipLattice(d, b, 1e100, true);
		RecipLattice.Insert( Empty(), SameDimensionZeroVector(b[0]) );
		std::vector<GeometryVector> KPoints;
		auto IterateFunction = [& KPoints, &d] (const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
		{
			bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
			for(DimensionType i=0; i<d; i++)
			{
				if(PeriodicShift[i]>0)
				{
					NonTrivial=true;
					break;
				}
				else if(PeriodicShift[i]<0)
				{
					NonTrivial=false;
					break;
				}
			}
			if(NonTrivial)
				KPoints.push_back(shift);
		};
		RecipLattice.IterateThroughNeighbors(SameDimensionZeroVector(b[0]), rK, IterateFunction);
		auto CompareFunction = [] (const GeometryVector & left, const GeometryVector & right) ->bool
		{
			return left.Modulus2()<right.Modulus2();
		};
		std::sort(KPoints.begin(), KPoints.end(), CompareFunction);

		if(this->pCCPot!=nullptr)
			delete this->pCCPot;
		ShiftedCCPotential Pot(d);
		double SumVk=0.0;
		for(auto iter=KPoints.begin(); iter!=KPoints.end(); iter++)
		{
			double temp=1.0-std::sqrt(iter->Modulus2())/rK;
			double Vk=temp*temp;
			Vk=1.0;
			Pot.constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, Vk));
			SumVk+=Vk;
		}

		size_t nbr=config.NumParticle();
		double v=config.PeriodicVolume();


		//double temp=config.PeriodicVolume()/config.NbrParticle()-1.0/TargetDensity;
		//this->pCCPot->Shift=temp*temp*0.01;


		Pot.Shift=nbr*(nbr-1)/v/2.0-nbr/v*SumVk;
		Pot.SetConfiguration(config);

		this->pCCPot=new ShiftedCCPotential_ForMC(Pot);

		for(DimensionType i=0; i<d; i++)
			this->CachedBasisVector[i]=config.GetBasisVector(i);
	}
};



#include <RandomSequentialAddition.h>
double MinDistance(const Configuration & a, size_t * pi, size_t * pj);
bool IsSaturated(const SpherePacking & pa, double radius);
//#include "KissFFT/kiss_fft.h"
//#include "KissFFT/kiss_fftnd.h"
void PlotPotential(double /*Volume*/)
{
}


class LinearBin
{
private:
	size_t NBins;
	double MinE, ERange;
public:
	LinearBin(double MaxE, double MinE, size_t NumBins) : NBins(NumBins), MinE(MinE), ERange(MaxE-MinE)
	{
	}
	size_t NumBins()
	{
		return NBins;
	}
	long GetBin(double E)
	{
		return std::floor( (E-MinE)/ERange*NBins );
	}
	double GetBinLowerBound(size_t NumBin)
	{
		return (double)(NumBin)/NBins*ERange+MinE;
	}
	double GetBinUpperBound(size_t NumBin)
	{
		return (double)(NumBin+1)/NBins*ERange+MinE;
	}
};


//move particles in conf so that its energy is within MinE and MaxE
//return whether successful or not
bool AdjustEnergy(Potential & pot, Configuration & conf, double MinE, double MaxE)
{
	size_t Count=0;
	RandomGenerator gen(12345);
	double d=1e-5;
	double mind=0.0, maxd= ::MaxDistance;
	for(;;)
	{
		Count++;
		if(Count>1000)
		{
			return false;
		}
		pot.SetConfiguration(conf);
		double E=pot.Energy();
		if(E>MinE && E<MaxE)
			break;
		else if(E>MaxE)
			RelaxStructure_NLOPT_Emin(conf, pot, 0.0, 0, 0.0, 10000, MaxE);
		else if(E<MinE)
		{
			//generate configuration with particle displacements
			Configuration conf2(conf);
			for(size_t i=0; i<conf2.NumParticle(); i++)
			{
				GeometryVector rel=conf2.GetRelativeCoordinates(i);
				for(DimensionType j=0; j<conf2.GetDimension(); j++)
					rel.x[j]+=(gen.RandomDouble()-0.5)*d;
				conf2.MoveParticle(i, rel);
			}
			//check energy
			pot.SetConfiguration(conf2);
			double newE=pot.Energy();
			if(newE>MinE && newE < MaxE)
			{
				conf=conf2;
			}
			else if(newE<MinE)
			{
				if(mind<d) mind=d;
				if(maxd== ::MaxDistance)
					d*=2;//no information of maximum displacement
				else
					d=(d+maxd)*0.5;
			}
			else if(newE>MaxE)
			{
				if(maxd>d) maxd=d;
				d=(d+mind)*0.5;
			}
		}
	}
	return true;
}
//set OutputPrefix to "No Output" to supress outputing files
std::vector<GeometryVector> CC_WLMC2(Potential & pot, const Configuration & StartConfig, std::string OutputPrefix, size_t cycle, size_t step, double MinEnergy, double MaxEnergy)
{
	DimensionType dim=StartConfig.GetDimension(); (void)dim;
	//Configuration * pc = GenerateConfig(NumP, std::pow(V, 1.0/dim), dim, gen, false);
	Configuration * pc = new Configuration(StartConfig);
	if(!AdjustEnergy(pot, *pc, MinEnergy, MaxEnergy))
		return std::vector<GeometryVector>();
	pot.SetConfiguration(*pc);

	//RelaxStructure(*pc, pot, pressure, 0.0);
	//RelaxStructure_NLOPT(*pc, pot, 0.0, 0, 0.0);
	pot.SetConfiguration(*pc);
	double StartEnergy = pot.Energy();
	std::cout<<'\n'<<"simulation start energy: "<<StartEnergy<<"\n";
	::Plot(OutputPrefix+std::string("InitConfig"), *pc);
	::Output(OutputPrefix+std::string("InitConfig"), *pc);

	AtomMove ma;
	ma.MinDistance=0.0;
	MCMove<Configuration, Potential *> & m=ma;


	LinearBin bin(MaxEnergy, MinEnergy, 500);
	auto GetEnthalpyFunc = [&] ( const Configuration & c ) -> double
	{
		pot.SetConfiguration(c);
		return pot.Energy();
	};
	WangLandauMonteCarlo<Configuration, Potential *, LinearBin> wlmc(*pc, 12345, GetEnthalpyFunc(*pc), bin);
	progress_display pd(cycle);
	for(size_t i=0; i<cycle; i++)
	{
		pd++;
		if(ma.AcceptRatio>ma.MyMinAcceptance-0.1)
			wlmc.SIncrease=5.0/(i+10);
		else
			wlmc.SIncrease=0.001;//increase s very cautiously when the system might not be equilibrated
		//At the beginning, the move classes might not have ideal sigma, leading to sharp peaks in histograms.

		wlmc.Move(step, m, &pot, GetEnthalpyFunc);
		if(OutputPrefix!="No Output")
		{
			std::stringstream s1, s2, s3;
			s1<<OutputPrefix<<"Conf"<<i;
			s2<<"H="<<GetEnthalpyFunc(*wlmc.pSys);
			::Plot(s1.str(), *wlmc.pSys, s2.str());
			::Output(s1.str(), *wlmc.pSys);
			s1<<"_Histogram";
			std::vector<GeometryVector> h=wlmc.GetHistogram();
			wlmc.ClearHistogram();
			PlotFunction_Grace(h, s1.str(), "H", "p(H)", "");
			s3<<OutputPrefix<<"Conf"<<i<<"_Entropy";
			h=wlmc.GetEntropy();
			PlotFunction_Grace(h, s3.str(), "H", "S(H)", "");
		}
	}
	delete pc;
	return wlmc.GetEntropy();
}
int CC_WLMC()
{
	std::string start, outputprefix;
	double chi, emin, emax;
	size_t l, cycle, step;
	std::cout<<"Starting structure=";
	std::cin>>start;
	Configuration c=ReadPos(start);
	std::cout<<"multiplicate in each direction by:";
	std::cin>>l;
	c=MultiplicateStructure(c, l);
	std::cout<<"Chi=";
	std::cin>>chi;
	ShiftedCCPotential * ppot = GeneratePotential(&c, chi, 1, 0, std::cin);
	std::cout<<"Final chi="<<chi<<'\n';
	std::cout<<"Prefix=";
	std::cin>>outputprefix;
	std::cout<<"E_max=";
	std::cin>>emax;
	std::cout<<"E_min=";
	std::cin>>emin;
	std::cout<<"cycle=";
	std::cin>>cycle;
	std::cout<<"step=";
	std::cin>>step;

	ShiftedCCPotential_ForMC pot(*ppot);
	std::vector<GeometryVector> temp2=CC_WLMC2(pot, c, outputprefix, cycle, step, emin, emax);
	PlotFunction_Grace(temp2, outputprefix+std::string("_FinalEntropy"), "E", "S(E)", "");

	//do linear regression
	double sxy=0.0, sx=0.0, sy=0.0, sx2=0.0;
	for(auto iter=temp2.begin(); iter!=temp2.end(); iter++)
	{
		double x=std::log(iter->x[0]);
		double y=iter->x[1];
		sxy+=x*y;
		sx+=x;
		sy+=y;
		sx2+=x*x;
	}
	double n=temp2.size();
	double s1=(sxy/n-sx*sy/n/n)/(sx2/n-sx*sx/n/n);
	std::cout<<"s1="<<s1<<'\n';
	delete ppot;
	return 0;
}
