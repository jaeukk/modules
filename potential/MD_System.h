#ifndef SYSTEM_INCLUDED
#define SYSTEM_INCLUDED

#include <functional>
#include <vector>
#include <cmath>
#include <float.h>

/* in $(core) */
#include <RandomGenerator.h>
#include <PeriodicCellList.h>
#include <etc.h>

//#include "Plots.h"

#include "Potential.h"





/** \brief Molecular-dynamics integrator for a system of point particles. */
class ParticleMolecularDynamics
{
private:
	//before using, potential must be set to use the current configuration.
	//get Acceleration of all particles
	/** \brief Compute the acceleration of all particles from the current configuration.
		\param potential  Potential providing the forces; must be set to the current configuration. */
	void _GetAcc(Potential & potential);

	/** \brief Compute the initial acceleration from the current configuration on the first integration step.
		The value constructors receive no Potential and so cannot fill Acceleration, but velocity-Verlet
		needs a(t0) for its first half-kick; this lazily computes it once (guarded by AccelerationInitialized).
		\param potential  Potential providing the forces. */
	void _EnsureInitialAcceleration(Potential & potential);
public:
	double TimeStep;
	double Time;//A Time Counter
	std::vector<double> Mass;
	Configuration Position;
	std::vector<GeometryVector> Velocity , Acceleration;
	size_t NumParticle;
	unsigned short Dimension;
	double MaxDeltaForce;//maximum delta F, user can query this variable to determine step size
	bool AccelerationInitialized;//true once Acceleration holds the force of the current Position (see _EnsureInitialAcceleration)

	//variable for Nose-Hoover algorithm
	double xi;
	double MassQ;

	////////////////////////////////////////////////////////
	/** \brief Initialize the Nose-Hoover thermostat.
		\param Q  Thermostat mass parameter (effective inertia of the heat bath). */
	void InitNoseHoover(double Q);

	/** \brief Construct from a configuration with a uniform particle mass.
		\param Position  Initial particle configuration.
		\param TimeStep  Integration time step.
		\param Mass      Mass assigned to every particle. */
	ParticleMolecularDynamics(const Configuration & Position, double TimeStep, double Mass);

	/** \brief Destructor. */
	~ParticleMolecularDynamics();

	/** \brief Construct by reading the MD state from a binary input stream.
		\param ifile  Input stream produced by WriteBinary. */
	ParticleMolecularDynamics(std::istream & ifile);

	/** \brief Write the MD state to a binary output stream.
		\param ofile  Output stream to write the state to. */
	void WriteBinary(std::ostream & ofile);

	/** \brief Return the potential energy of the current configuration.
		\param potential  Potential used to evaluate the energy.
		\return The potential energy. */
	double GetPotentialEnergy(Potential & potential);

	/** \brief microcanonical ensemble.
		\param repeat     Number of integration steps to perform.
		\param potential  Potential providing the forces. */
	void Evolve(size_t repeat, Potential & potential);//microcanonical ensemble
	/** \brief Evolve and then adjust the time step according to the observed per-step change in log total energy.
		\param repeat Number of integration steps to perform before adjustment.
		\param potential Potential providing the forces and energy.
		\param DeltaLogEnergyPerStep Target magnitude of the per-step change in log total energy. */
	void Evolve_AutoTimeStep(size_t repeat, Potential & potential, double DeltaLogEnergyPerStep);

	/** \brief Evolve while recording the maximum per-step change in force (into MaxDeltaForce).
		\param repeat     Number of integration steps to perform.
		\param potential  Potential providing the forces. */
	void Evolve_RecordDeltaForce(size_t repeat, Potential & potential);

	/** \brief Evolve in the canonical ensemble using the Nose-Hoover thermostat.
		\param repeat     Number of integration steps to perform.
		\param potential  Potential providing the forces.
		\param kT         Target temperature (in energy units). */
	void NoseHooverEvolve(size_t repeat, Potential & potential, double kT);

	/** \brief Evolve in the canonical ensemble using the Andersen thermostat (random velocity resets).
		\param repeat                Number of integration steps to perform.
		\param potential             Potential providing the forces.
		\param kT                    Target temperature (in energy units).
		\param ResetSpeedProbability  Per-step expected number of particle velocity resets.
		\param gen                   Random number generator. */
	void AndersonEvolve(size_t repeat, Potential & potential, double kT, double ResetSpeedProbability, RandomGenerator & gen);

	/** \brief Return the total kinetic energy of the system.
		\return The kinetic energy. */
	double GetKineticEnergy(void);

	/** \brief Draw particle velocities from the Maxwell-Boltzmann distribution at the given temperature.
		\param kT   Target temperature (in energy units).
		\param gen  Random number generator. */
	void SetRandomSpeed(double kT, RandomGenerator & gen);

	/** \brief Subtract the mean particle velocity so that the sum of velocities vanishes; for equal masses this also removes the total momentum. */
	void RemoveOverallMomentum(void);

	//search for the optimal time step size
	//variable TimeStep as well as this->TimeStep was adjusted so that (maximum delta force)/(average force)
	//between two steps is close to TargetDeltaForce
	//
	//This function DESTROYS equilibrium. Re-equilibrate after calling this function
	/** \brief Adjust the time step until the maximum componentwise acceleration change, normalized by the mean acceleration magnitude, is close to TargetDeltaForce; this destroys equilibrium.
		\param potential Potential providing the forces.
		\param[in,out] TimeStep Trial time step, updated in place and copied to this->TimeStep.
		\param TargetDeltaForce Target normalized acceleration-change ratio.
		\param kT Temperature used to draw trial velocities. */
	void FindOptimalTimeStep(Potential & potential, double & TimeStep, double TargetDeltaForce, double kT);
};

// Basis vector moves are not implemented. i.e. simulation box is not changed
/** \brief Molecular-dynamics simulated annealing driver for finding low-energy particle configurations. */
class ParticleMDAnnealing
{
public:
	ParticleMolecularDynamics * pSys;
	double alpha;
	double Tcut;
	double Tinit;
	double Pressure;
	RandomGenerator gen;
	Configuration * pList;
	std::ostream * Output;
	static bool StopSignal;
	bool SameLocalMinimaDetection;
	bool BoxMCMove;
	bool AutoTimeStep;
	bool Continue;//if true, read the last configuration from MDAnneal.ConfigPack and continue from there
					//also if true, will append to MDAnneal.ConfigPack rather than overwrite it.
	std::string ConfigPackName;//default: "MDAnneal"

	std::function<void(Configuration & List, Potential & pot, double Pressure, double MinDistance)> RelaxStructureFunction;

	//ParticleMDAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, double MinDistance=0, bool UseSortedList=false);//AlgorithmLevel : 0:Structure Optimization, 1:Monte Carlo, >1:Parallel Tempering with n copies
	/** \brief Construct molecular-dynamics annealing from an initial configuration.
		\param pInitConfig Initial configuration, copied into the internal MD system.
		\param RandomSeed Random-number seed.
		\param Pressure Pressure used in the enthalpy.
		\param pot Interaction potential.
		\param BoxMCMove If true, also attempt Monte Carlo cell moves; otherwise the cell remains fixed.
		\param output Optional logging stream.
		\param CoolingScheduleRescale Greater than one cools more slowly; less than one cools faster.
		\param UseSortedList Compatibility parameter; ignored because pInitConfig is copied with its existing representation. */
	ParticleMDAnnealing(const Configuration & pInitConfig, int RandomSeed, double Pressure, Potential & pot, bool BoxMCMove = false, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, bool UseSortedList=false);

	/** \brief Run the MD simulated-annealing procedure to find a low-energy configuration.
		\param pot       Potential defining the interaction.
		\param Callback  Function invoked once per temperature step (e.g. for progress reporting). */
	void Anneal(Potential & pot, std::function<void (void)> Callback=[](){});

	/** \brief Destructor; deletes the internal MD system. */
	~ParticleMDAnnealing()
	{
		//if(this->pList!=nullptr)
		//	delete this->pList;
		if(this->pSys!=nullptr)
			delete this->pSys;
	}
};

#endif
