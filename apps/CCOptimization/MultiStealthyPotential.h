/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com
*	Date	: November 2025 */

/** \file MultiStealthyPotential.h
 * \brief Header file to declare a class to compute multi-stealthy collective coordinate potential + soft-core repulsions */

#ifndef MultiStealthyPotential_H
#define MultiStealthyPotential_H

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
/* in $(pot) */
#include <Potential.h>
#include <CombinedEnergyForce.h>
#include <CollectiveCoordinatePotential.h>

#include <algorithm>
#include <memory>

struct Constraint_Fourier
{
	GeometryVector k;//absolute k, not relative to basis vectors
	double V;
	double rho_real, rho_imag;
	/** \brief Construct a Fourier-space constraint with a wavevector and energy weight.
	 *  \param k Absolute wavevector of the constraint.
	 *  \param V Fourier-space energy weight at this wavevector. */
	Constraint_Fourier(const GeometryVector & k, double V): k(k), V(V)
	{
	}
};

/* Constraint energies in Fourier space for a given subset of particles */
struct EnergyConstraints
{
	std::vector<size_t> idx ;
	bool RhoValid = false;
	std::vector<Constraint_Fourier> constraints;

	/** \brief Add a Fourier-space constraint at a wavevector for a species
	 *  \param k  A wavevector
	 *  \param v  Fourier-space energy, must be non-negative. */
	void AddConstraint(const GeometryVector & k, double v){
		this->constraints.emplace_back(k, v);
	}

	// /** \brief Print the defined constraints with the indices of species
	//  *  \param ofile  Destination output stream. */
	// void PrintConstraints(std::ostream & ofile = std::cout){
	// 	// particle indices
	// 	ofile << "indices:";
	// 	for (auto i = idx.begin(); i != idx.end(); i++)
	// 		ofile << "\t" << &i ;

	// 	ofile <<"\n";
	// 	// constraints
	// 	for (size_t i =0; i<this->constraints.size(); i++){
	// 		ofile << "@ k= " << this->constraints[i].k << ": v="<<this->constraints[i].V <<"\n";
	// 	}
	// 	ofile << "\n\n";
	// }
};

/* collective-coordinate potential in a species. */
struct species_potential{
	std::unique_ptr<ShiftedCCPotential> CCPotential;
	std::vector<size_t> idx ;   // particle i in CCPotential -> pConfig  =  particle idx[i] in the mother configuration.

	/** \brief Construct an empty species potential with a fresh collective-coordinate potential.
	 *  \param dim  Spatial dimension of the configuration. */
	species_potential (DimensionType dim): CCPotential(std::make_unique<ShiftedCCPotential>(dim)), idx() {}
	// non-copyable to preserve unique ownership of CCPotential
	/** \brief Deleted copy constructor (preserves unique ownership of CCPotential). */
	species_potential(const species_potential&) = delete;

	/** \brief Deleted copy assignment (preserves unique ownership of CCPotential). */
	species_potential& operator=(const species_potential&) = delete;

	/** \brief Defaulted move constructor. */
	species_potential(species_potential&&) = default;

	/** \brief Defaulted move assignment. */
	species_potential& operator=(species_potential&&) = default;

	/** \brief Print the defined constraints with the indices of species
	 *  \param ofile  Destination output stream. */
	void PrintConstraints(std::ostream & ofile = std::cout){
		// particle indices
		ofile << "indices:";
		for (auto i = idx.begin(); i != idx.end(); i++)
			ofile << "\t" << *i ;

		ofile <<"\n";
		// constraints
		for (size_t i =0; i<CCPotential->constraints.size(); i++){
			ofile << "@ k= " << CCPotential->constraints[i].k << ": v="<<CCPotential->constraints[i].V <<"\n";
		}
		ofile << "\n\n";
	}

};


struct spec_index{
	int species;
	size_t index;

	/** \brief Construct a species/index pair locating a particle within its species sub-configuration.
	 *  \param s    Species id (-1 if unspecified).
	 *  \param idx  Index of the particle within that species. */
	spec_index(int s, size_t idx): species(s), index(idx){}
};

class MultiStealthyPotential : public Potential, public CombinedEnergyForce
{
protected:
	std::unique_ptr<Configuration> pConfig;	//< packing of interest
	std::vector<double> radii; 		//< particle radii
	double rc;						//< the critical radius over which the contact potential become zero.
	bool WithRepulsion = true;
	bool SubconfigValid = false;
	//	bool RhoValid;
	//	double Volume;
	std::vector<spec_index> idx2species;

	/** \brief Work of one evaluation's concurrent FINUFFT path, prepared serially by PrepareSpecies. */
	struct SpeciesWork
	{
		/** \brief One FINUFFT particle block of one species. */
		struct Task
		{
			long species;	//< species index
			size_t block;	//< block index, or NoBlock for a species whose rho is already valid
			double cost;	//< predicted time (ShiftedCCPotential::FinufftBlockTime), only for the order
		};
		static constexpr size_t NoBlock = static_cast<size_t>(-1);	//< Task::block of a task without a Type 1 block
		std::vector<char> concurrent;	//< one flag per species, nonzero for the concurrent path
		std::vector<Task> rho;	//< phase-1 tasks: Type 1 blocks (or NoBlock), largest predicted time first
		std::vector<Task> force;	//< phase-2 tasks: Type 2 blocks of species with several, largest first
		std::vector<size_t> forceBlocks;	//< per species: Type 2 blocks (PrepareFinufftForces)
		std::vector<long> rhoLeft, forceLeft;	//< per species: phase-1 and phase-2 tasks not yet finished (shared)
		std::vector<char> ready;	//< per species: StartForces succeeded (written in phase 1, read in phase 2)
		long nextChunk = 0;	//< next contact chunk to take (shared by the team)
		long pendingRho = 0;	//< phase-1 tasks not yet finished (shared by the team)
	};

	// ponytail: single-thread cost of the contact walk per particle [microseconds], by dimension, measured on 2026-10-14 on
	// the development machine (fused walk, bench_finufft setup: packing fraction 0.2, cells of the cutoff). It only tells
	// the force gate how much other work fills idle threads (FinufftForceGate's otherWork); denser packings walk more
	// neighbors. Ceiling: another machine or packing fraction; upgrade: time ContactEnergyAndForce there.
	static constexpr double ContactCostParticle[3] = {0.05, 0.2, 1.1};

	/** \brief Prepare every species' FINUFFT plans (ShiftedCCPotential::PrepareFinufft and, with forces,
	 *  PrepareFinufftForces; serially, so no plan is made inside a parallel region), mark the species evaluated on
	 *  the concurrent path and list their block tasks. A species qualifies when its GetRho will run FINUFFT
	 *  (FinufftEligible): its rho is then split into Type 1 particle blocks (StartRho, block count from the inputs
	 *  only), and with forces its Type 2 forces into Type 2 blocks when the joint cost gate
	 *  (ShiftedCCPotential::FinufftForceGate over the concurrent species, with the contact walk as other work)
	 *  accepted them; each Type 1 block runs as one
	 *  task on one thread, and so does each Type 2 block of a species with several (RunSpeciesWork). Every other species (Direct, or FINUFFT falling back to the direct sum) keeps
	 *  its own evaluation, parallel over wavevectors with the full team. The marks and tasks depend only on
	 *  constraints, cells, particle numbers, settings and the team size, never on timing.
	 *  \param forces  Also decide and prepare the Type 2 force transforms (AllForce); false for Energy.
	 *  \return The flags and tasks; the rho of a concurrent species is started (StartRho(forces)). */
	SpeciesWork PrepareSpecies(bool forces);

	/** \brief The concurrent FINUFFT species of one evaluation; every thread of a parallel region of the team calls it.
	 *  Phase 1: the Type 1 block tasks, dynamically, largest first; the thread that finishes a species' last block
	 *  continues with that species on one thread (OnOneThread): FinishRho (rho is the sum of the blocks in block order,
	 *  or the direct fallback), the species energy if wanted, StartForces if forces, and a single Type 2 block with
	 *  FinishForces right away. Threads without a task take contact chunks, while phase-1 tasks remain if phase 2
	 *  follows, else until none is left. Barrier. Phase 2, only if some species has several Type 2 blocks: those
	 *  blocks, the thread finishing a species' last one running FinishForces, then the remaining contact chunks.
	 *  Barrier. Every result has its own slot and species state is handed over by seq_cst countdowns, so nothing
	 *  depends on the schedule; an evaluation without split Type 2 work has the one barrier of the unsplit schedule.
	 *  \param work                PrepareSpecies result (its shared counters are used).
	 *  \param speciesEnergy       Energy slot of each species to fill, or nullptr.
	 *  \param speciesForce        Type 2 force buffer of each species, or nullptr for no forces.
	 *  \param[out] bulk           Per species: nonzero if its Type 2 forces succeeded (set only with forces).
	 *  \param items               Number of contact items.
	 *  \param chunk               Contact items per chunk.
	 *  \param contact             Callable taking an item index (long), run once per item. */
	template <class Contact>
	void RunSpeciesWork(SpeciesWork & work, std::vector<double> * speciesEnergy,
		std::vector<std::vector<GeometryVector>> * speciesForce, std::vector<char> & bulk, long items, long chunk,
		const Contact & contact);

	/** \brief Add particle i's collective-coordinate force, scaled by u0 N_s/N, to its contact force, with the
	 *  arithmetic of Force(i): from its species' Type 2 forces when bulk, else from the species' direct sum.
	 *  \param i             Particle index.
	 *  \param contactForce  e0 times the contact force on particle i (ContactForce).
	 *  \param bulk          Per species: nonzero if speciesForce holds its forces.
	 *  \param speciesForce  Type 2 forces per species.
	 *  \return              The total force on particle i. */
	GeometryVector WithSpeciesForce(size_t i, const GeometryVector & contactForce, const std::vector<char> & bulk,
		const std::vector<std::vector<GeometryVector>> & speciesForce);

	/** \brief Soft-core contact force on one particle.
	 *  \param i  Particle index.
	 *  \return   e0 times the contact force on particle i (zero when repulsion is off). */
	GeometryVector ContactForce(size_t i) const;

	/** \brief Contact energy slot and contact force of one particle from one neighbor traversal: the same
	 *  neighbors, order and expressions as the contact sum of Energy (pairs with i < SourceParticle) and as
	 *  ContactForce (every SourceParticle != i), so both results are bitwise those of the separate walks.
	 *  \param i           Particle index.
	 *  \param[out] force  e0 times the contact force on particle i, as ContactForce(i).
	 *  \return            Contact energy of the pairs (i, j > i), before the factor e0 (zero when repulsion is off). */
	double ContactEnergyAndForce(size_t i, GeometryVector & force) const;

	/** \brief Total energy from per-species and per-particle slots, summed in index order.
	 *  \param speciesEnergy  Energy of each species potential.
	 *  \param contactEnergy  Contact energy slot of each particle (empty when repulsion is off).
	 *  \return               u0 times the sum of N_s/N times the species energies, plus e0 times the contact sum. */
	double TotalEnergy(const std::vector<double> & speciesEnergy, const std::vector<double> & contactEnergy) const;

public :
	double e0 = 1.0;	//< energy scale for contact potential
	double u0 = 1.0;	//< energy scale for collective-coordinate potential
	
	// constraints[i]: Fourier-space energy constraints for a species i.
	std::vector<species_potential> species_constraints;

	/** \brief Constructor
	 *  \param dimension  Spatial dimension of the configuration.
	 *  \param epsilon    Energy scale of the soft-core contact potential; repulsion is disabled when near zero. */
	MultiStealthyPotential(DimensionType dimension, double epsilon): Potential(dimension){
		this->pConfig=nullptr;
		this->rc = -1.0;
		this->ParallelNumThread=1;
		this->e0 = epsilon;
		this->WithRepulsion = (std::abs(this->e0) > 1e-10);
		if(this->WithRepulsion){
			std::cout<< "Soft-core repulsion is on\n";
		}
		else{
			std::cout<< "Soft-core repulsion is off\n";
		}
	}

	/** \brief Defaulted destructor. */
	~MultiStealthyPotential() override = default;

	/** \brief Deleted copy constructor (preserves unique ownership of pConfig). */
	MultiStealthyPotential(const MultiStealthyPotential&) = delete;

	/** \brief Deleted copy assignment (preserves unique ownership of pConfig). */
	MultiStealthyPotential& operator=(const MultiStealthyPotential&) = delete;

	/** \brief Defaulted move constructor. */
	MultiStealthyPotential(MultiStealthyPotential&&) = default;

	/** \brief Defaulted move assignment. */
	MultiStealthyPotential& operator=(MultiStealthyPotential&&) = default;

	/* -------------------------------------
	 * 	Methods to define configuration and potential
	 * -------------------------------------*/

	/** \brief Set the number of species in the system. */
	/** \brief Set the number of species in this configuration.
	 *  \param num  Number of species to allocate. */
	void SetNumSpecies(const size_t num){
		for (size_t i = 0; i < num ; i++){
			species_constraints.emplace_back(this->Dimension);
		}
		//std::cout << "species = "<< species_constraints.size() <<"\n\n";
	}

	/** \brief Determine the species of a particle
	 *  \param prt_idx  Global particle index in the full configuration.
	 *  \param species  Target species id to assign. */
	void SetSpeciesOfParticle(size_t prt_idx, size_t species);

	/** \brief Add a constraint to a species.
	 *  \param species  Species id to target.
	 *  \param k        Wavevector of the constraint.
	 *  \param V        Constraint coefficients (energy weights). */
	void AddConstraint(size_t species, const GeometryVector & k, const std::vector<double> V);

	/** \brief Select the collective-density backend of every species defined so far.
	 *  \param backend  Backend passed to ShiftedCCPotential::SetRhoBackend.
	 *  \throw std::invalid_argument if the backend is not compiled in (see USE_FINUFFT in the makefile). */
	void SetRhoBackend(ShiftedCCPotential::RhoBackend backend){
		for (auto & species : this->species_constraints)
			species.CCPotential->SetRhoBackend(backend);
	}

	/** \brief Update rho values for selected species.
	 *  \param species  Species id to evaluate. */
	void GetRho(size_t species);

	/** \brief Print constraints to an output stream.
	 *  \param ofile  Destination output stream. */
	void PrintConstraints(std::ostream & ofile = std::cout);

	/** Output encodings: readable .txt or versioned Multi-Stealthy Packing Binary .mspb. */
	enum class ExportFormat { Text, Binary };

	/** \brief Export the currently stored configuration as text or MSPB binary.
	 *  Columns are id, species, fractional coordinates (sx, sy, sz, then s3, ...),
	 *  and radius in Cartesian units. The header begins with "# "; the 2D data
	 *  columns match the sweep final.tsv files.
	 *  Particle order and actual assigned species/radii are preserved. Comment-prefixed
	 *  metadata includes dimension, particle count, basis_vector-0, basis_vector-1, ...,
	 *  coordinates relative, and radius_units cartesian to recover the periodic cell.
	 *  Call SetConfiguration(finalPosition) after optimization to refresh the snapshot.
	 *  Throws on missing configuration, inconsistent metadata, unassigned species,
	 *  or output failure. Does not change the destination stream's formatting.
	 *  \param ofile Destination stream; open binary files with std::ios::binary.
	 *  Caller owns opening/closing files and must check for errors on close.
	 *  \param format Text by default; Binary writes the layout in EXPORT_FORMAT.md.
	 *  Encoding is explicit, not inferred from a filename. */
	void ExportConfiguration(std::ostream & ofile, ExportFormat format = ExportFormat::Text) const;

	/** Return an in-memory snapshot with the stored cell, relative coordinates, and
	 *  individual radii. Refresh with SetConfiguration(finalPosition) first. */
	SpherePacking GetPackingSnapshot() const;

	/** \brief Set the working configuration.
	 *  \param Position  Full particle configuration. */
	virtual void SetConfiguration(const Configuration & Position);

	/** \brief Set particle radii explicitly.
	 *  \param radii  Radii list matching the particles. */
	void SetSecureRadii(const std::vector<double> & radii);

	/** \brief Size the contact cell list to the contact cutoff and prepare the neighbor search. */
	void PrepareContactCells();

	/** \brief Set the configuration from a sphere packing.
	 *  \param pck  Input sphere packing. */
	void SetPacking(const SpherePacking & pck);

	/** \brief Rebind every species potential, if stale, to an indexed view of *pConfig in species idx order. */
	void UpdateSubconfigs();

	/** \brief Get the contact cutoff radius.
	 *  \return  The contact cutoff radius rc. */
	double GetRc(void) const {return this->rc;};

	/* -------------------------------------
	 * 	Implementations of this class.
	 * -------------------------------------*/

	/** \brief Compute total potential energy.
	 *  Species that run FINUFFT (see PrepareSpecies) are evaluated inside one parallel region as Type 1 block tasks
	 *  while free threads take the contact sum (RunSpeciesWork); other species are evaluated first, one by one, with
	 *  the full team. Species and particles own result slots that are summed in index order, so the result
	 *  does not depend on thread scheduling.
	 *  \return  The total potential energy. */
	virtual double Energy();

	/** \brief Compute force on a particle.
	 *  \param result  Output force vector.
	 *  \param i       Particle index. */
	virtual void Force(GeometryVector & result, size_t i);

	/** \brief Compute forces on all particles.
	 *  Fills every species' rho first, so it is valid without a prior Energy call (MD order); species that
	 *  run FINUFFT are filled inside the same parallel region as the force loop, as block tasks, together with
	 *  their Type 2 forces where the joint cost gate accepted them (ShiftedCCPotential::FinufftForceGate), while
	 *  free threads take the contact forces (RunSpeciesWork); every other particle gets Force(i) arithmetic. See
	 *  GetLastForceBackendUsed of each species potential.
	 *  \param results  Output vector of per-particle forces. */
	virtual void AllForce(std::vector<GeometryVector> & results);

	/** \brief Compute the total energy and every force in one evaluation, bitwise equal to Energy() followed
	 *  by AllForce(results) (the optimizer's gradient callback, see CombinedEnergyForce). The species calls are
	 *  those of the two (FINUFFT species as block tasks, the others with the full team); the contact part walks each
	 *  particle's neighbors once (ContactEnergyAndForce) instead of twice. With FINUFFT species, threads free
	 *  from block tasks take the contact walk (RunSpeciesWork), and after the barrier each particle adds its species
	 *  force as AllForce does. Keeps no state beyond the call.
	 *  \param results  Output vector of per-particle forces.
	 *  \return         The total potential energy. */
	double EnergyAndAllForce(std::vector<GeometryVector> & results) override;
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure);

};


#endif 
