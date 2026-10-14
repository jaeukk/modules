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

class MultiStealthyPotential : public Potential
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
	std::vector<Configuration> subc;

	/** \brief Whether Energy and AllForce run the species concurrently, one thread per species.
	 *  True when any species selects the FINUFFT backend: its Type 1 call is single-threaded, so only
	 *  running species side by side uses the other threads. Direct keeps the serial species loop,
	 *  parallel over wavevectors, which balances better when S is not a multiple of the thread count.
	 *  \return true for the species-concurrent path. */
	bool SpeciesConcurrent() const;

	/** \brief Fill species s's rho and return its CC energy, running that species on the calling thread.
	 *  The species' ParallelNumThread is 1 during the call and restored afterwards; touches only species s.
	 *  \param s  Species id.
	 *  \return   The species' collective-coordinate energy. */
	double SpeciesEnergyOnOneThread(size_t s);


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

	/** \brief Set the configuration from a sphere packing.
	 *  \param pck  Input sphere packing. */
	void SetPacking(const SpherePacking & pck);

	/** \brief Rebuild species sub-configurations. */
	void UpdateSubconfigs();

	/** \brief Get the contact cutoff radius.
	 *  \return  The contact cutoff radius rc. */
	double GetRc(void) const {return this->rc;};

	/* -------------------------------------
	 * 	Implementations of this class.
	 * -------------------------------------*/

	/** \brief Compute total potential energy.
	 *  With FINUFFT (see SpeciesConcurrent) one parallel region runs every species on its own thread
	 *  while the remaining threads start on the contact sum. Species and particles own result slots that
	 *  are summed in index order, so the result does not depend on thread scheduling.
	 *  \return  The total potential energy. */
	virtual double Energy();

	/** \brief Compute force on a particle.
	 *  \param result  Output force vector.
	 *  \param i       Particle index. */
	virtual void Force(GeometryVector & result, size_t i);

	/** \brief Compute forces on all particles.
	 *  Fills every species' rho first, so it is valid without a prior Energy call (MD order); with
	 *  FINUFFT the species are filled concurrently inside the same parallel region as the force loop.
	 *  \param results  Output vector of per-particle forces. */
	virtual void AllForce(std::vector<GeometryVector> & results);
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure);

};


#endif 
