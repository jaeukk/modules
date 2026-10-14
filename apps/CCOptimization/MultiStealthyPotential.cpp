#include "MultiStealthyPotential.h"
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <climits>
#include <thread>

SpherePacking MultiStealthyPotential::GetPackingSnapshot() const
{
	if (!pConfig || radii.size() != pConfig->NumParticle())
		throw std::logic_error("GetPackingSnapshot: missing configuration or radii.");
	SpherePacking packing(*pConfig, 0.0);
	for (size_t i = 0; i < radii.size(); ++i) {
		if (!std::isfinite(radii[i]) || radii[i] <= 0.0)
			throw std::logic_error("GetPackingSnapshot: invalid radius.");
		packing.GetCharacteristics(i) = radii[i];
	}
	return packing;
}

void MultiStealthyPotential::ExportConfiguration(std::ostream &ofile, ExportFormat format) const
{
	if (format != ExportFormat::Text && format != ExportFormat::Binary)
		throw std::invalid_argument("ExportConfiguration: unknown format.");
	if (!pConfig)
		throw std::logic_error("ExportConfiguration: no configuration has been set.");
	const size_t count = pConfig->NumParticle();
	if (idx2species.size() != count || radii.size() != count)
		throw std::logic_error("ExportConfiguration: inconsistent particle metadata.");
	for (size_t i = 0; i < count; ++i) {
		if (idx2species[i].species < 0 ||
			static_cast<size_t>(idx2species[i].species) >= species_constraints.size())
			throw std::logic_error("ExportConfiguration: particle has no valid species.");
		if (!std::isfinite(radii[i]) || radii[i] <= 0.0)
			throw std::logic_error("ExportConfiguration: particle has an invalid radius.");
	}

	if (format == ExportFormat::Binary) {
		static_assert(CHAR_BIT == 8 && sizeof(double) == 8 &&
			std::numeric_limits<double>::is_iec559, "MSPB requires IEEE-754 binary64.");
		// Encode individual fields explicitly: no native struct padding or endianness.
		std::string bytes("MSPB\r\n\x1a\n", 8);
		auto integer = [&bytes](std::uint64_t value, unsigned int width) {
			for (unsigned int b = 0; b < width; ++b) {
				const unsigned char octet = static_cast<unsigned char>(value & 0xff);
				bytes.append(reinterpret_cast<const char *>(&octet), 1);
				value >>= 8;
			}
		};
		auto real = [&integer](double value) {
			if (!std::isfinite(value))
				throw std::logic_error("ExportConfiguration: nonfinite binary value.");
			std::uint64_t bits;
			std::memcpy(&bits, &value, sizeof(bits));
			integer(bits, 8);
		};
		const auto dimension = pConfig->GetDimension();
		integer(1, 4); // Format version.
		integer(dimension, 4);
		integer(count, 8);
		integer(species_constraints.size(), 8);
		integer(3, 4); // Bit 0: relative coordinates; bit 1: Cartesian radii.
		for (DimensionType d = 0; d < dimension; ++d) {
			const auto basis = pConfig->GetBasisVector(d);
			for (DimensionType axis = 0; axis < dimension; ++axis) real(basis.x[axis]);
		}
		for (size_t i = 0; i < count; ++i) {
			integer(i, 8);
			integer(static_cast<std::uint32_t>(idx2species[i].species), 4);
			const auto &s = pConfig->GetRelativeCoordinates(i);
			for (DimensionType d = 0; d < dimension; ++d) real(s.x[d]);
			real(radii[i]);
		}
		ofile.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
		if (!ofile) throw std::runtime_error("ExportConfiguration: output failed.");
		return;
	}

	// Format in a private buffer: preserve caller flags, precision, and locale.
	std::ostringstream table;
	table.imbue(std::locale::classic());
	table << std::setprecision(std::numeric_limits<double>::max_digits10);
	const auto dimension = pConfig->GetDimension();
	table << "# dimension " << static_cast<unsigned int>(dimension) << '\n';
	table << "# particles " << count << '\n';
	for (DimensionType d = 0; d < dimension; ++d) {
		const auto basis = pConfig->GetBasisVector(d);
		table << "# basis_vector-" << static_cast<unsigned int>(d);
		for (DimensionType axis = 0; axis < dimension; ++axis) {
			if (!std::isfinite(basis.x[axis]))
				throw std::logic_error("ExportConfiguration: nonfinite basis vector.");
			table << ' ' << basis.x[axis];
		}
		table << '\n';
	}
	table << "# coordinates relative\n# radius_units cartesian\n";
	table << "# id\tspecies";
	for (DimensionType d = 0; d < dimension; ++d) {
		if (d < 3) table << "\ts" << "xyz"[d];
		else table << "\ts" << static_cast<unsigned int>(d);
	}
	table << "\tradius\n";
	for (size_t i = 0; i < count; ++i) {
		const auto &s = pConfig->GetRelativeCoordinates(i);
		table << i << '\t' << idx2species[i].species;
		for (DimensionType d = 0; d < dimension; ++d) {
			if (!std::isfinite(s.x[d]))
				throw std::logic_error("ExportConfiguration: nonfinite coordinate.");
			table << '\t' << s.x[d];
		}
		table << '\t' << radii[i] << '\n';
	}
	const std::string text = table.str();
	ofile.write(text.data(), static_cast<std::streamsize>(text.size()));
	if (!ofile)
		throw std::runtime_error("ExportConfiguration: output failed.");
}

/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: November 2025 */

/** \file MultiStealthyPotential.cpp 
 * \brief c++ implementations of multi-stealthy potential */



/** \brief Assign a particle to a species and refresh bookkeeping indices.
 *  \param prt_idx Global particle index in the full configuration.
 *  \param species Target species id. */
void MultiStealthyPotential::SetSpeciesOfParticle(size_t prt_idx, size_t species)
{
	int oldspecies = this->idx2species[prt_idx].species;

	/* Already in this species: nothing to do. Appending unconditionally (as before) would
	   duplicate prt_idx in the species index list and push idx2species[].index past the end. */
	if (oldspecies == static_cast<int>(species))
		return;

	if (oldspecies >= 0){
		/* Remove prt_idx from its previously assigned species. */
		std::vector<size_t> & oldidx = this->species_constraints[oldspecies].idx;
		oldidx.erase(std::remove(oldidx.begin(), oldidx.end(), prt_idx), oldidx.end());
		/* The erase shifts every later particle in that species down by one, so their cached
		   within-species index is now stale; rebuild the old species' reverse indices. */
		for (size_t j = 0; j < oldidx.size(); j++)
			this->idx2species[oldidx[j]].index = j;
	}
	else {
		/* species were not specified so far. */
		std::cout << "species of particle " << prt_idx << " is initialized to " << species << "\n";
	}

	/* Append prt_idx to the new species and record its within-species index. */
	this->species_constraints[species].idx.emplace_back(prt_idx);
	this->idx2species[prt_idx].species = static_cast<int>(species);
	this->idx2species[prt_idx].index = this->species_constraints[species].idx.size()-1;

	this->SubconfigValid = false;
}



/** \brief Append a Fourier-space constraint to a given species.
 *  \param species Species id to target.
 *  \param k Wavevector for the constraint.
 *  \param V List of constraint coefficients. */
void MultiStealthyPotential::AddConstraint(size_t species, const GeometryVector &k, std::vector<double> V)
{
	this->species_constraints[species].CCPotential->AddConstraint(k, V);
	this->species_constraints[species].CCPotential->RhoValid = false;
}

/** \brief Compute the structure factor contributions for a species.
 *  \param s Species id to evaluate. */
void MultiStealthyPotential::GetRho(size_t s)
{
	assert(this->pConfig!=nullptr);
	assert(s < this->species_constraints.size());
	this->UpdateSubconfigs();
	this->species_constraints[s].CCPotential->GetRho();

	// double Volume=ptr_pot->pConfig->PeriodicVolume();

	// size_t end = this->species_constraints[s].constraints.size();
	// size_t num_in_species = this->species_constraints[s].idx.size();
	// #pragma omp parallel for num_threads(this->ParallelNumThread)
	// for (size_t i = 0; i < end; i++){
	// 	/* At a wavevector for the constraint*/
	// 	double tr=0, ti=0;//real and imag parts of result
		
	// 	for(size_t j = 0; j < num_in_species; j++)
	// 	{
	// 		/* repeat over particles in the same species */
	// 		size_t p_idx = this->species_constraints[s].idx[j];
	// 		GeometryVector x=pConfig->GetCartesianCoordinates(p_idx);
	// 		tr+=std::cos(this->species_constraints[s].constraints[j].k.Dot(x));
	// 		ti+=std::sin(this->species_constraints[s].constraints[j].k.Dot(x));
	// 	}
	// 	this->species_constraints[s].constraints[i].rho_real=tr;
	// 	this->species_constraints[s].constraints[i].rho_imag=ti;
	// }
	// this->species_constraints[s].RhoValid = true;
}

/** \brief Print all species constraints and indices.
 *  \param ofile Destination output stream. */
void MultiStealthyPotential::PrintConstraints(std::ostream &ofile)
{
	for (size_t s = 0; s < this->species_constraints.size(); s ++ ){
		ofile << "speciecs "<< s << "\n";
		ofile << "\t (particle index, radius)";
		for (size_t i = 0; i < this->species_constraints[s].idx.size(); i++){
			size_t idx = this->species_constraints[s].idx[i];
			ofile << "\t(" << idx << ",\t" << this->radii[idx] <<")\n"; //this->pConfig->GetCharacteristics(idx) <<")\n";
		}

		ofile << "\n";
		for (size_t i =0; i< this->species_constraints[s].CCPotential->constraints.size(); i++){
			ofile << "@ k= " << this->species_constraints[s].CCPotential->constraints[i].k << ": v=" << this->species_constraints[s].CCPotential->constraints[i].V <<"\n";
		}
		ofile << "----- \n\n";
	}
}

/** \brief Set the working configuration from particle positions.
 *  \param Position Full particle configuration. */
void MultiStealthyPotential::SetConfiguration(const Configuration &Position)
{
	if (this->rc < 0.0){
		// an empty configuration
		this->idx2species.resize(Position.NumParticle(), spec_index(-1,0));
		std::cout << "species of all particles are set to be -1 (unspecified)\n";
		std::cout << "the critical radius of the contact potential is set to be "<< this->rc <<std::endl;

		std::cout << "No input for particle radii\n";
		std::cout << "Particle radii are set to be 0.5 by default.\n";	
		this->rc = 1.0;
		if (this->radii.size() != Position.NumParticle()){
			//default radius 0.5 (so a pair diameter radii[i]+radii[j]=1.0=rc); storing -1.0 gave
			//negative pair diameters and thus nonsensical contact energies/forces.
			this->radii.resize(Position.NumParticle(), 0.5);
		}
	}
	
	// Invalidate the species views before their configuration is replaced.
	this-> SubconfigValid = false;
	pConfig = std::make_unique<Configuration>(Position, "a");
	this->PrepareContactCells();

	this-> UpdateSubconfigs();
	for (size_t s = 0; s < this->species_constraints.size(); s++){
		this->species_constraints[s].CCPotential->RhoValid = false;
	}
}

/** \brief Size the contact cell list to the contact cutoff and prepare the neighbor search.
 *  The copied configuration keeps the caller's cells, e.g. 10 per side for random packings built from
 *  GetUnitCubicBox(d, 0.1) and rescaled, which makes every contact query scan about 0.09 N particles
 *  (an O(N^2) contact sum). Cells of width max(rc, (V/N)^(1/d)) hold O(1) particles each, and the second
 *  term keeps the cell count at most about N when radii are tiny. */
void MultiStealthyPotential::PrepareContactCells()
{
	const double n = std::max<double>(1.0, double(this->pConfig->NumParticle()));
	const double spacing = std::pow(this->pConfig->PeriodicVolume() / n, 1.0 / this->Dimension);
	this->pConfig->SetCellSize(std::max(this->rc, spacing));
	this->pConfig->PrepareIterateThroughNeighbors(this->rc);
}

/** \brief Set particle radii for a configuration.
 *  \param radii Radii list matching particles. */
void MultiStealthyPotential::SetSecureRadii(const std::vector<double> &radii)
{
	assert(radii.size() == this->pConfig->NumParticle());
	for (double radius : radii)
		if (!std::isfinite(radius) || radius <= 0.0)
			throw std::invalid_argument("MultiStealthyPotential::SetSecureRadii: radii must be finite and positive.");
	this->radii = std::vector<double> (radii);
	this->rc = radii.empty() ? 0.0 : 2.01 * *std::max_element(radii.begin(), radii.end());
	this->PrepareContactCells();
}

/** \brief Initialize configuration and radii from a packing object.
 *  \param pck Input sphere packing. */
void MultiStealthyPotential::SetPacking(const SpherePacking &pck)
{
	if (this->rc < 0.0){
		// an empty configuration
		this->idx2species.resize(pck.NumParticle(),spec_index(-1,0));
		std::cout << "species of all particles are set to be -1 (unspecified)\n";
		std::cout << "the critical radius of the contact potential is set to be "<< this->rc <<std::endl;
	}
	this->radii.resize(pck.NumParticle(), -1.0);
	this->rc = 0.0;
	for(size_t i =0; i < pck.NumParticle(); i++){
		this->rc = std::max(this-> rc , pck.GetCharacteristics(i));
		this->radii[i] = pck.GetCharacteristics(i);
	}
	this->rc *= 2.01;
	// Invalidate the species views before their configuration is replaced.
	this-> SubconfigValid = false;
	pConfig = std::make_unique<Configuration>(pck, "a");
	this->PrepareContactCells();
	
	/* bind the species views. */
	this-> UpdateSubconfigs();
	for (size_t s = 0; s < this->species_constraints.size(); s++){
		this->species_constraints[s].CCPotential->RhoValid = false;
	}
}

/** \brief Rebind every species potential, if stale, to an indexed view of *pConfig.
 *  Replaces the former per-species sub-Configuration copies with the same particles and order, so
 *  results are unchanged. The views point into *pConfig: every replacement of pConfig must first
 *  clear SubconfigValid, and every read of a species potential must first call this function. */
void MultiStealthyPotential::UpdateSubconfigs()
{
	if (! this->SubconfigValid){
		for (size_t s = 0; s < this->species_constraints.size(); s ++ ){
			//propagate the thread count to each species potential; otherwise only the wrapper is
			//multi-threaded and every inner ShiftedCCPotential stays single-threaded.
			this->species_constraints[s].CCPotential->ParallelNumThread = this->ParallelNumThread;
			// The local particle order must match idx and idx2species, including reassignment.
			this->species_constraints[s].CCPotential->SetIndexedConfiguration(*this->pConfig, this->species_constraints[s].idx);
		}
		this->SubconfigValid = true;
	}
}

/** \brief Prepare every species' FINUFFT plans serially, mark the species evaluated concurrently and list their blocks.
 *  \param forces Also decide the force paths (FinufftForceGate over the concurrent species) and prepare them.
 *  \return Flags, block tasks (largest predicted time first; ties by species, then block) and countdowns. */
MultiStealthyPotential::SpeciesWork MultiStealthyPotential::PrepareSpecies(bool forces)
{
	const size_t numSpecies = this->species_constraints.size();
	SpeciesWork work;
	work.concurrent.assign(numSpecies, 0);
	std::vector<const ShiftedCCPotential *> pots;
	for (size_t s = 0; s < numSpecies; s++) {
		work.concurrent[s] = this->species_constraints[s].CCPotential->PrepareFinufft();
		if (work.concurrent[s])
			pots.push_back(this->species_constraints[s].CCPotential.get());
	}
	std::vector<size_t> forceBlocks(numSpecies, 0);
	if (forces) {
		// pots holds FINUFFT species only, so the dimension is 1-3 when the contact estimate is read.
		const double contact = this->WithRepulsion && !pots.empty() ? ContactCostParticle[this->Dimension - 1] * double(this->pConfig->NumParticle()) : 0.0;
		const std::vector<size_t> type2 = ShiftedCCPotential::FinufftForceGate(pots, this->ParallelNumThread, contact);
		for (size_t s = 0, j = 0; s < numSpecies; s++)
			forceBlocks[s] = this->species_constraints[s].CCPotential->PrepareFinufftForces(work.concurrent[s] ? type2[j++] : 0);
	}
	work.forceBlocks.assign(numSpecies, 0);
	work.rhoLeft.assign(numSpecies, 0);
	work.forceLeft.assign(numSpecies, 0);
	work.state.assign(numSpecies, 0);
	for (size_t s = 0; s < numSpecies; s++) {
		if (!work.concurrent[s])
			continue;
		ShiftedCCPotential & pot = *this->species_constraints[s].CCPotential;
		const long sp = static_cast<long>(s);
		const size_t rhoBlocks = pot.StartRho(forces);
		work.forceBlocks[s] = forceBlocks[s];
		// The species' last phase-1 task also runs a single Type 2 block; several go to phase 2.
		const double tail = forceBlocks[s] == 1 && rhoBlocks <= 1 ? pot.FinufftBlockTime(2) : 0.0;
		for (size_t b = 0; b < std::max<size_t>(rhoBlocks, 1); b++)
			work.rho.push_back({sp, rhoBlocks == 0 ? SpeciesWork::NoBlock : b, pot.FinufftBlockTime(1) + tail});
		work.rhoLeft[s] = static_cast<long>(std::max<size_t>(rhoBlocks, 1));
		if (forceBlocks[s] > 1) {
			for (size_t b = 0; b < forceBlocks[s]; b++)
				work.force.push_back({sp, b, pot.FinufftBlockTime(2)});
			work.forceLeft[s] = static_cast<long>(forceBlocks[s]);
		}
	}
	auto largestFirst = [](const SpeciesWork::Task & x, const SpeciesWork::Task & y) {
		if (x.cost != y.cost)
			return x.cost > y.cost;
		return x.species != y.species ? x.species < y.species : x.block < y.block;
	};
	std::sort(work.rho.begin(), work.rho.end(), largestFirst);
	std::sort(work.force.begin(), work.force.end(), largestFirst);
	return work;
}

/** \brief Run one call on a species potential with its ParallelNumThread set to 1, restored afterwards.
 *  A species evaluated on one thread of a parallel region then opens no nested team, and its energy partial
 *  sums use a single slot, so the result does not depend on the team size.
 *  \param pot  Species potential.
 *  \param call Callable taking the potential.
 *  \return What call returns. */
template <class Call>
static auto OnOneThread(ShiftedCCPotential & pot, Call call)
{
	unsigned long threads = 1;
	std::swap(threads, pot.ParallelNumThread);
	const auto result = call(pot);
	std::swap(threads, pot.ParallelNumThread);
	return result;
}

/** \brief The concurrent FINUFFT species of one evaluation, called by every thread of a parallel region.
 *  \param work          PrepareSpecies result; its shared counters are used.
 *  \param speciesEnergy Energy slot of each species to fill, or nullptr.
 *  \param speciesForce  Type 2 force buffer of each species, or nullptr for no forces.
 *  \param bulk          Output (with forces): nonzero for each species whose Type 2 forces succeeded.
 *  \param items         Number of contact items.
 *  \param chunk         Contact items per chunk.
 *  \param contact       Callable run once per contact item index. */
template <class Contact>
void MultiStealthyPotential::RunSpeciesWork(SpeciesWork & work, std::vector<double> * speciesEnergy,
	std::vector<std::vector<GeometryVector>> * speciesForce, std::vector<char> & bulk, long items, long chunk,
	const Contact & contact)
{
	auto pot = [this](long s) -> ShiftedCCPotential & { return *this->species_constraints[s].CCPotential; };
	/* One contact chunk from the shared counter; false when none is left. */
	auto chunkDone = [&]() {
		long next;
#pragma omp atomic capture
		next = work.nextChunk++;
		const long lo = next * chunk;
		for (long i = lo; i < std::min(items, lo + chunk); i++)
			contact(i);
		return lo < items;
	};
	const long numRho = static_cast<long>(work.rho.size()), numForce = static_cast<long>(work.force.size());
	/* After the last Type 1 task of species s, on this thread: rho, energy, Type 2 input, and a single Type 2 block;
	   species with several Type 2 blocks are released to the threads waiting for them (state 1, or 2 without Type 2). */
	auto finishSpecies = [&](long s) {
		const bool ready = OnOneThread(pot(s), [&](ShiftedCCPotential & p) {
			p.FinishRho();
			if (speciesEnergy)
				(*speciesEnergy)[s] = p.Energy();
			return speciesForce && p.StartForces((*speciesForce)[s]) > 0;
		});
		if (work.forceBlocks[s] > 1) {
			const int state = ready ? 1 : 2;
#pragma omp atomic write seq_cst
			work.state[s] = state;
		}
		else if (ready) {
			pot(s).RunForceBlock(0, (*speciesForce)[s]);
			bulk[s] = pot(s).FinishForces();
		}
	};
	/* One dynamic loop over the Type 1 tasks, then the Type 2 tasks; monotonic, so a thread runs its tasks in order and
	   holds no unstarted Type 1 task while it waits: a thread holding a Type 2 block of a species not yet released
	   takes contact chunks (or yields) meanwhile, and the release cannot wait on it. seq_cst countdowns and releases
	   hand the species state between threads. */
#pragma omp for schedule (monotonic: dynamic, 1) nowait
	for (long t = 0; t < numRho + numForce; t++) {
		if (t < numRho) {
			const long s = work.rho[t].species;
			if (work.rho[t].block != SpeciesWork::NoBlock)
				pot(s).RunRhoBlock(work.rho[t].block);
			long left;
#pragma omp atomic capture seq_cst
			left = --work.rhoLeft[s];
			if (left == 0)
				finishSpecies(s);
			continue;
		}
		const long s = work.force[t - numRho].species;
		int state;
		for (;;) {
#pragma omp atomic read seq_cst
			state = work.state[s];
			if (state != 0)
				break;
			if (!chunkDone())
				std::this_thread::yield();
		}
		if (state != 1)
			continue;
		pot(s).RunForceBlock(work.force[t - numRho].block, (*speciesForce)[s]);
		long left;
#pragma omp atomic capture seq_cst
		left = --work.forceLeft[s];
		if (left == 0)
			bulk[s] = pot(s).FinishForces();
	}
	while (chunkDone())
		;
#pragma omp barrier
}

/** \brief Total force on particle i from its contact force and its species force (Force(i) arithmetic).
 *  \param i            Particle index.
 *  \param contactForce e0 times the contact force on particle i.
 *  \param bulk         Per species: nonzero if speciesForce holds its forces.
 *  \param speciesForce Type 2 forces per species.
 *  \return The total force on particle i. */
GeometryVector MultiStealthyPotential::WithSpeciesForce(size_t i, const GeometryVector & contactForce, const std::vector<char> & bulk,
	const std::vector<std::vector<GeometryVector>> & speciesForce)
{
	const spec_index & at = this->idx2species[i];
	GeometryVector force_CC(static_cast<int> (this->Dimension));
	if (bulk[at.species])
		force_CC = speciesForce[at.species][at.index];
	else
		this->species_constraints[at.species].CCPotential->Force(force_CC, at.index);
	force_CC.MultiplyFrom(this->u0 * (double)this->species_constraints[at.species].idx.size()/(double)this->pConfig->NumParticle());
	return contactForce + force_CC;
}

/** \brief Compute total energy from stealthy and repulsive contributions.
 *  \return Total potential energy. */
double MultiStealthyPotential::Energy()
{
	this->UpdateSubconfigs();
	//every particle must belong to a species; an unassigned particle (species -1) is not in any
	//species index list and would be silently dropped from the collective-coordinate sum below.
	for (size_t i = 0; i < this->idx2species.size(); i++)
		if (this->idx2species[i].species < 0)
			throw std::runtime_error("MultiStealthyPotential::Energy: particle has no assigned species; assign every particle before evaluating the energy.");

	// Keep arithmetic order independent of worker assignment and reduction order.
	// Each species and each particle owns one slot; only the final, ordered sums combine slots.
	const long numSpecies = static_cast<long>(this->species_constraints.size());
	const long numContact = this->WithRepulsion ? static_cast<long>(this->pConfig->NumParticle()) : 0;
	std::vector<double> speciesEnergy(numSpecies, 0.0), contactEnergy(numContact, 0.0);
	SpeciesWork work = this->PrepareSpecies(false);
	/* Direct or FINUFFT fallback: stealthy potential contribution (CC) species by species, each parallel over wavevectors. */
	for (long s = 0; s < numSpecies; s++)
		if (!work.concurrent[s])
			speciesEnergy[s] = this->species_constraints[s].CCPotential->Energy();
	/* compute the contact potential contribution (CP) of particle i */
	auto contact = [this](size_t i) {
		double E = 0.0;
		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc,
			[&i, &E, this](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t SourceParticle)->void {
			/*Prevent double-counting */
			if (i < SourceParticle) {
				double r = sqrt(shift.Modulus2()) / (radii[i] + radii[SourceParticle]);
				if (r < 1.0){
					E += (1.0 - r)*(1.0 - r);
				}
			}
			}
		);
		return E;
	};
	if (work.rho.empty()) {
		/* No species runs concurrently: the contact sum as before, statically split over the team. */
		if (numContact > 0) {
#pragma omp parallel for schedule (static) num_threads(this->ParallelNumThread)
			for (long i = 0; i < numContact; i++)
				contactEnergy[i] = contact(i);
		}
	}
	else {
		/* FINUFFT: Type 1 block tasks, while free threads take the contact contribution (CP).
		   Keep the team at ParallelNumThread even when there are fewer tasks: a team size that differs from
		   AllForce's makes libgomp destroy and recreate threads on every evaluation. About four contact
		   chunks per thread let the threads freed by the block tasks share the rest, even for small N;
		   the cap of 64 keeps large N on fine-grained chunks (N/(4T) was 6-18% slower at N=1e4). */
		const long chunk = std::clamp(numContact / (4L * static_cast<long>(std::max<unsigned long>(this->ParallelNumThread, 1))), 1L, 64L);
		std::vector<char> bulk(numSpecies, 0);
#pragma omp parallel num_threads(this->ParallelNumThread)
		this->RunSpeciesWork(work, &speciesEnergy, nullptr, bulk, numContact, chunk,
			[&](long i) { contactEnergy[i] = contact(i); });
	}
	return this->TotalEnergy(speciesEnergy, contactEnergy);
}

/** \brief Total energy from per-species and per-particle slots, summed in index order.
 *  \param speciesEnergy Energy of each species potential.
 *  \param contactEnergy Contact energy slot of each particle (empty when repulsion is off).
 *  \return u0 times the sum of N_s/N times the species energies, plus e0 times the contact sum. */
double MultiStealthyPotential::TotalEnergy(const std::vector<double> &speciesEnergy, const std::vector<double> &contactEnergy) const
{
	double CC = 0.0, CP = 0.0;
	for (size_t s = 0; s < speciesEnergy.size(); s++)
		CC += ((double)this->species_constraints[s].idx.size()/(double)this->pConfig->NumParticle())*speciesEnergy[s];
	for (double energy : contactEnergy)
		CP += energy;
	return this->u0 * CC + this->e0 * CP;
}

/** \brief Evaluate force on a single particle.
 *  \param force Output force vector.
 *  \param i Particle index. */
void MultiStealthyPotential::Force(GeometryVector &force, size_t i)
{
	this->UpdateSubconfigs();
	// stealthy potential contribution
	GeometryVector force_CC(static_cast<int> (this->Dimension));
	//an unassigned particle keeps species -1; casting that to size_t would index species_constraints
	//out of bounds. Require every particle to have a species before force evaluation.
	if (this->idx2species[i].species < 0)
		throw std::runtime_error("MultiStealthyPotential::Force: particle has no assigned species; assign every particle before evaluating forces.");
	size_t s = (size_t) this->idx2species[i].species; // species
	size_t idx_in_s = this->idx2species[i].index;
	this->species_constraints[s].CCPotential -> Force(force_CC, idx_in_s);

	force_CC.MultiplyFrom(this->u0 * (double)this->species_constraints[s].idx.size()/(double)this->pConfig->NumParticle());

	// EnergyConstraints * ptr_s = & (this->species_constraints[s]);
	// if(ptr_s->RhoValid == false){
	// 	this->GetRho(s);
	// }
	
	// // iterate over wavevectors 
	// for (int q = 0; q < ptr_s -> constraints.size(); q++){
	// 	GeometryVector k = ptr_s->constraints[q].k;
	// 	GeometryVector r = this->pConfig->GetCartesianCoordinates(i);

	// 	force_CC = force_CC + k* ptr_s->constraints[q].V * (ptr_s->constraints[q].rho_real * sin(k.Dot(r))  +  ptr_s->constraints[q].rho_imag * cos(k.Dot(r)));
	// }
	// force_CC.MultiplyFrom(this->u0/this->Volume);

	// contact potential contribution (CP)
	force = this->ContactForce(i) + force_CC;
}

/** \brief Soft-core contact force on one particle.
 *  \param i Particle index.
 *  \return e0 times the contact force on particle i (zero when repulsion is off). */
GeometryVector MultiStealthyPotential::ContactForce(size_t i) const
{
	GeometryVector force_CP(static_cast<int> (this->Dimension));
	if (this->WithRepulsion){
		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc,
			[&i, &force_CP, this](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t SourceParticle)->void {
			/*Prevent self-interaction */
			if (i != SourceParticle) {
				double r = sqrt(shift.Modulus2());
				double sigma = this->radii[i] + this->radii[SourceParticle];

				if (r < 1e-15){
					/* two particles are completely overlapped.*/
					GeometryVector x(shift);
					x.x[0] = (i > SourceParticle)? -1.0 : 1.0 ;
					force_CP = force_CP - 2.0 / sigma * x;
				}
				else if ( r < sigma ){
					force_CP = force_CP + 2.0 / sigma * (1.0 / r - 1.0 / sigma) * (-1.*shift);
				}
			}
		});
	}
	force_CP.MultiplyFrom(this->e0);
	return force_CP;
}

/** \brief Contact energy slot and contact force of one particle from one neighbor traversal.
 *  The expressions are those of Energy's contact sum and of ContactForce, operation for operation.
 *  \param i Particle index.
 *  \param force Output: e0 times the contact force on particle i, as ContactForce(i).
 *  \return Contact energy of the pairs (i, j > i), before the factor e0 (zero when repulsion is off). */
double MultiStealthyPotential::ContactEnergyAndForce(size_t i, GeometryVector &force) const
{
	double E = 0.0;
	force = GeometryVector(static_cast<int> (this->Dimension));
	if (this->WithRepulsion){
		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc,
			[&i, &E, &force, this](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t SourceParticle)->void {
			/*Prevent self-interaction */
			if (i != SourceParticle) {
				double r = sqrt(shift.Modulus2());
				double sigma = this->radii[i] + this->radii[SourceParticle];

				/*Energy: prevent double-counting */
				if (i < SourceParticle) {
					double rs = r / sigma;
					if (rs < 1.0){
						E += (1.0 - rs)*(1.0 - rs);
					}
				}
				if (r < 1e-15){
					/* two particles are completely overlapped.*/
					GeometryVector x(shift);
					x.x[0] = (i > SourceParticle)? -1.0 : 1.0 ;
					force = force - 2.0 / sigma * x;
				}
				else if ( r < sigma ){
					force = force + 2.0 / sigma * (1.0 / r - 1.0 / sigma) * (-1.*shift);
				}
			}
		});
	}
	force.MultiplyFrom(this->e0);
	return E;
}

/** \brief Evaluate forces for all particles.
 *  Every species' FINUFFT plans are prepared serially first (PrepareSpecies). Species that do not run FINUFFT
 *  fill rho one by one, each parallel over wavevectors with the full team. Then one parallel region: species
 *  that run FINUFFT fill rho as Type 1 block tasks and, where the cost gate accepted, all their collective-coordinate
 *  forces as Type 2 block tasks, while free threads take the contact forces (RunSpeciesWork); after the barrier, the
 *  particle loop adds to each contact force the species force scaled by u0 N_s/N, from Type 2 or else from the direct
 *  sum of Force(i), so a declined or failed transform gives the previous direct results. Also valid without a prior
 *  Energy (MD order).
 *  \param results Output vector of forces. */
void MultiStealthyPotential::AllForce(std::vector<GeometryVector> &results)
{
	this->UpdateSubconfigs();
	for (size_t i = 0; i < this->idx2species.size(); i++)
		if (this->idx2species[i].species < 0)
			throw std::runtime_error("MultiStealthyPotential::AllForce: particle has no assigned species; assign every particle before evaluating forces.");
	// Plans are made here, serially; the parallel region below only executes them.
	SpeciesWork work = this->PrepareSpecies(true);
	const long numSpecies = static_cast<long>(this->species_constraints.size());
	for (long s = 0; s < numSpecies; s++)
		if (!work.concurrent[s])
			this->species_constraints[s].CCPotential->GetRho();
	long end = this->pConfig->NumParticle();
	results.resize(end);
	std::vector<char> bulk(numSpecies, 0);
	std::vector<std::vector<GeometryVector>> speciesForce(numSpecies);
	const long chunk = std::clamp(end / (4L * static_cast<long>(std::max<unsigned long>(this->ParallelNumThread, 1))), 1L, 64L);
	unsigned long tempNumThread = 1;
	std::swap(tempNumThread, this->ParallelNumThread);
	if (work.rho.empty()) {
#pragma omp parallel for num_threads(tempNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
	}
	else {
#pragma omp parallel num_threads(tempNumThread) default(shared)
		{
			this->RunSpeciesWork(work, nullptr, &speciesForce, bulk, end, chunk,
				[&](long i) { results[i] = this->ContactForce(i); });
#pragma omp for
			for (long i = 0; i<end; i++)
				results[i] = this->WithSpeciesForce(i, results[i], bulk, speciesForce);
		}
	}
	std::swap(tempNumThread, this->ParallelNumThread);
}

/** \brief Total energy and all forces in one evaluation, bitwise equal to Energy() then AllForce(results).
 *  Prepares the species once, as AllForce (PrepareSpecies(true); Energy's PrepareSpecies(false) is the same
 *  without the force setup, which leaves the Type 2 plans alone). Species that do not run FINUFFT get
 *  their energy (which fills rho) one by one with the full team, as in Energy. Then one parallel region: without
 *  FINUFFT species, one static particle loop gives each particle its contact slot and the force arithmetic of
 *  Force(i); with them, RunSpeciesWork runs their Type 1 blocks, energies and (gate permitting) rho-reusing Type 2
 *  blocks, as the two calls do, while free threads take the contact walk, and after the barrier the particle loop
 *  adds each species force as AllForce does. Slots are summed in index order (TotalEnergy).
 *  \param results Output vector of forces.
 *  \return Total potential energy. */
double MultiStealthyPotential::EnergyAndAllForce(std::vector<GeometryVector> &results)
{
	this->UpdateSubconfigs();
	for (size_t i = 0; i < this->idx2species.size(); i++)
		if (this->idx2species[i].species < 0)
			throw std::runtime_error("MultiStealthyPotential::EnergyAndAllForce: particle has no assigned species; assign every particle before evaluating.");
	// Plans are made here, serially; the parallel region below only executes them.
	SpeciesWork work = this->PrepareSpecies(true);
	const long numSpecies = static_cast<long>(this->species_constraints.size());
	const long end = static_cast<long>(this->pConfig->NumParticle());
	const long numContact = this->WithRepulsion ? end : 0;
	std::vector<double> speciesEnergy(numSpecies, 0.0), contactEnergy(numContact, 0.0);
	for (long s = 0; s < numSpecies; s++)
		if (!work.concurrent[s]) {
			speciesEnergy[s] = this->species_constraints[s].CCPotential->Energy();
			this->species_constraints[s].CCPotential->GetRho();	// as AllForce; a no-op once Energy filled rho
		}
	results.resize(end);
	std::vector<char> bulk(numSpecies, 0);
	std::vector<std::vector<GeometryVector>> speciesForce(numSpecies);
	/* contact slot of particle i, from the same walk as its contact force */
	auto contact = [&](long i, GeometryVector & force) {
		const double E = this->ContactEnergyAndForce(i, force);
		if (numContact > 0)
			contactEnergy[i] = E;
	};
	// As in Energy: team of ParallelNumThread; about four contact chunks per thread, at most 64 particles each.
	const long chunk = std::clamp(end / (4L * static_cast<long>(std::max<unsigned long>(this->ParallelNumThread, 1))), 1L, 64L);
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
	{
		if (work.rho.empty()) {
			// No FINUFFT species: rho is ready, so each particle is finished in one pass.
#pragma omp for schedule (static)
			for (long i = 0; i < end; i++) {
				GeometryVector contactForce;
				contact(i, contactForce);
				results[i] = this->WithSpeciesForce(i, contactForce, bulk, speciesForce);
			}
		}
		else {
			this->RunSpeciesWork(work, &speciesEnergy, &speciesForce, bulk, end, chunk,
				[&](long i) { contact(i, results[i]); });
#pragma omp for
			for (long i = 0; i < end; i++)
				results[i] = this->WithSpeciesForce(i, results[i], bulk, speciesForce);
		}
	}
	return this->TotalEnergy(speciesEnergy, contactEnergy);
}
