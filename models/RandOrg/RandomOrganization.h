#ifndef RANDOMORGANIZATION_H
#define RANDOMORGANIZATION_H

#include <vector>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <omp.h>

//Random organization : Bertheir's model

//----------------For square systems

/** \brief Gives a list of active particles.
 *  \param List       Particle configuration as a list of position vectors.
 *  \param result     Output list of active-particle indices.
 *  \param exclusionR Exclusion radius below which two particles are active.
 *  \param SystSize   Linear size of the (square) periodic system. */
void IsActive(std::vector<GeometryVector> & List, std::vector<int> & result, double exclusionR, double SystSize);


/** \brief Berthier's model: displace active particles by a random kick.
 *  \param List       Particle configuration as a list of position vectors (modified in place).
 *  \param exclusionR Exclusion radius below which two particles are active.
 *  \param delta      Amplitude of the random displacement applied to active particles.
 *  \param SystSize   Linear size of the (square) periodic system.
 *  \return           true if the system is in the active phase, false if passive. */
bool B_Model(std::vector<GeometryVector> & List, double exclusionR, double delta, double SystSize);

/** \brief Berthier's model using a caller-supplied random number generator.
 *  \param List Particle configuration as a list of position vectors, modified in place.
 *  \param exclusionR Separation below which two particles are active.
 *  \param delta Magnitude of the random displacement applied to each active particle.
 *  \param SystSize Linear size of the square periodic system.
 *  \param g Random number generator used to draw displacement directions.
 *  \return true if at least one particle was active, false otherwise. */
bool B_Model(std::vector<GeometryVector> & List, double exclusionR, double delta, double SystSize, RandomGenerator & g);



//------------------For general periodic systems
/** \brief Check active particles.
 *  \param List           particle configuration.
 *  \param ActivePrt_index clear it and write a list of particle indices of active particles.
 *  \param exclusionR     exclusion radius. */
void IsActive(const Configuration & List, std::vector<size_t> & ActivePrt_index, double exclusionR);

/** \brief Find overlapping particles in a polydisperse periodic system.
 *  Each particle characteristic is interpreted as its exclusion radius.
 *  \param List Particle configuration to inspect.
 *  \param ActivePrt_index Cleared and filled with the indices of active particles.
 *  \param exclusionR Radius used to set the neighbor-search bound; it must be at least the largest particle exclusion radius. */
void IsActive(const PeriodicCellList<double> &List, std::vector<size_t> & ActivePrt_index, double exclusionR);

/** \brief Run Bertheir's model.
 *  \param List            particle configuration (modified in place).
 *  \param delta           amplitude of the random displacement of active particles in the relative coordinates.
 *  \param activePrt_index move particles in activePrt_index by the given random displacement.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(Configuration &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);

/** \brief Run Berthier's model for polydisperse systems.
 *  \param List            Polydisperse particle configuration whose characteristic stores each particle's radius (modified in place).
 *  \param delta           Amplitude of the random displacement, scaled by each active particle's radius.
 *  \param activePrt_index Indices of the active particles to be displaced.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(PeriodicCellList<double> &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);

/** \brief This is the combination of the above two methods (detect active particles, then displace them).
 *  \param List            particle configuration (modified in place).
 *  \param exclusionR      radius of the exclusion region.
 *  \param delta           Amplitude of the random displacement of active particles.
 *  \param activePrt_index clear index and write a list of particle indices of active particles.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(Configuration & List, double exclusionR, double delta, std::vector<size_t> & activePrt_index, RandomGenerator &g);



#endif