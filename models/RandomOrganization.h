#ifndef RANDOMORGANIZATION_H
#define RANDOMORGANIZATION_H

#include <vector>
#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include "RandomGenerator.h"
#include <omp.h>

//Random organization : Bertheir's model

//----------------For square systems

/** \brief Gives a list of active particles.
 *  \param List       Particle configuration as a list of position vectors.
 *  \param result     Output list of active-particle indices.
 *  \param exclusionR Exclusion radius below which two particles are active.
 *  \param SystSize   Linear size of the (square) periodic system. */
void IsActive(std::vector<GeometryVector> & List, std::vector<int> & result, double exclusionR, double SystSize);


/** \brief Berthier's model: displace active particles by a random kick. (CLAUDE)
 *  \param List       Particle configuration as a list of position vectors (modified in place).
 *  \param exclusionR Exclusion radius below which two particles are active.
 *  \param delta      Amplitude of the random displacement applied to active particles.
 *  \param SystSize   Linear size of the (square) periodic system.
 *  \return           true if the system is in the active phase, false if passive. */
bool B_Model(std::vector<GeometryVector> & List, double exclusionR, double delta, double SystSize);



//------------------For general periodic systems
/** \brief Check active particles.
 *  \param List           particle configuration.
 *  \param ActivePrt_index clear it and write a list of particle indices of active particles.
 *  \param exclusionR     exclusion radius. */
void IsActive(const Configuration & List, std::vector<size_t> & ActivePrt_index, double exclusionR);

/** \brief Check active particles in a polydisperse periodic system. (CLAUDE)
 *  \param List            Particle configuration whose characteristic stores each particle's exclusion radius.
 *  \param ActivePrt_index Cleared and filled with the indices of active particles.
 *  \param exclusionR      Exclusion radius used to bound the neighbor search. */
void IsActive(const PeriodicCellList<double> &List, std::vector<size_t> & ActivePrt_index, double exclusionR);
/** \brief Run Bertheir's model.
 *  \param List            particle configuration (modified in place).
 *  \param delta           amplitude of the random displacement of active particles in the relative coordinates.
 *  \param activePrt_index move particles in activePrt_index by the given random displacement.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(Configuration &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);

/** \brief Run Berthier's model for polydisperse systems. (CLAUDE)
 *  \param List            Polydisperse particle configuration whose characteristic stores each particle's radius (modified in place).
 *  \param delta           Amplitude of the random displacement, scaled by each active particle's radius.
 *  \param activePrt_index Indices of the active particles to be displaced.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(PeriodicCellList<double> &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);

/** \brief This is the combination of the above two methods (detect active particles, then displace them). (CLAUDE)
 *  \param List            particle configuration (modified in place).
 *  \param exclusionR      radius of the exclusion region.
 *  \param delta           Amplitude of the random displacement of active particles.
 *  \param activePrt_index clear index and write a list of particle indices of active particles.
 *  \param g               Random number generator used to draw the displacement directions. */
void B_model(Configuration & List, double exclusionR, double delta, std::vector<size_t> & activePrt_index, RandomGenerator &g);



#endif