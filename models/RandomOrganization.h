#ifndef RANDOMORGANIZATION_H
#define RANDOMORGANIZATION_H

#include <vector>
#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include "RandomGenerator.h"
#include <omp.h>

//Random organization : Bertheir's model

//----------------For square systems

//Gives a list of active particles
void IsActive(std::vector<GeometryVector> & List, std::vector<int> & result, double exclusionR, double SystSize);


//Berthier's model
//true : active phase
//false : passive phase
bool B_Model(std::vector<GeometryVector> & List, double exclusionR, double delta, double SystSize);



//------------------For general periodic systems
// Check active particles
// List: particle configuration
// ActivePrt_index : clear it and write a list of particle indices of active particles
// exclusionR : exclusion radius
void IsActive(const Configuration & List, std::vector<size_t> & ActivePrt_index, double exclusionR);

void IsActive(const PeriodicCellList<double> &List, std::vector<size_t> & ActivePrt_index, double exclusionR);
//Run Bertheir's model : 
//delta : amplitude of the random displacement of active particles in the relative coordinates
//move particles in activePrt_index by the given random displacement.
void B_model(Configuration &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);

//polydisperse systems
void B_model(PeriodicCellList<double> &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g);
//returns the number of active particles
//List : particle configuration
//exclusion R : radius of the exclusion region
//activePrt_index : clear index and write a lits of particle indices of active particles.

//This is the combination of the above two methods
void B_model(Configuration & List, double exclusionR, double delta, std::vector<size_t> & activePrt_index, RandomGenerator &g);



#endif