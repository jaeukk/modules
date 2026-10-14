/** \file CombinedEnergyForce.h
	\brief Optional interface of a potential that computes its energy and every force in one call. */

#ifndef COMBINEDENERGYFORCE_INCLUDED
#define COMBINEDENERGYFORCE_INCLUDED

#include <GeometryVector.h>

#include <vector>

/** \brief Optional second base class of a Potential that computes the total energy and every force in one call,
	for example with one neighbor traversal instead of two. The structure optimizer (Objective in
	StructureOptimization.cpp) finds it with dynamic_cast when a gradient in the particle coordinates is needed;
	for any other potential it calls Energy() and then AllForce() as before. It is a separate class, not a virtual
	function of Potential, so that Potential's layout and the vtable of every class derived from it stay as they
	were: objects compiled before this interface existed keep working without a rebuild. */
class CombinedEnergyForce
{
public:
	/** \brief Compute the total energy and the force on every particle of the bound configuration. An implementation
		must return bitwise what Energy() followed by AllForce(forces) would return.
		\param[out] forces  Per-particle forces, as AllForce.
		\return The total potential energy, as Energy(). */
	virtual double EnergyAndAllForce(std::vector<GeometryVector> & forces) = 0;

	/** \brief Virtual destructor. */
	virtual ~CombinedEnergyForce() = default;
};

#endif
