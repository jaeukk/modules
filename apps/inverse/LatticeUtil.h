#ifndef INVERSE_LATTICEUTIL_H
#define INVERSE_LATTICEUTIL_H

/** \file LatticeUtil.h
 *  \brief The only two helpers the CC-inversion app needed from the legacy
 *  grab-bag header etc1.h, relocated here so etc1.{h,cpp} can be retired.
 *
 *  - GetLattice_config : build a Configuration of lattice points
 *                        x = sum_i n_i basis[i], n_i = 0 .. Nmax_i - 1.
 *  - Find_LeastCommonDivider : smallest prime factor of n (1 if n is prime).
 *
 *  Both depend only on the `core` module. They live in the app because (per the
 *  modularization rules) they are used by only this one consumer; if a future
 *  module needs them more broadly, promote them into core/etc.{h,cpp}.
 */

#include <cmath>
#include <GeometryVector.h>
#include <PeriodicCellList.h>

/** \brief Generate lattice points x = n_0*basis[0]+...+n_{d-1}*basis[d-1], n_i = 0..Nmax_i-1.
 *  \param basis  Lattice basis vectors; their count sets the spatial dimension.
 *  \param Nmax   Per-dimension counts: component i gives the number of lattice steps along basis[i].
 *  \return       A Configuration containing the generated lattice points. */
inline Configuration GetLattice_config(const std::vector<GeometryVector> &basis, const GeometryVector & Nmax) {
	int d = basis.size();
	GeometryVector * Basis;
	Basis = new GeometryVector[d];
	for (int i = 0; i < d; i++) {
		Basis[i] = Nmax.x[i]*basis[i];
	}
	Configuration result(d, Basis, Nmax.x[0]/100.0);
	auto N_narySystem = [&Nmax, d](int i)->GeometryVector {
		GeometryVector temp(static_cast<int> (d));
		for (int j = 0; j < d; j++) {
			temp.x[j] = (double)(i % (int)Nmax.x[j]) / Nmax.x[j];
			i = (int)floor(i / Nmax.x[j]);
		}
		return temp;
	};

	int NumPrt = 1;
	for (int i = 0; i < d; i++) { NumPrt *= Nmax.x[i]; }
	GeometryVector aa(static_cast<int>(d));
	for (int i = 0; i < NumPrt; i++) {
		aa = N_narySystem(i);
		result.Insert("a", aa);
	}

	delete[] Basis;
	return result;
}

/** \brief Same, with an isotropic Nmax in every dimension.
 *  \param basis  Lattice basis vectors; their count sets the spatial dimension.
 *  \param Nmax   Number of lattice steps along every basis direction.
 *  \return       A Configuration containing the generated lattice points. */
inline Configuration GetLattice_config(const std::vector<GeometryVector> &basis, int Nmax) {
	DimensionType d = basis.size();
	GeometryVector Nmax_(static_cast<int> (d));
	for (int i = 0; i < d; i++) {
		Nmax_.x[i] = Nmax;
	}
	return GetLattice_config(basis, Nmax_);
}

/** \brief Smallest prime factor of n; returns 1 if n is prime.
 *  \param n  Integer to factor.
 *  \return   The smallest prime factor of n, or 1 if n is prime. */
inline int Find_LeastCommonDivider(int n) {
	int max = (int)std::floor(std::sqrt((double)n));
	for (int i = 2; i <= max; i++) {
		if (n % i == 0)
			return i;
	}
	return 1;
}

#endif
