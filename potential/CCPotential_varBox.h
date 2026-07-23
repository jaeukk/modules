/**
*	Author	: Ge Zhang
*	Email	:  - 	
*	Date	: 2015 */

/** \file CCPotential_varBox.h
 * \brief Variable-cell (box-deformation) support for collective-coordinate potentials:
 *  the ShiftedCCPotential_varBoxable mixin and the varBox family (stealthy, Eq. (6) of
 *  paper-366, plus the S2 variants). Split out of CollectiveCoordinatePotential.h
 *  (2026-07-07); the only CC header needing GSL linear algebra (LU solves/inverts
 *  for the basis-vector derivatives dBi/dki and EnergyDerivativeToBasisVectors). */

#ifndef CCPOTENTIAL_VARBOX_INCLUDED
#define CCPOTENTIAL_VARBOX_INCLUDED

#include <gsl/gsl_linalg.h>

#include "CollectiveCoordinatePotential.h"


//this class allows box deformation for ShiftedCCPotential or its derivative classes, provides:
//k point list management
//user-specified tildeV(\mathbf k)
//d \mathbf b_i/da_{lm} and d \mathbf k_i/da_{lm}, where b_i is the ith reciprocal basis vector, k_i is the ith k point, 
// and a_{lm} is BasisVector[l].x[m]
//
//does NOT provide:
//management of ShiftedCCPotential::Shift
//implementation of Potential::EnergyDerivativeToBasisVectors 
/** \brief Mixin adding variable-cell (box-deformation) support to a collective-coordinate potential. */
template<typename Base>
class ShiftedCCPotential_varBoxable : public Base
{
protected:
	double rK;
	GeometryVector CachedBasisVector[::MaxDimension];
	GeometryVector CachedReciprocalBasisVector[::MaxDimension];
	bool SameBasisVector;
	//a k point is made of c[i]*b[i]
	/** \brief Integer coefficients expressing a k-point as a combination of reciprocal basis vectors. */
	struct kCoeff
	{
		double c[::MaxDimension];
		/** \brief Construct a zero coefficient set. */
		kCoeff()
		{
			for (int i = 0; i < ::MaxDimension; i++)
				c[i] = 0;
		}
	};
	std::vector<kCoeff> kCoeffs;
	/** \brief Weight tildeV(k) assigned to a wavevector; Eq. (2) in paper-366.
		\param k  Wavevector.
		\return The weight tildeV(k). */
	virtual double tildeVk(const GeometryVector & k) = 0;//Eq. (2) in paper-366
	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of tildeV(k) with respect to k. */
	virtual GeometryVector NablatildeVk(const GeometryVector & k) = 0;//gradient of tildeVk
	/** \brief d \mathbf b_i / da_{ lm }: derivative of reciprocal basis vectors w.r.t. a basis-vector component.
		\param l  Index of the basis vector being varied.
		\param m  Coordinate component of that basis vector being varied.
		\return Per-reciprocal-vector derivatives. */
	std::vector<GeometryVector> dBi(DimensionType l, DimensionType m)//d \mathbf b_i / da_{ lm }
	{
		std::vector<GeometryVector> result;
		//The system matrix (the basis-vector matrix) is identical for every reciprocal vector, so
		//allocate and LU-decompose it once and reuse the factorization for each right-hand side,
		//instead of the former per-vector alloc+decompose. (Loop indices renamed to avoid shadowing.)
		gsl_matrix * a = gsl_matrix_alloc(this->Dimension, this->Dimension);
		for (DimensionType row = 0; row < this->Dimension; row++)
			for (DimensionType col = 0; col < this->Dimension; col++)
				gsl_matrix_set(a, row, col, CachedBasisVector[row].x[col]);
		gsl_permutation * p = gsl_permutation_alloc(this->Dimension);
		int s;
		gsl_linalg_LU_decomp(a, p, &s);

		gsl_vector * x = gsl_vector_alloc(this->Dimension);
		gsl_vector * b = gsl_vector_alloc(this->Dimension);
		for (DimensionType i = 0; i < this->Dimension; i++)
		{
			gsl_vector_set_zero(b);
			gsl_vector_set(b, l, (-1.0)*CachedReciprocalBasisVector[i].x[m]);
			gsl_linalg_LU_solve(a, p, b, x);
			GeometryVector res(this->Dimension);
			for (DimensionType row = 0; row < this->Dimension; row++)
				res.x[row] = gsl_vector_get(x, row);
			result.push_back(res);
		}

		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_permutation_free(p);
		gsl_matrix_free(a);
		return result;
	}

	/** \brief d \mathbf k_i / da_{ lm }: derivative of every constrained k-point w.r.t. a basis-vector component.
		\param l  Index of the basis vector being varied.
		\param m  Coordinate component of that basis vector being varied.
		\return Per-constraint k-point derivatives. */
	std::vector<GeometryVector> dki(DimensionType l, DimensionType m)//d \mathbf k_i / da_{ lm }
	{
		std::vector<GeometryVector> result;
		std::vector<GeometryVector> db = dBi(l, m);
		for (size_t i = 0; i < this->constraints.size(); i++)
		{
			GeometryVector res(this->Dimension);
			for (DimensionType j = 0; j < this->Dimension; j++)
				res.AddFrom(kCoeffs[i].c[j] * db[j]);
			result.push_back(res);
		}
		return result;
	}
public:
	/** \brief Construct a variable-cell collective-coordinate potential.
		\param d   Spatial dimension of the system.
		\param rK  Cutoff radius in reciprocal space for selecting constrained k-points. */
	ShiftedCCPotential_varBoxable(DimensionType d, double rK) : Base(d), rK(rK)
	{
		this->SameBasisVector = false;
		this->Base::Shift = 0.0;
		for (DimensionType i = 0; i < MaxDimension; i++)
			CachedBasisVector[i] = GeometryVector(d);
	}

	/** \brief Bind the configuration and regenerate the constrained k-point list when the cell has changed.
		\param config  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & config)
	{
		//this function does not set Base::Shift
		DimensionType d = this->Potential::Dimension;
		//test if basis vector has never changed
		//if so, we don't need to re-generate the constrained k points
		for (DimensionType i = 0; i < d; i++)
		{
			if (this->CachedBasisVector[i] != config.GetBasisVector(i))
			{
				SameBasisVector = false;
				break;
			}
		}
		if (SameBasisVector)
		{
			this->Base::SetConfiguration(config);
			return;
		}

		GeometryVector b[::MaxDimension];
		for (DimensionType i = 0; i<d; i++)
			b[i] = config.GetReciprocalBasisVector(i);
		PeriodicCellList<Empty> RecipLattice(d, b, 1e100, true);
		RecipLattice.Insert(Empty(), SameDimensionZeroVector(b[0]));
		std::vector<GeometryVector> KPoints;
		std::vector<kCoeff> & temp = this->kCoeffs;
		temp.clear();
		auto IterateFunction = [&KPoints, &d, &temp](const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * PeriodicShift, const size_t /*Sourceparticle*/) ->void
		{
			bool NonTrivial = false;//is not vector zero and is not symmetrical to another vector
			for (DimensionType i = 0; i<d; i++)
			{
				if (PeriodicShift[i]>0)
				{
					NonTrivial = true;
					break;
				}
				else if (PeriodicShift[i]<0)
				{
					NonTrivial = false;
					break;
				}
			}
			if (NonTrivial)
			{
				KPoints.push_back(shift);
				kCoeff a;
				for (DimensionType i = 0; i < d; i++)
					a.c[i] = PeriodicShift[i];
				temp.push_back(a);
			}
		};
		//Must not refine basis vectors. Doing so will mess up correspondence between constraints and kCoeff
		//RecipLattice.TryRefineBasisVectors();
		RecipLattice.IterateThroughNeighbors(SameDimensionZeroVector(b[0]), rK, IterateFunction);

		this->Base::constraints.clear();
		for (auto iter = KPoints.begin(); iter != KPoints.end(); iter++)
		{
			double Vk = this->tildeVk(*iter);
			this->Base::constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, Vk));
		}
		this->Base::SetConfiguration(config);

		for (DimensionType i = 0; i < d; i++)
		{
			this->CachedBasisVector[i] = config.GetBasisVector(i);
			this->CachedReciprocalBasisVector[i] = config.GetReciprocalBasisVector(i);
		}
		//the cached basis now matches config, so a subsequent bind with the same cell can skip
		//regeneration (this flag was previously left false, defeating the k-point cache)
		this->SameBasisVector = true;
		//carried from old code, however, here one cannot sort them without changing order of kCoeffs
		//auto CompareFunction = [](const GeometryVector & left, const GeometryVector & right) ->bool
		//{
		//	return left.Modulus2()<right.Modulus2();
		//};
		//std::sort(KPoints.begin(), KPoints.end(), CompareFunction);
	}
};

//Potential in Eq. (6) of paper-366, need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);//Eq. (2) in paper-366
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
/** \brief Variable-cell stealthy collective-coordinate potential (Eq. (6) of paper-366). */
class ShiftedCCPotential_varBox : public ShiftedCCPotential_varBoxable<ShiftedCCPotential>
{
protected:
	size_t nbr;
public:
	/** \brief Construct a variable-cell stealthy potential.
		\param d   Spatial dimension of the system.
		\param rK  Cutoff radius in reciprocal space for selecting constrained k-points. */
	ShiftedCCPotential_varBox(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential>(d, rK), nbr(0)
	{}
	/** \brief Bind the configuration and recompute the energy shift when the cell or particle count changes.
		\param config  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector || nbr != config.NumParticle())
		{
			//we need to re-calculate this->Shift of basis vector of nbr has changed
			//calculate SumVk
			double SumVk = 0.0;
			for (auto iter = this->constraints.begin(); iter != this->constraints.end(); ++iter)
				SumVk += iter->V;

			nbr = config.NumParticle();
			double v = config.PeriodicVolume();

			this->Shift = nbr*(nbr - 1) / v / 2.0 - nbr / v*SumVk;
		}
	}

	/** \brief Energy derivative with respect to the cell basis vectors.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative w.r.t. component m of basis vector n.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		double vF = this->Volume;
		for (DimensionType n = 0; n<Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*(RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i] - nbr);
				temp /= vF;

				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure - phi / vF)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};

//// \phi = \sum_{\mathbf k} tildeV(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2, need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
//virtual double S0k(const GeometryVector & k);
//virtual GeometryVector NablaS0k(const GeometryVector & k);//gradient of S0k
/** \brief Variable-cell S2 collective-coordinate potential with a user-supplied target S0(k). */
class ShiftedCCPotential_S2_varBox : public ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>
{
protected:
	/** \brief Target structure factor S0(k) at a wavevector.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k) = 0;

	/** \brief gradient of S0k.
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k) = 0;//gradient of S0k
public:
	/** \brief Construct a variable-cell S2 potential.
		\param d   Spatial dimension of the system.
		\param rK  Cutoff radius in reciprocal space for selecting constrained k-points. */
	ShiftedCCPotential_S2_varBox(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>(d, rK)
	{}
	/** \brief Bind the configuration and recompute the target S0(k) when the cell changes.
		\param config  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector)
		{
			//we need to re-calculate S0 of basis vector has changed
			this->S0.resize(this->constraints.size());
			for (size_t i = 0; i < this->constraints.size(); i++)
				this->S0[i] = S0k(this->constraints[i].k);
		}
	}

	/** \brief Energy derivative with respect to the cell basis vectors.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative w.r.t. component m of basis vector n.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		for (DimensionType n = 0; n < Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
				{
					double ds = ((RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i]) / N - S0[i]);
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*ds*ds;
					temp += this->constraints[i].V*(-2.0)*ds*NablaS0k(this->constraints[i].k).Dot(dk[i]);
				}
				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		//double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};

//// \phi = \sum_{\mathbf k} tildeV(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2-V0, 
//where V0 is given by Eq. (28) of paper-309.
//need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
//virtual double S0k(const GeometryVector & k);
//virtual GeometryVector NablaS0k(const GeometryVector & k);//gradient of S0k
/** \brief Variable-cell S2 collective-coordinate potential subtracting the V0 baseline of Eq. (28) of paper-309. */
class ShiftedCCPotential_S2_varBox_noV0 : public ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>
{
protected:
	/** \brief Target structure factor S0(k) at a wavevector.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k) = 0;

	/** \brief gradient of S0k.
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k) = 0;//gradient of S0k
public:
	/** \brief Construct a variable-cell S2 (no-V0) potential.
		\param d   Spatial dimension of the system.
		\param rK  Cutoff radius in reciprocal space for selecting constrained k-points. */
	ShiftedCCPotential_S2_varBox_noV0(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>(d, rK)
	{}
	/** \brief Bind the configuration, recomputing S0(k) and the V0 shift when the cell changes.
		\param config  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector)
		{
			//we need to re-calculate S0 of basis vector has changed
			this->S0.resize(this->constraints.size());
			this->Shift = 0.0;
			for (size_t i = 0; i < this->constraints.size(); i++)
			{
				double s0 = S0k(this->constraints[i].k);
				this->S0[i] = s0;
				this->Shift -= (s0 - 1)*(s0 - 1)*this->constraints[i].V;
			}
		}
	}

	/** \brief Energy derivative with respect to the cell basis vectors.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative w.r.t. component m of basis vector n.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		for (DimensionType n = 0; n < Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
				{
					double ds = ((RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i]) / N - S0[i]);
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*(ds*ds-(S0[i]-1.0)*(S0[i]-1.0)) + this->constraints[i].V*(1.0-S0[i] - ds)*2.0*NablaS0k(this->constraints[i].k).Dot(dk[i]);
				}
				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		//double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};

#endif
