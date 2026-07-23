/**
*	Author	: Ge Zhang
*	Email	:  - 	
*	Date	: 2015 */

/** \file CCPotential_Kernels.h
 * \brief Catalog of concrete collective-coordinate kernels: tildeV(k) weight mixins
 *  (K2/KM2/KM3/KM4/Lemniscate), target-S0(k) mixins (S0K2/S0Kn/S0Klc/S0KM2/S0KMLog),
 *  and the named potential typedefs (K2Potential, K2K2Potential, ...).
 *  Split out of CollectiveCoordinatePotential.h (2026-07-07). */

#ifndef CCPOTENTIAL_KERNELS_INCLUDED
#define CCPOTENTIAL_KERNELS_INCLUDED

#include "CCPotential_varBox.h"

/** \brief Mixin providing the weight tildeV(k) = (1-k)^2 and its gradient. */
template<typename T>
class K2Template : public T
{
protected:
	/** \brief Eq. (2) in paper-366: weight tildeV(k) = (1-|k|)^2.
		\param k  Wavevector.
		\return The weight tildeV(k). */
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1.0 - kl;
		return temp*temp;
	}

	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of tildeV(k). */
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1.0 - kl;
		return ((-2.0)*temp / kl)*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	K2Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};

//V(k)=(1/k-1)^2
/** \brief Mixin providing the weight tildeV(k) = (1/k-1)^2 and its gradient. */
template<typename T> class KM2Template : public T
{
protected:
	/** \brief Eq. (2) in paper-366: weight tildeV(k) = (1/|k|-1)^2.
		\param k  Wavevector.
		\return The weight tildeV(k). */
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp;
	}

	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of tildeV(k). */
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-2.0)*temp / kl / kl / kl*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	KM2Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};
//V(k)=(1/k-1)^3
/** \brief Mixin providing the weight tildeV(k) = (1/k-1)^3 and its gradient. */
template<typename T> class KM3Template : public T
{
protected:
	/** \brief Eq. (2) in paper-366: weight tildeV(k) = (1/|k|-1)^3.
		\param k  Wavevector.
		\return The weight tildeV(k). */
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp*temp;
	}

	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of tildeV(k). */
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-3.0)*temp*temp / kl / kl / kl*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	KM3Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};
//V(k)=(1/k-1)^4
/** \brief Mixin providing the weight tildeV(k) = (1/k-1)^4 and its gradient. */
template<typename T> class KM4Template : public T
{
protected:
	/** \brief Eq. (2) in paper-366: weight tildeV(k) = (1/|k|-1)^4.
		\param k  Wavevector.
		\return The weight tildeV(k). */
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp*temp*temp;
	}

	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of tildeV(k). */
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-4.0)*temp*temp*temp / kl / kl / kl*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	KM4Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
}; 

/** \brief Mixin providing a lemniscate-shaped (2D) k-space weight tildeV(k) and its gradient. */
template<typename T>
class LemniscateVkTemplate : public T
{
protected:
	double Rc2;
	/** \brief Eq. (2) in paper-366: lemniscate-shaped k-space weight (2D only).
		\param k  Wavevector.
		\return The weight tildeV(k) (zero outside the lemniscate region). */
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double k2 = k.Modulus2();
		double x2my2 = k.x[0] * k.x[0] - k.x[1] * k.x[1];
		double ratio = k2*k2 / Rc2 / x2my2;
		if (ratio < 1 && x2my2>0)
		{
			double temp = ratio - 1;
			//return temp*temp;
			return temp*temp*temp*temp;
		}
		else
			return 0.0;
	}

	/** \brief gradient of tildeVk.
		\param k  Wavevector.
		\return The gradient of the lemniscate weight (zero outside the lemniscate region). */
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double k2 = k.Modulus2();
		double x2my2 = k.x[0] * k.x[0] - k.x[1] * k.x[1];
		double ratio = k2*k2 / Rc2 / x2my2;
		if (ratio < 1 && x2my2>0)
		{
			double temp = ratio - 1;
			double commonpart = 2 * temp / Rc2 / x2my2 / x2my2;
			//return GeometryVector(4 * k.x[0] * k2*x2my2 - 2 * k.x[0] * k2*k2, 4 * k.x[1] * k2*x2my2 + 2 * k.x[1] * k2*k2)*commonpart;
			return GeometryVector(4 * k.x[0] * k2*x2my2 - 2 * k.x[0] * k2*k2, 4 * k.x[1] * k2*x2my2 + 2 * k.x[1] * k2*k2)*(commonpart*2*temp*temp);
		}
		else
			return GeometryVector(0.0, 0.0);
	}
public:
	/** \brief Construct the mixin (2D).
		\param Rc  Lemniscate radius parameter. */
	LemniscateVkTemplate(double Rc) : T(2, Rc), Rc2(Rc*Rc)
	{}
};



typedef K2Template<ShiftedCCPotential_varBox> K2Potential;
typedef KM3Template<ShiftedCCPotential_varBox> KM3Potential;


//S0(k)=k^2
/** \brief Mixin providing the target S0(k) = k^2 and its gradient. */
template<typename T> class S0K2Template : public T
{
protected:
	/** \brief Target structure factor S0(k) = |k|^2.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k)
	{
		return k.Modulus2();
	}

	/** \brief Gradient of the target S0(k).
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return 2.0*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	S0K2Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};

//S0(k)=k^n
/** \brief Mixin providing the target S0(k) = k^n and its gradient. */
template<typename T> class S0KnTemplate : public T
{
protected:
	double n;
	/** \brief Target structure factor S0(k) = |k|^n.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k)
	{
		return std::pow(k.Modulus2(), 0.5*n);
	}

	/** \brief Gradient of the target S0(k).
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return n*std::pow(k.Modulus2(), 0.5*n-1)*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system.
		\param n  Exponent of the power-law target. */
	S0KnTemplate(DimensionType d, double n, double rK = 1.0) : T(d, rK), n(n)
	{}
};

//S0(k)=(k^n+1)/2
/** \brief Mixin providing the target S0(k) = (k^n+1)/2 and its gradient. */
template<typename T> class S0KlcTemplate : public T
{
protected:
	double n;
	/** \brief Target structure factor S0(k) = (|k|^n+1)/2.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k)
	{
		return 0.5*(std::pow(k.Modulus2(), 0.5*n)+1);
	}

	/** \brief Gradient of the target S0(k).
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return 0.5*n*std::pow(k.Modulus2(), 0.5*n - 1)*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system.
		\param n  Exponent of the power-law target. */
	S0KlcTemplate(DimensionType d, double n, double rK = 1.0) : T(d, rK), n(n)
	{}
};

/** \brief Mixin providing the target S0(k) = (1/k-1)^2 and its gradient. */
template<typename T> class S0KM2Template : public T
{
protected:
	/** \brief Target structure factor S0(k) = (1/|k|-1)^2.
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp;
	}

	/** \brief Gradient of the target S0(k).
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-2.0)*temp / kl / kl / kl*k;
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	S0KM2Template(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};


/** \brief Mixin providing the target S0(k) = -log(k) and its gradient. */
template<typename T> class S0KMLogTemplate : public T
{
protected:
	/** \brief Target structure factor S0(k) = -0.5*log(|k|^2).
		\param k  Wavevector.
		\return The target S0(k). */
	virtual double S0k(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl2 = k.Modulus2();
		return -0.5*std::log(kl2);
	}

	/** \brief Gradient of the target S0(k).
		\param k  Wavevector.
		\return The gradient of S0(k). */
	virtual GeometryVector NablaS0k(const GeometryVector & k)//gradient of tildeVk
	{
		return k*(-1.0/k.Modulus2());
	}
public:
	/** \brief Construct the mixin.
		\param d  Spatial dimension of the system. */
	S0KMLogTemplate(DimensionType d, double rK = 1.0) : T(d, rK)
	{}
};



typedef S0K2Template<K2Template<ShiftedCCPotential_S2_varBox> > K2K2Potential;
typedef S0K2Template<K2Template<ShiftedCCPotential_S2_varBox_noV0> > K2K2Potential_noV0;
typedef S0K2Template<KM3Template<ShiftedCCPotential_S2_varBox_noV0> > KM3K2Potential_noV0;


typedef S0K2Template<KM2Template<ShiftedCCPotential_S2_varBox> > KM2K2Potential;
typedef S0K2Template<KM3Template<ShiftedCCPotential_S2_varBox> > KM3K2Potential;
typedef S0K2Template<KM4Template<ShiftedCCPotential_S2_varBox> > KM4K2Potential;

typedef S0KM2Template<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KM2Potential;

typedef S0KMLogTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KMLogPotential;

typedef S0KnTemplate<KM2Template<ShiftedCCPotential_S2_varBox> > KM2KnPotential;
typedef S0KnTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KnPotential;
typedef S0KnTemplate<KM4Template<ShiftedCCPotential_S2_varBox> > KM4KnPotential;

typedef S0KlcTemplate<KM2Template<ShiftedCCPotential_S2_varBox> > KM2KlcPotential;
typedef S0KlcTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KlcPotential;
typedef S0KlcTemplate<KM4Template<ShiftedCCPotential_S2_varBox> > KM4KlcPotential;

#endif
