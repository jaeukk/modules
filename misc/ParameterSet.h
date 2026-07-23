#ifndef PARAMETERSET_INCLUDED
#define PARAMETERSET_INCLUDED

#include <iostream>
#include <Potential.h>
#include <RandomGenerator.h>

class ParameterSet
{
public:
	double * Parameters;
	/** \brief Construct an empty parameter set with no allocated parameters. */
	ParameterSet()
	{
		this->Parameters = nullptr;
	}

	/** \brief Copy-construct a parameter set, deep-copying the parameter array.
	 *  \param src  Parameter set to copy from. */
	ParameterSet(const ParameterSet & src)
	{
		this->Parameters = nullptr;
		if(src.Parameters==nullptr)
			return;
		size_t nbr=src.NumParameters();
		this->Parameters = new double [nbr];
		std::memcpy(this->Parameters, src.Parameters, sizeof(double)*nbr);
	}

	/** \brief Move-construct a parameter set, taking ownership of src's parameter array.
	 *  \param src  Parameter set to move from. */
	ParameterSet(ParameterSet && src)
	{
		this->Parameters = src.Parameters;
		src.Parameters = nullptr;
	}

	/** \brief Copy assignment is deleted; use clone() for polymorphic deep copies. */
	ParameterSet & operator=(const ParameterSet &) = delete;
	/** \brief Move assignment is deleted. */
	ParameterSet & operator=(ParameterSet &&) = delete;

	/** \brief Destroy the parameter set and release the parameter array. */
	virtual ~ParameterSet()
	{
		if(this->Parameters != nullptr)
			delete [] this->Parameters;
	}

	/** \brief Build the pair potential corresponding to the current parameters.
	 *  \return  Newly allocated potential owned by the caller. */
	virtual IsotropicPairPotential * GetPotential(void) const =0;

	/** \brief calculate the penalty from OptimizableParameters
	 *  \return  The penalty value. */
	virtual double Penalty(void) const =0;//calculate the penalty from OptimizableParameters
	/** \brief get range of each parameters
	 *  \param[out] LowerBound  Array receiving the lower bound of each parameter.
	 *  \param[out] UpperBound  Array receiving the upper bound of each parameter. */
	virtual void GetBound(double * LowerBound, double * UpperBound) const =0;//get range of each parameters
	/** \brief evolve to a more complex Parameterization, return true if succed, return false if can't evolve any more
	 *  \return  True on success, false if the parameterization cannot evolve further. */
	virtual bool Evolve(void) =0;//evolve to a more complex Parameterization, return true if succed, return false if can't evolve any more
	/** \brief for some non-continuous parameter, iterate through it. Intended to be used when it has never evolved
	 *  \return  True if a next iterate exists, false otherwise. */
	virtual bool Iterate(void) =0;//for some non-continuous parameter, iterate through it. Intended to be used when it has never evolved
	/** \brief Create a deep copy of this parameter set.
	 *  \return  Newly allocated clone owned by the caller. */
	virtual ParameterSet * clone() const =0;

	/** \brief Serialize the parameter set to an output stream.
	 *  \param out  Destination stream. */
	virtual void Write(std::ostream & out) const =0;

	/** \brief Deserialize the parameter set from an input stream.
	 *  \param in  Source stream. */
	virtual void Read(std::istream & in) =0;
	//virtual double RecommendedTolerence(void) const =0;
	//virtual double RecommendedMaxEnergy(void) const =0;//when optimizing specific volume, the smallest specific volume needs to be so small that this value is achieved
	/** \brief when r<MinDistance, the potential should have no interesting behavior, guaranteed to fix unchanged during the life of this object
	 *  \return  The minimum distance below which the potential is featureless. */
	virtual double MinDistance(void) const =0;//when r<MinDistance, the potential should have no interesting behavior, guaranteed to fix unchanged during the life of this object
	/** \brief the Rcut of the potential should not exceed this, guaranteed to fix unchanged during the life of this object
	 *  \return  The maximum allowed cutoff distance. */
	virtual double MaxDistance(void) const =0;//the Rcut of the potential should not exceed this, guaranteed to fix unchanged during the life of this object

	/** \brief Return the number of parameters in this set.
	 *  \return  Number of parameters. */
	virtual size_t NumParameters(void) const =0;

	/** \brief a signal that global optimization will be carried out */
	virtual void SignalGlobalOptimization(void)//a signal that global optimization will be carried out
	{
	}

	/** \brief vice versa */
	virtual void SignalLocalOptimization(void)//vice versa
	{
	}

	/** \brief when there is a well in the parameterization, it is guaranteed that (it's width)/(it's location) is greater than this number
	 *  \return  Lower bound on the relative size of the sharpest well. */
	virtual double SharpestRelativeSize(void) =0;//when there is a well in the parameterization, it is guaranteed that (it's width)/(it's location) is greater than this number
	/** \brief Randomize each parameter uniformly within its bound.
	 *  \param gen  Random number generator. */
	void Randomize(RandomGenerator & gen)
	{
		size_t n=this->NumParameters();
		double * lb=new double[n];
		double * ub=new double[n];
		this->GetBound(lb, ub);
		for(size_t i=0; i<n; i++)
			this->Parameters[i]=lb[i]+(ub[i]-lb[i])*gen.RandomDouble();

		delete [] ub;
		delete [] lb;
	}
};
	
class CombinedParameterization : public ParameterSet
{
private:
	ParameterSet * p1;
	ParameterSet * p2;
	double * pParameters1;
	double * pParameters2;
public:
	/** \brief Combine two parameter sets into one, taking ownership of both.
	 *  \param p1  First parameterization; this object takes ownership.
	 *  \param p2  Second parameterization; this object takes ownership. */
	CombinedParameterization(ParameterSet * p1, ParameterSet * p2) : ParameterSet()
	{
		this->p1=p1;
		this->p2=p2;
		this->Allocate();
	}

	/** \brief Copy construction is deleted; owning p1/p2 make shallow copies unsafe. Use clone(). */
	CombinedParameterization(const CombinedParameterization &) = delete;
	/** \brief Copy assignment is deleted; owning p1/p2 make shallow copies unsafe. Use clone(). */
	CombinedParameterization & operator=(const CombinedParameterization &) = delete;

	/** \brief Destroy the combined parameterization and its two sub-parameterizations. */
	~CombinedParameterization()
	{
		this->DeAllocate();
		delete this->p1;
		delete this->p2;
	}

	/** \brief Build the combined pair potential from both sub-parameterizations.
	 *  \return  Newly allocated combined potential owned by the caller. */
	virtual IsotropicPairPotential * GetPotential(void) const
	{
		IsotropicPairPotential * result = new CombinedIsotropicPairPotential(this->p1->GetPotential(), this->p2->GetPotential());
		return result;
	}

	/** \brief Return the summed penalty of both sub-parameterizations.
	 *  \return  Sum of the two penalties. */
	virtual double Penalty(void) const
	{
		double result = this->p1->Penalty()+this->p2->Penalty();
		return result;
	}

	/** \brief Fill the combined parameter bounds from both sub-parameterizations.
	 *  \param[out] LowerBound  Array receiving the lower bound of each parameter.
	 *  \param[out] UpperBound  Array receiving the upper bound of each parameter. */
	virtual void GetBound(double * LowerBound, double * UpperBound) const
	{
		size_t Num1=p1->NumParameters();
		this->p1->GetBound(LowerBound, UpperBound);
		this->p2->GetBound(LowerBound+Num1, UpperBound+Num1);
	}

	/** \brief Evolve the second sub-parameterization to a more complex form.
	 *  \return  True on success, false if it cannot evolve further. */
	virtual bool Evolve(void)
	{
		this->DeAllocate();
		double result=this->p2->Evolve();
		this->Allocate();
		return result;
	}

	/** \brief Iterate the second sub-parameterization through its discrete states.
	 *  \return  True if a next iterate exists, false otherwise. */
	virtual bool Iterate(void)
	{
		this->DeAllocate();
		double result=this->p2->Iterate();
		this->Allocate();
		return result;
	}

	/** \brief Serialize both sub-parameterizations to an output stream.
	 *  \param out  Destination stream. */
	virtual void Write(std::ostream & out) const
	{
		this->p1->Write(out);
		this->p2->Write(out);
	}

	/** \brief Deserialize both sub-parameterizations from an input stream.
	 *  \param in  Source stream. */
	virtual void Read(std::istream & in)
	{
		this->DeAllocate();
		this->p1->Read(in);
		this->p2->Read(in);
		this->Allocate();
	}

	/** \brief Return the smaller MinDistance of the two sub-parameterizations.
	 *  \return  Minimum of the two MinDistance values. */
	virtual double MinDistance(void) const
	{
		double m1=this->p1->MinDistance();
		double m2=this->p2->MinDistance();

		return m1<m2?m1:m2;
	}

	/** \brief Return the larger MaxDistance of the two sub-parameterizations.
	 *  \return  Maximum of the two MaxDistance values. */
	virtual double MaxDistance(void) const
	{
		double m1=this->p1->MaxDistance();
		double m2=this->p2->MaxDistance();

		return m1>m2?m1:m2;
	}

	/** \brief Return the total parameter count of both sub-parameterizations.
	 *  \return  Sum of the two parameter counts. */
	virtual size_t NumParameters(void) const
	{
		return this->p1->NumParameters()+this->p2->NumParameters();
	}

	/** \brief Return the smaller SharpestRelativeSize of the two sub-parameterizations.
	 *  \return  Minimum of the two SharpestRelativeSize values. */
	virtual double SharpestRelativeSize(void)
	{
		double s1=this->p1->SharpestRelativeSize();
		double s2=this->p2->SharpestRelativeSize();

		return s1<s2?s1:s2;
	}

	/** \brief Deep-copy this combined parameterization, cloning both sub-parameterizations.
	 *  \return  Newly allocated clone owned by the caller. */
	virtual ParameterSet * clone() const
	{
		return new CombinedParameterization(this->p1->clone(), this->p2->clone());
	}

	/** \brief Allocate the combined parameter array and alias each sub-parameterization into it. */
	void Allocate()
	{
		size_t Num1=p1->NumParameters();
		size_t Num2=p2->NumParameters();
		this->Parameters = new double [Num1+Num2];
		this->pParameters1 = this->Parameters;
		this->pParameters2 = this->Parameters+Num1;
		std::memcpy(pParameters1, p1->Parameters, Num1*sizeof(double));
		std::memcpy(pParameters2, p2->Parameters, Num2*sizeof(double));
		std::swap(this->pParameters1, this->p1->Parameters);
		std::swap(this->pParameters2, this->p2->Parameters);
	}

	/** \brief Flush aliased parameters back to each sub-parameterization and free the combined array. */
	void DeAllocate()
	{
		size_t Num1=p1->NumParameters();
		size_t Num2=p2->NumParameters();
		std::memcpy(pParameters1, p1->Parameters, Num1*sizeof(double));
		std::memcpy(pParameters2, p2->Parameters, Num2*sizeof(double));
		std::swap(this->pParameters1, this->p1->Parameters);
		std::swap(this->pParameters2, this->p2->Parameters);
		delete [] this->Parameters;
		this->Parameters=nullptr;
		this->pParameters1=nullptr;
		this->pParameters2=nullptr;
	}
};

#endif