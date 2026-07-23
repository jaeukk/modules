/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com	
*	Date	: Feb. 2019 */

/** \file RepulsiveCCPotential.h
 * \brief Header file to declare a class to compute "stealthy" collective coordinate potential + Lenard-Jones potential. */

#ifndef RepulsiveCCPotential_H
#define RepulsiveCCPotential_H

#include <memory>

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
/* in $(pot) */
#include <Potential.h>
#include <CollectiveCoordinatePotential.h>

class RepulsiveCCPotential : public Potential 
{
protected:
	const Configuration * pConfig;	//< A configuration of interest
	double rc;						//< The critical radius over which the potential becomes 0.
	const double e_critical = 4.0*(pow(2, -12/6) - pow(2, -6/6));//< LJ potential at the critical radius.
	bool WithRepulsion = true;
public:
	double e0,	//< energy scale for LJ
		sigma;	//< length scale for LJ
	// Owned collective-coordinate kernel (was a raw pointer that leaked; the
	// unique_ptr releases it in the destructor).
	std::unique_ptr<ShiftedCCPotential> CCPotential;

	/** \brief Construct the soft-core and collective-coordinate potential, selecting the collective-coordinate kernel from S0.
	 *  \param dim Spatial dimension.
	 *  \param epsilon Energy scale of the quadratic soft-core repulsion; repulsion is disabled when its magnitude is near zero.
	 *  \param sigma Soft-core contact distance.
	 *  \param S0 Kernel selector: zero creates ShiftedCCPotential and nonzero creates ShiftedCCPotential_S2; the value itself is not stored as a target, which must be supplied when constraints are added. */
	RepulsiveCCPotential(DimensionType dim, double epsilon, double sigma, double S0 = 0.0) : Potential (dim) {
		if (S0 == 0.0){
			CCPotential = std::make_unique<ShiftedCCPotential>(dim);
		}
		else{
			CCPotential = std::make_unique<ShiftedCCPotential_S2>(dim);
		}
		this->e0 = epsilon;	this->sigma = sigma;	this->rc = sigma;/*this->rc = pow(2, 1.0/6.0)*sigma;*/
		this->WithRepulsion = (std::abs(this->e0) > 1e-10);//(this->e0 > 1e-5) || (this->e0 < -1e-5);
		if(this->WithRepulsion){
			std::cout<< "Soft-core repulsion is on\n";
		}
		else{
			std::cout<< "Soft-core repulsion is off\n";
		}
	}

	/** \brief Destructor. */
	virtual ~RepulsiveCCPotential() {}

	/** \brief Compute the truncated LJ potential and collective-coordiante potentials for a given configuraiton. */
	virtual double Energy() {
		double LJ = 0.0L, CC = 0.0L;
		if (this->WithRepulsion){
			//this->pConfig->PrepareIterateThroughNeighbors(this->rc);
	#pragma omp parallel for schedule (guided, 1) num_threads(this->ParallelNumThread) reduction(+:LJ)
			for (size_t i = 0; i < pConfig->NumParticle(); i++) {
				double E = 0.0L;
					pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->sigma, 
						[&i, &E, this](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t SourceParticle)->void {
						/*Prevent double-counting */
						if (i < SourceParticle) { 
							double r = sqrt(shift.Modulus2()) / this->sigma;
							if (r < 1.0){
								E += this->e0*(1.0 - r)*(1.0 - r);
							}
						}
						}
					);
//	#pragma omp atomic
				LJ += E;
			}
		}
		CC = this->CCPotential->Energy();
		if(Verbosity > 3){
			double E = LJ + CC;
			std::cout << E << "\t(soft-core = " << 100.*LJ/E <<" %, \t CC = " << 100.*CC/E <<" %)"<<std::endl;
			// std::cout << "LJ = " << LJ <<std::endl;
			// std::cout << "CC = " << CC <<std::endl;
		}
		return LJ + CC;
	}

	/** \brief Accumulate only the soft-core (contact) force on particle i.
	 *  \param[in,out] force_LJ  Contact-force accumulator; the contact contribution is added to it.
	 *  \param[in] i             Particle index. */
	void ForceLJ(GeometryVector & force_LJ, size_t i) {
		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->sigma,
			[&i, &force_LJ, this](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t SourceParticle)->void {
			/*Prevent self-interaction */
			if (i != SourceParticle) {
				double r = sqrt(shift.Modulus2());
				if (r < 1e-15){
					/* two particles are overlapped.*/
					GeometryVector x(shift);
					x.x[0] = (i > SourceParticle)? -1.0 : 1.0 ;
					force_LJ = force_LJ - 2.0*this->e0 / this->sigma * x;
				}
				else{
					force_LJ = force_LJ - 2.0*this->e0 / this->sigma * (1.0 / r - 1.0 / this->sigma) * shift;
				}
				//force_LJ = force_LJ + 2.0*this->e0 / this->sigma * (1.0 / sqrt(shift.Modulus2()) - 1.0 / this->sigma) * shift;
			}
		});
	}

	/** \brief Compute the force exerted at particle i
	 *  \param[out] force	Net force in a GeometryVector.
	 *  \param[in] i			Particle index. */
	virtual void Force(GeometryVector & force, size_t i) {
		GeometryVector force_LJ(static_cast<int> (this->Dimension));

		if (this->WithRepulsion){
			this->ForceLJ(force_LJ, i);
		}
		GeometryVector force_cc;
		this->CCPotential->Force(force_cc, i);
		force = force_LJ + force_cc;
	}

	/** \brief Compute the force exerted at all particles.
	 *  \param[out] result	result[i] is the force exerted at particle i. */
	virtual void AllForce(std::vector<GeometryVector> & result) {
		this->CCPotential->AllForce(result);
			if (this->WithRepulsion){
			//this->pConfig->PrepareIterateThroughNeighbors(this->rc);
	#pragma omp parallel for schedule (guided) num_threads(this->ParallelNumThread)
			for (size_t i = 0; i < pConfig->NumParticle(); i++) {
				//result already holds the collective-coordinate force from CCPotential->AllForce;
				//add ONLY the soft-core contact force here (calling Force would double-count CC).
				GeometryVector force_LJ(static_cast<int>(this->Dimension));
				this->ForceLJ(force_LJ, i);
				result[i] = result[i] + force_LJ;
			}
		}
	}

	/** \brief Define a configuration to compute energy and forces.
	 *  \param c	A configuration of interest. */
	virtual void SetConfiguration(const Configuration & c) {
		this->pConfig = &c;
		if(this->WithRepulsion){
			this->pConfig->PrepareIterateThroughNeighbors(this->sigma);
		}
		this->CCPotential->SetConfiguration(c);
	}

};


#endif