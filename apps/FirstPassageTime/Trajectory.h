/** \file Trajectory.h
  * \brief One random-walker's state and its log serialization. */

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <GeometryVector.h>
#include <string>
#include <sstream>
#include <iostream>

/** \struct A class to store simple data for trajectories. */
struct Trajectory{
	DimensionType dimension = 1;
	double sigma_host = 0.0, sigma_particle = 0.0;
	double R2_host = 0.0, R2_particle = 0.0;
	double R2_s = 0.0, tau_s = 0.0;	/* data near the interface*/
	double tau_total = 0.0;

	GeometryVector pos_curr;	/* test particle position in the relative coordinates. */
	GeometryVector delta_pos;   /* a displacement from the initial position in the Cartesian */

	/** \brief Default-construct an empty trajectory with all accumulators set to zero. */
	Trajectory(){
		sigma_host = 0.0; sigma_particle = 0.0;
		R2_host = 0.0;	R2_particle = 0.0;
		R2_s = 0.0; tau_s = 0.0;
	}

	/** \brief Construct a trajectory for the given dimension and phase conductivities.
	 *  \param d       Spatial dimension of the trajectory.
	 *  \param sigma1  Conductivity of the host phase.
	 *  \param sigma2  Conductivity of the particle phase. */
	Trajectory(DimensionType d, double sigma1, double sigma2){
		dimension = d;
		sigma_host = sigma1;
		sigma_particle = sigma2;

		pos_curr = GeometryVector(static_cast<int> (d));
		delta_pos = GeometryVector(static_cast<int> (d));
	}

	/** \brief Copy-construct a trajectory from an existing one.
	 *  \param src  Trajectory to copy. */
	Trajectory(const Trajectory & src);

	/** \brief Defaulted copy-assignment; member-wise copy of all value-type members.
	 *  \param src  Trajectory to copy from.
	 *  \return     Reference to this trajectory. */
	Trajectory & operator=(const Trajectory & src) = default;

	/** \brief Read a logged trajectory record, replacing the positions and accumulating its stored jump statistics into the current trajectory.
	 * \param ifile Input stream positioned at the trajectory record. */
	void ReadLog(std::istream & ifile);

	/** \brief Return the squared net displacement of the test particle.
	 *  \return  Squared modulus of the total displacement from the initial position. */
	inline double GetX2() const{
		return delta_pos.Modulus2();
	}

	/** \brief Accumulate a squared host-phase jump and the associated elapsed time.
	 *  \param R2_h          Squared jump length taken in the host phase.
	 *  \param tau_assigned  Time to add directly when the host phase is insulating (sigma_host == 0). */
	inline void Update_R2_host(double R2_h, double tau_assigned = 0.0){
		
		R2_host += R2_h;
		if (sigma_host == 0.0){
			tau_total += tau_assigned;
		}
		else{
			tau_total += R2_h / (2.0*this->dimension * sigma_host);
		}
	}

	/** \brief Accumulate a squared particle-phase jump and the associated elapsed time.
	 *  \param R2_p          Squared jump length taken in the particle phase.
	 *  \param tau_assigned  Time to add directly when the particle phase is insulating (sigma_particle == 0). */
	inline void Update_R2_particle(double R2_p, double tau_assigned = 0.0){
		R2_particle += R2_p;
		if (sigma_particle == 0.0){
			tau_total += tau_assigned;
		}
		else{
			tau_total += R2_p / (2.0*this->dimension * sigma_particle);
		}
	}

	/** \brief Accumulate a squared near-interface jump and its mean first-passage time.
	 *  \param R2_surface   Squared jump length taken near the interface.
	 *  \param tau_surface  Mean first-passage time elapsed during the near-interface jump. */
	inline void Update_R2_s(double R2_surface, double tau_surface){
		this->R2_s += R2_surface;
		this->tau_s += tau_surface;
		this->tau_total += tau_surface;
	}

	/** \brief Return the total elapsed (first-passage) time accumulated along the trajectory.
	 *  \return  Total accumulated time. */
	inline double GetTimeTotal() const{
		return this->tau_total;
	}
};

/** \brief Write a trajectory's state to an output stream in the logging format.
 *  \param os  Output stream to write to.
 *  \param t   Trajectory to serialize.
 *  \return    Reference to the output stream. */
std::ostream & operator << (std::ostream & os, const Trajectory & t);

#endif
