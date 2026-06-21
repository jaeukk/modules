/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: April 2021 */


/** \file FirstPassageTime.h
  * \brief Header file for a FirstPassageTime simulation
  * for two-phase packings.
  * */

#ifndef FIRSTPASSAGETIME_H
#define FIRSTPASSAGETIME_H

#include <Dispersions.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>
#include <RandomGenerator.h>
#include <etc.h>

#include <omp.h>

//TODO: Superconducting phase
//TODO: Insulting phase

/** \struct A class to store simple data for trajectories. */
struct Trajectory{
	DimensionType dimension = 1;
	double sigma_host = 0.0, sigma_particle = 0.0;
	double R2_host = 0.0, R2_particle = 0.0;
	double R2_s = 0.0, tau_s = 0.0;	/* data near the interface*/
	double tau_total = 0.0;

	GeometryVector pos_curr;	/* test particle position in the relative coordinates. */
	GeometryVector delta_pos;   /* a displacement from the initial position in the Cartesian */

	/** \brief Default-construct an empty trajectory with all accumulators set to zero. (CLAUDE) */
	Trajectory(){
		sigma_host = 0.0; sigma_particle = 0.0;
		R2_host = 0.0;	R2_particle = 0.0;
		R2_s = 0.0; tau_s = 0.0;
	}

	/** \brief Construct a trajectory for the given dimension and phase conductivities. (CLAUDE)
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

	/** \brief Copy-construct a trajectory from an existing one. (CLAUDE)
	 *  \param src  Trajectory to copy. */
	Trajectory(const Trajectory & src);

	/** \brief Read a trajectory's state from an input stream produced by the logging format. (CLAUDE)
	 *  \param ifile  Input stream positioned at the trajectory record. */
	void ReadLog(std::istream & ifile);

	/** \brief Return the squared net displacement of the test particle. (CLAUDE)
	 *  \return  Squared modulus of the total displacement from the initial position. */
	inline double GetX2() const{
		return delta_pos.Modulus2();
	}

	/** \brief Accumulate a squared host-phase jump and the associated elapsed time. (CLAUDE)
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

	/** \brief Accumulate a squared particle-phase jump and the associated elapsed time. (CLAUDE)
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

	/** \brief Accumulate a squared near-interface jump and its mean first-passage time. (CLAUDE)
	 *  \param R2_surface   Squared jump length taken near the interface.
	 *  \param tau_surface  Mean first-passage time elapsed during the near-interface jump. */
	inline void Update_R2_s(double R2_surface, double tau_surface){
		this->R2_s += R2_surface;
		this->tau_s += tau_surface;
		this->tau_total += tau_surface;
	}

	/** \brief Return the total elapsed (first-passage) time accumulated along the trajectory. (CLAUDE)
	 *  \return  Total accumulated time. */
	inline double GetTimeTotal() const{
		return this->tau_total;
	}
};

/** \brief Write a trajectory's state to an output stream in the logging format. (CLAUDE)
 *  \param os  Output stream to write to.
 *  \param t   Trajectory to serialize.
 *  \return    Reference to the output stream. */
std::ostream & operator << (std::ostream & os, const Trajectory & t);

struct FirstPassageTimeParameters{
	double sigma_host, sigma_particle;
	double alpha; 
	
	size_t dimension;
	size_t num_test_prts;
	size_t num_threads;

	/** \brief Construct the simulation parameters and derive the conductivity ratio alpha. (CLAUDE)
	 *  \param sigma_host      Conductivity of the host phase.
	 *  \param sigma_particle  Conductivity of the particle phase.
	 *  \param Dimension       Spatial dimension of the simulation. */
	FirstPassageTimeParameters(double sigma_host, double sigma_particle, size_t Dimension){
		this->sigma_host = sigma_host; this->sigma_particle = sigma_particle;
		this->dimension = Dimension;
		
		if(sigma_host*sigma_particle > 0.0){
			this->alpha = sigma_particle / sigma_host;}
		else{
			this->alpha = 0.0;
		}
	}

	/** \brief Reset the trajectory if its parameters differ from the current object; otherwise keep it. (CLAUDE)
	 *  \param traj  Trajectory to validate and possibly re-initialize.
	 *  \return      True if the trajectory matched and can be continued, false if it was reset. */
	bool ContinueTrajectory(Trajectory & traj);

	/** \brief Write the simulation parameters and a set of trajectories to a ".trajlog" file. (CLAUDE)
	 *  \param name   Base name of the log file (".trajlog" is appended).
	 *  \param trajs  Trajectories to serialize. */
	void WriteLogs(const std::string & name, const std::vector<Trajectory> & trajs);
	/** \brief Read trajectories from a ".trajlog" file, checking that the header matches the current object. (CLAUDE)
	 *  \param name   Base name of the log file (".trajlog" is appended).
	 *  \param trajs  Trajectory vector to populate from the file.
	 *  \return       True if the file existed and its parameters were consistent, false otherwise. */
	bool ReadLogs(const std::string & name, std::vector<Trajectory> & trajs);

/** \brief The probability that a test particle goes to the host phase.
 * The probability that a test particle in phase 1 (\theta = \pi/2) or phase 2 (\theta = 3\pi/2) hits \partial \Omega_1 in 2D for the first time.
 * It corresponds to p_1(r, \theta = \pi/2) or p_1(r, \theta = 3\pi/2), respectively, in the paper (phase 1 = host phase).
 * \param x = r/R	r = The distance from the flat interface, R = The radius of an imaginary sphere.
 * \param pt_in_phase1 = True if the particle is in the particle phase.
 * \return Probability that a test particle first hit \partial\Omega_1 without hitting \partial\Omega_2.*/
	inline double p1_2D(double x, bool pt_in_phase1){
		if(this->alpha > 0.0){
			/* General case */
			if(x < 1e-10){
				return 1.0/(1.0+alpha);
			}
			else{
				if (pt_in_phase1){
					return 1.0/(1.0+alpha) * (1.0 + 4.0*alpha / pi * atan(x));
				}
				else{
					return 1.0/(1.0+alpha) * (1.0 - 4.0 / pi * atan(x));
				}
			}
		}
		else{
			/* if one phase is insulating, 
			 zero probability to cross the interface. */
			if(pt_in_phase1)
				return 1.0;
			else
				return 0.0;
		}
	}

	/** \brief The mean time that a test particle in phase 1 (\theta = \pi/2) or phase 2 (\theta = 3\pi/2) hits \partial \Omega in 2D for the first time.
	 * \param r	The distance from the flat interface,
	 * \param R The radius of an imaginary sphere*
	 * \param pt_in_phase1 = True if the particle is in the particle phase.
	 * \return The mean time that a test particle first hit the boundary of the imaginary disk.*/
	inline double tau_s_2D(double r, double R, bool pt_in_phase1){
		double x = r/R;
		if(this->alpha > 0.0){
			double tau_1 = R*R/(4.0*sigma_host);
			if(x < 1e-10){
				return tau_1 * 2.0 / (1.0 + alpha);
			}
			else{
				double rR_Rr = x + 1.0/x;
				if (pt_in_phase1){
					double corr = (1.0 - alpha * x*x +(alpha-1.0)/pi * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );
					return tau_1* 2.0/(1.0+alpha)* corr;
				}
				else{
					double corr = (1.0 - x*x/alpha -(alpha-1.0)/pi/alpha * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );
					return tau_1* 2.0/(1.0+alpha)* corr;
				}
			}
		}
		else{
			/* one phase is insulating. */
			if ( x < 1e-10){
				return R*R/ 2.0 /(this->sigma_host + this->sigma_particle);
			}
			else{
				double rR_Rr = x + 1.0/x;
				double corr = (1.0 -1.0/pi * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );

				return R*R/ 2.0 /(this->sigma_host + this->sigma_particle)* corr;
			}
		}
	}


	/* 2D version for packings of polygonal cells. */


	/** \brief Run a fixed number of walk-on-spheres steps for a single test particle, logging debug info to file. (CLAUDE)
	 *  \param config      Polygonal-cell packing the test particle moves through.
	 *  \param num_steps   Number of jump steps to perform.
	 *  \param trajectory  Trajectory to advance and accumulate statistics into.
	 *  \param rng         Random number generator. */
	void PerformSimulation2D_test(const Dispersion & config, size_t num_steps, Trajectory & trajectory, RandomGenerator & rng);

	/** \brief Advance a single test particle until the accumulated first-passage time reaches tau_max. (CLAUDE)
	 *  \param config      Polygonal-cell packing the test particle moves through.
	 *  \param tau_max     Additional first-passage time to simulate.
	 *  \param trajectory  Trajectory to advance and accumulate statistics into.
	 *  \param rng         Random number generator. */
	void PerformSimulation2D_time(const Dispersion & config, double tau_max, Trajectory & trajectory, RandomGenerator & rng);

	/** \brief Run the first-passage-time simulation over many configurations in parallel and collect ensemble statistics. (CLAUDE)
	 *  \param GetConfiguration  Functor returning the polygonal-cell packing for a given configuration index.
	 *  \param num_configs       Number of configurations to process.
	 *  \param time_max          Target simulation time per test particle.
	 *  \param results           Output per-configuration statistics (mean time, mean X^2, and their standard errors).
	 *  \param log_filename      Base name for per-configuration log files; empty disables logging. */
	void PerformSimulation2D_parallel(std::function<const Dispersion(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	/** \brief Simulate a single configuration and record ensemble statistics at successive time slices. (CLAUDE)
	 *  \param config         Polygonal-cell packing to simulate.
	 *  \param time_max        Total simulation time.
	 *  \param time_interval   Time increment between recorded slices.
	 *  \param results         Output per-time-slice statistics (mean time, mean X^2, and their standard errors).
	 *  \param log_filename    Base name for the log file; empty disables logging. */
	void PerformSimulation2D_TimeEvolution(const Dispersion & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	/* 2D version for disk(sphere) packings. */

	/** \brief Run a fixed number of walk-on-spheres steps for a single test particle, logging debug info to file. (CLAUDE)
	 *  \param config      Disk (sphere) packing the test particle moves through.
	 *  \param num_steps   Number of jump steps to perform.
	 *  \param trajectory  Trajectory to advance and accumulate statistics into.
	 *  \param rng         Random number generator. */
	void PerformSimulation2D_test(const SpherePacking & config, size_t num_steps, Trajectory & trajectory, RandomGenerator & rng);

	/** \brief Advance a single test particle until the accumulated first-passage time reaches tau_max. (CLAUDE)
	 *  \param config      Disk (sphere) packing the test particle moves through.
	 *  \param tau_max     Additional first-passage time to simulate.
	 *  \param trajectory  Trajectory to advance and accumulate statistics into.
	 *  \param rng         Random number generator. */
	void PerformSimulation2D_time(const SpherePacking & config, double tau_max, Trajectory & trajectory, RandomGenerator & rng);

	/** \brief Run the first-passage-time simulation over many configurations in parallel and collect ensemble statistics. (CLAUDE)
	 *  \param GetConfiguration  Functor returning the disk (sphere) packing for a given configuration index.
	 *  \param num_configs       Number of configurations to process.
	 *  \param time_max          Target simulation time per test particle.
	 *  \param results           Output per-configuration statistics (mean time, mean X^2, and their standard errors).
	 *  \param log_filename      Base name for per-configuration log files; empty disables logging. */
	void PerformSimulation2D_parallel(std::function<const SpherePacking(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	/** \brief Simulate a single configuration and record ensemble statistics at successive time slices. (CLAUDE)
	 *  \param config         Disk (sphere) packing to simulate.
	 *  \param time_max        Total simulation time.
	 *  \param time_interval   Time increment between recorded slices.
	 *  \param results         Output per-time-slice statistics (mean time, mean X^2, and their standard errors).
	 *  \param log_filename    Base name for the log file; empty disables logging. */
	void PerformSimulation2D_TimeEvolution(const SpherePacking & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename = "");

protected:
	/** \brief Perform a single walk-on-spheres jump for a test particle in a polygonal-cell packing. (CLAUDE)
	 *  \param config       Polygonal-cell packing the test particle moves through.
	 *  \param trajectory   Trajectory to advance and accumulate statistics into.
	 *  \param rng          Random number generator.
	 *  \param debug_info   Output buffer filled with per-step diagnostics when debug is true.
	 *  \param debug        Whether to populate debug_info for this step. */
	void JumpOneStep_2D(const Dispersion & config, Trajectory & trajectory, RandomGenerator & rng, std::vector<double> & debug_info, bool debug = false);

	/** \brief Perform a single walk-on-spheres jump for a test particle in a disk (sphere) packing. (CLAUDE)
	 *  \param config       Disk (sphere) packing the test particle moves through.
	 *  \param trajectory   Trajectory to advance and accumulate statistics into.
	 *  \param rng          Random number generator.
	 *  \param debug_info   Output buffer filled with per-step diagnostics when debug is true.
	 *  \param debug        Whether to populate debug_info for this step. */
	void JumpOneStep_2D(const SpherePacking & config, Trajectory & trajectory, RandomGenerator & rng, std::vector<double> & debug_info, bool debug = false);
};

#endif