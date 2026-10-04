#ifndef __TAU_N_H__
#define __TAU_N_H__

/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: Feburary 2024 */

/** \file tau_n.h
 * \brief The tau_n class: tau_n order metric of 1D point patterns computed directly
 * from a ConfigPack, without storing g_n (which is memory-intensive).
 *
 * Unlike gn, positions are binned on a global grid of width dx (not relative to each
 * particle), and g_n at a tuple of relative bin offsets is estimated from differences of
 * global bin indices.  tau_n(L) accumulates [1 - g_n]^2 over tuples with max offset L.*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>

// uncomment to disable assert()
//#define NDEBUG 
#include <cassert>

// in $(core) 
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>

/** \brief Sparse occupancy of a periodic 1D grid: bins[bin index] = number of points. */
struct nlist
{
	std::map<int,int> bins;	///< occupied bins only (bin index -> count)
	int period = 0;			///< number of bins in the periodic box

	/** \brief Default constructor (period 0; unusable until reassigned). */
	nlist(){};

	/** \brief Empty grid with a given period.
	 *  \param num_bins  Number of bins in the periodic box. */
	nlist(int num_bins) {
		this->period = num_bins;
		this->bins = std::map<int,int>();
	}

	/** \brief Set the count of a bin (index taken modulo period).
	 *  \param key  Bin index (any integer).
	 *  \param val  New count. */
	inline void update(int key, int val){
		this->bins[this->modulo(key)] = val;
	}

	/** \brief Count of a bin (index taken modulo period); 0 if empty.
	 *  \param key  Bin index (any integer).
	 *  \return Number of points in the bin. */
	inline int get(int key) const {
		auto iter = this->bins.find(this->modulo(key));
		if (iter == this->bins.end())
			return 0;
		else
			return iter->second;
	}

	/** \brief Non-negative remainder of k modulo period.
	 *  \param k  Bin index.
	 *  \return k mod period in [0, period). */
	inline int modulo(int k) const {
		return ((k %= this->period) < 0) ? k+this->period : k;
	}
};

/** \brief tau_n order metric of an ensemble of 1D point patterns, iterating over bin
 *  tuples (resumable from a binary checkpoint) instead of storing g_n. */
struct tau_n
{
	size_t n = 0;				///< order n (2, 3, or 4)
	std::string bin_name;		///< checkpoint file name (<prefix>.tau_<n>)
	std::string savename;		///< prefix of the output text files
	size_t num_threads = 1;		///< OpenMP threads

	size_t Nc = 0;				///< number of configurations

	double Lmax = 0.;			///< largest offset L [length]
	double dx = 0.;				///< bin width [length]

	int num_half = 0;			///< floor(Lmax/dx)
	size_t len_row = 0;			///< 2*num_half+1
	size_t tot_num_tuples = 0;	///< total number of offset tuples (x_21, ..., x_n1) to visit

	/** \brief Default constructor (empty). */
	tau_n() {};

	/** \brief Interactive constructor: reads n, the ConfigPack, Nc, threads, output and
	 *  checkpoint prefixes, Lmax and dx from ifile (Lmax/dx are ignored when resuming).
	 *  \param ofile  Output stream for prompts.
	 *  \param ifile  Input stream with the parameters. */
	tau_n(std::ostream & ofile, std::istream & ifile);

	/** \brief Write the checkpoint (progress index, grid, partial sums) to bin_name. */
	void WriteBinary() const;

	/** \brief Load the checkpoint from bin_name if it exists.
	 *  \return true if a checkpoint was loaded. */
	bool ReadBinary();

	/** \brief Write tau_n(L) and its error to <savename>-tau<n>. */
	void Write_tau_n();

	/** \brief Convert a row-major tuple index to bin offsets.
	 *  \param idx       Tuple index.
	 *  \param[out] pos  Offsets in units of dx (size n-1). */
	inline void idx2pos(int idx, std::vector<int> & pos){
		pos.resize(this->n-1);
		int val = idx;
		for(int i = this->n-2; i >= 0; i--){
			//pos[this->n -2-i] = (val % len_row) - num_half;
			if (i == 0)
				pos.at(i) = (val % len_row);
			else
				pos.at(i) = (val % len_row) - num_half;
			val /= len_row;
		}	
	}

	/** \brief Accumulate tau_n over offset tuples until done or the time limit is near,
	 *  then save the checkpoint and write tau_n.
	 *  \param ofile      Output stream for progress messages.
	 *  \param TimeLimit  Wall-clock deadline (absolute time_t) [s].
	 *  \param output_gn  Also write g_n per tuple to <savename>_g<n>.txt. */
	void Run(std::ostream & ofile, time_t TimeLimit, bool output_gn=false);

private:
	size_t curr_idx = 0;					///< next tuple index to visit
	std::string cp_name;					///< input ConfigPack name
	std::vector<GeometryVector> result;		///< result[L] = (L, sum of [1-g_n]^2, error sum)
	double normalization_factor = 0;		///< 1/dx^{n-1}
	ConfigurationPack CP;					///< input ensemble

	/** \brief Unnormalized g_n of one configuration at given bin offsets, from its occupancy grid.
	 *  \param nlist_            Occupancy grid of the configuration.
	 *  \param idx_relative_pos  Bin offsets i_2, ..., i_n.
	 *  \param idx_overlapped    Overlap indices: -1 = distinct bin, k = same bin as point k.
	 *  \return Sum over x_1 of the product of effective counts, divided by N. */
	double gn_from_nlist (const nlist & nlist_, const std::vector<int> & idx_relative_pos, const std::vector<int> idx_overlapped);
};

#endif // __TAU_N_H__