#ifndef TWOPHASEDISPERSION_INCLUDED
#define TWOPHASEDISPERSION_INCLUDED

/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	March	2021 */

/** \file TwoPhaseDispersion.h
*	\brief Header file to define a two-phase medium 
*	under periodic boundary conditions and finite-difference method. 
*	The TwoPhaseDispersion class is separated from "FFH.h" file.
*/


#include <vector>
#include <omp.h>
#include <fstream>

#include "Tensors.h"
#include <PeriodicCellList.h>
#include <RandomGenerator.h>


/* A class to store a sphere packing and discretize it. */
struct TwoPhaseDispersion {
	SpherePacking packing;	/*< dispersions with the retangular simulation box */
	std::vector<size_t> prt_indicator;		/* prt_indicator[i] = index of the pixel
											*	belongs to the particle phase.	*/
	std::vector<double> weights;	/* The weights obtained by a subpixel smoothing method.*/
	bool update_indicator = false;
											//	TwoPhaseDispersion();
	bool use_weights = false;	//True, if the weights variable is defined.

	/** \brief Construct an empty two-phase dispersion (empty packing, zero resolution). (CLAUDE) */
	TwoPhaseDispersion() :TwoPhaseDispersion(SpherePacking(), 0) {};
	/** \brief Construct a dispersion from a sphere packing with a uniform per-side resolution. (CLAUDE)
	 *  \param packing  The sphere packing defining the inclusion phase and simulation box.
	 *  \param num      Number of voxels along each side of the grid. */
	TwoPhaseDispersion(const SpherePacking & packing, size_t num);
	/** \brief Copy constructor. (CLAUDE)
	 *  \param src  The dispersion to copy. */
	TwoPhaseDispersion(const TwoPhaseDispersion & src) : packing(src.packing), prt_indicator(src.prt_indicator),
		num_side_voxels(src.num_side_voxels), num_total_voxels(src.num_total_voxels) {
		this->update_indicator = src.update_indicator;
	};
	/** \brief Construct a dispersion by reading a discretized array from a file. (CLAUDE)
	 *  \param filename  Prefix of the .DisArray file to read. */
	TwoPhaseDispersion(const std::string & filename): TwoPhaseDispersion()
	{
		this->Read(filename);
	}
	//~TwoPhaseDispersion() {};

	/** \brief Set the grid resolution from a per-dimension list of voxel counts. (CLAUDE)
	 *  \param num_voxels  Number of voxels along each side. */
	void SetResolutions(const std::vector<size_t> & num_voxels);
	/** \brief Set an isotropic grid resolution (same voxel count on every side). (CLAUDE)
	 *  \param num_side_voxels  Number of voxels along each side. */
	void SetResolutions(const size_t num_side_voxels) {
		this->SetResolutions(std::vector<size_t>(this->packing.GetDimension(), num_side_voxels));
	};
	/** \brief Adjust cell sizes to be optimal for FFTW (not yet implemented). (CLAUDE) */
	void OptimizeSize();
	/** \brief Identify the voxels whose centers lie inside the particle phase. (CLAUDE) */
	void ComputeIndicator();

	/** \brief Return the phase-indicator weight of the voxel at the given index. (CLAUDE)
	 *  \param index  Index into the prt_indicator / weights array.
	 *  \return       The subpixel weight if weights are defined, otherwise 1. */
	double WeightedPhaseIndicator(size_t index) const;
	/** \brief A function to convert a dD component to an index of 1D array (for fftw).
	 *  \param[in] component	=[n1,n2,...nd].
	 *  \return	size_t	Corresponding index.  */
	inline size_t Components2Index(const std::vector<size_t> & component) const;
	/** \brief A function to convert an index of 1D array (for fftw) to a dD component.
	 *  \param[in] index	of 1D array (fftw_complex).
	 *  \param[out]	components	the resulting dD component. */
	inline void Index2Components(const size_t index, std::vector<size_t> & components) const;
	/** \brief Return the per-side voxel resolutions of the grid. (CLAUDE) */
	const std::vector<size_t> & GetResolutions() const { return this->num_side_voxels; };
	/* Save a text file to draw a scatter plot. 
	 * Date format:
	 * # dimension
	 * # resolutions
	 * # basis...
	 * # components \t ratio of the inclusion phase.
	 */
	void Print(const std::string & prefix) const;
	/** \brief Read a discretized two-phase array (basis, resolution, voxel weights) from a .DisArray file. (CLAUDE)
	 *  \param prefix  Prefix of the .DisArray file (without extension). */
	void Read(const std::string & prefix);

	/** \brief Build the rank-1 tensor field of DFT wavevectors for this grid. (CLAUDE)
	 *  \return  A tensor field whose i-th component holds the i-th wavevector coordinate at each pixel. */
	TensorField GetKs() const;
private:
	std::vector<size_t> num_side_voxels;	/*< the number of pixels in each side */
	size_t num_total_voxels;				/*< total number of voxels */

};



inline size_t TwoPhaseDispersion::Components2Index(const std::vector<size_t>& component) const
{	/* raw-major format !! */
	if (component.size() == 3) {
		return component[2] + this->num_side_voxels[2] * (component[1] + this->num_side_voxels[1] * component[0]);
	}
	else if (component.size() == 2) {
		return component[1] + this->num_side_voxels[1] * component[0];
	}
	return size_t();
}

inline void TwoPhaseDispersion::Index2Components(const size_t index, std::vector<size_t> & components) const
{
	components = std::vector<size_t>(this->num_side_voxels.size());
	size_t value = index;
	for (size_t i = components.size() - 1; i > 0; i--) {
		components[i] = value % (this->num_side_voxels[i]);
		value = (value - components[i]) / this->num_side_voxels[i];
	}
	components[0] = value;
}


#endif