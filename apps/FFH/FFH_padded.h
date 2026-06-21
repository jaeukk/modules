#ifndef FFH2_INCLUDED
#define FFH2_INCLUDED

#include "FFH.h"

/* a class to test the truncated Green's function... */
struct DynamicDielectricProblem_padded : public StaticDielectricProblem {
	
	/** \brief Construct an empty zero-padded dynamic dielectric problem. (CLAUDE) */
	DynamicDielectricProblem_padded();
	/** \brief Construct a zero-padded dynamic dielectric problem from a geometry and two permittivities. (CLAUDE)
	 *  \param src       The two-phase dispersion (geometry).
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase. */
	DynamicDielectricProblem_padded(const TwoPhaseDispersion & src,
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	/** \brief Destructor; frees the direct-space propagator. (CLAUDE) */
	~DynamicDielectricProblem_padded() {
		Propagator_direct.~TensorField();
	}
	/** \brief Define the problem geometry on a zero-padded grid and set the wavevector. (CLAUDE)
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Constant tensor for the matrix phase.
	 *  \param inclusion_phase  Constant tensor for the inclusion phase.
	 *  \param k1               Wavevector of the field in the reference phase. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1);
	/** \brief Set the permittivities and the incident wavevector. (CLAUDE)
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase.
	 *  \param k1        Wavevector of the field in the reference phase. */
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2, const GeometryVector & k1);
	/** \brief Set the incident wavevector and shift the DFT wavevectors accordingly. (CLAUDE)
	 *  \param k1  Wavevector of the field in the reference phase. */
	void SetWavevector(const GeometryVector & k1);

	/** \brief Prepare the solver, requiring the applied field to be perpendicular to k1. (CLAUDE)
	 *  \param applied_electric_field  The applied (constant) electric field. */
	virtual void Prepare(const ConstantTensor & applied_electric_field);

	/** \brief Crop a direct-space field from the padded grid back to the physical grid. (CLAUDE)
	 *  \param[out] out  The cropped field on the physical grid.
	 *  \param in        The field on the padded grid. */
	void CropTensorField_direct(TensorField & out, const TensorField & in) const;
	/** \brief Crop a Fourier-space field from the padded grid back to the physical grid. (CLAUDE)
	 *  \param[out] out  The cropped field on the physical grid.
	 *  \param in        The field on the padded grid. */
	void CropTensorField_Fourier(TensorField & out, const TensorField & in) const;
protected:
	GeometryVector k1;	/* wavevector of the electric field in the reference phase!! */
	TensorField Propagator_direct;

	/** \brief Compute the displacement field D in Fourier space (cropped to the physical grid). (CLAUDE)
	 *  \param[out] out  The resulting displacement field. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;
	/** \brief Compute the displacement field D in direct space, cropped and multiplied by e^{ik1*x}. (CLAUDE)
	 *  \param[out] out  The resulting displacement field. */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const {
		TensorField aux;
		StaticDielectricProblem::ComputeDivergenceFreeField_Direct(aux);
		CropTensorField_direct(out, aux);
		ShiftFieldByk1(out);
	};
	/** \brief Compute the electric field E in Fourier space (cropped to the physical grid). (CLAUDE)
	 *  \param[out] out  The resulting electric field. */
	virtual void ComputeField_Fourier(TensorField & out) const;
	/** \brief Compute the electric field E in direct space, cropped and multiplied by e^{ik1*x}. (CLAUDE)
	 *  \param[out] out  The resulting electric field. */
	virtual void ComputeField_Direct(TensorField & out) const {
		TensorField aux;
		StaticDielectricProblem::ComputeField_Direct(aux);
		CropTensorField_direct(out, aux);
		ShiftFieldByk1(out);
	};

	/** \brief Fourier transform of the Green's function.
	 *  G_{ij}(q+k_1) of the free-space Green's function . */
	virtual void ComputePropagator_weak();
	//virtual void  ComputeL_weak();

	/** \brief Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	//virtual void ComputeL_strong();


	/** \brief Evaluate the direct-space H-function (double-gradient Green's function). (CLAUDE)
	 *  \param position  Direct-space position.
	 *  \param k1        Wavevector of the field in the reference phase.
	 *  \param d         Spatial dimension.
	 *  \param dx        Grid spacing (used for the self-interaction at the origin).
	 *  \return          The rank-2 H-tensor at the position. */
	inline ConstantTensor Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx);

private:
	/** \brief Shift the tensor field in the direct space by k1.
	 *  \param[in,out] out  The tensor field multiplied by the e^{ik1*x} phase factor. */
	void ShiftFieldByk1(TensorField & out) const;
};

#endif