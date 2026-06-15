#ifndef FFH2_INCLUDED
#define FFH2_INCLUDED

#include "FFH.h"

/* a class to test the truncated Green's function... */
struct DynamicDielectricProblem_padded : public StaticDielectricProblem {
	
	DynamicDielectricProblem_padded();
	DynamicDielectricProblem_padded(const TwoPhaseDispersion & src,
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	~DynamicDielectricProblem_padded() {
		Propagator_direct.~TensorField();
	}
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1);
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2, const GeometryVector & k1);
	void SetWavevector(const GeometryVector & k1);

	virtual void Prepare(const ConstantTensor & applied_electric_field);

	void CropTensorField_direct(TensorField & out, const TensorField & in) const;
	void CropTensorField_Fourier(TensorField & out, const TensorField & in) const;
protected:
	GeometryVector k1;	/* wavevector of the electric field in the reference phase!! */
	TensorField Propagator_direct;

	/* Compute the displacement field D. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;
	/* It should multiply e^{ik1*x} */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const {
		TensorField aux;
		StaticDielectricProblem::ComputeDivergenceFreeField_Direct(aux);
		CropTensorField_direct(out, aux);
		ShiftFieldByk1(out);
	};
	/* Compute the electric field E. */
	virtual void ComputeField_Fourier(TensorField & out) const;
	virtual void ComputeField_Direct(TensorField & out) const {
		TensorField aux;
		StaticDielectricProblem::ComputeField_Direct(aux);
		CropTensorField_direct(out, aux);
		ShiftFieldByk1(out);
	};

	/* Fourier transform of the Green's function.
	 * G_{ij}(q+k_1) of the free-space Green's function . */
	virtual void ComputePropagator_weak();
	//virtual void  ComputeL_weak();

	/* Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	//virtual void ComputeL_strong();


	inline ConstantTensor Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx);

private:
	/* Shift the tensor field in the direct space by k1.*/
	void ShiftFieldByk1(TensorField & out) const;
};

#endif