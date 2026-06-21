#ifndef TENSOR_INCLUDED
#define TENSOR_INCLUDED

/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	November	2019 */

/** \file Tensors.h
*	\brief Header file to define classes for tensors and tensor fields, 
* and assocaited operations. Each component is complex-valued.
* ConstantTensor and TensorFields are designed to interface Tensor and fftw_complex in FFTW3. */

/*=====================================
	Notations:

	In this code, all tensors are linearized to be flexible for change rank and dimensions.
	In other words, the component indices (i0 i1 i2 i3) of a rank-4 tensor T in d dimensions
	corresponds to the (i0 + d*(i1 + d*(i2 + d*(i3))))-th element of the object T.
	In the rest of this code, for brevity, we call the former the component of T ("T_comp");
	the latter as the index of T ("T_idx").

	For a tensor field object T, each element T[i] stands for a scalar field 
	in the i-th "index" of the T tensor field.
	Each element of T then stores the complex-valued field information 
	in a linearized form in the same manner.
	
  =====================================*/
#include <complex>	/* this should be compiled earlier than fftw3!! */
#include <cmath>
#include <fftw3.h>	/* library for computing the discerete fast Fourier transform. */
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <assert.h>


const size_t size_fftw = sizeof(fftw_complex);


/* \brief A class to define constant tensors.
 * Each components are complex-valued. */
struct ConstantTensor {
	std::vector<std::complex<double>> value;
	size_t d, rank, total_comp;

	/** \brief Construct an empty constant tensor (dimension 0, rank 0). (CLAUDE) */
	ConstantTensor();
	/** \brief Construct a zero-initialized constant tensor of given dimension and rank. (CLAUDE)
	 *  \param dimension  Spatial dimension d.
	 *  \param rank       Tensor rank (number of indices). */
	ConstantTensor(size_t dimension, size_t rank);
	/** \brief Copy constructor. (CLAUDE)
	 *  \param src  The constant tensor to copy. */
	ConstantTensor(const ConstantTensor & src);

	/** \brief Mutable access to the i-th linearized component. (CLAUDE)
	 *  \param i  Linearized component index.
	 *  \return   Reference to the component (or component 0 if out of range). */
	inline std::complex<double> & operator[](size_t i) {
		if (i < total_comp)
			return value[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return value[0];
		}
	}

	/** \brief Read-only access to the i-th linearized component. (CLAUDE)
	 *  \param i  Linearized component index.
	 *  \return   The component value (or 0 if out of range). */
	inline const std::complex<double> operator[](size_t i) const {
		if (i < total_comp)
			return value[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return std::complex<double>(0, 0);
		}
	}
	/** \brief Copy-assignment operator. (CLAUDE)
	 *  \param rhs  The constant tensor to assign from.
	 *  \return     Reference to this tensor. */
	ConstantTensor & operator = (const ConstantTensor & rhs);

};

/* operators */
/** \brief Component-wise sum of two constant tensors. (CLAUDE)
 *  \param left   Left operand.
 *  \param right  Right operand.
 *  \return       The element-wise sum. */
ConstantTensor operator + (const ConstantTensor & left, const ConstantTensor & right);
/** \brief Component-wise difference of two constant tensors. (CLAUDE)
 *  \param left   Left operand.
 *  \param right  Right operand.
 *  \return       The element-wise difference (left - right). */
ConstantTensor operator - (const ConstantTensor & left, const ConstantTensor & right);
/** \brief Scalar multiplication of a constant tensor. (CLAUDE)
 *  \param left   Complex scalar factor.
 *  \param right  The constant tensor.
 *  \return       The scaled tensor. */
ConstantTensor operator * (const std::complex<double> & left, const ConstantTensor & right);

/** \brief Return the hydrostatic projection tensor.
 *  \param dimension  Spatial dimension d.
 *  \return           The rank-4 hydrostatic projection tensor. */
ConstantTensor Lambda_h(size_t dimension);

/** \brief Return the shear projection tensor.
 *  \param dimension  Spatial dimension d.
 *  \return           The rank-4 shear projection tensor. */
ConstantTensor Lambda_s(size_t dimension);

/** \brief Return the rank-2 identity tensor in d dimensions. (CLAUDE)
 *  \param dimension  Spatial dimension d.
 *  \return           The identity tensor. */
ConstantTensor Identity(size_t dimension);

/** \brief Return the trace of a rank-2 constant tensor. (CLAUDE)
 *  \param src  A rank-2 constant tensor.
 *  \return     The sum of its diagonal components. */
inline std::complex<double> Trace(const ConstantTensor & src) {
	std::complex<double> out(0.0L, 0.0L);
	if (src.rank == 2) {
		for (int i = 0; i < src.d; i++) {
			out += src[i*(1 + src.d)];
		}
	}
	else
		std::cerr << "src should be a rank-2 tensor\n";
	
	return out;
}
/* \brief A class to define Tensor field.
 * the field*/
struct TensorField {
	std::vector<fftw_complex *> field;
	size_t d, rank, total_comp;
	std::vector<size_t> num_points;	/* store size of fields */
	size_t num_total_pixels;

	/** \brief Construct an empty tensor field. (CLAUDE) */
	TensorField();
	/** \brief Construct a zero-initialized tensor field over a grid. (CLAUDE)
	 *  \param dimension   Spatial dimension d.
	 *  \param rank        Tensor rank (number of indices).
	 *  \param num_points  Number of grid points along each dimension. */
	TensorField(size_t dimension, size_t rank, const std::vector<size_t> & num_points);
	/** \brief Copy constructor (deep-copies the fftw buffers). (CLAUDE)
	 *  \param src  The tensor field to copy. */
	TensorField(const TensorField & src);
	/** \brief Destructor; frees the allocated fftw buffers. (CLAUDE) */
	~TensorField();

	/** \brief Store the element-wise sum of two tensor fields into out. (CLAUDE)
	 *  \param out  Output tensor field receiving A + B.
	 *  \param A    First operand.
	 *  \param B    Second operand. */
	void Plus(TensorField & out, const TensorField & A, const TensorField & B);

	/** \brief directly access to the component of tensor field with components.
	 *  i=(1,2) -> index = 2*d + 1.
	 *  \param i  Per-index component vector.
	 *  \return   Pointer to the scalar field of that component. */
	inline fftw_complex * comp(const std::vector<size_t> &i);
	/** \brief directly access to the component of tensor field with components (const overload).
	 *  i=(1,2) -> index = 2*d + 1.
	 *  \param i  Per-index component vector.
	 *  \return   Pointer to the scalar field of that component. */
	inline const fftw_complex * comp(const std::vector<size_t> &i) const;

	/** \brief Add another tensor field element-wise in place. (CLAUDE)
	 *  \param rhs  Tensor field to add.
	 *  \return     Reference to this tensor field. */
	TensorField & operator += (const TensorField &rhs);
	/** \brief Subtract another tensor field element-wise in place. (CLAUDE)
	 *  \param rhs  Tensor field to subtract.
	 *  \return     Reference to this tensor field. */
	TensorField & operator -= (const TensorField & rhs);
	/** \brief Copy-assignment (reallocates buffers if sizes differ). (CLAUDE)
	 *  \param src  Tensor field to assign from.
	 *  \return     Reference to this tensor field. */
	TensorField & operator = (const TensorField & src);
	/** \brief Multiply every component by a complex scalar in place. (CLAUDE)
	 *  \param rhs  Complex scalar factor.
	 *  \return     Reference to this tensor field. */
	TensorField & operator *= (const std::complex<double> rhs);

	/** \brief The tensor field of the i-th component.
	 *  \param i  Linearized component index.
	 *  \return   Pointer to the i-th scalar field (or component 0 if out of range). */
	inline fftw_complex* operator[](size_t i) {
		if (i < total_comp)
			return field[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return field[0];
		}
	}

	/** \brief The tensor field of the i-th component (const overload). (CLAUDE)
	 *  \param i  Linearized component index.
	 *  \return   Pointer to the i-th scalar field (or component 0 if out of range). */
	inline const fftw_complex* operator[](size_t i) const {
		if (i < total_comp)
			return field[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return field[0];
		}
	}
	/** \brief return the position index corresponding to the comp.
	 *  \param comp  Per-dimension grid coordinates.
	 *  \return      The linearized pixel index. */
	inline size_t at(const std::vector<size_t>& comp) const
	{
		assert(comp.size() == d);
		size_t out = comp[0];
		for (int i = 1; i < d; i++) {
			out = comp[i] + this->num_points[i] * out;
		}

		return out;
	}

	/** \brief Convert a linearized pixel index to per-dimension grid coordinates. (CLAUDE)
	 *  \param index  Linearized pixel index.
	 *  \param[out] out  The resulting grid coordinates. */
	inline void GetPositionComponent(size_t index, std::vector<size_t> & out) const {
		assert(index < this->num_total_pixels);
		out = std::vector<size_t>(this->num_points);
		size_t value = index;
		for (size_t i = out.size() - 1; i > 0; i--) {
			out[i] = value % (this->num_points[i]);
			value = (value - out[i]) / this->num_points[i];
		}
		out[0] = value;
	}
	/** \brief Write tensor-field components to text files.
	 *  \param prefix     Output filename prefix.
	 *  \param component  Component to print; -1 = all.
	 *  \param part       0 = real, 1 = imaginary, 2 = both. */
	void Print(const std::string & prefix, int component = -1, int part = 0) const;

	/** \brief Print tensor fields with the origin placed at the center of the plot.
	 *  It will be useful to visualize tensorfield in the Fourier domain.
	 *  \param prefix     Output filename prefix.
	 *  \param component  Component to print; -1 = all.
	 *  \param part       0 = real, 1 = imaginary, 2 = both. */
	void Print_shifted_origin(const std::string & prefix, int component = -1, int part = 0) const;

	/** \brief Convert a linearized component index to its per-index components. (CLAUDE)
	 *  \param index       Linearized component index.
	 *  \param[out] components  The resulting rank-length component vector. */
	inline void GetComponent(size_t index, std::vector<size_t> & components) const {
		components = std::vector<size_t>(rank, 0);
		size_t value = index;
		if (index < total_comp) {
			for (int i = 0; i < rank; i++) {
				components[i] = value % d;
				value -= components[i];	value /= d;
			}
		}
		else
			std::cerr << "TensorField::index exceeds the range\n";
	}
};

/* 3 basic operations of fftw_complex */
/** \brief Scalar multiplication of a tensor field. (CLAUDE)
 *  \param left   Complex scalar factor.
 *  \param right  The tensor field.
 *  \return       The scaled tensor field. */
inline TensorField operator*(const std::complex<double>& left, const TensorField & right)
{
	TensorField out(right);
	out *= left;
	return out;
}

/** \brief Element-wise sum of two tensor fields. (CLAUDE)
 *  \param left   Left operand.
 *  \param right  Right operand.
 *  \return       The element-wise sum. */
inline TensorField operator+(const TensorField & left, const TensorField & right)
{
	TensorField out(left);
	out += right;
	return out;
}

/** \brief Element-wise difference of two tensor fields. (CLAUDE)
 *  \param left   Left operand.
 *  \param right  Right operand.
 *  \return       The element-wise difference (left - right). */
inline TensorField operator-(const TensorField & left, const TensorField & right)
{
	TensorField out(left);
	out -= right;
	return out;
}

/** \brief Per-pixel inner product of two tensor fields A and B. (CLAUDE)
 *  \param[out] out  Resulting tensor field of the contracted product.
 *  \param A    First tensor field.
 *  \param B    Second tensor field.
 *  \param num_contraction  Number of indices contracted between A and B. */
void InnerProduct(TensorField & out, const TensorField & A, const TensorField &B, size_t num_contraction);
/** \brief Per-pixel inner product of A with the complex conjugate of B. (CLAUDE)
 *  \param A    First tensor field.
 *  \param B    Second tensor field (conjugated in the contraction).
 *  \param num_contraction  Number of indices contracted between A and B.
 *  \return     Reference to the resulting tensor field. */
TensorField & InnerProduct_conjugate(const TensorField & A, const TensorField &B, size_t num_contraction);

/** \brief Stream-extract a complex value into an fftw_complex. (CLAUDE)
 *  \param input  Input stream.
 *  \param rhs    The fftw_complex to fill.
 *  \return       The input stream. */
std::istream & operator >> (std::istream & input, const fftw_complex & rhs);
/** \brief Stream-insert an fftw_complex value. (CLAUDE)
 *  \param input  Output stream.
 *  \param rhs    The fftw_complex to print.
 *  \return       The output stream. */
std::ostream & operator << (std::ostream & input, const fftw_complex & rhs);

/** \brief Compute the spatial mean of a tensor field. (CLAUDE)
 *  \param src  The tensor field.
 *  \return     The constant tensor of per-component averages. */
ConstantTensor MeanField(const TensorField & src);

/** \brief Compute the inner product of \sum_{i} A[i]B[i][position] at a given position and component.
 *  \param[in] position	The place where the inner product is computed.
 *  \param[in] component The component of the resulting tensor.
 *  \param[in] A	a constant tensor.
 *  \param[in] B	a tensor field
 *  \param[in] num_comp	The number of indices that are contracted.
 *  \return the complex-valued component. */
std::complex<double> InnerProduct(size_t position, size_t component, const ConstantTensor & A, const TensorField & B, size_t num_comp);

/** \brief Compute the inner product between a constant tensor A and a tensor field B at each position. (CLAUDE)
 *  \param[out] out  Resulting tensor field.
 *  \param A    The constant tensor.
 *  \param B    The tensor field.
 *  \param num_contraction  Number of indices contracted. */
void InnerProduct(TensorField & out, const ConstantTensor & A, const TensorField & B, size_t num_contraction);

/** \brief Compute the inner product of two constant tensors A and B.
 *  \param[in] A	a constant tensor.
 *  \param[in] B	a tensor field
 *  \param[in] num_comp	The number of indices that are contracted.
 *  \return A resulting tensor. */
ConstantTensor InnerProduct(const ConstantTensor & A, const ConstantTensor & B, size_t num_comp);




inline fftw_complex * TensorField::comp(const std::vector<size_t>& i)
{	/*i[0] + d*i[1] + d^2 *i[2] + ... */
	size_t index = 0;
	for (size_t r = this->rank - 1; r >= 0; r--)
		index = this->d*index + i[r];

	//std::cout << index << "\t";
	return &field[index][0];
}

inline const fftw_complex * TensorField::comp(const std::vector<size_t>& i) const
{
	size_t index = 0;
	for (size_t r = this->rank - 1; r >= 0; r--)
		index = this->d*index + i[r];

	//std::cout << index << "--\t";
	return &field[index][0];
}




#endif