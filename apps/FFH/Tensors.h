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

	ConstantTensor();
	ConstantTensor(size_t dimension, size_t rank);
	ConstantTensor(const ConstantTensor & src);

	inline std::complex<double> & operator[](size_t i) {
		if (i < total_comp)
			return value[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return value[0];
		}
	}

	inline const std::complex<double> operator[](size_t i) const {
		if (i < total_comp)
			return value[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return std::complex<double>(0, 0);
		}
	}
	ConstantTensor & operator = (const ConstantTensor & rhs);

};

/* operators */
ConstantTensor operator + (const ConstantTensor & left, const ConstantTensor & right);
ConstantTensor operator - (const ConstantTensor & left, const ConstantTensor & right);
ConstantTensor operator * (const std::complex<double> & left, const ConstantTensor & right);

/* @return the hydrostatic projection tensor. */
ConstantTensor Lambda_h(size_t dimension);

/* @return the shear projection tensor. */
ConstantTensor Lambda_s(size_t dimension);

ConstantTensor Identity(size_t dimension);

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

	TensorField();
	TensorField(size_t dimension, size_t rank, const std::vector<size_t> & num_points);
	TensorField(const TensorField & src);
	~TensorField();

	void Plus(TensorField & out, const TensorField & A, const TensorField & B);

	/* directly access to the component of tensor field with components.
	* i=(1,2) -> index = 2*d + 1. */
	inline fftw_complex * comp(const std::vector<size_t> &i);
	inline const fftw_complex * comp(const std::vector<size_t> &i) const;

	TensorField & operator += (const TensorField &rhs);
	TensorField & operator -= (const TensorField & rhs);
	TensorField & operator = (const TensorField & src);
	TensorField & operator *= (const std::complex<double> rhs);

	/* @return The tensor field of the i-th component. */
	inline fftw_complex* operator[](size_t i) {
		if (i < total_comp)
			return field[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return field[0];
		}
	}

	inline const fftw_complex* operator[](size_t i) const {
		if (i < total_comp)
			return field[i];
		else {
			std::cerr << "ConstantTensor::[] cannot access\n";
			return field[0];
		}
	}
	/* return the position index corresponding to the comp. */
	inline size_t at(const std::vector<size_t>& comp) const
	{
		assert(comp.size() == d);
		size_t out = comp[0];
		for (int i = 1; i < d; i++) {
			out = comp[i] + this->num_points[i] * out;
		}

		return out;
	}

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
	/*
	component -1 = all,
	part	0 = real, 1 = imaginary, 2 = both. */
	void Print(const std::string & prefix, int component = -1, int part = 0) const;
	
	/* Print tensor fields with the origin placed at the center of the plot.
	 It will be useful to visualize tensorfield in the Fourier domain. */
	void Print_shifted_origin(const std::string & prefix, int component = -1, int part = 0) const;
	
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
inline TensorField operator*(const std::complex<double>& left, const TensorField & right)
{
	TensorField out(right);
	out *= left;
	return out;
}

inline TensorField operator+(const TensorField & left, const TensorField & right)
{
	TensorField out(left);
	out += right;
	return out;
}

inline TensorField operator-(const TensorField & left, const TensorField & right)
{
	TensorField out(left);
	out -= right;
	return out;
}

void InnerProduct(TensorField & out, const TensorField & A, const TensorField &B, size_t num_contraction);
TensorField & InnerProduct_conjugate(const TensorField & A, const TensorField &B, size_t num_contraction);

std::istream & operator >> (std::istream & input, const fftw_complex & rhs);
std::ostream & operator << (std::ostream & input, const fftw_complex & rhs);

ConstantTensor MeanField(const TensorField & src);

/* Compute the inner product of \sum_{i} A[i]B[i][position] at a given position and component.
* @param[in] position	The place where the inner product is computed.
* @param[in] component The component of the resulting tensor.
* @param[in] A	a constant tensor.
* @param[in] B	a tensor field
* @param[in] num_contraction	The number of indices that are contracted.
* @return the complex-valued component. */
std::complex<double> InnerProduct(size_t position, size_t component, const ConstantTensor & A, const TensorField & B, size_t num_comp);

/* Compute the inner product between a constant tensor A and a tensor field B at each position. */
void InnerProduct(TensorField & out, const ConstantTensor & A, const TensorField & B, size_t num_contraction);

/* Compute the inner product of two constant tensors A and B.
* @return A resulting tensor.
* @param[in] A	a constant tensor.
* @param[in] B	a tensor field
* @param[in] num_contraction	The number of indices that are contracted.*/
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