/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	November	2019 */

/** \file Tensors.cpp
*	\brief Function implementations of ConstantTensor and TensorField. */

#include "Tensors.h"

/* Implementation of ConstantTensor!! */
ConstantTensor::ConstantTensor() : d(0), rank(0), value(std::vector < std::complex<double> >()) {}

ConstantTensor::ConstantTensor(size_t dimension, size_t rank) : d(dimension), rank(rank), total_comp((size_t)pow(dimension, rank))
{
	value = std::vector<std::complex<double>>(this->total_comp, std::complex<double>(0.0L, 0.0L));
}

ConstantTensor::ConstantTensor(const ConstantTensor & src) : ConstantTensor(src.d, src.rank)
{
	value = std::vector<std::complex<double>>(src.value);
}

ConstantTensor & ConstantTensor::operator=(const ConstantTensor & rhs)
{
	this->d = rhs.d;	this->rank = rhs.rank; this->total_comp = rhs.total_comp;
	this->value.resize(total_comp);
	for (int i = 0; i < total_comp; i++)
		this->value[i] = rhs.value[i];

	return *this;
}



ConstantTensor operator+(const ConstantTensor & left, const ConstantTensor & right)
{
	assert(left.d == right.d && left.rank == right.rank);
	ConstantTensor result(left);
	for (size_t i = 0; i < result.total_comp; i++)
		result[i] += right[i];

	return result;
}

ConstantTensor operator-(const ConstantTensor & left, const ConstantTensor & right)
{
	assert(left.d == right.d && left.rank == right.rank);
	ConstantTensor result(left);
	for (size_t i = 0; i < result.total_comp; i++)
		result[i] -= right[i];

	return result;
}

ConstantTensor operator*(const std::complex<double>& left, const ConstantTensor & right)
{
	ConstantTensor result(right);
	for (size_t i = 0; i < result.total_comp; i++)
		result[i] *= left;

	return result;
}

ConstantTensor Lambda_h(size_t dimension) {
	ConstantTensor out(dimension, 4);
	for (size_t i = 0; i < dimension; i++)
		for (size_t j = 0; j < dimension; j++)
			out[i*(1 + dimension) + dimension * dimension*(1 + dimension)*j] = std::complex<double>(1.0L, 0.0L);
	out = 1.0 / (double)dimension * out;
	return out;
}

ConstantTensor Lambda_s(size_t dimension) {
	ConstantTensor out(dimension, 4);
	for (size_t i = 0; i < dimension; i++)
		for (size_t j = 0; j < dimension; j++)
			for (size_t k = 0; k < dimension; k++)
				for (size_t l = 0; l < dimension; l++) {
					if (i == k && j == l)
						out[i + dimension * (j + dimension * (k + dimension * l))] += std::complex<double>(1.0L, 0.0L);
					if (i == l && j == k)
						out[i + dimension * (j + dimension * (k + dimension * l))] += std::complex<double>(1.0L, 0.0L);
				}

	return 0.5*out - Lambda_h(dimension);
}

ConstantTensor Identity(size_t dimension)
{
	ConstantTensor out(dimension, 2);
	for (int i = 0; i < out.d; i++) {
		out[i*(1 + dimension)] = std::complex<double>(1.0L, 0.0);
	}
	return out;
}

std::complex<double> InnerProduct(size_t position, size_t component, const ConstantTensor & A, const TensorField & B, size_t num_contraction)
{
	std::complex<double> result(0.0L, 0.0L);
	if (A.d == B.d && A.rank >= num_contraction &&
		B.rank >= num_contraction) {

		/* component = A_front_index + d^(rank of A - contracted ranks) *(B_back_index)  */
		size_t A_fixed_index = A.total_comp, dummy_index = 1;
		for (int i = 0; i < num_contraction; i++) {
			A_fixed_index /= A.d;
			dummy_index *= A.d;
		}

		size_t A_front_index = component % A_fixed_index,
			B_back_index = (size_t)floor(component / A_fixed_index);

		for (size_t dummy = 0; dummy < dummy_index; dummy++) {
			size_t index_of_A = A_front_index + A_fixed_index * dummy,
				index_of_B = dummy + dummy_index * B_back_index;
			result += (A[index_of_A] *
				std::complex<double>(B.field[index_of_B][position][0], B.field[index_of_B][position][1]));
		}
		return result;
	}
	else {
		std::cerr << "Inner product is not available.\n";
	}
}

void InnerProduct(TensorField & out, const ConstantTensor & A, const TensorField & B, size_t num_contraction)
{
	if (A.d == B.d && A.rank >= num_contraction && B.rank >= num_contraction) {
		
		out = TensorField(A.d, A.rank + B.rank - 2 * num_contraction, B.num_points);
		size_t A_free_index = (size_t)pow(A.d, A.rank - num_contraction),
			B_free_index = (size_t)pow(B.d, B.rank - num_contraction),
			dummy_index = (size_t)pow(out.d, num_contraction);

		/* compute inner product
		result_{i_0,i_1...,i_{n-num_contraction}, j_{num_contraction},...j_m}
		= \sum _{i_{n-num_contraction + 1 ... i_n}} A_{i_0 ... i_n} B_{j_0 ... j_m}
		*/
#pragma omp parallel for schedule(guided)
		for (int x = 0; x < out.num_total_pixels; x++) {
			for (size_t a_free = 0; a_free < A_free_index; a_free++) {
				for (size_t b_free = 0; b_free < B_free_index; b_free++) {
					size_t index_of_result = a_free + A_free_index * b_free;
					//					std::cout << index_of_result <<"--\n"; 
					for (size_t dummy = 0; dummy < dummy_index; dummy++) {
						size_t index_of_A = a_free + A_free_index * dummy,
							index_of_B = dummy + dummy_index * b_free;
						//						std::cout << index_of_A << ", " << index_of_B << "\n";
						std::complex<double> result = A[index_of_A] * std::complex<double>(B[index_of_B][x][0], B[index_of_B][x][1]);
						/* real part */
#pragma omp atomic
						out[index_of_result][x][0] += result.real();
						/* imaginary part */
#pragma omp atomic
						out[index_of_result][x][1] += result.imag();
					}
					//					std::cout << "\n\n";
				}
			}

		}
	}
	else {
		std::cerr << "Inner product is not available.\n";
		out = TensorField();
	}
}

ConstantTensor InnerProduct(const ConstantTensor & A, const ConstantTensor & B, size_t num_contraction)
{
	ConstantTensor out;
	if (A.d == B.d && A.rank >= num_contraction && B.rank >= num_contraction) {

		out = ConstantTensor(A.d, A.rank + B.rank - 2 * num_contraction);
		size_t A_free_index = (size_t)pow(A.d, A.rank - num_contraction),
			B_free_index = (size_t)pow(B.d, B.rank - num_contraction),
			dummy_index = (size_t)pow(out.d, num_contraction);

		/* compute inner product
		result_{i_0,i_1...,i_{n-num_contraction}, j_{num_contraction},...j_m}
		= \sum _{i_{n-num_contraction + 1 ... i_n}} A_{i_0 ... i_n} B_{j_0 ... j_m}
		*/
#pragma omp parallel for schedule(guided)
		for (int a_free = 0; a_free < A_free_index; a_free++) {
			for (size_t b_free = 0; b_free < B_free_index; b_free++) {
				size_t index_of_result = a_free + A_free_index * b_free;
				for (size_t dummy = 0; dummy < dummy_index; dummy++) {
					size_t index_of_A = a_free + A_free_index * dummy,
						index_of_B = dummy + dummy_index * b_free;
#pragma omp critical
					{
						out[index_of_result] += (A[index_of_A] * B[index_of_B]);
					}
				}
			}
		}
	}
	else {
		std::cerr << "Inner product is not available.\n";
		out = ConstantTensor();
	}
	return out;
}


/* Function implementations of TensorField!! */

TensorField::TensorField() : field(std::vector<fftw_complex *>()),
d(0), rank(0), num_points(std::vector<size_t>(0)), num_total_pixels(0), total_comp(0)
{}

TensorField::TensorField(size_t dimension, size_t rank, const std::vector<size_t>& num_points) :
	d(dimension), rank(rank), num_points(num_points), total_comp((size_t)pow(dimension, rank))
{
	if (this->d == num_points.size()) {
		this->num_total_pixels = 1;
		for (int i = 0; i < this->d; i++)
			this->num_total_pixels *= this->num_points[i];

		this->field.resize(this->total_comp);
		for (int i = 0; i < total_comp; i++) {
			//if (this->field[i] != nullptr)
			//	fftw_free(this->field[i]);

			this->field[i] = (fftw_complex*)fftw_malloc(this->num_total_pixels* size_fftw);
			//this->field[i] = (fftw_complex*)fftw_alloc_complex(this->num_total_pixels);
			for (int j = 0; j < this->num_total_pixels; j++) {
				this->field[i][j][0] = 0.0L;	this->field[i][j][1] = 0.0L;
			}
		}
	}
	else {
		std::cout << "inconsistent dimensions\n";
		return;
	}
}

TensorField::TensorField(const TensorField & src) : TensorField(src.d, src.rank, src.num_points)
{
#pragma omp parallel for
	for (int index = 0; index < total_comp; index++) {
		for (int x = 0; x < num_total_pixels; x++) {
			this->field[index][x][0] = src.field[index][x][0];
			this->field[index][x][1] = src.field[index][x][1];
			//&(this->field[index][x]) = *(src.field[index][x]);
			//&(this->field[index][x]+1) = *(src.field[index][x]+1);
		}
	}
}

TensorField::~TensorField()
{

	if (field.size() > 0) {
		//	std::cout << "deallocate " << field[0] << "\n";
		for (int i = 0; i < this->total_comp; i++) {
			if (field[i] != nullptr)
				fftw_free(field[i]);
		}
		field.clear();
	}

	//	std::for_each(field.begin(), field.end(), fftw_free())
	//	this->field.clear();
	//	for (auto iter = this->field.begin(); iter!= this->field.end(); iter++) {
	//		fftw_free(*iter);
	//	}
}

void TensorField::Plus(TensorField & out, const TensorField & A, const TensorField & B)
{
	assert(A.num_total_pixels == B.num_total_pixels);
	assert(A.total_comp == B.total_comp);

	out = TensorField(A);
	for (int index = 0; index < out.total_comp; index++) {
		for (int x = 0; x < out.num_total_pixels; x++) {
			out.field[index][x][0] += B.field[index][x][0];
			out.field[index][x][1] += B.field[index][x][1];
		}
	}

}

TensorField & TensorField::operator+=(const TensorField & rhs)
{
	if (this->num_total_pixels == rhs.num_total_pixels && this->total_comp == rhs.total_comp) {
#pragma omp parallel for 
		for (int i = 0; i < this->total_comp; i++) {
			for (int x = 0; x < this->num_total_pixels; x++) {
#pragma omp atomic
				this->field[i][x][0] += rhs.field[i][x][0];
#pragma omp atomic
				this->field[i][x][1] += rhs.field[i][x][1];
			}
		}
	}
	else {
		std::cerr << "Wrong dimension!\n";
	}
	return *this;
}

TensorField & TensorField::operator-=(const TensorField & rhs)
{
	if (this->num_total_pixels == rhs.num_total_pixels && this->total_comp == rhs.total_comp) {
#pragma omp parallel for 
		for (int i = 0; i < this->total_comp; i++) {
			for (int x = 0; x < this->num_total_pixels; x++) {
#pragma omp atomic
				this->field[i][x][0] -= rhs.field[i][x][0];
#pragma omp atomic
				this->field[i][x][1] -= rhs.field[i][x][1];
			}
		}
	}
	else {
		std::cerr << "Wrong dimension!\n";
	}
	return *this;
}

TensorField & TensorField::operator = (const TensorField & src) {
	/* sizes are inconsistent => reallocate memory */
	if (d != src.d || rank != src.rank || num_total_pixels != src.num_total_pixels) {
		d = src.d;	rank = src.rank; num_points = std::vector<size_t>(src.num_points);
		total_comp = (size_t)pow(d, rank);		num_total_pixels = src.num_total_pixels;
		for (int i = 0; i < this->field.size(); i++)
			if (field[i] != nullptr)
				fftw_free(field[i]);

		this->field.resize(total_comp);
		for (int index = 0; index < total_comp; index++) {
			this->field[index]=(fftw_complex*)fftw_malloc(this->num_total_pixels* size_fftw);
		}

	}



#pragma omp parallel for
	for (int index = 0; index < total_comp; index++) {
		for (int x = 0; x < this->num_total_pixels; x++) {
			this->field[index][x][0] = src.field[index][x][0];
			this->field[index][x][1] = src.field[index][x][1];
		}
	}
	return *this;
}
//TensorField & TensorField::operator*=(double factor)
//{
//#pragma omp parallel for
//	for (int index = 0; index < total_comp; index++) {
//		for (int x = 0; x < this->num_total_pixels; x++) {
//			this->field[index][x][0] *= factor;
//			this->field[index][x][1] *= factor;
//		}
//	}
//	return *this;
//}

TensorField & TensorField::operator*=(const std::complex<double> rhs)
{
#pragma omp parallel for
	for (int index = 0; index < total_comp; index++) {
		for (int x = 0; x < this->num_total_pixels; x++) {
			std::complex<double> val(this->field[index][x][0], this->field[index][x][1]);
			val = val * rhs;
			this->field[index][x][0] = val.real();
			this->field[index][x][1] = val.imag();
		}
	}
	return *this;
}


void TensorField::Print(const std::string & prefix, int component, int part) const
{
	//TODO: full implementations
	size_t comp_init=0, comp_fin=0;
	size_t part_init, part_fin, part_curr = part;
	if (component < 0) {
		comp_init = 0;	comp_fin = this->total_comp;
	}
	else if (component < this->total_comp) {
		comp_init = component; comp_fin = component + 1;
	}
	else {
		std::cerr << "Print():: out of range\n";
		return;
	}
	for (size_t i = comp_init; i < comp_fin; i++) {
		std::string component, name;
		{
			std::vector<size_t> comp;	this->GetComponent(i, comp);
			for(int j=0; j<rank; j++)
				component += std::to_string(comp[j]);
			if (rank == 0)
				component = std::to_string(0);
		}

		if (part == 0)
			name = prefix + "__" + component + "_Re.txt";
		else if(part==1)
			name = prefix + "__" + component + "_Im.txt";
	
		std::fstream output(name, std::fstream::out);
		//TODO:: combine these cases....
		if (d == 1) {
			for (size_t x = 0; x < this->num_points[0]; x++) {
				output << field[i][x][part_curr] << "\t";
			}
		}
		else if (d == 2) {
			std::vector<size_t> comp({ 0,0 });
			for (size_t y = 0; y < num_points[1]; y++) {
				for (size_t x = 0; x < num_points[0]; x++) {
					comp[0] = x; comp[1] = y;
					output << field[i][y + num_points[1] * x][part_curr] << "\t";
				}
				output << "\n";
			}
		}
		else if (d == 3) {
			std::cout << "Not implemented.\n";
		}

		output.close();
	}
}

void TensorField::Print_shifted_origin(const std::string & prefix, int component, int part) const
{
	size_t comp_init = 0, comp_fin = 0;
	size_t part_init, part_fin, part_curr = part;
	if (component < 0) {
		comp_init = 0;	comp_fin = this->total_comp;
	}
	else if (component < this->total_comp) {
		comp_init = component; comp_fin = component + 1;
	}
	else {
		std::cerr << "Print():: out of range\n";
		return;
	}

	for (size_t i = comp_init; i < comp_fin; i++) {
		std::string component, name;
		{
			std::vector<size_t> comp;	this->GetComponent(i, comp);
			for (int j = 0; j < rank; j++)
				component += std::to_string(comp[j]);
			if (rank == 0)
				component = std::to_string(0);
		}

		if (part == 0)
			name = prefix + "_shifted_" + component + "_Re.txt";
		else if (part == 1)
			name = prefix + "_shifted_" + component + "_Im.txt";

		std::fstream output(name, std::fstream::out);		//TODO:: combine these cases....
		if (d == 1) {
			int x1 = 0;
			for (int x = 0; x < this->num_points[0]; x++) {
				x1 = (int)floor(this->num_points[0] / 2.0) + x;
				x1 = (x1 > this->num_points[0]) ? x1 - this->num_points[0] : x1;
				output << field[i][x1][part_curr] << "\t";
			}
		}
		else if (d == 2) {
			int x1=0, x2=0;
			std::vector<size_t> comp({ 0,0 });
			for (size_t y = 0; y < num_points[1]; y++) {
				x2 = (int)floor(this->num_points[1] / 2.0)+1 + y;
				if (x2 >= this->num_points[1])
					x2 -= this->num_points[1];

				for (size_t x = 0; x < num_points[0]; x++) {
					x1 = (int)floor(this->num_points[0] / 2.0)+1 + x;
					if(x1 >= this->num_points[0])
						x1 -= this->num_points[0];
					
					output << field[i][x2 + num_points[0] * x1][part_curr] << "\t";
				}
				output << "\n";
			}
		}
		else if (d == 3) {
			std::cout << "Not implemented.\n";
		}

		output.close();
	}
}





void InnerProduct(TensorField & result, const TensorField & A, const TensorField & B, size_t num_contraction)
{
	if (A.d == B.d && A.num_total_pixels == B.num_total_pixels
		&& A.rank >= num_contraction && B.rank >= num_contraction) {

		result = TensorField(A.d, A.rank + B.rank - 2 * num_contraction, A.num_points);
		size_t A_free_index = (size_t)pow(A.d, A.rank - num_contraction),
			B_free_index = (size_t)pow(B.d, B.rank - num_contraction),
			dummy_index = (size_t)pow(result.d, num_contraction);

		/* compute inner product
		result_{i_0,i_1...,i_{n-num_contraction}, j_{num_contraction},...j_m}
		= \sum _{i_{n-num_contraction + 1 ... i_n}} A_{i_0 ... i_n} B_{j_0 ... j_m}
		*/
#pragma omp parallel for schedule(guided)
		for (int x = 0; x < result.num_total_pixels; x++) {
			for (size_t a_free = 0; a_free < A_free_index; a_free++) {
				for (size_t b_free = 0; b_free < B_free_index; b_free++) {
					size_t index_of_result = a_free + A_free_index * b_free;
					//					std::cout << index_of_result <<"--\n"; 
					for (size_t dummy = 0; dummy < dummy_index; dummy++) {
						size_t index_of_A = a_free + A_free_index * dummy,
							index_of_B = dummy + dummy_index * b_free;
						std::complex<double> val = std::complex<double>(A.field[index_of_A][x][0], A.field[index_of_A][x][1])*std::complex<double>(B.field[index_of_B][x][0], B.field[index_of_B][x][1]);

						/* real part */
#pragma omp atomic
						result.field[index_of_result][x][0] += val.real();
						/* imaginary part */
#pragma omp atomic
						result.field[index_of_result][x][1] += val.imag();

					}
					//					std::cout << "\n\n";
				}
			}

		}
	}
	else {
		std::cerr << "Inner product is not available.\n";
		result = TensorField();
	}

}

ConstantTensor MeanField(const TensorField & src)
{
	ConstantTensor out(src.d, src.rank);
#pragma omp parallel for
	for (int i = 0; i < src.total_comp; i++) {
		std::complex<double> sum(0.0L, 0.0L);
		for (int x = 0; x < src.num_total_pixels; x++) {
			sum += std::complex<double>(src[i][x][0], src[i][x][1]);
		}
#pragma omp critical
		{
			out[i] = sum;
		}
		
	}
	
	return 1.0 / (double)src.num_total_pixels * out;
}
