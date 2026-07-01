/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Created	:	2015? 
 *	Comments: Jaeuk Kim (November 2025).
 */

/** \file GeometryVector.h
 *	\brief Header file to define vectors in Euclidean space and their operatios. 
 *	Some type-related warnings are resolved. */


#ifndef GEOMETRYVECTOR_INCLUDED
#define GEOMETRYVECTOR_INCLUDED

#include "etc.h"

#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>
#include "RandomGenerator.h"

#ifndef NDEBUG
#define GEOMETRYVECTOR_RECORDDIMENSION
#endif


//In debug mode, this class also records the dimension of the GeometryVector.
/** \brief Fixed-size vector in Euclidean space supporting basic arithmetic. */
class GeometryVector
{
public:
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	DimensionType Dimension;
#endif

	/** \brief Record the current dimension (debug builds only).
	 *  \param d  Dimension to record. */
	void SetDimension(DimensionType d)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=d;
#endif
	}

	double x[ ::MaxDimension];

	/** \brief Default-construct a zero vector of unspecified dimension. */
	GeometryVector()
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=0;
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	//do not change this "int" to "DimensionType"!
	//otherwize GeometryVector(2) would be ambiguous.

	/** \brief Construct a zero vector with an explicit dimension.
	 *  \param dimension  Dimension to assign to the vector. */
	GeometryVector(int dimension)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=dimension;
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0.0L;
	}

	/** \brief Construct a 1D vector.
	 *  \param xx  First (x) coordinate. */
	GeometryVector(double xx)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=1;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		for(int i=1; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 2D vector.
	 *  \param xx  First (x) coordinate.
	 *  \param yy  Second (y) coordinate. */
	GeometryVector(double xx, double yy)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=2;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		for(int i=2; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 3D vector.
	 *  \param xx  First (x) coordinate.
	 *  \param yy  Second (y) coordinate.
	 *  \param zz  Third (z) coordinate. */
	GeometryVector(double xx, double yy, double zz)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=3;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		this->x[2]=zz;
		for(int i=3; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 4D vector.
	 *  \param xx  First (x) coordinate.
	 *  \param yy  Second (y) coordinate.
	 *  \param zz  Third (z) coordinate.
	 *  \param aa  Fourth coordinate. */
	GeometryVector(double xx, double yy, double zz, double aa)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=4;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		this->x[2]=zz;
		this->x[3]=aa;
		for(int i=4; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Copy-construct from another vector.
	 *  \param src  Source vector to copy from. */
	GeometryVector(const GeometryVector & src)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(src.Dimension<= ::MaxDimension);
		this->Dimension=src.Dimension;
#endif
		//std::memcpy(this->x, src.x, ::MaxDimension*sizeof(double));
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=src.x[i];
	}

	/** \brief Add components from another vector in-place.
	 *  \param right  Vector whose components are added to this one. */
	void AddFrom(const GeometryVector & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]+=right.x[i];
	}

	/** \brief Subtract components from another vector in-place.
	 *  \param right  Vector whose components are subtracted from this one. */
	void MinusFrom(const GeometryVector & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]-=right.x[i];
	}

	/** \brief Scale all components by a constant.
	 *  \param right  Scalar multiplier applied to every component. */
	void MultiplyFrom(const double & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(this->Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]*=right;
	}

	/** \brief Inner product with another vector.
	 *  \param right  The other operand of the inner product.
	 *  \return       The scalar inner product. */
	double Dot(const GeometryVector & right) const		// equivalent to the inner product
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		double result=0;
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		for(int i=0; i<this->Dimension; i++)
#else
		for(int i=0; i< ::MaxDimension; i++)
#endif
				result+=this->x[i]*right.x[i];
		return result;
	}
	//cross product

	/** \brief 3D cross product.
	 *  \param right  The right-hand operand of the cross product.
	 *  \return       The 3D vector orthogonal to both operands. */
	GeometryVector Cross(const GeometryVector & right) const
	{
		GeometryVector result((int)3);
	//#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension == 3);
		assert(this->Dimension == 3);
	//#endif

		result.x[0] = this->x[1] * right.x[2] - this->x[2] * right.x[1];
		result.x[1] = this->x[2] * right.x[0] - this->x[0] * right.x[2];
		result.x[2] = this->x[0] * right.x[1] - this->x[1] * right.x[0];

		return result;
	}


	/** \brief Squared Euclidean norm. */
	double Modulus2(void) const
	{
		return this->Dot(*this);
	}

	/** \brief Exact component-wise equality test.
	 *  \param right  Vector to compare against.
	 *  \return       True if all compared components are exactly equal. */
	bool IsEqual(const GeometryVector & right) const
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
		for(int i=0; i<right.Dimension; i++)
			if(this->x[i]!=right.x[i])
				return false;
#else
		for(int i=0; i< ::MaxDimension; i++)
			if(this->x[i]!=right.x[i])
				return false;
#endif
		return true;
	}

	/** \brief Stream out the first dim coordinates separated by tabs.
	 *  \param os   Output stream to write to.
	 *  \param dim  Number of leading coordinates to write. */
	void OutputCoordinate(std::ostream & os, DimensionType dim) const
	{
		for(int i=0; i< dim; i++)
			os<<this->x[i]<<" \t";
	}

	/** \brief Read dim coordinates from a stream.
	 *  \param is   Input stream to read from.
	 *  \param dim  Number of leading coordinates to read. */
	void InputCoordinate(std::istream & is, DimensionType dim)
	{
		for(int i=0; i< dim; i++)
			is>>this->x[i];
	}

	/** \brief Write coordinates in binary form.
	 *  \param ofile      Output stream to write the binary data to.
	 *  \param dimension  Number of coordinates to write. */
	void WriteBinary(std::ostream & ofile, DimensionType dimension) const;

	/** \brief Read coordinates from binary form.
	 *  \param ifile      Input stream to read the binary data from.
	 *  \param dimension  Number of coordinates to read. */
	void ReadBinary(std::istream & ifile, DimensionType dimension);

	friend inline GeometryVector SameDimensionZeroVector(const GeometryVector & src);
};
/** \brief Create another GeometryVector with same dimension as src, but the value is zero.
 *  \param src  Vector whose dimension is matched.
 *  \return     A zero vector of the same dimension as src. */
//create another GeometryVector with same dimension as src, but the value is zero
inline GeometryVector SameDimensionZeroVector(const GeometryVector & src)
{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	return GeometryVector(src.Dimension);
#else
	return GeometryVector(0);
#endif
}



/** \brief Stream out all MaxDimension coordinates followed by a newline.  
 *  \param os  Output stream to write to.
 *  \param a   Vector to write.
 *  \return    The output stream. */
inline std::ostream & operator << (std::ostream & os, const GeometryVector & a)
{
	//os<<a.Dimension<<" \t";
	a.OutputCoordinate(os, ::MaxDimension);
	os<<'\n';
	return os;
}

/** \brief Read all MaxDimension coordinates from a stream.
 *  \param is  Input stream to read from.
 *  \param a   Vector to receive the read coordinates.
 *  \return    The input stream. */
inline std::istream & operator >> (std::istream & is, GeometryVector & a)
{
	//is>>a.Dimension;
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	assert(a.Dimension<= ::MaxDimension);
#endif
	a.InputCoordinate(is, ::MaxDimension);
	return is;
}

/** \brief Component-wise vector addition.
 *  \param left   Left-hand operand.
 *  \param right  Right-hand operand.
 *  \return       The component-wise sum. */
inline GeometryVector operator + (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.AddFrom(right);
	return result;
}

/** \brief Component-wise vector subtraction.
 *  \param left   Left-hand operand (minuend).
 *  \param right  Right-hand operand (subtrahend).
 *  \return       The component-wise difference. */
inline GeometryVector operator - (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.MinusFrom(right);
	return result;
}

/** \brief Scale a vector by a scalar on the right.
 *  \param left   Vector to scale.
 *  \param right  Scalar multiplier.
 *  \return       The scaled vector. */
inline GeometryVector operator * (const GeometryVector & left, const double & right)
{
	GeometryVector result(left);
	result.MultiplyFrom(right);
	return result;
}

/** \brief Scale a vector by a scalar on the left.
 *  \param left   Scalar multiplier.
 *  \param right  Vector to scale.
 *  \return       The scaled vector. */
inline GeometryVector operator * (const double & left, const GeometryVector & right)
{
	GeometryVector result(right);
	result.MultiplyFrom(left);
	return result;
}

/** \brief Exact component-wise equality comparison.
 *  \param left   Left-hand operand.
 *  \param right  Right-hand operand.
 *  \return       True if the vectors are exactly equal. */
inline bool operator == (const GeometryVector & left, const GeometryVector & right)
{
	return left.IsEqual(right);
}

/** \brief Exact component-wise inequality comparison.
 *  \param left   Left-hand operand.
 *  \param right  Right-hand operand.
 *  \return       True if the vectors differ in any component. */
inline bool operator != (const GeometryVector & left, const GeometryVector & right)
{
	return !left.IsEqual(right);
}

// RandomVector(...) -> renamed RandomVectorInCube(...) and moved to RandomGenerator.h
// (declared there; reachable here because GeometryVector.h includes it).
// The old generator-less 2-arg overload was retired (it default-seeded a fresh
// RandomGenerator on every call, so consecutive calls replayed the same draw).

// RandomUnitVector(...) and GetRandomVector(...) moved to RandomGenerator.h
// (declared there; reachable here because GeometryVector.h includes it).
/** \brief Volume of the parallelepiped spanned by dimension basis vectors.
 *  \param vecs       Array of basis vectors (length dimension).
 *  \param dimension  Number of vectors and space dimension.
 *  \return           Absolute value of the determinant, or 0 if dimension is 0. */
double Volume (const GeometryVector * vecs, DimensionType dimension);

/** \brief Volume of the simplex spanned by d basis vectors (parallelepiped volume / d!).
 *  \param vecs  Array of basis vectors (length d).
 *  \param d     Number of vectors and space dimension.
 *  \return      The simplex volume. */
 double SimplexVolume (const GeometryVector * vecs, DimensionType d);

 /** \brief Rotate basis vectors into lower-triangular form with descending lengths.
 *  Sorts so that l1>=l2>=l3..., then QR-rotates so the 1st basis vector is
 *  (c11,0,0,...), the 2nd is (c21,c22,0,...), etc. (volume is preserved).
 *  \param bas  Basis vectors to standardize in place (length d).
 *  \param d    Space dimension. */
void StandardizeOrientation(GeometryVector * bas, DimensionType d);

/** \brief Write a list of vectors to a stream, with optional comment header fields.
 *  \param result  Vectors to write, one per line (all MaxDimension coordinates).
 *  \param ofile   Output stream.
 *  \param fields  Optional field names printed as a leading '#'-commented header. */
void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields = std::vector<std::string>());

/** \brief Write a list of vectors to file OutFilePrefix.txt, with optional header fields.
 *  \param result         Vectors to write.
 *  \param OutFilePrefix  Output file path prefix (".txt" is appended).
 *  \param fields         Optional field names printed as a leading header.*/
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields = std::vector<std::string>());

/** \brief Read a list of vectors from a stream, ignoring '#'-commented lines.
 *  \param result      Output container, cleared then filled with the read vectors.
 *  \param ifile       Input stream.
 *  \param NumColumns  Number of columns (coordinates) to read per line.*/
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);

/** \brief Read a list of vectors from file InFilePrefix.txt, ignoring commented lines.
 *  \param result        Output container, cleared then filled with the read vectors.
 *  \param InFilePrefix  Input file path prefix (".txt" is appended).
 *  \param NumColumns    Number of columns (coordinates) to read per line. */
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);


//another implementation of vector
//the difference between GeometryVector and BigVector is that BigVectors can have dimensions exceeding MaxDimension
//BigVector can be used to represent vectors in configuration space
/** \brief Heap-allocated vector of arbitrary dimension (may exceed MaxDimension).*/
class BigVector
{
public:
	double * x;
	size_t dim;
	/** \brief Default-construct an empty (zero-dimensional) vector. */
	BigVector()
	{
		this->dim=0;
		this->x=nullptr;
	}

	/** \brief Construct a zero vector of the given dimension.
	 *  \param dimension  Number of components to allocate. */
	BigVector(size_t dimension)
	{
		this->dim=dimension;
		this->x=new double[this->dim];
		for(size_t i=0; i<dimension; i++)
			this->x[i]=0.0;
	}

	/** \brief Copy-construct from another vector (deep copy).
	 *  \param src  Source vector to copy from. */
	BigVector(const BigVector & src)
	{
		this->dim=src.dim;
		this->x=new double[this->dim];
		std::memcpy(this->x, src.x, sizeof(double)*this->dim);
	}

	/** \brief Move-construct from another vector, leaving the source empty.
	 *  \param src  Source vector whose storage is taken over (left empty). */
	BigVector(BigVector && src)
	{
		this->dim=src.dim;
		this->x=src.x;
		src.dim=0;
		src.x=nullptr;
	}

	/** \brief Destructor; releases the heap-allocated storage. */
	~BigVector()
	{
		if(this->x!=nullptr)
			delete [] this->x;
	}

	/** \brief Add components from another vector in-place.
	 *  \param right  Vector whose components are added to this one. */
	void AddFrom(const BigVector & right)
	{
		assert(right.dim==this->dim);
		for(size_t i=0; i<this->dim; i++)
			this->x[i]+=right.x[i];
	}

	/** \brief Subtract components from another vector in-place.
	 *  \param right  Vector whose components are subtracted from this one. */
	void MinusFrom(const BigVector & right)
	{
		assert(right.dim==this->dim);
		for(size_t i=0; i<this->dim; i++)
			this->x[i]-=right.x[i];
	}

	/** \brief Scale all components by a constant.
	 *  \param right  Scalar multiplier applied to every component. */
	void MultiplyFrom(const double & right)
	{
		for(size_t i=0; i<this->dim; i++)
			this->x[i]*=right;
	}

	/** \brief Inner product with another vector.
	 *  \param src  The other operand of the inner product.
	 *  \return     The scalar inner product. */
	double Dot(const BigVector & src) const
	{
		assert(src.dim==this->dim);
		double result=0.0;
		for(size_t i=0; i<this->dim; i++)
			result+=this->x[i]*src.x[i];
		return result;
	}

	/** \brief Squared Euclidean norm.*/
	double Modulus2(void) const
	{
		return this->Dot(*this);
	}

	/** \brief Stream out all dim coordinates separated by tabs.
	 *  \param os  Output stream to write to. */
	void OutputCoordinate(std::ostream & os) const
	{
		for(size_t i=0; i< dim; i++)
			os<<this->x[i]<<" \t";
	}

	/** \brief Read dim coordinates from a stream.
	 *  \param is  Input stream to read from. */
	void InputCoordinate(std::istream & is)
	{
		for(size_t i=0; i< dim; i++)
			is>>this->x[i];
	}
};
/** \brief Component-wise addition of two BigVectors.
 *  \param left   Left-hand operand.
 *  \param right  Right-hand operand.
 *  \return       The component-wise sum. */
inline BigVector operator + (const BigVector & left, const BigVector & right)
{
	BigVector result(left);
	result.AddFrom(right);
	return result;
}

/** \brief Component-wise subtraction of two BigVectors.
 *  \param left   Left-hand operand (minuend).
 *  \param right  Right-hand operand (subtrahend).
 *  \return       The component-wise difference. */
inline BigVector operator - (const BigVector & left, const BigVector & right)
{
	BigVector result(left);
	result.MinusFrom(right);
	return result;
}

/** \brief Scale a BigVector by a scalar on the right. 
 *  \param left   Vector to scale.
 *  \param right  Scalar multiplier.
 *  \return       The scaled vector. */
inline BigVector operator * (const BigVector & left, const double & right)
{
	BigVector result(left);
	result.MultiplyFrom(right);
	return result;
}

/** \brief Scale a BigVector by a scalar on the left.
 *  \param left   Scalar multiplier.
 *  \param right  Vector to scale.
 *  \return       The scaled vector. */
inline BigVector operator * (const double & left, const BigVector & right)
{
	return  right*left;
}


#endif
