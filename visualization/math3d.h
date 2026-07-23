// Math3d.h
// Header file for the Math3d library. The C-Runtime has math.h, this file and the
// accompanying math.c are meant to suppliment math.h by adding geometry/math routines
// useful for graphics, simulation, and physics applications (3D stuff).
// Richard S. Wright Jr.
#ifndef _MATH3D_LIBRARY__
#define _MATH3D_LIBRARY__

#include <math.h>
#include <memory.h>

///////////////////////////////////////////////////////////////////////////////
// Data structures and containers
// Much thought went into how these are declared. Many libraries declare these
// as structures with x, y, z data members. However structure alignment issues
// could limit the portability of code based on such structures, or the binary
// compatibility of data files (more likely) that contain such structures across
// compilers/platforms. Arrays are always tightly packed, and are more efficient 
// for moving blocks of data around (usually).
typedef float	M3DVector3f[3];		// Vector of three floats (x, y, z)
typedef double	M3DVector3d[3];		// Vector of three doubles (x, y, z)

typedef float M3DVector4f[4];		// Lesser used... Do we really need these?
typedef double M3DVector4d[4];		// Yes, occasionaly

typedef float M3DVector2f[2];		// 3D points = 3D Vectors, but we need a 
typedef double M3DVector2d[2];		// 2D representations sometimes... (x,y) order



// 3x3 matrix - column major. X vector is 0, 1, 2, etc.
//		0	3	6	
//		1	4	7
//		2	5	8
typedef float M3DMatrix33f[9];		// A 3 x 3 matrix, column major (floats) - OpenGL Style
typedef double M3DMatrix33d[9];		// A 3 x 3 matrix, column major (doubles) - OpenGL Style


// 4x4 matrix - column major. X vector is 0, 1, 2, etc.
//	0	4	8	12
//	1	5	9	13
//	2	6	10	14
//	3	7	11	15
typedef float M3DMatrix44f[16];		// A 4 X 4 matrix, column major (floats) - OpenGL style
typedef double M3DMatrix44d[16];	// A 4 x 4 matrix, column major (doubles) - OpenGL style


///////////////////////////////////////////////////////////////////////////////
// Useful constants
#define M3D_PI (3.14159265358979323846)
#define M3D_2PI (2.0 * M3D_PI)
#define M3D_PI_DIV_180 (0.017453292519943296)
#define M3D_INV_PI_DIV_180 (57.2957795130823229)


///////////////////////////////////////////////////////////////////////////////
// Useful shortcuts and macros
// Radians are king... but we need a way to swap back and forth
#define m3dDegToRad(x)	((x)*M3D_PI_DIV_180)
#define m3dRadToDeg(x)	((x)*M3D_INV_PI_DIV_180)

// Hour angles
#define m3dHrToDeg(x)	((x) * (1.0 / 15.0))
#define m3dHrToRad(x)	m3dDegToRad(m3dHrToDeg(x))

#define m3dDegToHr(x)	((x) * 15.0))
#define m3dRadToHr(x)	m3dDegToHr(m3dRadToDeg(x))


// Returns the same number if it is a power of
// two. Returns a larger integer if it is not a
// power of two. The larger integer is the next
// highest power of two.
/** \brief Return iValue if it is a power of two, else the next higher power of two.
 *  \param iValue  Value to round up to a power of two.
 *  \return        The smallest power of two that is >= iValue. */
inline unsigned int m3dIsPOW2(unsigned int iValue)
    {
    unsigned int nPow2 = 1;

    while(iValue > nPow2)
        {
        unsigned int next = (nPow2 << 1);
        if(next <= nPow2)   // left shift overflowed (wrapped to 0): no representable power of two >= iValue
            return 0;
        nPow2 = next;
        }

    return nPow2;
    }


///////////////////////////////////////////////////////////////////////////////
// Inline accessor functions for people who just can't count to 3 - Vectors
#define	m3dGetVectorX(v) (v[0])
#define m3dGetVectorY(v) (v[1])
#define m3dGetVectorZ(v) (v[2])
#define m3dGetVectorW(v) (v[3])

#define m3dSetVectorX(v, x)	((v)[0] = (x))
#define m3dSetVectorY(v, y)	((v)[1] = (y))
#define m3dSetVectorZ(v, z)	((v)[2] = (z))
#define m3dSetVectorW(v, w)	((v)[3] = (w))

///////////////////////////////////////////////////////////////////////////////
// Inline vector functions
// Load Vector with (x, y, z, w).
/** \brief Load a 2D vector with the components (x, y).
 *  \param v  Destination 2D vector.
 *  \param x  X component.
 *  \param y  Y component. */
inline void m3dLoadVector2(M3DVector2f v, float x, float y)
    { v[0] = x; v[1] = y; }
/** \brief Load a 2D vector with the components (x, y).
 *  \param v  Destination 2D vector.
 *  \param x  X component.
 *  \param y  Y component. */
inline void m3dLoadVector2(M3DVector2d v, float x, float y)
    { v[0] = x; v[1] = y; }
/** \brief Load a 3D vector with the components (x, y, z).
 *  \param v  Destination 3D vector.
 *  \param x  X component.
 *  \param y  Y component.
 *  \param z  Z component. */
inline void m3dLoadVector3(M3DVector3f v, float x, float y, float z)
	{ v[0] = x; v[1] = y; v[2] = z; }
/** \brief Load a 3D vector with the components (x, y, z).
 *  \param v  Destination 3D vector.
 *  \param x  X component.
 *  \param y  Y component.
 *  \param z  Z component. */
inline void m3dLoadVector3(M3DVector3d v, double x, double y, double z)
	{ v[0] = x; v[1] = y; v[2] = z; }
/** \brief Load a 4D vector with the components (x, y, z, w).
 *  \param v  Destination 4D vector.
 *  \param x  X component.
 *  \param y  Y component.
 *  \param z  Z component.
 *  \param w  W component. */
inline void m3dLoadVector4(M3DVector4f v, float x, float y, float z, float w)
	{ v[0] = x; v[1] = y; v[2] = z; v[3] = w;}
/** \brief Load a 4D vector with the components (x, y, z, w).
 *  \param v  Destination 4D vector.
 *  \param x  X component.
 *  \param y  Y component.
 *  \param z  Z component.
 *  \param w  W component. */
inline void m3dLoadVector4(M3DVector4d v, double x, double y, double z, double w)
	{ v[0] = x; v[1] = y; v[2] = z; v[3] = w;}


////////////////////////////////////////////////////////////////////////////////
// Copy vector src into vector dst
/** \brief Copy 2D vector src into vector dst.
 *  \param dst  Destination 2D vector.
 *  \param src  Source 2D vector. */
inline void	m3dCopyVector2(M3DVector2f dst, const M3DVector2f src) { memcpy(dst, src, sizeof(M3DVector2f)); }
/** \brief Copy 2D vector src into vector dst.
 *  \param dst  Destination 2D vector.
 *  \param src  Source 2D vector. */
inline void	m3dCopyVector2(M3DVector2d dst, const M3DVector2d src) { memcpy(dst, src, sizeof(M3DVector2d)); }

/** \brief Copy 3D vector src into vector dst.
 *  \param dst  Destination 3D vector.
 *  \param src  Source 3D vector. */
inline void	m3dCopyVector3(M3DVector3f dst, const M3DVector3f src) { memcpy(dst, src, sizeof(M3DVector3f)); }
/** \brief Copy 3D vector src into vector dst.
 *  \param dst  Destination 3D vector.
 *  \param src  Source 3D vector. */
inline void	m3dCopyVector3(M3DVector3d dst, const M3DVector3d src) { memcpy(dst, src, sizeof(M3DVector3d)); }

/** \brief Copy 4D vector src into vector dst.
 *  \param dst  Destination 4D vector.
 *  \param src  Source 4D vector. */
inline void	m3dCopyVector4(M3DVector4f dst, const M3DVector4f src) { memcpy(dst, src, sizeof(M3DVector4f)); }
/** \brief Copy 4D vector src into vector dst.
 *  \param dst  Destination 4D vector.
 *  \param src  Source 4D vector. */
inline void	m3dCopyVector4(M3DVector4d dst, const M3DVector4d src) { memcpy(dst, src, sizeof(M3DVector4d)); }


////////////////////////////////////////////////////////////////////////////////
// Add Vectors (r, a, b) r = a + b
/** \brief Add two 2D vectors: r = a + b.
 *  \param r  Destination 2D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors2(M3DVector2f r, const M3DVector2f a, const M3DVector2f b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1];  }
/** \brief Add two 2D vectors: r = a + b.
 *  \param r  Destination 2D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors2(M3DVector2d r, const M3DVector2d a, const M3DVector2d b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1];  }

/** \brief Add two 3D vectors: r = a + b.
 *  \param r  Destination 3D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors3(M3DVector3f r, const M3DVector3f a, const M3DVector3f b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; }
/** \brief Add two 3D vectors: r = a + b.
 *  \param r  Destination 3D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors3(M3DVector3d r, const M3DVector3d a, const M3DVector3d b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; }

/** \brief Add two 4D vectors: r = a + b.
 *  \param r  Destination 4D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors4(M3DVector4f r, const M3DVector4f a, const M3DVector4f b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; r[3] = a[3] + b[3]; }
/** \brief Add two 4D vectors: r = a + b.
 *  \param r  Destination 4D vector receiving the sum.
 *  \param a  First addend.
 *  \param b  Second addend. */
inline void m3dAddVectors4(M3DVector4d r, const M3DVector4d a, const M3DVector4d b)
	{ r[0] = a[0] + b[0];	r[1] = a[1] + b[1]; r[2] = a[2] + b[2]; r[3] = a[3] + b[3]; }

////////////////////////////////////////////////////////////////////////////////
// Subtract Vectors (r, a, b) r = a - b
/** \brief Subtract two 2D vectors: r = a - b.
 *  \param r  Destination 2D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors2(M3DVector2f r, const M3DVector2f a, const M3DVector2f b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1];  }
/** \brief Subtract two 2D vectors: r = a - b.
 *  \param r  Destination 2D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors2(M3DVector2d r, const M3DVector2d a, const M3DVector2d b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; }

/** \brief Subtract two 3D vectors: r = a - b.
 *  \param r  Destination 3D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors3(M3DVector3f r, const M3DVector3f a, const M3DVector3f b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2]; }
/** \brief Subtract two 3D vectors: r = a - b.
 *  \param r  Destination 3D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors3(M3DVector3d r, const M3DVector3d a, const M3DVector3d b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2]; }

/** \brief Subtract two 4D vectors: r = a - b.
 *  \param r  Destination 4D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors4(M3DVector4f r, const M3DVector4f a, const M3DVector4f b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2]; r[3] = a[3] - b[3]; }
/** \brief Subtract two 4D vectors: r = a - b.
 *  \param r  Destination 4D vector receiving the difference.
 *  \param a  Minuend.
 *  \param b  Subtrahend. */
inline void m3dSubtractVectors4(M3DVector4d r, const M3DVector4d a, const M3DVector4d b)
	{ r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2]; r[3] = a[3] - b[3]; }



///////////////////////////////////////////////////////////////////////////////////////
// Scale Vectors (in place)
/** \brief Scale a 2D vector in place by a scalar.
 *  \param v      2D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector2(M3DVector2f v, float scale)
	{ v[0] *= scale; v[1] *= scale; }
/** \brief Scale a 2D vector in place by a scalar.
 *  \param v      2D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector2(M3DVector2d v, double scale)
	{ v[0] *= scale; v[1] *= scale; }

/** \brief Scale a 3D vector in place by a scalar.
 *  \param v      3D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector3(M3DVector3f v, float scale)
	{ v[0] *= scale; v[1] *= scale; v[2] *= scale; }
/** \brief Scale a 3D vector in place by a scalar.
 *  \param v      3D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector3(M3DVector3d v, double scale)
	{ v[0] *= scale; v[1] *= scale; v[2] *= scale; }

/** \brief Scale a 4D vector in place by a scalar.
 *  \param v      4D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector4(M3DVector4f v, float scale)
	{ v[0] *= scale; v[1] *= scale; v[2] *= scale; v[3] *= scale; }
/** \brief Scale a 4D vector in place by a scalar.
 *  \param v      4D vector to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleVector4(M3DVector4d v, double scale)
	{ v[0] *= scale; v[1] *= scale; v[2] *= scale; v[3] *= scale; }


//////////////////////////////////////////////////////////////////////////////////////
// Cross Product
// u x v = result
// We only need one version for floats, and one version for doubles. A 3 component
// vector fits in a 4 component vector. If  M3DVector4d or M3DVector4f are passed
// we will be OK because 4th component is not used.
/** \brief Cross product: result = u x v.
 *  \param result  Destination 3D vector receiving u x v.
 *  \param u       First (left) operand.
 *  \param v       Second (right) operand. */
inline void m3dCrossProduct(M3DVector3f result, const M3DVector3f u, const M3DVector3f v)
	{
	result[0] = u[1]*v[2] - v[1]*u[2];
	result[1] = -u[0]*v[2] + v[0]*u[2];
	result[2] = u[0]*v[1] - v[0]*u[1];
	}

/** \brief Cross product: result = u x v.
 *  \param result  Destination 3D vector receiving u x v.
 *  \param u       First (left) operand.
 *  \param v       Second (right) operand. */
inline void m3dCrossProduct(M3DVector3d result, const M3DVector3d u, const M3DVector3d v)
	{
	result[0] = u[1]*v[2] - v[1]*u[2];
	result[1] = -u[0]*v[2] + v[0]*u[2];
	result[2] = u[0]*v[1] - v[0]*u[1];
	}

//////////////////////////////////////////////////////////////////////////////////////
// Dot Product, only for three component vectors
// return u dot v
/** \brief Dot product of two 3D vectors.
 *  \param u  First operand.
 *  \param v  Second operand.
 *  \return   u dot v. */
inline float m3dDotProduct(const M3DVector3f u, const M3DVector3f v)
	{ return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

/** \brief Dot product of two 3D vectors.
 *  \param u  First operand.
 *  \param v  Second operand.
 *  \return   u dot v. */
inline double m3dDotProduct(const M3DVector3d u, const M3DVector3d v)
	{ return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

//////////////////////////////////////////////////////////////////////////////////////
// Angle between vectors, only for three component vectors. Angle is in radians...
/** \brief Angle (in radians) between two 3D vectors.
 *  \param u  First vector.
 *  \param v  Second vector.
 *  \return   Angle between u and v in radians. */
inline float m3dGetAngleBetweenVectors(const M3DVector3f u, const M3DVector3f v)
    {
    double lu = sqrt(double(m3dDotProduct(u, u)));
    double lv = sqrt(double(m3dDotProduct(v, v)));
    if(lu == 0.0 || lv == 0.0)   // undefined angle for a zero-length vector
        return 0.0f;
    double c = double(m3dDotProduct(u, v)) / (lu * lv);
    if(c > 1.0) c = 1.0;         // clamp against round-off so acos never yields NaN
    else if(c < -1.0) c = -1.0;
    return float(acos(c));
    }

/** \brief Angle (in radians) between two 3D vectors.
 *  \param u  First vector.
 *  \param v  Second vector.
 *  \return   Angle between u and v in radians. */
inline double m3dGetAngleBetweenVectors(const M3DVector3d u, const M3DVector3d v)
    {
    double lu = sqrt(m3dDotProduct(u, u));
    double lv = sqrt(m3dDotProduct(v, v));
    if(lu == 0.0 || lv == 0.0)   // undefined angle for a zero-length vector
        return 0.0;
    double c = m3dDotProduct(u, v) / (lu * lv);
    if(c > 1.0) c = 1.0;         // clamp against round-off so acos never yields NaN
    else if(c < -1.0) c = -1.0;
    return acos(c);
    }

//////////////////////////////////////////////////////////////////////////////////////
// Get Square of a vectors length
// Only for three component vectors
/** \brief Squared length (magnitude) of a 3D vector.
 *  \param u  3D vector.
 *  \return   The squared length of u. */
inline float m3dGetVectorLengthSquared(const M3DVector3f u)
	{ return (u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]); }

/** \brief Squared length (magnitude) of a 3D vector.
 *  \param u  3D vector.
 *  \return   The squared length of u. */
inline double m3dGetVectorLengthSquared(const M3DVector3d u)
	{ return (u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]); }

//////////////////////////////////////////////////////////////////////////////////////
// Get lenght of vector
// Only for three component vectors.
/** \brief Length (magnitude) of a 3D vector.
 *  \param u  3D vector.
 *  \return   The length of u. */
inline float m3dGetVectorLength(const M3DVector3f u)
	{ return float(sqrt(double(m3dGetVectorLengthSquared(u)))); }

/** \brief Length (magnitude) of a 3D vector.
 *  \param u  3D vector.
 *  \return   The length of u. */
inline double m3dGetVectorLength(const M3DVector3d u)
	{ return sqrt(m3dGetVectorLengthSquared(u)); }

//////////////////////////////////////////////////////////////////////////////////////
// Normalize a vector
// Scale a vector to unit length. Easy, just scale the vector by it's length
/** \brief Normalize a 3D vector to unit length in place.
 *  \param u  3D vector to normalize in place. */
inline void m3dNormalizeVector(M3DVector3f u)
	{ m3dScaleVector3(u, 1.0f / m3dGetVectorLength(u)); }

/** \brief Normalize a 3D vector to unit length in place.
 *  \param u  3D vector to normalize in place. */
inline void m3dNormalizeVector(M3DVector3d u)
	{ m3dScaleVector3(u, 1.0 / m3dGetVectorLength(u)); }


//////////////////////////////////////////////////////////////////////////////////////
// Get the distance between two points. The distance between two points is just
// the magnitude of the difference between two vectors
// Located in math.cpp
/** \brief Squared distance between two 3D points.
 *  \param u  First point.
 *  \param v  Second point.
 *  \return   The squared distance between u and v. */
float m3dGetDistanceSquared(const M3DVector3f u, const M3DVector3f v);

/** \brief Squared distance between two 3D points.
 *  \param u  First point.
 *  \param v  Second point.
 *  \return   The squared distance between u and v. */
double m3dGetDistanceSquared(const M3DVector3d u, const M3DVector3d v);

/** \brief Distance between two 3D points.
 *  \param u  First point.
 *  \param v  Second point.
 *  \return   The distance between u and v. */
inline double m3dGetDistance(const M3DVector3d u, const M3DVector3d v)
{ return sqrt(m3dGetDistanceSquared(u, v)); }

/** \brief Distance between two 3D points.
 *  \param u  First point.
 *  \param v  Second point.
 *  \return   The distance between u and v. */
inline float m3dGetDistance(const M3DVector3f u, const M3DVector3f v)
{ return float(sqrt(m3dGetDistanceSquared(u, v))); }

/** \brief Squared magnitude (length) of a 3D vector. (CLAUDE)
 *  \param u  3D vector.
 *  \return   The squared magnitude of u. */
inline float m3dGetMagnitudeSquared(const M3DVector3f u) { return u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; }
/** \brief Squared magnitude (length) of a 3D vector. (CLAUDE)
 *  \param u  3D vector.
 *  \return   The squared magnitude of u. */
inline double m3dGetMagnitudeSquared(const M3DVector3d u) { return u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; }

/** \brief Magnitude (length) of a 3D vector. (CLAUDE)
 *  \param u  3D vector.
 *  \return   The magnitude of u. */
inline float m3dGetMagnitude(const M3DVector3f u) { return float(sqrt(m3dGetMagnitudeSquared(u))); }
/** \brief Magnitude (length) of a 3D vector. (CLAUDE)
 *  \param u  3D vector.
 *  \return   The magnitude of u. */
inline double m3dGetMagnitude(const M3DVector3d u) { return sqrt(m3dGetMagnitudeSquared(u)); }

	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrix functions
// Both floating point and double precision 3x3 and 4x4 matricies are supported.
// No support is included for arbitrarily dimensioned matricies on purpose, since
// the 3x3 and 4x4 matrix routines are the most common for the purposes of this
// library. Matrices are column major, like OpenGL matrices.
// Unlike the vector functions, some of these are going to have to not be inlined,
// although many will be.


// Copy Matrix
// Brain-dead memcpy
/** \brief Copy a 3x3 matrix src into dst.
 *  \param dst  Destination 3x3 matrix.
 *  \param src  Source 3x3 matrix. */
inline void m3dCopyMatrix33(M3DMatrix33f dst, const M3DMatrix33f src)
	{ memcpy(dst, src, sizeof(M3DMatrix33f)); }

/** \brief Copy a 3x3 matrix src into dst.
 *  \param dst  Destination 3x3 matrix.
 *  \param src  Source 3x3 matrix. */
inline void m3dCopyMatrix33(M3DMatrix33d dst, const M3DMatrix33d src)
	{ memcpy(dst, src, sizeof(M3DMatrix33d)); }

/** \brief Copy a 4x4 matrix src into dst.
 *  \param dst  Destination 4x4 matrix.
 *  \param src  Source 4x4 matrix. */
inline void m3dCopyMatrix44(M3DMatrix44f dst, const M3DMatrix44f src)
	{ memcpy(dst, src, sizeof(M3DMatrix44f)); }

/** \brief Copy a 4x4 matrix src into dst.
 *  \param dst  Destination 4x4 matrix.
 *  \param src  Source 4x4 matrix. */
inline void m3dCopyMatrix44(M3DMatrix44d dst, const M3DMatrix44d src)
	{ memcpy(dst, src, sizeof(M3DMatrix44d)); }

// LoadIdentity
// Implemented in Math3d.cpp
/** \brief Load a 3x3 matrix with the identity matrix.
 *  \param m  3x3 matrix to overwrite with the identity. */
void m3dLoadIdentity33(M3DMatrix33f m);

/** \brief Load a 3x3 matrix with the identity matrix.
 *  \param m  3x3 matrix to overwrite with the identity. */
void m3dLoadIdentity33(M3DMatrix33d m);

/** \brief Load a 4x4 matrix with the identity matrix.
 *  \param m  4x4 matrix to overwrite with the identity. */
void m3dLoadIdentity44(M3DMatrix44f m);

/** \brief Load a 4x4 matrix with the identity matrix.
 *  \param m  4x4 matrix to overwrite with the identity. */
void m3dLoadIdentity44(M3DMatrix44d m);

/////////////////////////////////////////////////////////////////////////////
// Get/Set Column.
/** \brief Read one column of a 3x3 matrix into a 3D vector.
 *  \param dst     Destination 3D vector receiving the column.
 *  \param src     Source 3x3 matrix.
 *  \param column  Index of the column to read (0-2). */
inline void m3dGetMatrixColumn33(M3DVector3f dst, const M3DMatrix33f src, int column)
	{ memcpy(dst, src + (3 * column), sizeof(float) * 3); }

/** \brief Read one column of a 3x3 matrix into a 3D vector.
 *  \param dst     Destination 3D vector receiving the column.
 *  \param src     Source 3x3 matrix.
 *  \param column  Index of the column to read (0-2). */
inline void m3dGetMatrixColumn33(M3DVector3d dst, const M3DMatrix33d src, int column)
	{ memcpy(dst, src + (3 * column), sizeof(double) * 3); }

/** \brief Write a 3D vector into one column of a 3x3 matrix.
 *  \param dst     Destination 3x3 matrix.
 *  \param src     Source 3D vector to write into the column.
 *  \param column  Index of the column to write (0-2). */
inline void m3dSetMatrixColumn33(M3DMatrix33f dst, const M3DVector3f src, int column)
	{ memcpy(dst + (3 * column), src, sizeof(float) * 3); }

/** \brief Write a 3D vector into one column of a 3x3 matrix.
 *  \param dst     Destination 3x3 matrix.
 *  \param src     Source 3D vector to write into the column.
 *  \param column  Index of the column to write (0-2). */
inline void m3dSetMatrixColumn33(M3DMatrix33d dst, const M3DVector3d src, int column)
	{ memcpy(dst + (3 * column), src, sizeof(double) * 3); }

/** \brief Read one column of a 4x4 matrix into a 4D vector.
 *  \param dst     Destination 4D vector receiving the column.
 *  \param src     Source 4x4 matrix.
 *  \param column  Index of the column to read (0-3). */
inline void m3dGetMatrixColumn44(M3DVector4f dst, const M3DMatrix44f src, int column)
	{ memcpy(dst, src + (4 * column), sizeof(float) * 4); }

/** \brief Read one column of a 4x4 matrix into a 4D vector.
 *  \param dst     Destination 4D vector receiving the column.
 *  \param src     Source 4x4 matrix.
 *  \param column  Index of the column to read (0-3). */
inline void m3dGetMatrixColumn44(M3DVector4d dst, const M3DMatrix44d src, int column)
	{ memcpy(dst, src + (4 * column), sizeof(double) * 4); }

/** \brief Write a 4D vector into one column of a 4x4 matrix.
 *  \param dst     Destination 4x4 matrix.
 *  \param src     Source 4D vector to write into the column.
 *  \param column  Index of the column to write (0-3). */
inline void m3dSetMatrixColumn44(M3DMatrix44f dst, const M3DVector4f src, int column)
	{ memcpy(dst + (4 * column), src, sizeof(float) * 4); }

/** \brief Write a 4D vector into one column of a 4x4 matrix.
 *  \param dst     Destination 4x4 matrix.
 *  \param src     Source 4D vector to write into the column.
 *  \param column  Index of the column to write (0-3). */
inline void m3dSetMatrixColumn44(M3DMatrix44d dst, const M3DVector4d src, int column)
	{ memcpy(dst + (4 * column), src, sizeof(double) * 4); }


// Get/Set row purposely omitted... use the functions below.
// I don't think row vectors are useful for column major ordering...
// If I'm wrong, add them later.


//////////////////////////////////////////////////////////////////////////////
// Get/Set RowCol - Remember column major ordering...
// Provides for element addressing
/** \brief Set the (row, col) element of a 3x3 matrix.
 *  \param m      3x3 matrix to modify.
 *  \param row    Row index (0-2).
 *  \param col    Column index (0-2).
 *  \param value  Value to store. */
inline void m3dSetMatrixRowCol33(M3DMatrix33f m, int row, int col, float value)
	{ m[(col * 3) + row] = value; }

/** \brief Get the (row, col) element of a 3x3 matrix.
 *  \param m    3x3 matrix to read.
 *  \param row  Row index (0-2).
 *  \param col  Column index (0-2).
 *  \return     The element at (row, col). */
inline float m3dGetMatrixRowCol33(const M3DMatrix33f m, int row, int col)
	{ return m[(col * 3) + row]; }

/** \brief Set the (row, col) element of a 3x3 matrix.
 *  \param m      3x3 matrix to modify.
 *  \param row    Row index (0-2).
 *  \param col    Column index (0-2).
 *  \param value  Value to store. */
inline void m3dSetMatrixRowCol33(M3DMatrix33d m, int row, int col, double value)
	{ m[(col * 3) + row] = value; }

/** \brief Get the (row, col) element of a 3x3 matrix.
 *  \param m    3x3 matrix to read.
 *  \param row  Row index (0-2).
 *  \param col  Column index (0-2).
 *  \return     The element at (row, col). */
inline double m3dGetMatrixRowCol33(const M3DMatrix33d m, int row, int col)
	{ return m[(col * 3) + row]; }

/** \brief Set the (row, col) element of a 4x4 matrix.
 *  \param m      4x4 matrix to modify.
 *  \param row    Row index (0-3).
 *  \param col    Column index (0-3).
 *  \param value  Value to store. */
inline void m3dSetMatrixRowCol44(M3DMatrix44f m, int row, int col, float value)
	{ m[(col * 4) + row] = value; }

/** \brief Get the (row, col) element of a 4x4 matrix.
 *  \param m    4x4 matrix to read.
 *  \param row  Row index (0-3).
 *  \param col  Column index (0-3).
 *  \return     The element at (row, col). */
inline float m3dGetMatrixRowCol44(const M3DMatrix44f m, int row, int col)
	{ return m[(col * 4) + row]; }

/** \brief Set the (row, col) element of a 4x4 matrix.
 *  \param m      4x4 matrix to modify.
 *  \param row    Row index (0-3).
 *  \param col    Column index (0-3).
 *  \param value  Value to store. */
inline void m3dSetMatrixRowCol44(M3DMatrix44d m, int row, int col, double value)
	{ m[(col * 4) + row] = value; }

/** \brief Get the (row, col) element of a 4x4 matrix.
 *  \param m    4x4 matrix to read.
 *  \param row  Row index (0-3).
 *  \param col  Column index (0-3).
 *  \return     The element at (row, col). */
inline double m3dGetMatrixRowCol44(const M3DMatrix44d m, int row, int col)
	{ return m[(col * 4) + row]; }


///////////////////////////////////////////////////////////////////////////////
// Extract a rotation matrix from a 4x4 matrix
// Extracts the rotation matrix (3x3) from a 4x4 matrix
/** \brief Extract the 3x3 rotation matrix from a 4x4 matrix.
 *  \param dst  Destination 3x3 rotation matrix.
 *  \param src  Source 4x4 matrix. */
inline void m3dExtractRotation(M3DMatrix33f dst, const M3DMatrix44f src)
	{
	memcpy(dst, src, sizeof(float) * 3); // X column
	memcpy(dst + 3, src + 4, sizeof(float) * 3); // Y column
	memcpy(dst + 6, src + 8, sizeof(float) * 3); // Z column
	}

// Ditto above, but for doubles
/** \brief Extract the 3x3 rotation matrix from a 4x4 matrix.
 *  \param dst  Destination 3x3 rotation matrix.
 *  \param src  Source 4x4 matrix. */
inline void m3dExtractRotation(M3DMatrix33d dst, const M3DMatrix44d src)
	{
	memcpy(dst, src, sizeof(double) * 3); // X column
	memcpy(dst + 3, src + 4, sizeof(double) * 3); // Y column
	memcpy(dst + 6, src + 8, sizeof(double) * 3); // Z column
	}

// Inject Rotation (3x3) into a full 4x4 matrix...
/** \brief Inject a 3x3 rotation matrix into the upper-left of a 4x4 matrix.
 *  \param dst  Destination 4x4 matrix.
 *  \param src  Source 3x3 rotation matrix. */
inline void m3dInjectRotation(M3DMatrix44f dst, const M3DMatrix33f src)
	{
	memcpy(dst, src, sizeof(float) * 4);
	memcpy(dst + 4, src + 4, sizeof(float) * 4);
	memcpy(dst + 8, src + 8, sizeof(float) * 4);
	}

// Ditto above for doubles
/** \brief Inject a 3x3 rotation matrix into the upper-left of a 4x4 matrix.
 *  \param dst  Destination 4x4 matrix.
 *  \param src  Source 3x3 rotation matrix. */
inline void m3dInjectRotation(M3DMatrix44d dst, const M3DMatrix33d src)
	{
	memcpy(dst, src, sizeof(double) * 4);
	memcpy(dst + 4, src + 4, sizeof(double) * 4);
	memcpy(dst + 8, src + 8, sizeof(double) * 4);
	}


////////////////////////////////////////////////////////////////////////////////
// MultMatrix
// Implemented in Math.cpp
/** \brief Multiply two 4x4 matrices: product = a * b. (CLAUDE)
 *  \param product  Destination 4x4 matrix receiving the product.
 *  \param a        Left 4x4 matrix.
 *  \param b        Right 4x4 matrix. */
void m3dMatrixMultiply44(M3DMatrix44f product, const M3DMatrix44f a, const M3DMatrix44f b);

/** \brief Multiply two 4x4 matrices: product = a * b. (CLAUDE)
 *  \param product  Destination 4x4 matrix receiving the product.
 *  \param a        Left 4x4 matrix.
 *  \param b        Right 4x4 matrix. */
void m3dMatrixMultiply44(M3DMatrix44d product, const M3DMatrix44d a, const M3DMatrix44d b);

/** \brief Multiply two 3x3 matrices: product = a * b. (CLAUDE)
 *  \param product  Destination 3x3 matrix receiving the product.
 *  \param a        Left 3x3 matrix.
 *  \param b        Right 3x3 matrix. */
void m3dMatrixMultiply33(M3DMatrix33f product, const M3DMatrix33f a, const M3DMatrix33f b);

/** \brief Multiply two 3x3 matrices: product = a * b. (CLAUDE)
 *  \param product  Destination 3x3 matrix receiving the product.
 *  \param a        Left 3x3 matrix.
 *  \param b        Right 3x3 matrix. */
void m3dMatrixMultiply33(M3DMatrix33d product, const M3DMatrix33d a, const M3DMatrix33d b);


// Transform - Does rotation and translation via a 4x4 matrix. Transforms
// a point or vector.
// By-the-way __inline means I'm asking the compiler to do a cost/benefit analysis. If 
// these are used frequently, they may not be inlined to save memory. I'm experimenting
// with this....
/** \brief Transform a 3D point/vector by a 4x4 matrix (rotation and translation). (CLAUDE)
 *  \param vOut  Destination 3D vector receiving the transformed result.
 *  \param v     Source 3D point or vector.
 *  \param m     4x4 transformation matrix. */
__inline void m3dTransformVector3(M3DVector3f vOut, const M3DVector3f v, const M3DMatrix44f m)
    {
    vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12];// * v[3];
    vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13];// * v[3];
    vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14];// * v[3];
	//vOut[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }

// Ditto above, but for doubles
/** \brief Transform a 3D point/vector by a 4x4 matrix (rotation and translation). (CLAUDE)
 *  \param vOut  Destination 3D vector receiving the transformed result.
 *  \param v     Source 3D point or vector.
 *  \param m     4x4 transformation matrix. */
__inline void m3dTransformVector3(M3DVector3d vOut, const M3DVector3d v, const M3DMatrix44d m)
    {
    vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12];// * v[3];	 
    vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13];// * v[3];	
    vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14];// * v[3];	
	//vOut[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }

/** \brief Transform a 4D point/vector by a 4x4 matrix (full homogeneous transform). (CLAUDE)
 *  \param vOut  Destination 4D vector receiving the transformed result.
 *  \param v     Source 4D point or vector.
 *  \param m     4x4 transformation matrix. */
__inline void m3dTransformVector4(M3DVector4f vOut, const M3DVector4f v, const M3DMatrix44f m)
    {
    vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12] * v[3];
    vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13] * v[3];	
    vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3];	
	vOut[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }

// Ditto above, but for doubles
/** \brief Transform a 4D point/vector by a 4x4 matrix (full homogeneous transform). (CLAUDE)
 *  \param vOut  Destination 4D vector receiving the transformed result.
 *  \param v     Source 4D point or vector.
 *  \param m     4x4 transformation matrix. */
__inline void m3dTransformVector4(M3DVector4d vOut, const M3DVector4d v, const M3DMatrix44d m)
    {
    vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12] * v[3];	 
    vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13] * v[3];	
    vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3];	
	vOut[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }



// Just do the rotation, not the translation... this is usually done with a 3x3
// Matrix.
/** \brief Rotate a 3D vector by a 3x3 matrix (rotation only, no translation).
 *  \param vOut  Destination 3D vector receiving the rotated result.
 *  \param p     Source 3D vector.
 *  \param m     3x3 rotation matrix. */
__inline void m3dRotateVector(M3DVector3f vOut, const M3DVector3f p, const M3DMatrix33f m)
	{
    vOut[0] = m[0] * p[0] + m[3] * p[1] + m[6] * p[2];	
    vOut[1] = m[1] * p[0] + m[4] * p[1] + m[7] * p[2];	
    vOut[2] = m[2] * p[0] + m[5] * p[1] + m[8] * p[2];	
	}

// Ditto above, but for doubles
/** \brief Rotate a 3D vector by a 3x3 matrix (rotation only, no translation).
 *  \param vOut  Destination 3D vector receiving the rotated result.
 *  \param p     Source 3D vector.
 *  \param m     3x3 rotation matrix. */
__inline void m3dRotateVector(M3DVector3d vOut, const M3DVector3d p, const M3DMatrix33d m)
	{
    vOut[0] = m[0] * p[0] + m[3] * p[1] + m[6] * p[2];	
    vOut[1] = m[1] * p[0] + m[4] * p[1] + m[7] * p[2];	
    vOut[2] = m[2] * p[0] + m[5] * p[1] + m[8] * p[2];	
	}


// Scale a matrix (I don't beleive in Scaling matricies ;-)
// Yes, it's faster to loop backwards... These could be 
// unrolled... but eh... if you find this is a bottleneck,
// then you should unroll it yourself
/** \brief Scale all elements of a 3x3 matrix in place by a scalar. (CLAUDE)
 *  \param m      3x3 matrix to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleMatrix33(M3DMatrix33f m, float scale)
{ for(int i = 8; i >=0; i--) m[i] *= scale; }

/** \brief Scale all elements of a 3x3 matrix in place by a scalar. (CLAUDE)
 *  \param m      3x3 matrix to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleMatrix33(M3DMatrix33d m, double scale)
{ for(int i = 8; i >=0; i--) m[i] *= scale; }

/** \brief Scale all elements of a 4x4 matrix in place by a scalar. (CLAUDE)
 *  \param m      4x4 matrix to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleMatrix44(M3DMatrix44f m, float scale)
{ for(int i = 15; i >=0; i--) m[i] *= scale; }

/** \brief Scale all elements of a 4x4 matrix in place by a scalar. (CLAUDE)
 *  \param m      4x4 matrix to scale in place.
 *  \param scale  Scalar multiplier. */
inline void m3dScaleMatrix44(M3DMatrix44d m, double scale)
{ for(int i = 15; i >=0; i--) m[i] *= scale; }


// Create a Rotation matrix
// Implemented in math.cpp
/** \brief Build a 3x3 rotation matrix for a rotation by angle (radians) about axis (x, y, z). (CLAUDE)
 *  \param m      Destination 3x3 matrix.
 *  \param angle  Rotation angle in radians.
 *  \param x      X component of the rotation axis.
 *  \param y      Y component of the rotation axis.
 *  \param z      Z component of the rotation axis. */
void m3dRotationMatrix33(M3DMatrix33f m, float angle, float x, float y, float z);

/** \brief Build a 3x3 rotation matrix for a rotation by angle (radians) about axis (x, y, z). (CLAUDE)
 *  \param m      Destination 3x3 matrix.
 *  \param angle  Rotation angle in radians.
 *  \param x      X component of the rotation axis.
 *  \param y      Y component of the rotation axis.
 *  \param z      Z component of the rotation axis. */
void m3dRotationMatrix33(M3DMatrix33d m, double angle, double x, double y, double z);

/** \brief Build a 4x4 rotation matrix for a rotation by angle (radians) about axis (x, y, z). (CLAUDE)
 *  \param m      Destination 4x4 matrix.
 *  \param angle  Rotation angle in radians.
 *  \param x      X component of the rotation axis.
 *  \param y      Y component of the rotation axis.
 *  \param z      Z component of the rotation axis. */
void m3dRotationMatrix44(M3DMatrix44f m, float angle, float x, float y, float z);

/** \brief Build a 4x4 rotation matrix for a rotation by angle (radians) about axis (x, y, z). (CLAUDE)
 *  \param m      Destination 4x4 matrix.
 *  \param angle  Rotation angle in radians.
 *  \param x      X component of the rotation axis.
 *  \param y      Y component of the rotation axis.
 *  \param z      Z component of the rotation axis. */
void m3dRotationMatrix44(M3DMatrix44d m, double angle, double x, double y, double z);

// Create a Translation matrix. Only 4x4 matrices have translation components
/** \brief Build a 4x4 translation matrix for the offset (x, y, z).
 *  \param m  Destination 4x4 matrix.
 *  \param x  Translation along X.
 *  \param y  Translation along Y.
 *  \param z  Translation along Z. */
inline void m3dTranslationMatrix44(M3DMatrix44f m, float x, float y, float z)
{ m3dLoadIdentity44(m); m[12] = x; m[13] = y; m[14] = z; }

/** \brief Build a 4x4 translation matrix for the offset (x, y, z).
 *  \param m  Destination 4x4 matrix.
 *  \param x  Translation along X.
 *  \param y  Translation along Y.
 *  \param z  Translation along Z. */
inline void m3dTranslationMatrix44(M3DMatrix44d m, double x, double y, double z)
{ m3dLoadIdentity44(m); m[12] = x; m[13] = y; m[14] = z; }


// Translate matrix. Only 4x4 matrices supported
/** \brief Add a translation (x, y, z) to a 4x4 matrix in place.
 *  \param m  4x4 matrix to modify in place.
 *  \param x  Translation to add along X.
 *  \param y  Translation to add along Y.
 *  \param z  Translation to add along Z. */
inline void m3dTranslateMatrix44(M3DMatrix44f m, float x, float y, float z)
{ m[12] += x; m[13] += y; m[14] += z; }

/** \brief Add a translation (x, y, z) to a 4x4 matrix in place.
 *  \param m  4x4 matrix to modify in place.
 *  \param x  Translation to add along X.
 *  \param y  Translation to add along Y.
 *  \param z  Translation to add along Z. */
inline void m3dTranslateMatrix44(M3DMatrix44d m, double x, double y, double z)
{ m[12] += x; m[13] += y; m[14] += z; }


// Scale matrix. Only 4x4 matrices supported
/** \brief Apply a non-uniform scale (x, y, z) to a 4x4 matrix in place.
 *  \param m  4x4 matrix to modify in place.
 *  \param x  Scale factor along X.
 *  \param y  Scale factor along Y.
 *  \param z  Scale factor along Z. */
inline void m3dScaleMatrix44(M3DMatrix44f m, float x, float y, float z)
{ m[0] *= x; m[5] *= y; m[10] *= z; }

/** \brief Apply a non-uniform scale (x, y, z) to a 4x4 matrix in place.
 *  \param m  4x4 matrix to modify in place.
 *  \param x  Scale factor along X.
 *  \param y  Scale factor along Y.
 *  \param z  Scale factor along Z. */
inline void m3dScaleMatrix44(M3DMatrix44d m, double x, double y, double z)
{ m[0] *= x; m[5] *= y; m[10] *= z; }


// Transpose/Invert - Only 4x4 matricies supported
#define TRANSPOSE44(dst, src)            \
{                                        \
    for (int j = 0; j < 4; j++)          \
    {                                    \
        for (int i = 0; i < 4; i++)      \
        {                                \
            dst[(j*4)+i] = src[(i*4)+j]; \
        }                                \
    }                                    \
}

/** \brief Transpose a 4x4 matrix: dst = transpose(src). (CLAUDE)
 *  \param dst  Destination 4x4 matrix receiving the transpose.
 *  \param src  Source 4x4 matrix. */
inline void m3dTransposeMatrix44(M3DMatrix44f dst, const M3DMatrix44f src)
{ TRANSPOSE44(dst, src); }
/** \brief Transpose a 4x4 matrix: dst = transpose(src). (CLAUDE)
 *  \param dst  Destination 4x4 matrix receiving the transpose.
 *  \param src  Source 4x4 matrix. */
inline void m3dTransposeMatrix44(M3DMatrix44d dst, const M3DMatrix44d src)
{ TRANSPOSE44(dst, src); }
/** \brief Invert a 4x4 matrix: dst = inverse(src). (CLAUDE)
 *  \param dst  Destination 4x4 matrix receiving the inverse.
 *  \param src  Source 4x4 matrix.
 *  \return     True on success, false if src is singular. */
bool m3dInvertMatrix44(M3DMatrix44f dst, const M3DMatrix44f src);

/** \brief Invert a 4x4 matrix: dst = inverse(src). (CLAUDE)
 *  \param dst  Destination 4x4 matrix receiving the inverse.
 *  \param src  Source 4x4 matrix.
 *  \return     True on success, false if src is singular. */
bool m3dInvertMatrix44(M3DMatrix44d dst, const M3DMatrix44d src);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Other Miscellaneous functions

// Find a normal from three points
// Implemented in math3d.cpp
/** \brief Find the normal vector of a triangle defined by three points (counter-clockwise winding).
 *  \param result  Destination 3D vector receiving the normal.
 *  \param point1  First triangle vertex.
 *  \param point2  Second triangle vertex.
 *  \param point3  Third triangle vertex. */
void m3dFindNormal(M3DVector3f result, const M3DVector3f point1, const M3DVector3f point2,
							const M3DVector3f point3);

/** \brief Find the normal vector of a triangle defined by three points (counter-clockwise winding).
 *  \param result  Destination 3D vector receiving the normal.
 *  \param point1  First triangle vertex.
 *  \param point2  Second triangle vertex.
 *  \param point3  Third triangle vertex. */
void m3dFindNormal(M3DVector3d result, const M3DVector3d point1, const M3DVector3d point2,
							const M3DVector3d point3);



// Calculates the signed distance of a point to a plane
/** \brief Signed distance from a point to a plane.
 *  \param point  3D point to test.
 *  \param plane  Plane equation (A, B, C, D) where Ax+By+Cz+D = 0.
 *  \return       The signed distance from the point to the plane. */
inline float m3dGetDistanceToPlane(const M3DVector3f point, const M3DVector4f plane)
           { return point[0]*plane[0] + point[1]*plane[1] + point[2]*plane[2] + plane[3]; }

/** \brief Signed distance from a point to a plane.
 *  \param point  3D point to test.
 *  \param plane  Plane equation (A, B, C, D) where Ax+By+Cz+D = 0.
 *  \return       The signed distance from the point to the plane. */
inline double m3dGetDistanceToPlane(const M3DVector3d point, const M3DVector4d plane)
           { return point[0]*plane[0] + point[1]*plane[1] + point[2]*plane[2] + plane[3]; }


// Get plane equation from three points and a normal
/** \brief Compute the plane equation coefficients (A, B, C, D) from three points.
 *  \param planeEq  Destination 4D vector receiving the plane coefficients (A, B, C, D).
 *  \param p1       First point on the plane.
 *  \param p2       Second point on the plane.
 *  \param p3       Third point on the plane. */
void m3dGetPlaneEquation(M3DVector4f planeEq, const M3DVector3f p1, const M3DVector3f p2, const M3DVector3f p3);

/** \brief Compute the plane equation coefficients (A, B, C, D) from three points.
 *  \param planeEq  Destination 4D vector receiving the plane coefficients (A, B, C, D).
 *  \param p1       First point on the plane.
 *  \param p2       Second point on the plane.
 *  \param p3       Third point on the plane. */
void m3dGetPlaneEquation(M3DVector4d planeEq, const M3DVector3d p1, const M3DVector3d p2, const M3DVector3d p3);

// Determine if a ray intersects a sphere
/** \brief Test whether a ray intersects a sphere and return the distance to the hit.
 *  \param point         Origin point of the ray.
 *  \param ray           Direction of the ray.
 *  \param sphereCenter  Center of the sphere.
 *  \param sphereRadius  Radius of the sphere.
 *  \return              Distance to the intersection; < 0 if no intersection, 0 if tangent. */
double m3dRaySphereTest(const M3DVector3d point, const M3DVector3d ray, const M3DVector3d sphereCenter, double sphereRadius);

/** \brief Test whether a ray intersects a sphere and return the distance to the hit.
 *  \param point         Origin point of the ray.
 *  \param ray           Direction of the ray.
 *  \param sphereCenter  Center of the sphere.
 *  \param sphereRadius  Radius of the sphere.
 *  \return              Distance to the intersection; < 0 if no intersection, 0 if tangent. */
float m3dRaySphereTest(const M3DVector3f point, const M3DVector3f ray, const M3DVector3f sphereCenter, float sphereRadius);

// Etc. etc.

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Faster (and more robust) replacements for gluProject
/** \brief Project a 3D point to 2D window coordinates (discarding Z), like gluProject. (CLAUDE)
 *  \param vPointOut    Destination 2D window coordinates.
 *  \param mModelView   Model-view matrix.
 *  \param mProjection  Projection matrix.
 *  \param iViewPort    Viewport as (x, y, width, height).
 *  \param vPointIn     Source 3D point in object space. */
void m3dProjectXY( M3DVector2f vPointOut, const M3DMatrix44f mModelView, const M3DMatrix44f mProjection, const int iViewPort[4], const M3DVector3f vPointIn);

/** \brief Project a 3D point to 3D window coordinates (keeping Z), like gluProject. (CLAUDE)
 *  \param vPointOut    Destination 3D window coordinates.
 *  \param mModelView   Model-view matrix.
 *  \param mProjection  Projection matrix.
 *  \param iViewPort    Viewport as (x, y, width, height).
 *  \param vPointIn     Source 3D point in object space. */
void m3dProjectXYZ(M3DVector3f vPointOut, const M3DMatrix44f mModelView, const M3DMatrix44f mProjection, const int iViewPort[4], const M3DVector3f vPointIn);



//////////////////////////////////////////////////////////////////////////////////////////////////
// This function does a three dimensional Catmull-Rom "spline" interpolation between p1 and p2
/** \brief Catmull-Rom spline interpolation between vP1 and vP2 in 3D. (CLAUDE)
 *  \param vOut  Destination 3D vector receiving the interpolated point.
 *  \param vP0   Control point before the segment.
 *  \param vP1   Start point of the interpolated segment.
 *  \param vP2   End point of the interpolated segment.
 *  \param vP3   Control point after the segment.
 *  \param t     Interpolation parameter in [0, 1]. */
void m3dCatmullRom(M3DVector3f vOut, M3DVector3f vP0, M3DVector3f vP1, M3DVector3f vP2, M3DVector3f vP3, float t);

/** \brief Catmull-Rom spline interpolation between vP1 and vP2 in 3D. (CLAUDE)
 *  \param vOut  Destination 3D vector receiving the interpolated point.
 *  \param vP0   Control point before the segment.
 *  \param vP1   Start point of the interpolated segment.
 *  \param vP2   End point of the interpolated segment.
 *  \param vP3   Control point after the segment.
 *  \param t     Interpolation parameter in [0, 1]. */
void m3dCatmullRom(M3DVector3d vOut, M3DVector3d vP0, M3DVector3d vP1, M3DVector3d vP2, M3DVector3d vP3, double t);

//////////////////////////////////////////////////////////////////////////////////////////////////
// Compare floats and doubles...
/** \brief Test whether two floats are within a given epsilon of each other.
 *  \param fCandidate  First value.
 *  \param fCompare    Second value.
 *  \param fEpsilon    Maximum allowed absolute difference.
 *  \return            True if |fCandidate - fCompare| < fEpsilon. */
inline bool m3dCloseEnough(float fCandidate, float fCompare, float fEpsilon)
    {
    return (fabs(fCandidate - fCompare) < fEpsilon);
    }

/** \brief Test whether two doubles are within a given epsilon of each other.
 *  \param dCandidate  First value.
 *  \param dCompare    Second value.
 *  \param dEpsilon    Maximum allowed absolute difference.
 *  \return            True if |dCandidate - dCompare| < dEpsilon. */
inline bool m3dCloseEnough(double dCandidate, double dCompare, double dEpsilon)
    {
    return (fabs(dCandidate - dCompare) < dEpsilon);
    }

////////////////////////////////////////////////////////////////////////////
// Used for normal mapping. Finds the tangent bases for a triangle...
// Only a floating point implementation is provided.
/** \brief Calculate the tangent basis vector for a triangle (used in normal mapping).
 *  \param pvTriangle   Three triangle vertices.
 *  \param pvTexCoords  Three texture coordinates corresponding to the vertices.
 *  \param N            Surface normal of the triangle.
 *  \param vTangent     Destination 3D vector receiving the tangent basis. */
void m3dCalculateTangentBasis(const M3DVector3f pvTriangle[3], const M3DVector2f pvTexCoords[3], const M3DVector3f N, M3DVector3f vTangent);

////////////////////////////////////////////////////////////////////////////
// Smoothly step between 0 and 1 between edge1 and edge 2
/** \brief Smooth Hermite interpolation of x from 0 to 1 between edge1 and edge2.
 *  \param edge1  Lower edge of the interpolation interval.
 *  \param edge2  Upper edge of the interpolation interval.
 *  \param x      Input value.
 *  \return       Smoothly interpolated value in [0, 1]. */
double m3dSmoothStep(double edge1, double edge2, double x);

/** \brief Smooth Hermite interpolation of x from 0 to 1 between edge1 and edge2.
 *  \param edge1  Lower edge of the interpolation interval.
 *  \param edge2  Upper edge of the interpolation interval.
 *  \param x      Input value.
 *  \return       Smoothly interpolated value in [0, 1]. */
float m3dSmoothStep(float edge1, float edge2, float x);

/////////////////////////////////////////////////////////////////////////////
// Planar shadow Matrix
/** \brief Build a 4x4 projection matrix that flattens geometry onto a plane to cast a shadow.
 *  \param proj      Destination 4x4 shadow projection matrix.
 *  \param planeEq   Plane equation (A, B, C, D) of the shadow-receiving plane.
 *  \param vLightPos Position of the light source. */
void m3dMakePlanarShadowMatrix(M3DMatrix44d proj, const M3DVector4d planeEq, const M3DVector3d vLightPos);

/** \brief Build a 4x4 projection matrix that flattens geometry onto a plane to cast a shadow.
 *  \param proj      Destination 4x4 shadow projection matrix.
 *  \param planeEq   Plane equation (A, B, C, D) of the shadow-receiving plane.
 *  \param vLightPos Position of the light source. */
void m3dMakePlanarShadowMatrix(M3DMatrix44f proj, const M3DVector4f planeEq, const M3DVector3f vLightPos);

/** \brief Find the point on a ray closest to a given point in space.
 *  \param vPointOnRay   Destination 3D vector receiving the closest point on the ray.
 *  \param vRayOrigin    Origin of the ray.
 *  \param vUnitRayDir   Unit direction vector of the ray.
 *  \param vPointInSpace The point in space.
 *  \return              The squared distance between the ray and the point. */
double m3dClosestPointOnRay(M3DVector3d vPointOnRay, const M3DVector3d vRayOrigin, const M3DVector3d vUnitRayDir,
							const M3DVector3d vPointInSpace);

/** \brief Find the point on a ray closest to a given point in space.
 *  \param vPointOnRay   Destination 3D vector receiving the closest point on the ray.
 *  \param vRayOrigin    Origin of the ray.
 *  \param vUnitRayDir   Unit direction vector of the ray.
 *  \param vPointInSpace The point in space.
 *  \return              The squared distance between the ray and the point. */
float m3dClosestPointOnRay(M3DVector3f vPointOnRay, const M3DVector3f vRayOrigin, const M3DVector3f vUnitRayDir,
							const M3DVector3f vPointInSpace);

#endif