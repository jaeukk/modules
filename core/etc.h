#ifndef ETC_INCLUDED
#define ETC_INCLUDED


const double LengthPrecision=1e-4;
const double EnergyPrecision=1e-4;
const double MaxDistance = 1e10;
const double MaxEnergy=1e10;
const double pi=3.1415926535897932384626433832795028841971693993751058209;
const double E = 2.71828182845904523536;
typedef unsigned short DimensionType;
const DimensionType MaxDimension = 4;
const bool UseCorrectedPotential = false;

extern bool AlsoWriteEPS;//when write bmp files, also write eps

//extern unsigned long MCThreadsPerConfiguration;
extern unsigned long MCNumParallelConfigurations;
extern unsigned long PTParallelConfigurationsPerSystem;
extern unsigned long OptimizationNumThreadsPerTrial;
extern unsigned long OptimizationNumThreads;

extern double TwoBodyDistance_MaxLength;

extern double PlotPhononFrequencies_MaxFrequencySquared_Override;

#ifdef USE_PNG
const char * const FigureFormat = ".png";
#else
const char * const FigureFormat = ".bmp";
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

/** \brief Diagnostic log stream.  Backed by a lazily-opened "DetailedLog.txt":
 *  the file is created only when something is actually written, so apps that
 *  never log leave no DetailedLog.txt behind.*/
extern std::ostream logfile;

extern size_t Verbosity;

/** \brief Empty placeholder type used as a trivial per-particle payload.  */
struct Empty
{
	/** \brief Default constructor.*/
	Empty()
	{}
};

extern time_t ProgramStart;

//implemented in : structure optimization, 
extern time_t TimeLimit;


/** \brief Coefficient C_n such that an n-ball of radius R has volume C_n*R^n. 
 *  \param n  Space dimension.
 *  \return   The dimension-dependent unit-ball volume coefficient. */
double HyperSphere_C(DimensionType n);

/** \brief Surface area of an n-dimensional sphere of radius R.
 *  \param n  Space dimension.
 *  \param R  Sphere radius.
 *  \return   Surface area of the (n-1)-sphere. */
double HyperSphere_SurfaceArea(DimensionType n, double R);

/** \brief Volume of an n-dimensional ball of radius R.
 *  \param n  Space dimension.
 *  \param R  Ball radius.
 *  \return   Volume of the n-ball. */
double HyperSphere_Volume(DimensionType n, double R);

/** \brief scaled intersection volume of two spheres.
 *  \param d  Space dimension.
 *  \param x  Reduced separation r/(2R).
 *  \return   Scaled intersection volume. */
double alpha(DimensionType d, double x);


/** \brief compute the radius of spheres from the reduced covered fraction of spheres
 *  \param fraction    Reduced covered fraction.
 *  \param NumDensity  Number density of spheres.
 *  \param d           Space dimension.
 *  \return            Sphere radius. */
double R_from_ReducedFraction(double fraction, double NumDensity, DimensionType d);

/** \brief compute the reduced fraction from the given info
 *  \param r           Sphere radius.
 *  \param NumDensity  Number density of spheres.
 *  \param d           Space dimension.
 *  \return            Reduced covered fraction. */
double ReducedFraction_from_R(double r, double NumDensity, DimensionType d);

#include <boost/math/special_functions/gamma.hpp>
/** \brief Volume of a d-dimensional superdisk/superball of deformation exponent p and half-edge L. 
 *  \param d  Space dimension (2 or 3).
 *  \param p  Deformation exponent.
 *  \param L  Half edge length.
 *  \return   Superdisk/superball volume. */
double Superdisk_vol(DimensionType d, double p, double L);

/** \brief Surface area of a d-dimensional superdisk/superball.
 *  \param d  Space dimension (2 or 3).
 *  \param p  Deformation exponent.
 *  \param L  Half edge length.
 *  \return   Superdisk/superball surface area. */
double Superdisk_surface_area(DimensionType d, double p, double L);

/** \brief Text-mode progress bar that prints tics to a stream as a count advances. */
class progress_display
{
public:
	/** \brief Construct the bar and draw its scale for the expected total count. 
	 *  \param expected_count  Total count corresponding to a full bar.
	 *  \param os              Output stream to draw on.
	 *  \param s1              Leading string before the scale.
	 *  \param s2              Leading string before the ruler.
	 *  \param s3              Leading string before the tics. */
	explicit progress_display( unsigned long expected_count,
		std::ostream & os = (Verbosity>1) ? std::cout : logfile ,
		const std::string & s1 = "\n", //leading strings
		const std::string & s2 = "",
		const std::string & s3 = "" )
		// os is hint; implementation may ignore, particularly in embedded systems
		: m_os(os), m_s1(s1), m_s2(s2), m_s3(s3) { restart(expected_count); }

	/** \brief Redraw the scale and reset the count for a new expected total. 
	 *  \param expected_count  Total count corresponding to a full bar. */
	void restart( unsigned long expected_count )
		//  Effects: display appropriate scale
		//  Postconditions: count()==0, expected_count()==expected_count
	{
		_count = _next_tic_count = _tic = 0;
		_expected_count = expected_count;

		m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
			<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
			<< std::endl  // endl implies flush, which ensures display
			<< m_s3;
		if ( !_expected_count ) _expected_count = 1;  // prevent divide by zero
	} // restart

	/** \brief Advance the count by increment, drawing tics as thresholds are crossed.
	 *  \param increment  Amount to add to the current count.
	 *  \return           The new count. */
	unsigned long  operator+=( unsigned long increment )
		//  Effects: Display appropriate progress tic if needed.
		//  Postconditions: count()== original count() + increment
		//  Returns: count().
	{
		if ( (_count += increment) >= _next_tic_count ) { display_tic(); }
		return _count;
	}

	/** \brief Pre-increment the count by one.
	 *  \return  The new count. */
	unsigned long  operator++()           { return operator+=( 1 ); }
	/** \brief Post-increment the count by one.
	 *  \return  The count before incrementing. */
	unsigned long  operator++(int)           { return (operator+=( 1 )-1); }
	/** \brief Current count.
	 *  \return  The current count. */
	unsigned long  count() const          { return _count; }
	/** \brief Expected total count.
	 *  \return  The expected count for a full bar. */
	unsigned long  expected_count() const { return _expected_count; }

private:
	std::ostream &     m_os;  // may not be present in all imps
	const std::string  m_s1;  // string is more general, safer than 
	const std::string  m_s2;  //  const char *, and efficiency or size are
	const std::string  m_s3;  //  not issues

	unsigned long _count, _expected_count, _next_tic_count;
	unsigned int  _tic;
	/** \brief Print the tics required to match the current count and schedule the next.*/
	void display_tic()
	{
		// use of floating point ensures that both large and small counts
		// work correctly.  static_cast<>() is also used several places
		// to suppress spurious compiler warnings. 
		unsigned int tics_needed =
			static_cast<unsigned int>(
			(static_cast<double>(_count)/_expected_count)*50.0 );
		do { m_os << '*' << std::flush; } while ( ++_tic < tics_needed );
		_next_tic_count = 
			static_cast<unsigned long>((_tic/50.0)*_expected_count);
		if ( _count == _expected_count ) {
			if ( _tic < 51 ) m_os << '*';
			m_os << std::endl;
		}
	} // display_tic
};

#include <vector>
/** \brief do a linear fit y=c1*x1+c2*x2+...+cn*xn
 *  \param y       input y[i] is y of the ith data point.
 *  \param x       x[j][i] is xj of the ith data point.
 *  \param pchisq  Optional output for the chi-squared of the fit.
 *  \param pR2     Optional output for the coefficient of determination.
 *  \return        a vector of <c1, c2, ..., cn>. */
std::vector<double> MultiVariableLinearFit(const std::vector<double> & y, const std::vector< std::vector<double> > & x, double * pchisq = nullptr, double * pR2 = nullptr);

/** \brief In windows, the return value of this function changes 10,000,000 every second
 *  \return  A high-resolution clock value since the epoch. */
long GetPreciseClock(void);

/** \brief Compute the intersection volume of two spheres.
 *	\param d	Space dimension
 *	\param r1	Radius of sphere1.
 *	\param r2	Radius of sphere2.
 *	\param r	Separation of centers of two spheres.
 *	\return		Intersection volume. */
inline double v2(DimensionType d, double r1, double r2, double r) {
	if (r1 > 0 && r2 > 0 && r >= 0) {
		double R2 = std::max(r1, r2);
		double R1 = std::min(r1, r2);
		double dr2 = R2*R2 - R1*R1;

		if (r <= (R2 - R1) || r<1e-15) {
			return HyperSphere_Volume(d, R1);
		}
		else if (r < sqrt(dr2)) {
			return HyperSphere_Volume(d, R1)*(1.0 - 0.5* alpha(d, (dr2 / r - r) / (2.0*R1))) + 0.5*HyperSphere_Volume(d, R2)*alpha(d, (dr2 / r + r) / (2 * R2));
		}
		else if (r < R1 + R2) {
			return 0.5*(HyperSphere_Volume(d, R1)*alpha(d, (r - dr2 / r) / (2.0*R1)) + HyperSphere_Volume(d, R2)*alpha(d, (r + dr2 / r) / (2.0*R2)));
		}
		else {
			return 0;
		}

	}
	else if(r1 ==0.0 || r2 ==0.0)	{
		return 0.0;
	}
	else {
		std::cout << "etc1::v2; error in input" << std::endl;
		return 0;
	}

}

#endif
