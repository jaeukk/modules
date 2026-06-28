#include "etc.h"
#include <cmath>

/** \brief A filebuf that opens "DetailedLog.txt" only on the first write, so
 *  applications that never log anything do not leave behind an empty
 *  DetailedLog.txt in their working directory. */
namespace {
	class LazyLogBuf : public std::filebuf
	{
		bool opened = false;
		void ensure_open()
		{
			if (!opened) {
				opened = true;	// set first: a failed open should not retry on every char
				this->open("DetailedLog.txt", std::fstream::app | std::fstream::out);
			}
		}
	protected:
		int_type overflow(int_type c = traits_type::eof()) override
		{
			ensure_open();
			return std::filebuf::overflow(c);
		}
		std::streamsize xsputn(const char * s, std::streamsize n) override
		{
			ensure_open();
			return std::filebuf::xsputn(s, n);
		}
	};
	LazyLogBuf logbuf;	// constructed before logfile (same TU, declared first)
}

std::ostream logfile(&logbuf);

//unsigned long MCThreadsPerConfiguration=0;
unsigned long MCNumParallelConfigurations=0;
unsigned long PTParallelConfigurationsPerSystem=8;
unsigned long OptimizationNumThreadsPerTrial=0;
unsigned long OptimizationNumThreads=7;

bool AlsoWriteEPS = true;//false;

double TwoBodyDistance_MaxLength=4;

double PlotPhononFrequencies_MaxFrequencySquared_Override=0.0;
//size_t Verbosity = 2;
time_t ProgramStart;
time_t TimeLimit;

#include <gsl/gsl_sf_gamma.h>

double HyperSphere_C(DimensionType n)
{
	switch (n) {
	case 1:
		return 2;
	case 2:
		return pi;
	case 3:
		return 4.1887902047863909846;
	case 4:
		return 4.9348022005446793094;
	default:
		return std::pow(::pi, n / 2.0) / gsl_sf_gamma(n / 2.0 + 1);
	}
}
double HyperSphere_SurfaceArea(DimensionType n, double R)
{
	return n*HyperSphere_C(n)*std::pow(R, static_cast<double>(n)-1.0);
}
double HyperSphere_Volume(DimensionType n, double R)
{
	return HyperSphere_C(n)*std::pow(R, static_cast<double>(n));
}

double alpha(DimensionType d, double x) { //x= r/2*R
	if (x > 1.0) {
		return 0;
	}
	else if (x < 0) {
		std::cout << "etc::alpha : sth wrong in x \n";
		return 0;
	}
	else{
		switch (d) {
		case 1:
			return 1.0 - x;
		case 2:
			return 2.0 / pi *(acos(x) - x*sqrt(1.0 - x*x));
		case 3:
			return 1.0 - 1.5*x + 0.5*x*x*x;
		case 4:
			return 2.0 / pi*(acos(x) - (5.0 *x - 2.0 *x*x*x)/3.0 *sqrt(1-x*x));
		case 5:
			return 1.0 - 1.875*x + 1.25 *x*x*x - 0.875 *x*x*x*x*x;
		default:
			std::cout << "etc::alpha : dimension out of the range - please add more generalized definition \n";
			return 0;
		}
	}
}


double R_from_ReducedFraction(double fraction, double NumDensity, DimensionType d) {
	return pow(fraction / (HyperSphere_Volume(d,1.0)*NumDensity), 1.0 / (double)d);
}
double ReducedFraction_from_R(double r, double NumDensity, DimensionType d) {
	return NumDensity*HyperSphere_Volume(d, r);
}

double Superdisk_vol(DimensionType d, double p, double L) {
	double vc = 0;
	double p_ = 1.0 / p;
	switch (d) {
	case 2:
		vc = p_*sqrt(pi) / pow(2, p_ - 1.0)*boost::math::tgamma(p_ / 2.0) / boost::math::tgamma((p_ + 1.0) / 2.0);
		break;
	case 3: //http://gisaxs.com/index.php/Form_Factor:Superball
		vc = 2.0 / 9.0 *p_ *pow(boost::math::tgamma(0.5*p_)*boost::math::tgamma(p_) / boost::math::tgamma(1.5*p_), 2);
		break;
	default:
		std::cout << "Superdisk_area : wrong dimension";
	}
	return pow(L, d)*vc;
}
double Superdisk_surface_area(DimensionType d, double p, double L) {
	return (double)d*Superdisk_vol(d,p,L) / L;
}


#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
//do a linear fit y=c1*x1+c2*x2+...+cn*xn
//input y[i] is y of the ith data point
//      x[j][i] is xj of the ith data point
//return a vector of <c1, c2, ..., cn>
std::vector<double> MultiVariableLinearFit(const std::vector<double> & y, const std::vector< std::vector<double> > & x, double * pchisq, double * pR2)
{
	size_t ndata = y.size();
	size_t nvar = x.size();
	gsl_multifit_linear_workspace * pwor = gsl_multifit_linear_alloc(ndata, nvar);
	gsl_matrix * px = gsl_matrix_alloc(ndata, nvar);
	gsl_matrix * pcov = gsl_matrix_alloc(nvar, nvar);
	gsl_vector * py = gsl_vector_alloc(ndata);
	gsl_vector * pc = gsl_vector_alloc(nvar);

	for (size_t i = 0; i < ndata; i++)
		*gsl_vector_ptr(py, i) = y[i];
	for (size_t i = 0; i < ndata; i++)
		for (size_t j = 0; j < nvar; j++)
			*gsl_matrix_ptr(px, i, j) = x[j][i];

	double chisq;
	gsl_multifit_linear(px, py, pc, pcov, &chisq, pwor);

	std::vector<double> result(nvar, 0.0);
	for (size_t i = 0; i < nvar; i++)
		result[i] = gsl_vector_get(pc, i);

	if (pchisq != nullptr)
		(*pchisq) = chisq;
	if (pR2 != nullptr)
		(*pR2) = 1.0 - chisq / gsl_stats_tss(py->data, py->stride, py->size);

	gsl_vector_free(pc);
	gsl_vector_free(py);
	gsl_matrix_free(pcov);
	gsl_matrix_free(px);
	gsl_multifit_linear_free(pwor);

	return result;
}

#include <chrono>
long GetPreciseClock(void)
{
	std::chrono::system_clock sc;
	return sc.now().time_since_epoch().count();
}

