/** \file test_finufft_type1.cpp
 *  \brief Self-contained checks of the optional FINUFFT Type 1 backend of ShiftedCCPotential::GetRho.
 *
 *  Built three ways by run_tests.sh (see README.md):
 *   - direct-only: no FINUFFT macro, headers, or library; FinufftType1 must be unavailable;
 *   - -DCCPOTENTIAL_USE_FINUFFT: every lattice case must really run FINUFFT and match Direct;
 *   - -DCCPOTENTIAL_USE_FINUFFT -DHAVETWOFIXEDPARTICLES: same, with the +2 real offset.
 *  Exits with status 0 only if every check passes. */

#include "CollectiveCoordinatePotential.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
typedef ShiftedCCPotential::RhoBackend Backend;
typedef std::vector<long> Mode;

const unsigned long TestThreads=2;	//task-local cap on OpenMP and FINUFFT threads
const bool HaveFinufft=ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1);
#ifdef HAVETWOFIXEDPARTICLES
const double FixedParticles=2.0;
#else
const double FixedParticles=0.0;
#endif
size_t FinufftEvaluations=0;	//number of compared evaluations that really used FINUFFT

/** \brief Counts checks and failures and prints each failing check. */
struct Checker
{
	size_t checks=0;
	size_t failures=0;

	/** \brief Record one boolean check.
	 *  \param ok    Whether the check passed.
	 *  \param what  Description printed on failure. */
	void Check(bool ok, const std::string & what)
	{
		checks++;
		if(!ok)
		{
			failures++;
			std::cout<<"FAIL: "<<what<<'\n';
		}
	}

	/** \brief Check |a-b| <= atol + rtol*scale (and that both are finite).
	 *  \param a      First value.
	 *  \param b      Second value.
	 *  \param scale  Magnitude the relative tolerance is scaled by.
	 *  \param rtol   Relative tolerance.
	 *  \param atol   Absolute tolerance.
	 *  \param what   Description printed on failure.
	 *  \param index  Element index printed on failure.
	 *  \return |a-b|. */
	double Close(double a, double b, double scale, double rtol, double atol, const std::string & what, size_t index=0)
	{
		double err=std::abs(a-b);
		bool ok=std::isfinite(a) && std::isfinite(b) && err<=atol+rtol*scale;
		if(ok)
			this->Check(true, what);
		else
		{
			std::ostringstream os;
			os.precision(17);
			os<<what<<" [index "<<index<<"]: a="<<a<<" b="<<b<<" |a-b|="<<err<<" allowed="<<atol+rtol*scale;
			this->Check(false, os.str());
		}
		return err;
	}
};

/** \brief Maximum absolute deviations found by CompareBackends. */
struct Errors
{
	double rho=0.0;
	double energy=0.0;
	double force=0.0;
};

/** \brief Build a configuration with uniformly random fractional coordinates.
 *  \param dim    Spatial dimension.
 *  \param rows   Basis vectors, one row per vector.
 *  \param n      Number of particles.
 *  \param rng    Random engine (advanced by this call).
 *  \return The configuration. */
Configuration MakeConfig(DimensionType dim, const std::vector<std::vector<double> > & rows, size_t n, std::mt19937_64 & rng)
{
	std::vector<GeometryVector> basis(dim, GeometryVector(static_cast<int>(dim)));
	double maxLength=0.0;
	for(DimensionType d=0; d<dim; d++)
	{
		for(DimensionType c=0; c<dim; c++)
			basis[d].x[c]=rows[d][c];
		maxLength=std::max(maxLength, std::sqrt(basis[d].Modulus2()));
	}
	Configuration config(dim, basis.data(), 10.0*maxLength);
	std::uniform_real_distribution<double> u(0.0, 1.0);
	for(size_t j=0; j<n; j++)
	{
		GeometryVector s(static_cast<int>(dim));
		for(DimensionType d=0; d<dim; d++)
			s.x[d]=u(rng);
		config.Insert("A", s);
	}
	return config;
}

/** \brief Cartesian reciprocal-lattice vector k = sum_d m_d b_d of a configuration.
 *  \param config  Configuration providing the reciprocal basis.
 *  \param m       Integer indices.
 *  \return The wavevector. */
GeometryVector LatticeK(const Configuration & config, const Mode & m)
{
	GeometryVector k(static_cast<int>(config.GetDimension()));
	for(DimensionType d=0; d<config.GetDimension(); d++)
		k.AddFrom(static_cast<double>(m[d])*config.GetReciprocalBasisVector(d));
	return k;
}

/** \brief All integer modes in [-K, K]^dim (includes zero, negative and positive modes).
 *  \param dim  Spatial dimension.
 *  \param K    Largest absolute index.
 *  \return The modes in lexicographic order. */
std::vector<Mode> AllModes(DimensionType dim, long K)
{
	std::vector<Mode> result;
	Mode m(dim, -K);
	while(true)
	{
		result.push_back(m);
		DimensionType d=0;
		for(; d<dim; d++)
		{
			if(m[d]<K)
			{
				m[d]++;
				break;
			}
			m[d]=-K;
		}
		if(d==dim)
			break;
	}
	return result;
}

/** \brief Shuffle modes and insert duplicated entries at random positions.
 *  \param modes       Modes to scramble.
 *  \param duplicates  Number of duplicate entries to insert.
 *  \param rng         Random engine.
 *  \return The scrambled mode list. */
std::vector<Mode> Scramble(std::vector<Mode> modes, size_t duplicates, std::mt19937_64 & rng)
{
	std::shuffle(modes.begin(), modes.end(), rng);
	for(size_t t=0; t<duplicates && !modes.empty(); t++)
	{
		Mode copy=modes[std::uniform_int_distribution<size_t>(0, modes.size()-1)(rng)];
		modes.insert(modes.begin()+std::uniform_int_distribution<size_t>(0, modes.size())(rng), copy);
	}
	return modes;
}

/** \brief Build a direct-backend potential with random weights V in [vmin, vmax].
 *  \param config  Configuration providing the reciprocal basis.
 *  \param modes   Integer modes of the constraints.
 *  \param rng     Random engine.
 *  \param vmin    Smallest weight.
 *  \param vmax    Largest weight.
 *  \return The potential (unbound). */
ShiftedCCPotential MakePotential(const Configuration & config, const std::vector<Mode> & modes, std::mt19937_64 & rng, double vmin=-0.5, double vmax=1.5)
{
	ShiftedCCPotential pot(config.GetDimension());
	pot.ParallelNumThread=TestThreads;
	std::uniform_real_distribution<double> u(vmin, vmax);
	for(const Mode & m : modes)
		pot.constraints.emplace_back(LatticeK(config, m), u(rng));
	return pot;
}

/** \brief Magnitudes used to scale relative tolerances. */
struct Scales
{
	double rho, energy, force;
};

/** \brief Tolerance scales: rho ~ N, energy ~ sum|V| N^2/vol, force ~ 4 sum|V| N kmax/vol.
 *  \param pot     Potential with constraints.
 *  \param config  Bound configuration.
 *  \return The scales. */
Scales MakeScales(const ShiftedCCPotential & pot, const Configuration & config)
{
	double n=config.NumParticle()+FixedParticles;
	double sumV=0.0, kmax=0.0;
	for(const auto & c : pot.constraints)
	{
		sumV+=std::abs(c.V);
		kmax=std::max(kmax, std::sqrt(c.k.Modulus2()));
	}
	double vol=config.PeriodicVolume();
	return Scales{std::max(n, 1.0), sumV*n*n/vol, 4.0*sumV*n*kmax/vol};
}

/** \brief Compare rho, Energy, Force and AllForce of a Direct copy and a FINUFFT-selected copy of proto.
 *  When the FINUFFT copy falls back, the results must be bitwise identical to Direct.
 *  \param ck             Checker.
 *  \param proto          Potential providing constraints and Shift.
 *  \param config         Configuration to bind.
 *  \param label          Case name.
 *  \param expectFinufft  Whether FINUFFT must actually be used (only enforced if compiled in).
 *  \param eps            FINUFFT tolerance.
 *  \param maxModes       Mode cap for the FINUFFT copy (0 keeps the default).
 *  \param rtol           Relative tolerance for FINUFFT results.
 *  \return Maximum deviations. */
Errors CompareBackends(Checker & ck, const ShiftedCCPotential & proto, const Configuration & config, const std::string & label, bool expectFinufft, double eps=1e-12, std::int64_t maxModes=0, double rtol=1e-9)
{
	Errors e;
	ShiftedCCPotential direct(proto), fast(proto);
	direct.ParallelNumThread=TestThreads;
	fast.ParallelNumThread=TestThreads;
	direct.SetRhoBackend(Backend::Direct);
	if(HaveFinufft)
		fast.SetRhoBackend(Backend::FinufftType1);
	fast.SetFinufftTolerance(eps);
	if(maxModes>0)
		fast.SetFinufftMaxModes(maxModes);
	direct.SetConfiguration(config);
	fast.SetConfiguration(config);
	direct.GetRho();
	fast.GetRho();

	const bool usedFinufft=fast.GetLastRhoBackendUsed()==Backend::FinufftType1;
	ck.Check(direct.GetLastRhoBackendUsed()==Backend::Direct, label+": reference used Direct");
	if(expectFinufft && HaveFinufft)
	{
		ck.Check(usedFinufft, label+": FINUFFT was actually used");
		ck.Check(fast.GetLastFinufftStatus()==0, label+": FINUFFT status 0 (got "+std::to_string(fast.GetLastFinufftStatus())+")");
	}
	else
	{
		ck.Check(!usedFinufft, label+": fell back to Direct");
		ck.Check(fast.GetLastFinufftStatus()==ShiftedCCPotential::FinufftNotCalled, label+": rejected before calling FINUFFT (status "+std::to_string(fast.GetLastFinufftStatus())+")");
	}
	if(usedFinufft)
		FinufftEvaluations++;

	const double r=usedFinufft ? rtol : 0.0;
	const double a=usedFinufft ? 1e-10 : 0.0;
	const Scales s=MakeScales(proto, config);
	const size_t nk=proto.constraints.size();
	ck.Check(fast.RhoReal.size()==nk && fast.RhoImag.size()==nk, label+": rho sizes");
	ck.Check(fast.RhoValid, label+": RhoValid after GetRho");
	for(size_t i=0; i<nk && i<fast.RhoReal.size(); i++)
	{
		e.rho=std::max(e.rho, ck.Close(fast.RhoReal[i], direct.RhoReal[i], s.rho, r, a, label+": Re rho", i));
		e.rho=std::max(e.rho, ck.Close(fast.RhoImag[i], direct.RhoImag[i], s.rho, r, a, label+": Im rho", i));
	}
	e.energy=ck.Close(fast.Energy(), direct.Energy(), s.energy, r, a, label+": Energy");

	std::vector<GeometryVector> fd, ff;
	direct.AllForce(fd);
	fast.AllForce(ff);
	ck.Check(fd.size()==config.NumParticle() && ff.size()==fd.size(), label+": AllForce sizes");
	for(size_t j=0; j<ff.size() && j<fd.size(); j++)
		for(DimensionType d=0; d<config.GetDimension(); d++)
			e.force=std::max(e.force, ck.Close(ff[j].x[d], fd[j].x[d], s.force, r, a, label+": AllForce", j));
	//the FINUFFT rho backend also selects Type 2 forces; its fallback matches Direct bitwise (checked above)
	const bool type2=fast.GetLastForceBackendUsed()==ShiftedCCPotential::ForceBackend::FinufftType2;
	ck.Check(type2==usedFinufft, label+": AllForce used Type 2 exactly when rho used FINUFFT");

	//single-particle Force is always direct: bitwise equal to a direct AllForce, within tolerance of Type 2
	const size_t np=config.NumParticle();
	std::vector<size_t> picks;
	if(np>0)
		picks={0, np/2, np-1};
	for(size_t j : picks)
	{
		GeometryVector f;
		fast.Force(f, j);
		for(DimensionType d=0; d<config.GetDimension(); d++)
			ck.Close(f.x[d], ff[j].x[d], s.force, type2 ? r : 0.0, type2 ? a : 0.0, label+": Force vs AllForce", j);
	}
	return e;
}

/** \brief Default selection, availability, setter validation, and disabled-backend behavior. */
void TestBackendSelection(Checker & ck)
{
	std::mt19937_64 rng(101);
	Configuration config=MakeConfig(1, {{2.0}}, 8, rng);
	ShiftedCCPotential pot=MakePotential(config, AllModes(1, 3), rng);
	pot.SetConfiguration(config);

	ck.Check(ShiftedCCPotential::IsRhoBackendAvailable(Backend::Direct), "Direct is always available");
#ifdef CCPOTENTIAL_USE_FINUFFT
	ck.Check(HaveFinufft, "FinufftType1 available with CCPOTENTIAL_USE_FINUFFT");
#else
	ck.Check(!HaveFinufft, "FinufftType1 unavailable without CCPOTENTIAL_USE_FINUFFT");
#endif
	ck.Check(pot.GetRhoBackend()==Backend::Direct, "default backend is Direct");
	ck.Check(pot.GetFinufftTolerance()==1e-12, "default tolerance 1e-12");
	ck.Check(pot.GetFinufftMaxModes()==(std::int64_t(1)<<22), "default max modes 2^22");
	ck.Check(pot.GetLastFinufftStatus()==ShiftedCCPotential::FinufftNotCalled, "status FinufftNotCalled before GetRho");
	pot.GetRho();
	ck.Check(pot.RhoValid && pot.GetLastRhoBackendUsed()==Backend::Direct, "default GetRho uses Direct");

	bool threw=false;
	try
	{
		pot.SetRhoBackend(Backend::FinufftType1);
	}
	catch(const std::invalid_argument &)
	{
		threw=true;
	}
	if(HaveFinufft)
	{
		ck.Check(!threw, "selecting FinufftType1 succeeds when compiled in");
		ck.Check(!pot.RhoValid, "selecting a backend invalidates rho");
		ck.Check(pot.GetRhoBackend()==Backend::FinufftType1, "FinufftType1 selected");
	}
	else
	{
		ck.Check(threw, "selecting disabled FinufftType1 throws invalid_argument");
		ck.Check(pot.GetRhoBackend()==Backend::Direct, "disabled selection leaves Direct");
		ck.Check(pot.RhoValid, "disabled selection leaves the cache untouched");
		pot.GetRho();
		ck.Check(pot.GetLastRhoBackendUsed()==Backend::Direct, "disabled build evaluates with Direct");
	}

	pot.GetRho();
	const Backend before=pot.GetRhoBackend();
	threw=false;
	try
	{
		pot.SetRhoBackend(static_cast<Backend>(42));
	}
	catch(const std::invalid_argument &)
	{
		threw=true;
	}
	ck.Check(threw && pot.GetRhoBackend()==before && pot.RhoValid, "unknown backend value throws and changes nothing");

	const double bad[]={0.0, -1e-9, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()};
	for(double eps : bad)
	{
		threw=false;
		try
		{
			pot.SetFinufftTolerance(eps);
		}
		catch(const std::invalid_argument &)
		{
			threw=true;
		}
		ck.Check(threw && pot.GetFinufftTolerance()==1e-12 && pot.RhoValid, "invalid tolerance throws and changes nothing");
	}
	threw=false;
	try
	{
		pot.SetFinufftMaxModes(0);
	}
	catch(const std::invalid_argument &)
	{
		threw=true;
	}
	ck.Check(threw && pot.GetFinufftMaxModes()==(std::int64_t(1)<<22) && pot.RhoValid, "invalid max modes throws and changes nothing");

	pot.SetFinufftTolerance(1e-10);
	ck.Check(!pot.RhoValid && pot.GetFinufftTolerance()==1e-10, "tolerance change invalidates rho");
	pot.GetRho();
	pot.SetFinufftMaxModes(1000);
	ck.Check(!pot.RhoValid && pot.GetFinufftMaxModes()==1000, "max-modes change invalidates rho");
}

/** \brief Analytic sign and phase checks: rho(k) = sum_j exp(+i k.r_j). */
void TestSign(Checker & ck)
{
	std::vector<Backend> backends{Backend::Direct};
	if(HaveFinufft)
		backends.push_back(Backend::FinufftType1);
	for(Backend b : backends)
	{
		const std::string name=(b==Backend::Direct) ? "Direct" : "FinufftType1";
		{
			GeometryVector basis[1]={GeometryVector(static_cast<int>(1))};
			basis[0].x[0]=1.0;
			Configuration config(1, basis, 10.0);
			GeometryVector s(static_cast<int>(1));
			s.x[0]=0.25;
			config.Insert("A", s);
			ShiftedCCPotential pot(1);
			pot.constraints.emplace_back(LatticeK(config, Mode{1}), 1.0);
			pot.constraints.emplace_back(LatticeK(config, Mode{-1}), 1.0);
			pot.constraints.emplace_back(LatticeK(config, Mode{0}), 1.0);
			pot.SetRhoBackend(b);
			pot.SetConfiguration(config);
			pot.GetRho();
			ck.Check(pot.GetLastRhoBackendUsed()==b, "sign 1D: "+name+" was used");
			ck.Close(pot.RhoReal[0], 0.0+FixedParticles, 1.0, 0.0, 1e-10, "sign 1D m=+1 Re, "+name);
			ck.Close(pot.RhoImag[0], 1.0, 1.0, 0.0, 1e-10, "sign 1D m=+1 Im=+1, "+name);
			ck.Close(pot.RhoImag[1], -1.0, 1.0, 0.0, 1e-10, "sign 1D m=-1 Im=-1, "+name);
			ck.Close(pot.RhoReal[2], 1.0+FixedParticles, 1.0, 0.0, 1e-10, "sign 1D m=0 Re=N, "+name);
		}
		{
			std::vector<GeometryVector> basis(2, GeometryVector(static_cast<int>(2)));
			basis[0].x[0]=1.0;
			basis[1].x[0]=0.4;
			basis[1].x[1]=0.9;
			Configuration config(2, basis.data(), 10.0);
			GeometryVector s(static_cast<int>(2));
			s.x[0]=0.25;
			s.x[1]=0.125;
			config.Insert("A", s);
			ShiftedCCPotential pot(2);
			pot.constraints.emplace_back(LatticeK(config, Mode{1, 1}), 1.0);
			pot.constraints.emplace_back(LatticeK(config, Mode{-2, 3}), 1.0);
			pot.SetRhoBackend(b);
			pot.SetConfiguration(config);
			pot.GetRho();
			ck.Check(pot.GetLastRhoBackendUsed()==b, "sign 2D skew: "+name+" was used");
			const double p0=2.0*::pi*(0.25+0.125), p1=2.0*::pi*(-2*0.25+3*0.125);
			ck.Close(pot.RhoReal[0], std::cos(p0)+FixedParticles, 1.0, 0.0, 1e-10, "sign 2D (1,1) Re, "+name);
			ck.Close(pot.RhoImag[0], std::sin(p0), 1.0, 0.0, 1e-10, "sign 2D (1,1) Im, "+name);
			ck.Close(pot.RhoReal[1], std::cos(p1)+FixedParticles, 1.0, 0.0, 1e-10, "sign 2D (-2,3) Re, "+name);
			ck.Close(pot.RhoImag[1], std::sin(p1), 1.0, 0.0, 1e-10, "sign 2D (-2,3) Im, "+name);
		}
	}
}

/** \brief Empty constraint list and empty particle list, with the FINUFFT backend selected when available. */
void TestEmpty(Checker & ck)
{
	const Backend b=HaveFinufft ? Backend::FinufftType1 : Backend::Direct;
	std::mt19937_64 rng(202);
	{
		Configuration config=MakeConfig(2, {{1.0, 0.0}, {0.3, 1.1}}, 12, rng);
		ShiftedCCPotential pot(2);
		pot.ParallelNumThread=TestThreads;
		pot.Shift=0.7;
		pot.SetRhoBackend(b);
		pot.SetConfiguration(config);
		pot.GetRho();
		ck.Check(pot.RhoValid && pot.RhoReal.empty() && pot.RhoImag.empty(), "no constraints: empty rho, valid");
		ck.Check(pot.GetLastRhoBackendUsed()==Backend::Direct, "no constraints: reported Direct");
		ck.Check(pot.Energy()==0.7, "no constraints: Energy == Shift");
		std::vector<GeometryVector> f;
		pot.AllForce(f);
		bool zero=f.size()==12;
		for(const auto & v : f)
			zero=zero && v.x[0]==0.0 && v.x[1]==0.0;
		ck.Check(zero, "no constraints: zero forces");
	}
	{
		Configuration config=MakeConfig(3, {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}, 0, rng);
		ShiftedCCPotential pot=MakePotential(config, AllModes(3, 1), rng, 0.5, 1.5);
		pot.SetRhoBackend(b);
		pot.SetConfiguration(config);
		pot.GetRho();
		ck.Check(pot.RhoValid && pot.RhoReal.size()==27, "no particles: rho sized, valid");
		ck.Check(pot.GetLastRhoBackendUsed()==Backend::Direct, "no particles: reported Direct");
		bool ok=true;
		double expectedE=0.0;
		for(size_t i=0; i<pot.RhoReal.size(); i++)
		{
			ok=ok && pot.RhoReal[i]==FixedParticles && pot.RhoImag[i]==0.0;
			expectedE+=pot.constraints[i].V*FixedParticles*FixedParticles;
		}
		ck.Check(ok, "no particles: rho == fixed-particle offset");
		expectedE/=config.PeriodicVolume();
		ck.Close(pot.Energy(), expectedE, std::max(std::abs(expectedE), 1.0), 1e-12, 1e-14, "no particles: Energy");
		std::vector<GeometryVector> f;
		pot.AllForce(f);
		ck.Check(f.empty(), "no particles: AllForce empty");
	}
}

/** \brief Direct vs FINUFFT on full shuffled mode sets with duplicates, dims 1-3, orthogonal and skew cells. */
void TestLatticeCases(Checker & ck)
{
	struct Case
	{
		std::string label;
		DimensionType dim;
		std::vector<std::vector<double> > basis;
		size_t n;
		long K;
	};
	const std::vector<Case> cases={
		{"1D", 1, {{2.7}}, 64, 25},
		{"2D orthogonal", 2, {{1.3, 0.0}, {0.0, 0.8}}, 100, 7},
		{"2D skew", 2, {{1.3, 0.0}, {0.45, 0.9}}, 100, 7},
		{"2D rotated skew", 2, {{1.0, 0.3}, {-0.4, 1.1}}, 100, 6},
		{"3D orthogonal", 3, {{1.0, 0.0, 0.0}, {0.0, 1.2, 0.0}, {0.0, 0.0, 0.9}}, 60, 3},
		{"3D skew", 3, {{1.1, 0.0, 0.0}, {0.3, 0.95, 0.1}, {0.2, -0.25, 1.2}}, 60, 3},
	};
	std::mt19937_64 rng(303);
	for(const Case & c : cases)
	{
		Configuration config=MakeConfig(c.dim, c.basis, c.n, rng);
		std::vector<Mode> modes=Scramble(AllModes(c.dim, c.K), 7, rng);
		ShiftedCCPotential proto=MakePotential(config, modes, rng);
		proto.Shift=0.25;
		CompareBackends(ck, proto, config, c.label, true);
	}

	//looser tolerance still runs FINUFFT and agrees at that level
	Configuration config=MakeConfig(2, {{1.3, 0.0}, {0.45, 0.9}}, 80, rng);
	ShiftedCCPotential proto=MakePotential(config, Scramble(AllModes(2, 5), 3, rng), rng);
	CompareBackends(ck, proto, config, "2D skew eps=1e-6", true, 1e-6, 0, 1e-4);
}

/** \brief Cases that must fall back to the direct sum (bitwise identical results). */
void TestFallbacks(Checker & ck)
{
	std::mt19937_64 rng(404);
	{
		Configuration config=MakeConfig(2, {{1.3, 0.0}, {0.45, 0.9}}, 50, rng);
		ShiftedCCPotential proto=MakePotential(config, AllModes(2, 3), rng);
		proto.constraints.emplace_back(1.37*config.GetReciprocalBasisVector(0), 1.0);
		CompareBackends(ck, proto, config, "fallback: non-lattice k (2D)", false);
	}
	{
		Configuration config=MakeConfig(3, {{1.1, 0.0, 0.0}, {0.3, 0.95, 0.1}, {0.2, -0.25, 1.2}}, 30, rng);
		ShiftedCCPotential proto=MakePotential(config, AllModes(3, 2), rng);
		GeometryVector k=LatticeK(config, Mode{1, -1, 2})+1e-7*config.GetReciprocalBasisVector(1);
		proto.constraints.emplace_back(k, 1.0);
		CompareBackends(ck, proto, config, "fallback: near-lattice k off by 1e-7 (3D)", false);
	}
	{
		//index 1e6 off the lattice by 1e-8: |n|<2^30, grid 2e6+1 <= 2^22 and <= 64*(Nk+N) with N=32000,
		//so only the phase budget rejects it (the former max(0.1 eps, 256 DBL_EPSILON (1+|n|)) accepted it)
		Configuration config=MakeConfig(1, {{1.0}}, 32000, rng);
		ShiftedCCPotential proto=MakePotential(config, {Mode{1}, Mode{-3}}, rng);
		proto.constraints.emplace_back(LatticeK(config, Mode{1000000})+1e-8*config.GetReciprocalBasisVector(0), 1.0);
		CompareBackends(ck, proto, config, "fallback: high-index near-lattice k (1D, n=1e6+1e-8)", false);
	}
	{
		Configuration config=MakeConfig(2, {{1.0, 0.0}, {0.0, 1.0}}, 10, rng);
		//indices small enough for the phase budget; grid 31^2=961 > 64*(3+10) trips the sparsity guard
		ShiftedCCPotential proto=MakePotential(config, {Mode{15, 15}, Mode{1, 0}, Mode{-2, 1}}, rng);
		CompareBackends(ck, proto, config, "fallback: pathological sparse grid", false);
	}
	{
		Configuration config=MakeConfig(1, {{1.0}}, 10, rng);
		ShiftedCCPotential proto=MakePotential(config, {Mode{1}, Mode{long(1)<<40}}, rng);
		CompareBackends(ck, proto, config, "fallback: index beyond 2^30", false);
	}
	{
		Configuration config=MakeConfig(2, {{1.0, 0.0}, {0.2, 1.0}}, 100, rng);
		ShiftedCCPotential proto=MakePotential(config, AllModes(2, 5), rng);
		CompareBackends(ck, proto, config, "fallback: grid 121 > max modes 120", false, 1e-12, 120);
		CompareBackends(ck, proto, config, "grid 121 == max modes 121 runs FINUFFT", true, 1e-12, 121);
	}
	{
		Configuration config=MakeConfig(4, {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}}, 6, rng);
		ShiftedCCPotential proto=MakePotential(config, AllModes(4, 1), rng);
		CompareBackends(ck, proto, config, "fallback: unsupported dimension 4", false);
	}
	{
		//eps far below double roundoff: the phase budget rejects every mode before FINUFFT is called
		Configuration config=MakeConfig(2, {{1.0, 0.0}, {0.2, 1.0}}, 40, rng);
		ShiftedCCPotential proto=MakePotential(config, AllModes(2, 3), rng);
		CompareBackends(ck, proto, config, "fallback: eps=1e-20 below the phase budget", false, 1e-20);
	}
}

/** \brief Evaluate pot in place (it must already be bound to config) and compare with a fresh Direct copy.
 *  \param ck        Checker.
 *  \param pot       Potential evaluated in place.
 *  \param config    Configuration pot is bound to.
 *  \param expected  Backend pot must report as used.
 *  \param label     Case name. */
void CheckAgainstFreshDirect(Checker & ck, ShiftedCCPotential & pot, const Configuration & config, Backend expected, const std::string & label)
{
	pot.GetRho();
	ck.Check(pot.RhoValid && pot.GetLastRhoBackendUsed()==expected, label+": expected backend used");
	if(expected==Backend::FinufftType1)
		FinufftEvaluations++;
	ShiftedCCPotential ref(pot);
	ref.SetRhoBackend(Backend::Direct);
	ref.SetConfiguration(config);
	ref.GetRho();
	const Scales s=MakeScales(pot, config);
	ck.Check(pot.RhoReal.size()==ref.RhoReal.size(), label+": rho size");
	for(size_t i=0; i<ref.RhoReal.size() && i<pot.RhoReal.size(); i++)
	{
		ck.Close(pot.RhoReal[i], ref.RhoReal[i], s.rho, 1e-9, 1e-10, label+": Re rho", i);
		ck.Close(pot.RhoImag[i], ref.RhoImag[i], s.rho, 1e-9, 1e-10, label+": Im rho", i);
	}
	ck.Close(pot.Energy(), ref.Energy(), s.energy, 1e-9, 1e-10, label+": Energy");
}

/** \brief Cache invalidation on backend/tolerance/configuration changes, and no stale mode mapping after
 *  constraints are mutated. */
void TestCacheInvalidation(Checker & ck)
{
	std::mt19937_64 rng(505);
	Configuration config=MakeConfig(2, {{1.2, 0.0}, {0.35, 1.0}}, 60, rng);
	ShiftedCCPotential proto=MakePotential(config, Scramble(AllModes(2, 4), 2, rng), rng);
	ShiftedCCPotential direct(proto), pot(proto);
	direct.SetConfiguration(config);
	direct.GetRho();
	const Scales s=MakeScales(proto, config);
	const Backend fastBackend=HaveFinufft ? Backend::FinufftType1 : Backend::Direct;

	pot.SetRhoBackend(fastBackend);
	pot.SetConfiguration(config);
	pot.GetRho();
	ck.Check(pot.RhoValid && pot.GetLastRhoBackendUsed()==fastBackend, "cache: first evaluation backend");

	pot.SetRhoBackend(Backend::Direct);
	ck.Check(!pot.RhoValid, "cache: SetRhoBackend(Direct) invalidates");
	pot.GetRho();
	ck.Check(pot.GetLastRhoBackendUsed()==Backend::Direct && pot.RhoReal==direct.RhoReal && pot.RhoImag==direct.RhoImag, "cache: Direct recompute equals reference exactly");

	pot.SetRhoBackend(fastBackend);
	pot.GetRho();
	std::vector<double> tightReal=pot.RhoReal;
	pot.SetFinufftTolerance(1e-6);
	ck.Check(!pot.RhoValid, "cache: SetFinufftTolerance invalidates");
	pot.GetRho();
	ck.Check(pot.GetLastRhoBackendUsed()==fastBackend, "cache: eps=1e-6 evaluation backend");
	for(size_t i=0; i<pot.RhoReal.size(); i++)
		ck.Close(pot.RhoReal[i], direct.RhoReal[i], s.rho, 1e-4, 1e-10, "cache: eps=1e-6 rho", i);
	pot.SetFinufftTolerance(1e-12);
	pot.GetRho();
	for(size_t i=0; i<pot.RhoReal.size(); i++)
		ck.Close(pot.RhoReal[i], tightReal[i], s.rho, 1e-9, 1e-10, "cache: back to eps=1e-12", i);

	//mutate constraints of the same, already evaluated object: reorder, change one k, append a
	//larger mode; the mode mapping must be recomputed rather than reused
	std::reverse(pot.constraints.begin(), pot.constraints.end());
	pot.constraints[0].k=LatticeK(config, Mode{-3, 2});
	pot.constraints.emplace_back(LatticeK(config, Mode{6, -5}), 0.8);
	pot.RhoValid=false;
	CheckAgainstFreshDirect(ck, pot, config, fastBackend, "cache: mutated constraints");

	//moving a particle and rebinding invalidates
	Configuration moved(config);
	GeometryVector sNew=moved.GetRelativeCoordinates(3);
	sNew.x[0]+=0.123;
	moved.MoveParticle(3, sNew);
	pot.SetConfiguration(moved);
	ck.Check(!pot.RhoValid, "cache: SetConfiguration invalidates");
	CheckAgainstFreshDirect(ck, pot, moved, fastBackend, "cache: moved particle");
}

/** \brief Central finite-difference check of Force against Energy for each available backend. */
void TestFiniteDifference(Checker & ck)
{
	std::mt19937_64 rng(606);
	const Configuration base=MakeConfig(2, {{1.3, 0.0}, {0.45, 0.9}}, 30, rng);
	const ShiftedCCPotential proto=MakePotential(base, Scramble(AllModes(2, 4), 2, rng), rng, 0.2, 1.2);
	std::vector<Backend> backends{Backend::Direct};
	if(HaveFinufft)
		backends.push_back(Backend::FinufftType1);
	const double h=5e-5;
	for(Backend b : backends)
	{
		const std::string name=(b==Backend::Direct) ? "Direct" : "FinufftType1";
		Configuration c(base);
		ShiftedCCPotential pot(proto);
		pot.SetRhoBackend(b);
		pot.SetConfiguration(c);
		const Scales s=MakeScales(pot, c);
		for(size_t i : {size_t(0), size_t(17)})
		{
			pot.SetConfiguration(c);
			GeometryVector F;
			pot.Force(F, i);
			const GeometryVector s0=c.GetRelativeCoordinates(i);
			const GeometryVector r0=c.GetCartesianCoordinates(i);
			for(DimensionType d=0; d<2; d++)
			{
				double E[2];
				for(int side=0; side<2; side++)
				{
					GeometryVector r=r0;
					r.x[d]+=(side==0 ? h : -h);
					c.MoveParticle(i, c.CartesianCoord2RelativeCoord(r));
					pot.SetConfiguration(c);
					E[side]=pot.Energy();
					ck.Check(pot.GetLastRhoBackendUsed()==b, "finite difference: "+name+" was used");
				}
				c.MoveParticle(i, s0);
				const double fd=-(E[0]-E[1])/(2.0*h);
				ck.Close(F.x[d], fd, s.force, 1e-5, 1e-8, "finite-difference force, "+name, i);
			}
		}
	}
}

/** \brief Derived classes: weighted and multi-config overrides keep their direct sums; S2 inherits the backend. */
void TestDerived(Checker & ck)
{
	std::mt19937_64 rng(707);
	const Backend fastBackend=HaveFinufft ? Backend::FinufftType1 : Backend::Direct;
	Configuration config=MakeConfig(2, {{1.3, 0.0}, {0.45, 0.9}}, 40, rng);
	std::vector<Mode> modes=Scramble(AllModes(2, 4), 2, rng);
	ShiftedCCPotential base=MakePotential(config, modes, rng, 0.2, 1.2);
	base.SetRhoBackend(fastBackend);
	const size_t nk=base.constraints.size(), np=config.NumParticle();

	{
		std::vector<double> weights(np);
		std::uniform_real_distribution<double> u(0.5, 2.0);
		double sumW=0.0;
		for(double & w : weights)
		{
			w=u(rng);
			sumW+=w;
		}
		ShiftedCCPotential_withWeight weighted(base, weights);
		ck.Check(weighted.GetRhoBackend()==fastBackend, "withWeight: setting is inherited (and ignored by its GetRho)");
		weighted.SetConfiguration(config);
		weighted.GetRho();
		base.SetConfiguration(config);
		base.GetRho();
		double maxDiffUnweighted=0.0;
		for(size_t i=0; i<nk; i++)
		{
			double tr=FixedParticles, ti=0.0;
			for(size_t j=0; j<np; j++)
			{
				double phase=base.constraints[i].k.Dot(config.GetCartesianCoordinates(j));
				tr+=std::cos(phase)*weights[j];
				ti+=std::sin(phase)*weights[j];
			}
			ck.Close(weighted.RhoReal[i], tr, sumW, 1e-13, 1e-13, "withWeight: Re rho is the weighted direct sum", i);
			ck.Close(weighted.RhoImag[i], ti, sumW, 1e-13, 1e-13, "withWeight: Im rho is the weighted direct sum", i);
			maxDiffUnweighted=std::max(maxDiffUnweighted, std::abs(weighted.RhoReal[i]-base.RhoReal[i]));
		}
		ck.Check(maxDiffUnweighted>1e-3, "withWeight: result differs from the unweighted rho");
	}
	{
		const size_t nsub=2, npSub=np/nsub;
		ShiftedCCPotential_S2_MultiConfig multi(2, nsub);
		multi.ParallelNumThread=TestThreads;
		multi.constraints=base.constraints;
		multi.S0.assign(nk, 0.1);
		multi.SetRhoBackend(fastBackend);
		multi.SetConfiguration(config);
		multi.GetRho();
		ck.Check(multi.RhoReal.size()==nk*nsub, "MultiConfig: rho holds one block per sub-configuration");
		for(size_t nc=0; nc<nsub; nc++)
			for(size_t i=0; i<nk; i++)
			{
				double tr=0.0, ti=0.0;
				for(size_t j=nc*npSub; j<(nc+1)*npSub; j++)
				{
					double phase=multi.constraints[i].k.Dot(config.GetCartesianCoordinates(j));
					tr+=std::cos(phase);
					ti+=std::sin(phase);
				}
				ck.Close(multi.RhoReal[nc*nk+i], tr, double(npSub), 1e-13, 1e-13, "MultiConfig: Re rho per sub-configuration", i);
				ck.Close(multi.RhoImag[nc*nk+i], ti, double(npSub), 1e-13, 1e-13, "MultiConfig: Im rho per sub-configuration", i);
			}
	}
	{
		ShiftedCCPotential_S2 s2(2);
		s2.ParallelNumThread=TestThreads;
		std::uniform_real_distribution<double> u(0.0, 2.0);
		for(const auto & c : base.constraints)
			s2.AddConstraint(c.k, {c.V, u(rng)});
		ShiftedCCPotential_S2 s2direct(s2);
		s2.SetRhoBackend(fastBackend);
		s2.SetConfiguration(config);
		s2direct.SetConfiguration(config);
		const double e1=s2.Energy(), e0=s2direct.Energy();
		ck.Check(s2.GetLastRhoBackendUsed()==fastBackend, "S2: inherits the selected backend");
		if(fastBackend==Backend::FinufftType1)
			FinufftEvaluations++;
		double sumV=0.0, kmax=0.0, n=np+FixedParticles;
		for(const auto & c : s2.constraints)
		{
			sumV+=std::abs(c.V);
			kmax=std::max(kmax, std::sqrt(c.k.Modulus2()));
		}
		ck.Close(e1, e0, sumV*(n+2.0)*(n+2.0), 1e-9, 1e-10, "S2: Energy");
		for(size_t j : {size_t(0), size_t(np-1)})
		{
			GeometryVector f1, f0;
			s2.Force(f1, j);
			s2direct.Force(f0, j);
			for(DimensionType d=0; d<2; d++)
				ck.Close(f1.x[d], f0.x[d], 8.0*sumV*(n+2.0)*kmax, 1e-9, 1e-10, "S2: Force", j);
		}
	}
}

/** \brief Reproducible N=500 2D stealthy-like smoke comparison (not an optimizer benchmark).
 *  Uses the 299 shortest half-plane modes (one per +/-k pair), chi = 299/(2*(500-1)) ~ 0.3. */
void TestSmokeN500(Checker & ck)
{
	std::mt19937_64 rng(20260927);
	const size_t n=500, pairs=299;
	const double L=std::sqrt(double(n));
	Configuration config=MakeConfig(2, {{L, 0.0}, {0.0, L}}, n, rng);
	std::vector<Mode> modes;
	const long R=16;
	for(long m1=-R; m1<=R; m1++)
		for(long m0=-R; m0<=R; m0++)
			if(m0*m0+m1*m1<=R*R && (m1>0 || (m1==0 && m0>0)))
				modes.push_back(Mode{m0, m1});
	//total order (|m|^2, m1, m0) makes the cut at the last shell deterministic
	std::sort(modes.begin(), modes.end(), [](const Mode & a, const Mode & b)
	{
		const long na=a[0]*a[0]+a[1]*a[1], nb=b[0]*b[0]+b[1]*b[1];
		if(na!=nb)
			return na<nb;
		return a[1]!=b[1] ? a[1]<b[1] : a[0]<b[0];
	});
	modes.resize(pairs);
	ShiftedCCPotential proto=MakePotential(config, modes, rng);
	for(auto & c : proto.constraints)
		c.V=1.0;
	const Errors e=CompareBackends(ck, proto, config, "smoke N=500 2D", true);

	ShiftedCCPotential direct(proto), fast(proto);
	if(HaveFinufft)
		fast.SetRhoBackend(Backend::FinufftType1);
	direct.SetConfiguration(config);
	fast.SetConfiguration(config);
	auto t0=std::chrono::steady_clock::now();
	direct.GetRho();
	auto t1=std::chrono::steady_clock::now();
	fast.GetRho();
	auto t2=std::chrono::steady_clock::now();
	std::cout.precision(3);
	std::cout<<"  smoke N=500 2D: Nk="<<modes.size()<<" chi="<<modes.size()/(2.0*(n-1))
		<<" max|drho|="<<e.rho<<" |dE|="<<e.energy<<" max|dF|="<<e.force
		<<" GetRho direct="<<std::chrono::duration<double, std::milli>(t1-t0).count()<<" ms"
		<<" selected="<<std::chrono::duration<double, std::milli>(t2-t1).count()<<" ms (informational, "<<TestThreads<<" threads)\n";
}
}

/** \brief Run all checks.
 *  \return 0 if every check passed, 1 otherwise. */
int main()
{
	std::cout<<"FINUFFT backend compiled in: "<<(HaveFinufft ? "yes" : "no")
		<<"; HAVETWOFIXEDPARTICLES: "<<(FixedParticles!=0.0 ? "yes" : "no")<<'\n';
	Checker ck;
	TestBackendSelection(ck);
	TestSign(ck);
	TestEmpty(ck);
	TestLatticeCases(ck);
	TestFallbacks(ck);
	TestCacheInvalidation(ck);
	TestFiniteDifference(ck);
	TestDerived(ck);
	TestSmokeN500(ck);
	if(HaveFinufft)
		ck.Check(FinufftEvaluations>=10, "FINUFFT was really exercised (evaluations="+std::to_string(FinufftEvaluations)+")");
	std::cout<<"FINUFFT evaluations compared: "<<FinufftEvaluations<<'\n';
	std::cout<<(ck.failures==0 ? "PASS" : "FAIL")<<": "<<ck.checks-ck.failures<<"/"<<ck.checks<<" checks passed\n";
	return ck.failures==0 ? 0 : 1;
}
