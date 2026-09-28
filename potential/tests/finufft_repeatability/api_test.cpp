/** \file api_test.cpp
 * \brief API regression for ShiftedCCPotential::SetFinufftNumThreads / GetFinufftNumThreads.
 * run_repeatability.py builds this file twice: with -DCCPOTENTIAL_USE_FINUFFT (linked with
 * finufft_observer.o, -Wl,--wrap=finufft2d1, and libfinufft), and without the macro, the
 * FINUFFT include paths, or the library. Prints one JSON line; exits nonzero on the first failure.
 */
#include "MultiStealthyPotential.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CCPOTENTIAL_USE_FINUFFT
#include "finufft_observer.h"
#endif

namespace {
using Backend=ShiftedCCPotential::RhoBackend;
constexpr size_t N=60;
const double L=std::sqrt(double(N));
size_t checks=0;

/** \brief Count a check and throw if it fails.
 *  \param value    Condition that must hold.
 *  \param message  Description reported on failure. */
void check(bool value,const std::string& message) {
    ++checks;
    if(!value)throw std::runtime_error(message);
}

/** \brief Discard std::cout while noisy library setters run. */
struct Quiet {
    std::ostringstream sink;
    std::streambuf* old;
    /** \brief Start redirecting std::cout. */
    Quiet():old(std::cout.rdbuf(sink.rdbuf())) {}
    /** \brief Restore std::cout. */
    ~Quiet() {std::cout.rdbuf(old);}
};

/** \brief Build a fixed-seed random 2D periodic configuration.
 *  \return N uniformly random particles in a square box of side L. */
Configuration configuration() {
    GeometryVector basis[2]={GeometryVector(2),GeometryVector(2)};
    basis[0].x[0]=L;basis[1].x[1]=L;
    Configuration c(2,basis,1.0);
    std::mt19937_64 rng(20260927);
    std::uniform_real_distribution<double> u(0,1);
    for(size_t i=0;i<N;++i){GeometryVector s(2);s.x[0]=u(rng);s.x[1]=u(rng);c.Insert("a",s);}
    return c;
}

/** \brief Half-plane reciprocal-lattice wavevectors with |n_x|,|n_y| <= 3.
 *  \return Wavevectors 2 pi (n_x, n_y) / L. */
std::vector<GeometryVector> wavevectors() {
    std::vector<GeometryVector> result;
    for(int y=0;y<=3;++y)for(int x=-3;x<=3;++x)
        if(y>0||x>0){GeometryVector k(2);k.x[0]=2*::pi*x/L;k.x[1]=2*::pi*y/L;result.push_back(k);}
    return result;
}

/** \brief Create a potential with unit-weight constraints, bound to a configuration.
 *  \param c        Configuration to bind; must outlive the potential's evaluations.
 *  \param backend  Backend to select.
 *  \return The configured potential. */
ShiftedCCPotential potential(const Configuration& c,Backend backend) {
    ShiftedCCPotential p(2);
    for(const auto& k:wavevectors())p.AddConstraint(k,{1.0});
    p.SetConfiguration(c);
    p.SetRhoBackend(backend);
    return p;
}

/** \brief Option and cache state that an invalid setter call must leave unchanged. */
struct State {
    int threads;
    Backend backend;
    double tolerance;
    std::int64_t maxModes;
    unsigned long parallel;
    bool valid;
    std::vector<double> re,im;
    /** \brief Capture the state of a potential.
     *  \param p  Potential to inspect. */
    explicit State(const ShiftedCCPotential& p):threads(p.GetFinufftNumThreads()),backend(p.GetRhoBackend()),
        tolerance(p.GetFinufftTolerance()),maxModes(p.GetFinufftMaxModes()),parallel(p.ParallelNumThread),
        valid(p.RhoValid),re(p.RhoReal),im(p.RhoImag) {}
    /** \brief Exact comparison of every captured field.
     *  \param o  State to compare with.
     *  \return true if all fields, including every cached density value, are equal. */
    bool operator==(const State& o) const {
        return threads==o.threads&&backend==o.backend&&tolerance==o.tolerance&&maxModes==o.maxModes
            &&parallel==o.parallel&&valid==o.valid&&re==o.re&&im==o.im;
    }
};

/** \brief Largest absolute difference between the cached densities of two potentials.
 *  \param a  First potential with a valid cache.
 *  \param b  Second potential with a valid cache and the same constraints.
 *  \return max over modes of |Re a - Re b| and |Im a - Im b|. */
double rhoDifference(const ShiftedCCPotential& a,const ShiftedCCPotential& b) {
    double result=0;
    for(size_t i=0;i<a.RhoReal.size();++i)
        result=std::max({result,std::abs(a.RhoReal[i]-b.RhoReal[i]),std::abs(a.RhoImag[i]-b.RhoImag[i])});
    return result;
}

/** \brief Default construction gives one FINUFFT thread, the Direct backend, and no cache. */
void testDefaults() {
    ShiftedCCPotential p(2);
    check(p.GetFinufftNumThreads()==1,"default FINUFFT thread count is not 1");
    check(p.GetRhoBackend()==Backend::Direct,"default backend is not Direct");
    check(!p.RhoValid,"default-constructed cache is valid");
    check(p.GetLastFinufftStatus()==ShiftedCCPotential::FinufftNotCalled,"default FINUFFT status");
    ShiftedCCPotential_S2 s2(2);
    check(s2.GetFinufftNumThreads()==1,"ShiftedCCPotential_S2 default FINUFFT thread count is not 1");
    Quiet quiet;
    MultiStealthyPotential m(2,1);
    m.SetNumSpecies(2);
    for(const auto& s:m.species_constraints){
        check(s.CCPotential->GetFinufftNumThreads()==1,"species potential default FINUFFT thread count is not 1");
        check(s.CCPotential->GetRhoBackend()==Backend::Direct,"species potential default backend is not Direct");
    }
}

/** \brief Valid values are stored and returned by the getter. */
void testGetter() {
    ShiftedCCPotential p(2);
    for(int value:{2,4,7,1,INT_MAX,1}){
        p.SetFinufftNumThreads(value);
        check(p.GetFinufftNumThreads()==value,"getter does not return value "+std::to_string(value));
    }
}

/** \brief Invalid values throw std::invalid_argument and leave options and cache unchanged.
 *  \param c        Configuration to bind.
 *  \param backend  Backend used to fill the cache. */
void testInvalid(const Configuration& c,Backend backend) {
    auto p=potential(c,backend);
    p.SetFinufftNumThreads(3);
    p.Energy();
    check(p.RhoValid,"cache not valid after Energy");
    const State before(p);
    for(int bad:{0,-1,INT_MIN}){
        bool threw=false;
        try{p.SetFinufftNumThreads(bad);}catch(const std::invalid_argument&){threw=true;}
        check(threw,"value "+std::to_string(bad)+" did not throw std::invalid_argument");
        check(State(p)==before,"value "+std::to_string(bad)+" modified option or cache state");
    }
#ifdef CCPOTENTIAL_USE_FINUFFT
    const auto calls=finufft_observer::Calls();
    p.Energy();
    check(finufft_observer::Calls()==calls,"invalid setter caused a FINUFFT recomputation");
#endif
    ShiftedCCPotential fresh(2);
    try{fresh.SetFinufftNumThreads(0);}catch(const std::invalid_argument&){}
    check(fresh.GetFinufftNumThreads()==1&&!fresh.RhoValid,"invalid value modified a fresh object");
}

/** \brief A valid setter call invalidates the cache, including when the value is unchanged.
 *  \param c        Configuration to bind.
 *  \param backend  Backend used to fill the cache. */
void testInvalidation(const Configuration& c,Backend backend) {
    auto p=potential(c,backend);
    p.Energy();
    check(p.RhoValid,"cache not valid after Energy");
    p.SetFinufftNumThreads(2);
    check(!p.RhoValid,"changed valid value did not invalidate the cache");
    p.Energy();
    check(p.RhoValid,"cache not valid after second Energy");
    p.SetFinufftNumThreads(2);
    check(!p.RhoValid,"unchanged valid value did not invalidate the cache");
#ifdef CCPOTENTIAL_USE_FINUFFT
    const auto calls=finufft_observer::Calls();
    p.Energy();
    const size_t expected=backend==Backend::FinufftType1 ? 1 : 0;
    check(finufft_observer::Calls()==calls+expected,"unexpected FINUFFT call count after invalidation");
#endif
}

/** \brief Copies, assignment, clone, and converting constructors preserve the option. */
void testCopies() {
    ShiftedCCPotential p(2);
    for(const auto& k:wavevectors())p.AddConstraint(k,{1.0});
    p.SetFinufftNumThreads(3);
    ShiftedCCPotential copy(p);
    check(copy.GetFinufftNumThreads()==3,"copy constructor");
    ShiftedCCPotential assigned(2);
    assigned=p;
    check(assigned.GetFinufftNumThreads()==3,"copy assignment");
    std::unique_ptr<ShiftedCCPotential> cloned(p.clone());
    check(cloned->GetFinufftNumThreads()==3,"clone");
    LogCCPotential log(p);
    check(log.GetFinufftNumThreads()==3,"LogCCPotential converting constructor");
    std::unique_ptr<LogCCPotential> logClone(log.clone());
    check(logClone->GetFinufftNumThreads()==3,"LogCCPotential clone");
    ShiftedCCPotential_ForMC mc(p);
    check(mc.GetFinufftNumThreads()==3,"ShiftedCCPotential_ForMC converting constructor");
    ShiftedCCPotential_withWeight weighted(p,std::vector<double>(N,1.0));
    std::unique_ptr<ShiftedCCPotential_withWeight> weightedClone(weighted.clone());
    check(weighted.GetFinufftNumThreads()==3&&weightedClone->GetFinufftNumThreads()==3,"ShiftedCCPotential_withWeight copies");
    copy.SetFinufftNumThreads(5);
    check(p.GetFinufftNumThreads()==3&&copy.GetFinufftNumThreads()==5,"copies share option state");
}

/** \brief Changing ParallelNumThread, directly or through MultiStealthyPotential, leaves the option unchanged.
 *  \param c  Configuration to bind. */
void testParallelIndependence(const Configuration& c) {
    auto p=potential(c,Backend::Direct);
    p.ParallelNumThread=4;
    check(p.GetFinufftNumThreads()==1,"ParallelNumThread=4 changed the default FINUFFT option");
    p.Energy();
    check(p.GetFinufftNumThreads()==1,"Energy with ParallelNumThread=4 changed the FINUFFT option");
    p.SetFinufftNumThreads(3);
    p.ParallelNumThread=2;
    check(p.GetFinufftNumThreads()==3,"ParallelNumThread=2 changed an explicit FINUFFT option");
    const Backend backend=ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1) ? Backend::FinufftType1 : Backend::Direct;
    std::unique_ptr<MultiStealthyPotential> m;
    {
        Quiet quiet;
        m=std::make_unique<MultiStealthyPotential>(2,1);
        m->ParallelNumThread=4;
        m->SetConfiguration(c);
        m->SetNumSpecies(2);
        m->SetSecureRadii(std::vector<double>(N,0.1));
        for(size_t i=0;i<N;++i)m->SetSpeciesOfParticle(i,i%2);
        for(size_t s=0;s<2;++s){
            for(const auto& k:wavevectors())m->AddConstraint(s,k,{1.0});
            m->species_constraints[s].CCPotential->SetRhoBackend(backend);
        }
        m->SetConfiguration(c);
    }
#ifdef CCPOTENTIAL_USE_FINUFFT
    const auto calls=finufft_observer::Calls(),ones=finufft_observer::CallsWithNthreads(1);
#endif
    const double energy=m->Energy();
    check(std::isfinite(energy),"MultiStealthyPotential energy is not finite");
    for(const auto& s:m->species_constraints){
        check(s.CCPotential->ParallelNumThread==4,"ParallelNumThread was not propagated to species potentials");
        check(s.CCPotential->GetFinufftNumThreads()==1,"propagated ParallelNumThread changed the FINUFFT option");
        check(s.CCPotential->GetLastRhoBackendUsed()==backend,"species potential used an unexpected backend");
    }
#ifdef CCPOTENTIAL_USE_FINUFFT
    const auto made=finufft_observer::Calls()-calls;
    check(made>=2,"expected at least one FINUFFT call per species");
    check(finufft_observer::CallsWithNthreads(1)-ones==made,"species FINUFFT calls did not request nthreads=1");
#endif
}

#ifdef CCPOTENTIAL_USE_FINUFFT
/** \brief The option reaches finufft2d1 as opts.nthreads (1 by default, explicit 4), independent of ParallelNumThread.
 *  \param c  Configuration to bind. */
void testForwarding(const Configuration& c) {
    auto direct=potential(c,Backend::Direct);
    auto fast=potential(c,Backend::FinufftType1);
    const auto calls=finufft_observer::Calls();
    const auto ones=finufft_observer::CallsWithNthreads(1),fours=finufft_observer::CallsWithNthreads(4);
    direct.GetRho();
    check(finufft_observer::Calls()==calls,"Direct GetRho called FINUFFT");

    fast.GetRho();
    check(finufft_observer::Calls()==calls+1&&finufft_observer::CallsWithNthreads(1)==ones+1,"default option did not reach FINUFFT as nthreads=1");
    check(finufft_observer::LastNthreads()==1,"last observed nthreads is not 1");
    check(fast.GetLastRhoBackendUsed()==Backend::FinufftType1&&fast.GetLastFinufftStatus()==0,"FINUFFT fell back (default option)");
    check(rhoDifference(direct,fast)<1e-8,"FINUFFT density differs from direct (default option)");

    fast.SetFinufftNumThreads(4);
    fast.GetRho();
    check(finufft_observer::Calls()==calls+2&&finufft_observer::CallsWithNthreads(4)==fours+1,"explicit option 4 did not reach FINUFFT");
    check(finufft_observer::LastNthreads()==4,"last observed nthreads is not 4");
    check(fast.GetLastRhoBackendUsed()==Backend::FinufftType1&&fast.GetLastFinufftStatus()==0,"FINUFFT fell back (option 4)");
    check(rhoDifference(direct,fast)<1e-8,"FINUFFT density differs from direct (option 4)");

    fast.SetFinufftNumThreads(1);
    fast.ParallelNumThread=4;
    fast.GetRho();
    check(finufft_observer::Calls()==calls+3&&finufft_observer::LastNthreads()==1,"ParallelNumThread=4 was forwarded to FINUFFT");

    direct.SetFinufftNumThreads(4);
    direct.GetRho();
    check(finufft_observer::Calls()==calls+3,"Direct backend called FINUFFT");
    check(direct.GetLastFinufftStatus()==ShiftedCCPotential::FinufftNotCalled,"Direct backend reported a FINUFFT status");
}
#else
/** \brief Without the macro the option is stored but unused, and FINUFFT cannot be selected.
 *  \param c  Configuration to bind. */
void testMacroDisabled(const Configuration& c) {
    check(!ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1),"FinufftType1 reported available");
    auto p=potential(c,Backend::Direct);
    p.SetFinufftNumThreads(4);
    bool threw=false;
    try{p.SetRhoBackend(Backend::FinufftType1);}catch(const std::invalid_argument&){threw=true;}
    check(threw,"selecting FinufftType1 did not throw");
    check(p.GetRhoBackend()==Backend::Direct&&p.GetFinufftNumThreads()==4,"failed backend selection modified state");
    check(std::isfinite(p.Energy()),"direct energy is not finite");
    check(p.GetLastRhoBackendUsed()==Backend::Direct&&p.GetLastFinufftStatus()==ShiftedCCPotential::FinufftNotCalled,"direct evaluation state");
}
#endif
}

/** \brief Run every API check and print a one-line JSON receipt.
 *  \return 0 if all checks pass, 1 otherwise. */
int main() try {
#ifdef CCPOTENTIAL_USE_FINUFFT
    const bool macro=true;
#else
    const bool macro=false;
#endif
    const Configuration c=configuration();
    testDefaults();
    testGetter();
    testInvalid(c,Backend::Direct);
    testInvalidation(c,Backend::Direct);
    testCopies();
    testParallelIndependence(c);
#ifdef CCPOTENTIAL_USE_FINUFFT
    testInvalid(c,Backend::FinufftType1);
    testInvalidation(c,Backend::FinufftType1);
    testForwarding(c);
    const size_t observed=finufft_observer::Calls();
#else
    testMacroDisabled(c);
    const size_t observed=0;
#endif
    std::cout<<"{\"finufft_macro\":"<<(macro?"true":"false")<<",\"checks\":"<<checks
             <<",\"observed_finufft_calls\":"<<observed<<",\"passed\":true}\n";
    return 0;
} catch(const std::exception& e) {
    std::cout<<"{\"passed\":false,\"checks\":"<<checks<<",\"error\":\""<<e.what()<<"\"}\n";
    std::cerr<<e.what()<<'\n';
    return 1;
}
