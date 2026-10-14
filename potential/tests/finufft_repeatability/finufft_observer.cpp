/** \file finufft_observer.cpp
 * \brief Observation-only link wrapper for the FINUFFT guru interface; link with
 * -Wl,--wrap=finufft_makeplan,--wrap=finufft_setpts,--wrap=finufft_execute.
 * ShiftedCCPotential keeps one persistent plan per object and runs one finufft_execute per
 * transform, so a "call" is one finufft_execute; its thread count is the opts->nthreads that
 * the plan's finufft_makeplan received. Every argument, including the caller's options
 * pointer, is forwarded unchanged to the real functions. The wrapper only records the
 * transform count, the requested nthreads of each transform's plan, null options, and nonzero
 * return statuses of all three functions. At program exit it writes one receipt line to
 * stderr, which run_repeatability.py parses:
 *   FINUFFT_OBSERVER calls=C null_options=Z nonzero_status=S nthreads=N1:C1,N2:C2
 * (nthreads=none when no call was made).
 */
#include "finufft_observer.h"

#include <finufft.h>
#include <complex>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>

namespace {
/** \brief Observation state; a namespace-scope object so that the receipt prints even with zero calls. */
struct Observation {
    std::mutex lock;
    std::size_t calls=0,nullOptions=0,nonzeroStatus=0;
    int lastNthreads=0;
    std::map<int,std::size_t> histogram;
    std::map<finufft_plan,int> planNthreads;  //< opts->nthreads per plan; -1 for null options

    /** \brief Write the receipt line at program exit. */
    ~Observation() {
        std::cerr<<"FINUFFT_OBSERVER calls="<<calls<<" null_options="<<nullOptions
                 <<" nonzero_status="<<nonzeroStatus<<" nthreads=";
        if(histogram.empty())std::cerr<<"none";
        bool first=true;
        for(const auto& entry:histogram){
            if(!first)std::cerr<<',';
            first=false;
            std::cerr<<entry.first<<':'<<entry.second;
        }
        std::cerr<<'\n';
    }
} observed;
}

extern "C" int __real_finufft_makeplan(int,int,const int64_t*,int,int,double,finufft_plan*,const finufft_opts*);
extern "C" int __real_finufft_setpts(finufft_plan,int64_t,const double*,const double*,const double*,
    int64_t,const double*,const double*,const double*);
extern "C" int __real_finufft_execute(finufft_plan,std::complex<double>*,std::complex<double>*);

/** \brief Count a nonzero status returned by a real FINUFFT function.
 *  \param status  Returned status.
 *  \return status, unchanged. */
static int recordStatus(int status) {
    if(status!=0){
        std::lock_guard<std::mutex> guard(observed.lock);
        ++observed.nonzeroStatus;
    }
    return status;
}

/** \brief Call the real finufft_makeplan with identical arguments and remember the plan's requested nthreads.
 *  \return The status returned by the real finufft_makeplan. */
extern "C" int __wrap_finufft_makeplan(int type,int dim,const int64_t* nmodes,int iflag,int ntrans,double tol,
    finufft_plan* plan,const finufft_opts* opts) {
    const int status=__real_finufft_makeplan(type,dim,nmodes,iflag,ntrans,tol,plan,opts);
    if(*plan){
        std::lock_guard<std::mutex> guard(observed.lock);
        observed.planNthreads[*plan]=opts ? opts->nthreads : -1;
    }
    return recordStatus(status);
}

/** \brief Call the real finufft_setpts with identical arguments, recording a nonzero status.
 *  \return The status returned by the real finufft_setpts. */
extern "C" int __wrap_finufft_setpts(finufft_plan plan,int64_t m,const double* x,const double* y,const double* z,
    int64_t n,const double* s,const double* t,const double* u) {
    return recordStatus(__real_finufft_setpts(plan,m,x,y,z,n,s,t,u));
}

/** \brief Record one transform with its plan's requested nthreads, then call the real finufft_execute.
 *  A plan not made through the wrapper is counted as having null options (fails closed).
 *  \return The status returned by the real finufft_execute. */
extern "C" int __wrap_finufft_execute(finufft_plan plan,std::complex<double>* weights,std::complex<double>* result) {
    {
        std::lock_guard<std::mutex> guard(observed.lock);
        ++observed.calls;
        const auto found=observed.planNthreads.find(plan);
        const int nthreads=found==observed.planNthreads.end() ? -1 : found->second;
        if(nthreads>=0){
            observed.lastNthreads=nthreads;
            ++observed.histogram[nthreads];
        }else{
            observed.lastNthreads=0;
            ++observed.nullOptions;
        }
    }
    return recordStatus(__real_finufft_execute(plan,weights,result));
}

namespace finufft_observer {
std::size_t Calls() {
    std::lock_guard<std::mutex> guard(observed.lock);
    return observed.calls;
}

std::size_t CallsWithNthreads(int nthreads) {
    std::lock_guard<std::mutex> guard(observed.lock);
    auto found=observed.histogram.find(nthreads);
    return found==observed.histogram.end() ? 0 : found->second;
}

int LastNthreads() {
    std::lock_guard<std::mutex> guard(observed.lock);
    return observed.lastNthreads;
}
}
