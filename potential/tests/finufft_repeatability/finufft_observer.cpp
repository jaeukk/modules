/** \file finufft_observer.cpp
 * \brief Observation-only link wrapper for finufft2d1; link with -Wl,--wrap=finufft2d1.
 * Every argument, including the caller's options pointer, is forwarded unchanged to the real
 * function. The wrapper only records the call count, the opts->nthreads each call requested,
 * null options, and nonzero return statuses. At program exit it writes one receipt line to
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

extern "C" int __real_finufft2d1(int64_t,const double*,const double*,
    const std::complex<double>*,int,double,int64_t,int64_t,std::complex<double>*,const finufft_opts*);

/** \brief Record the requested options, then call the real finufft2d1 with identical arguments.
 *  \return The status returned by the real finufft2d1. */
extern "C" int __wrap_finufft2d1(int64_t nj,const double* x,const double* y,
    const std::complex<double>* c,int sign,double eps,int64_t nx,int64_t ny,
    std::complex<double>* out,const finufft_opts* opts) {
    {
        std::lock_guard<std::mutex> guard(observed.lock);
        ++observed.calls;
        if(opts){
            observed.lastNthreads=opts->nthreads;
            ++observed.histogram[opts->nthreads];
        }else{
            observed.lastNthreads=0;
            ++observed.nullOptions;
        }
    }
    const int status=__real_finufft2d1(nj,x,y,c,sign,eps,nx,ny,out,opts);
    if(status!=0){
        std::lock_guard<std::mutex> guard(observed.lock);
        ++observed.nonzeroStatus;
    }
    return status;
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
