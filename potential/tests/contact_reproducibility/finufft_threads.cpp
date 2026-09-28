/** \file finufft_threads.cpp
 * \brief Diagnostic link wrapper changing only FINUFFT's internal thread count.
 * Link with --wrap=finufft2d1. The main harness and all other OpenMP teams remain
 * unchanged; no global OpenMP thread limit is lowered by this wrapper.
 */
#include <finufft.h>
#include <iostream>

#ifndef CC_TEST_FINUFFT_THREADS
#define CC_TEST_FINUFFT_THREADS 1
#endif
namespace {
/** \brief Record each requested thread count and return status for the receipt. */
struct ThreadAudit {
    size_t calls=0,unexpectedIncoming=0,nonzeroStatus=0;
    ~ThreadAudit() {
        std::cerr<<"FINUFFT_THREAD_CONTROL nthreads="<<CC_TEST_FINUFFT_THREADS
                 <<" calls="<<calls<<" unexpected_incoming="<<unexpectedIncoming
                 <<" nonzero_status="<<nonzeroStatus<<'\n';
    }
} audit;
}

extern "C" int __real_finufft2d1(int64_t,const double*,const double*,
    const std::complex<double>*,int,double,int64_t,int64_t,std::complex<double>*,const finufft_opts*);

/** \brief Copy options, override only nthreads, and call the real FINUFFT function. */
extern "C" int __wrap_finufft2d1(int64_t nj,const double* x,const double* y,
    const std::complex<double>* c,int sign,double eps,int64_t nx,int64_t ny,
    std::complex<double>* out,const finufft_opts* incoming) {
    finufft_opts opts;
    if(incoming)opts=*incoming;else finufft_default_opts(&opts);
    ++audit.calls;
    if(!incoming||incoming->nthreads!=4)++audit.unexpectedIncoming;
    opts.nthreads=CC_TEST_FINUFFT_THREADS;
    const int status=__real_finufft2d1(nj,x,y,c,sign,eps,nx,ny,out,&opts);
    if(status!=0)++audit.nonzeroStatus;
    return status;
}
