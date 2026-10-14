/** \file finufft_threads.cpp
 * \brief Diagnostic link wrapper changing only FINUFFT's internal thread count.
 * Link with -Wl,--wrap=finufft2d1,--wrap=finufft_makeplan,--wrap=finufft_setpts,--wrap=finufft_execute.
 * Objects built before persistent plans call finufft2d1; current ShiftedCCPotential objects
 * make one plan (finufft_makeplan, where the options arrive) and run one finufft_execute per
 * transform. Both paths are covered: a call is one finufft2d1 or one finufft_execute. The main
 * harness and all other OpenMP teams remain unchanged; no global OpenMP thread limit is lowered
 * by this wrapper.
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

extern "C" int __real_finufft_makeplan(int,int,const int64_t*,int,int,double,finufft_plan*,const finufft_opts*);
extern "C" int __real_finufft_setpts(finufft_plan,int64_t,const double*,const double*,const double*,
    int64_t,const double*,const double*,const double*);
extern "C" int __real_finufft_execute(finufft_plan,std::complex<double>*,std::complex<double>*);

/** \brief Copy options, check the incoming request, override only nthreads, and make the real plan.
 *  \return The status returned by the real finufft_makeplan. */
extern "C" int __wrap_finufft_makeplan(int type,int dim,const int64_t* nmodes,int iflag,int ntrans,double tol,
    finufft_plan* plan,const finufft_opts* incoming) {
    finufft_opts opts;
    if(incoming)opts=*incoming;else finufft_default_opts(&opts);
    if(!incoming||incoming->nthreads!=4)++audit.unexpectedIncoming;
    opts.nthreads=CC_TEST_FINUFFT_THREADS;
    const int status=__real_finufft_makeplan(type,dim,nmodes,iflag,ntrans,tol,plan,&opts);
    if(status!=0)++audit.nonzeroStatus;
    return status;
}

/** \brief Call the real finufft_setpts unchanged, recording a nonzero status.
 *  \return The status returned by the real finufft_setpts. */
extern "C" int __wrap_finufft_setpts(finufft_plan plan,int64_t m,const double* x,const double* y,const double* z,
    int64_t n,const double* s,const double* t,const double* u) {
    const int status=__real_finufft_setpts(plan,m,x,y,z,n,s,t,u);
    if(status!=0)++audit.nonzeroStatus;
    return status;
}

/** \brief Count one transform of a plan made through the wrapper, then run it unchanged.
 *  \return The status returned by the real finufft_execute. */
extern "C" int __wrap_finufft_execute(finufft_plan plan,std::complex<double>* weights,std::complex<double>* result) {
    ++audit.calls;
    const int status=__real_finufft_execute(plan,weights,result);
    if(status!=0)++audit.nonzeroStatus;
    return status;
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
