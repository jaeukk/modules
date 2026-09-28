/** \file finufft_observer.h
 * \brief Query interface of the observation-only finufft2d1 link wrapper (finufft_observer.cpp).
 * Only meaningful in programs linked with finufft_observer.o and -Wl,--wrap=finufft2d1.
 */
#ifndef FINUFFT_OBSERVER_H
#define FINUFFT_OBSERVER_H

#include <cstddef>

namespace finufft_observer {
/** \brief Return the number of finufft2d1 calls observed so far.
 *  \return Total call count since program start. */
std::size_t Calls();

/** \brief Return the number of observed calls whose options requested a given thread count.
 *  \param nthreads  Value of finufft_opts::nthreads to count.
 *  \return Number of calls with opts->nthreads == nthreads (null options are not counted). */
std::size_t CallsWithNthreads(int nthreads);

/** \brief Return opts->nthreads of the most recent call.
 *  \return The last requested thread count, or 0 before any call or after a call with null options. */
int LastNthreads();
}

#endif
