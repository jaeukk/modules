/** \file finufft_observer.h
 * \brief Query interface of the observation-only FINUFFT guru-interface link wrapper (finufft_observer.cpp).
 * Only meaningful in programs linked with finufft_observer.o and
 * -Wl,--wrap=finufft_makeplan,--wrap=finufft_setpts,--wrap=finufft_execute. A call is one finufft_execute.
 */
#ifndef FINUFFT_OBSERVER_H
#define FINUFFT_OBSERVER_H

#include <cstddef>

namespace finufft_observer {
/** \brief Return the number of FINUFFT transforms (finufft_execute calls) observed so far.
 *  \return Total call count since program start. */
std::size_t Calls();

/** \brief Return the number of observed calls whose options requested a given thread count.
 *  \param nthreads  Value of finufft_opts::nthreads to count.
 *  \return Number of calls whose plan was made with opts->nthreads == nthreads (null options are not counted). */
std::size_t CallsWithNthreads(int nthreads);

/** \brief Return the opts->nthreads that the plan of the most recent call was made with.
 *  \return The last requested thread count, or 0 before any call or after a call with null options. */
int LastNthreads();
}

#endif
