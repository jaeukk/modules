/** \file Computation.cpp
 *  \brief Translation unit for the Computation task hierarchy.
 *
 *  In the legacy single-file design every Computation method body was written
 *  inline inside its class, so the definitions live in Computation.h.  This
 *  .cpp exists as the module's compilation unit (and the natural home for any
 *  future out-of-line method bodies); including the header is enough to type-
 *  check the whole hierarchy on its own.
 */

#include "Computation.h"
