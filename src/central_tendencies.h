// central_tendencies.h
#ifndef CENTRAL_TENDENCIES_H
#define CENTRAL_TENDENCIES_H

#include <Rcpp.h>


/// Compute the "center of gravity" statistic used by NNS.
///
/// @param xSEXP Input vector supplied from R.
/// @param discrete Whether to coerce the result to the discrete analogue used
///   by the package's discrete workflow.
/// @return A scalar SEXP containing the estimated center of gravity.
SEXP NNS_gravity_cpp(SEXP xSEXP, bool discrete);

/// Compute the mode (or modal class) depending on the supplied flags.
///
/// @param xSEXP Input vector supplied from R.
/// @param discrete Treat data as discrete values.
/// @param multi Return the multi-modal result when requested from the R layer.
/// @return A scalar or vector SEXP mirroring the behaviour of the R-facing
///   wrapper.
SEXP NNS_mode_cpp(SEXP xSEXP, bool discrete, bool multi);

#endif // CENTRAL_TENDENCIES_H