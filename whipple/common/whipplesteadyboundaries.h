#ifndef WHIPPLESTEADYBOUNDARIES_H
#define WHIPPLESTEADYBOUNDARIES_H
#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include "whipple.h"
#include "gslVecUtils.h"

// Functions for computation of static equilibrium boundary curve
// F:  [kinematic constraint; coefficient of g in steady turning eom]
static int static_f(const gsl_vector * x, void * params, gsl_vector * f);
static int static_df(const gsl_vector * x, void * params, gsl_matrix * J);
static int static_fdf(const gsl_vector * x, void * params,
                      gsl_vector * f, gsl_matrix * J);
// staticEq calculates the lean and pitch values of static equilibrium
// for each value of steer specified.
static int staticEq(gsl_vector * lean, gsl_vector * pitch,
                    const gsl_vector * steer, Whipple * bike);

// Functions for computation of configuration limit curves
// F:  [kinematic constraint; partial of kinematic constraint w.r.t. pitch]
static int cfglim_f(const gsl_vector * x, void * params, gsl_vector * f);
static int cfglim_df(const gsl_vector * x, void * params, gsl_matrix * J);
static int cfglim_fdf(const gsl_vector * x, void * params,
                      gsl_vector * f, gsl_matrix * J);

// staticEq calculates the lean and pitch values of static equilibrium
// for each value of steer specified.
static void cfglim(gsl_vector * lean_max, gsl_vector * pitch_max,
                   gsl_vector * lean_min, gsl_vector * pitch_min,
                   const gsl_vector * steer, Whipple * bike);

// F:  [kinematic constraint; coefficient of u5^2 in steady turning eom]
static int inf_f(const gsl_vector * x, void * params, gsl_vector * f);
static int inf_df(const gsl_vector * x, void * params, gsl_matrix * J);
static int inf_fdf(const gsl_vector * x, void * params,
                   gsl_vector * f, gsl_matrix * J);

// Give the steer discretized in the range of [0, pi], solve for the infinite
// speed steady turning lean.  Initial guess for the lean and pitch must be
// supplied, along with the index of the steer vector that those initial guess
// correspond to.  Ideally, these quantities are provide from the index that is
// returned by StaticEq(), along with the corresponding lean and pitch values
// at that index.
static int infspeed(gsl_vector * lean, gsl_vector * pitch,
                    double lean_ig, double pitch_ig, int ig_index,
                    const gsl_vector * steer, Whipple * bike);

// Convenience functions used when solving the various nonlinear systems

// If gsl_multiroot_fdfsolver_iterate() returns non-zero, then it either
// encountered a nan or inf in the function or its derivative, or for some
// reason no progress is being made.  Currently this calls exit(1), and I have
// not encountered these error codes in practice (thankfully)
static void iterateError(int status, const char * routine);

// Increases ftol by an order of magnitude and subtracts 1 from i so that loop
// iteration is repeated using a greater error tolerance.  Also display a
// warning to stderr.
static void increaseftol(double * ftol, int * i,
                         const char * routine, double steer);

#endif
