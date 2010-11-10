#ifndef WHIPPLESTEADYBOUNDARIES_H
#define WHIPPLESTEADYBOUNDARIES_H
#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include "whipple.h"
#include "gslVecUtils.h"

// Functions for computation of static equilibrium boundary curve
static int static_f(const gsl_vector * x, void * params, gsl_vector * f);
static int static_df(const gsl_vector * x, void * params, gsl_matrix * J);
static int static_fdf(const gsl_vector * x, void * params, gsl_vector * f,
               gsl_matrix * J);
static int staticEq(gsl_vector * lean, gsl_vector * pitch,
                    const gsl_vector * steer, Whipple * bike);


// Functions for computation of configuration limit curves,
// i.e., the minimum and maximum lean curves
static int cfglim_f(const gsl_vector * x, void * params, gsl_vector * f);
static int cfglim_df(const gsl_vector * x, void * params, gsl_matrix * J);
static int cfglim_fdf(const gsl_vector * x, void * params, gsl_vector * f,
                   gsl_matrix * J);

static void cfglim(gsl_vector * lean_max, gsl_vector * pitch_max,
                gsl_vector * lean_min, gsl_vector * pitch_min,
                const gsl_vector * steer, Whipple * bike);

// Convenience functions
static void iterateError(int status, const char * routine);
static void increaseftol(double * ftol, int * i,
                         const char * routine, double steer);
#endif
