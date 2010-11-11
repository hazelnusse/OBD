/* whipplesteadycalcs.h
 * 
 * Copyright (C) 2010 Dale Lukas Peterson
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef WHIPPLESTEADYBOUNDARIES_H
#define WHIPPLESTEADYBOUNDARIES_H
#define GSL_RANGE_CHECK_OFF
#include <string>
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
static void infspeed(gsl_vector * lean, gsl_vector * pitch,
                    double lean_ig, double pitch_ig, int ig_index,
                    const gsl_vector * steer, Whipple * bike);

// For all values of lean and steer within the region of feasible steady turns,
// calculate u5^2, then compute eigenvalues and eigenvectors for both forward
// and backwards speeds.
//static void steadyStability(gsl_vector * lean_zero, gsl_vector * pitch_zero,
//                            gsl_vector * lean_inf,  gsl_vector * pitch_inf,
//                            gsl_vector * steer,
//                       int ig_index, steadyOpts_t * options);


// Convenience functions used when solving the various nonlinear systems

// If gsl_multiroot_fdfsolver_iterate() returns non-zero, then it either
// encountered a nan or inf in the function or its derivative, or for some
// reason no progress is being made.  Currently this calls exit(1), and I have
// not encountered these error codes in practice (thankfully)
static void iterateError(int status, const char * routine);

// Increases ftol by an order of magnitude and subtracts 1 from i so that loop
// iteration is repeated using a greater error tolerance.  Also display a
// warning to stderr.
static void increaseftol(double * ftol, int * i, int iter_max,
                         const char * routine, double steer);

#endif
