/* gslVecUtils.cpp
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
#ifndef GSLVECUTILS_H
#define GSLVECUTILS_H
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE

#include "gslVecUtils.h"

gsl_vector * linspaceN(double start, double stop, size_t N)
{
  gsl_vector * vec;
  if (N < 1)
    return NULL;
  else if (N == 1) {
    vec = gsl_vector_alloc(N);
    gsl_vector_set(vec, 0, start);
  } else {
    vec = gsl_vector_alloc(N);
    double delta = (stop - start) /  (double) (N - 1);

    for (size_t i = 0; i < N; ++i)
      gsl_vector_set(vec, i, start + i*delta);
  }

  return vec;
} // linspace()

gsl_vector * linspace(double start, double stop, double delta)
{
  if (delta == 0.0)
    return NULL;
  size_t N = floor(fabs((start - stop) / delta)) + 1;
  gsl_vector * vec = gsl_vector_alloc(N);
  if (N == 1) {
    gsl_vector_set(vec, 0, start);
  } else {
    if (start < stop) 
      for (size_t i = 0; i < N; ++i)
        gsl_vector_set(vec, i, start + i*delta);
    else
      for (size_t i = 0; i < N; ++i)
        gsl_vector_set(vec, i, start - i*delta);
  }

  return vec;
}

// Returns 1 if all eigenvalues have negative real part, otherwise returns 0
int stable(gsl_vector_complex * evals)
{
  for (size_t i = 0; i < evals->size; ++i)
    if (!(GSL_REAL(gsl_vector_complex_get(evals, i)) < 0.0))
      return 0;
  return 1;
} // stable()

gsl_vector * zeros(size_t N)
{
  gsl_vector * vec = gsl_vector_alloc(N);

  for (size_t i = 0; i < N; ++i)
    gsl_vector_set(vec, i, 0.0);

  return vec;
} // zeros()

gsl_vector * ones(size_t N)
{
  gsl_vector * vec = gsl_vector_alloc(N);

  for (size_t i = 0; i < N; ++i)
    gsl_vector_set(vec, i, 1.0);

  return vec;
} // ones()
#endif
