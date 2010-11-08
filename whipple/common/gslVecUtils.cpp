/* gslVecUtils.cpp
 * 
 * Copyright (C) 2010 Dale Lukas Peterson
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef GSLVECUTILS_H
#define GSLVECUTILS_H

#include "gslVecUtils.h"

gsl_vector * linspaceN(double start, double stop, int N)
{
  gsl_vector * vec;
  if (N < 1)
    return NULL;
  else if (N == 1) {
    vec = gsl_vector_alloc(N);
    vec->data[0] = start;
  } else {
    int i;
    vec = gsl_vector_alloc(N);
    double delta = (stop - start) /  (double) (N - 1);

    for (i = 0; i < N; ++i)
      vec->data[i] = start + i*delta;
  }

  return vec;
} // linspace()

gsl_vector * linspace(double start, double stop, double delta)
{
  if (delta == 0.0)
    return NULL;
  int N = floor(fabs((start - stop) / delta)) + 1;
  gsl_vector * vec = gsl_vector_alloc(N);
  if (N == 1) {
    vec->data[0] = start;
  } else {
    int i;
    if (start < stop) 
      for (i = 0; i < N; ++i)
        vec->data[i] = start + i*delta;
    else
      for (i = 0; i < N; ++i)
        vec->data[i] = start - i*delta;
  }

  return vec;
}

// Returns 1 if all eigenvalues have negative real part, otherwise returns 0
int stable(gsl_vector_complex * evals)
{
  unsigned int i;
  for (i = 0; i < evals->size; ++i)
    if (!(GSL_REAL(gsl_vector_complex_get(evals, i)) < 0.0))
      return 0;
  return 1;
} // stable()

gsl_vector * zeros(unsigned long N)
{
  unsigned long i;
  gsl_vector * vec = gsl_vector_alloc(N);

  for (i = 0; i < N; ++i)
    vec->data[i] = 0.0;

  return vec;
} // zeros()

gsl_vector * ones(unsigned long N)
{
  unsigned long i;
  gsl_vector * vec = gsl_vector_alloc(N);

  for (i = 0; i < N; ++i)
    vec->data[i] = 1.0;

  return vec;
} // ones()
#endif
