/*
 * =====================================================================================
 *
 *       Filename:  gsl_vector_utils.c
 *
 *    Description:  Functions for convenient day to day vector operations
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California Davis
 *
 * =====================================================================================
 */

#ifndef GSLVECUTILS_H
#define GSLVECUTILS_H

#include "gslVecUtils.h"
#include <cmath>

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
