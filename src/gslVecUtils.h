/*
 * =====================================================================================
 *
 *       Filename:  gsl_vector_utils.h
 *
 *    Description:  Headers for convenient day to day vector operations
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California Davis
 *
 * =====================================================================================
 */

#include <gsl/gsl_vector.h>

gsl_vector * linspaceN(double start, double stop, int N);

gsl_vector * linspace(double start, double stop, double delta);

int stable(gsl_vector_complex * evals);

gsl_vector * zeros(unsigned long N);

gsl_vector * ones(unsigned long N);
