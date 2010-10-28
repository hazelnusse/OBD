#ifndef WHIPPLEUTILS_H
#define WHIPPLEUTILS_H

//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv.h>
#include <fstream>

//#include <math.h>
#include <string.h>
//#include <getopt.h>
#include "whipple.h"

// Convert Meijaard Bike Parameters to my model parameters
void convertParameters(WhippleParams *bout, const MJWhippleParams * bin);

// Initialize all parameters to those presented in benchmark paper
void setBenchmarkParameters(MJWhippleParams * bike);

// Read a text file to get the Meijaard parameters. If a parameter is not
// specified in the text file, the corresponding benchmark parameter is used.
void readMJWhippleParams(MJWhippleParams *mjbike, const char *filename);

// Read a text file to get the native parameters. If a parameter is not
// specified in the text file, the corresponding benchmark parameter is used.
void readWhippleParams(WhippleParams * bike, const char * filename);

// Read the state variables from file set state of Whipple object
void readState(double * state, const char * filename);

// Read the integration parameters from a text file
//void readIntegrationParams(BikeParams *bike, const char *filename);
#endif
