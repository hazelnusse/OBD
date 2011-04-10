/* whippleutils.h
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
#ifndef WHIPPLEUTILS_H
#define WHIPPLEUTILS_H

#include <fstream>
#include <cstring>
#include <gsl/gsl_odeiv.h>
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
void readState(double * bike, const char * filename);

// Read the integration parameters from a text file
//void readIntegrationParams(BikeParams *bike, const char *filename);
#endif
