/* whipplesteady.cpp
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
#include <iostream>
#include <string>
#include <getopt.h>
#include "whippleutils.h"
#include "whipple.h"

// Forward declaration
void processOptions(int argc, char ** argv, steadyOpts_t * options, Whipple * bike);
void parseisoargs(gsl_vector * v, char * optarg);

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();
  steadyOpts_t * opt = new steadyOpts_t;
  // set default options
  opt->outfolder[0] = '\0'; // Put results in current directory
  opt->all = false;         // Only calculate boundary curves by default
  opt->N = 1001;            // # of points to mesh steer in [0, pi]
  
  // Process command line options
  processOptions(argc, argv, opt, bb);

  // Make all the calculations requested
  bb->steadyCalcs(opt);

  cout << "Steady turning analysis completed.  Steady turning output written to "
       << opt->outfolder << '\n';
  // Free memory, close files
  delete bb;
  delete opt;
  return 0;
} // main

void processOptions(int argc, char ** argv, 
                    steadyOpts_t * opt, Whipple * bike)
{
  int c, option_index;
  bool verbose_flag = false;
  while (1) {
    static struct option long_options[] = {
      {"help",          no_argument,       0, 'h'},
      {"all",           no_argument,       0, 'a'},
      {"steerpoints",   required_argument, 0, 'N'},
      {"iso_v",         required_argument, 0, 's'},
      {"iso_t",         required_argument, 0, 't'},
      {"iso_mew",       required_argument, 0, 'f'},
      {"bmparams",      required_argument, 0, 'm'},
      {"params",        required_argument, 0, 'p'},
      {"output",        required_argument, 0, 'o'},
      {"verbose",       no_argument,       0, 'v'},
      {0, 0, 0, 0}};

    c = getopt_long(argc, argv, "ha:N:s:t:f:m:p:o:v", long_options, &option_index);

    if (c == -1) //Detect the end of the options.
      break;
    
    if (c == 'h') {
      cout << 
"usage: " << argv[0] << " [OPTION]\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n\n"
"  -h, --help                     display this help and exit\n"
"  -m, --mjparams=pfile           Meijaard bike parameters\n"
"  -p, --params=pfile             native bike model parameters\n"
"  -N, --steerpoints=N            number of points to mesh steer into\n"
"  -a, --all                      all feasible steady turning quantities\n"
"  -s, --iso_v=1,2,3,...          specify velocity level curves\n"
"  -t, --iso_t=1,2,3,...          specify torque level curves\n"
"  -f, --iso_mew=1,2,3,...        specify friction coefficient level curves\n"
"  -o, --output=folder            specify to store output data\n"
"  -v, --verbose                  display bike parameters used\n\n";
      exit(0);
    } else if (c == 'm') {
      MJWhippleParams * mjbike = new MJWhippleParams;
      WhippleParams * b = new WhippleParams;
      readMJWhippleParams(mjbike, optarg);
      convertParameters(b, mjbike);
      bike->setParameters(b);
      bike->evalConstants();
      bike->eoms();
      bike->computeOutputs();
      bike->steadyEqns();
      delete mjbike;
      delete b;
    } else if (c == 'p') {
      WhippleParams * b = new WhippleParams;
      readWhippleParams(b, optarg);
      bike->setParameters(b);
      bike->evalConstants();
      bike->eoms();
      bike->computeOutputs();
      bike->steadyEqns();
      delete b;
    } else if (c == 'a')
      opt->all = true;
    else if (c == 'N')
      opt->N = atoi(optarg);
    else if (c == 'o')
      strcpy(opt->outfolder, optarg);
    else if (c == 't')  // torque level curves specified
      parseisoargs(opt->iso_t, optarg);
    else if (c == 'f')  // friction coefficient level curves
      parseisoargs(opt->iso_mew, optarg);
    else if (c == 's') // velocity level curves
      parseisoargs(opt->iso_v, optarg);
    else if (c == 'v')  // verbose flag
      verbose_flag = true;
    else {
      cout << "Invalid option.\n";
      abort();
    }
  } // while()
  
  if (verbose_flag) {
    bike->printParameters();
    bike->printState();
  }
} // processOptions()
      
void parseisoargs(gsl_vector * v, char * optarg)
{
    double * iso = new double[30]; // up to 30 level curves
    char * l = strtok(optarg, ",");
    size_t i;
    for (i = 0; l; l = strtok(NULL, ","), ++i)
      iso[i] = atof(l);
    v = gsl_vector_alloc(i);
    for (i = 0; i < v->size; ++i)
      gsl_vector_set(v, i, iso[i]);
    delete [] iso;
} // parseisoargs()
