/* whipplesim.cpp
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
#include "OBDConfig.h"

// Forward declaration
void processOptions(int argc, char ** argv, char * outfolder, Whipple * bike);

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();
  char outfolder[512] = "";
  string filename;

  // Process command line options
  processOptions(argc, argv, outfolder, bb);

  // Write parameters
  filename = outfolder; filename += "simulation_parameters.txt";
  bb->writeParameters(filename.c_str());
  // Write data record file.
  filename = outfolder; filename += "sim_record.py";
  bb->writeSimRecord_dt(filename.c_str());
  // Open data file
  filename = outfolder; filename += "simulation.dat";
  ofstream OutputFile(filename.c_str(), ios::binary);

  // Write initial conditions
  OutputFile << bb;
  double tj, state[10] = {bb->q0, bb->q1, bb->q3, bb->q4, bb->q5, bb->q6,
                          bb->q7, bb->u1, bb->u3, bb->u5};

  // Perform numerical integration
  for (int j = 1; j < bb->fps*bb->tf + 1; ++j) {
    tj = ((double) j) / ((double) bb->fps);
    while (bb->t < tj)
      bb->evolve(tj, state);
    bb->computeOutputs(); // Compute outputs only at desired time points
    OutputFile << bb;     // Write quantities to file
  } // for j

  // Free memory, close files
  delete bb;
  OutputFile.close();
  cout << "Simulation completed.  Simulation output written to "
       << outfolder << endl;
  return 0;
} // main

void processOptions(int argc, char ** argv, char * outfolder, Whipple * bike)
{
  int c, option_index;
  bool verbose_flag = false;
  while (1) {
    static struct option long_options[] = {
      {"help",          no_argument,       0, 'h'},
      {"mjparams",      required_argument, 0, 'm'},
      {"params",        required_argument, 0, 'p'},
      {"state",         required_argument, 0, 's'},
      {"pitch_ig",      required_argument, 0, 'g'},
      {"output",        required_argument, 0, 'o'},
      {"tf",            required_argument, 0, 't'},
      {"fps",           required_argument, 0, 'f'},
      {"verbose",       no_argument,       0, 'v'},
      {0, 0, 0, 0}};

    c = getopt_long(argc, argv, "hm:p:s:g:o:t:f:v", long_options, &option_index);

    if (c == -1) //Detect the end of the options.
      break;

    if (c == 'h') {
      cout <<
argv[0] << " Version " << OBD_VERSION_MAJOR << "." << OBD_VERSION_MINOR <<
" commit " << OBD_VERSION_COMMIT << "\n"
"usage: " << argv[0] << " [OPTION]\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n\n"
"  -h, --help                         display this help and exit.\n"
"  -m, --mjparams=pfile               Meijaard bike parameters\n"
"  -p, --params=pfile                 native bike model parameters\n"
"  -g, --pitch_ig=FP_NUMBER           initial guess for the pitch root finder\n"
"  -o, --output=folder                write eigenvalues to folder\n"
"  -s, --state=statefile              file to specify initial simulation conditions\n"
"  -t, --tf=FP_NUMBER                 simulation time\n"
"  -f, --fps=FP_NUMBER                output data time interval\n"
"  -v, --verbose                      simulation time\n\n";
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
      delete mjbike;
      delete b;
    } else if (c == 'p') {
      WhippleParams * b = new WhippleParams;
      readWhippleParams(b, optarg);
      bike->setParameters(b);
      bike->evalConstants();
      bike->eoms();
      bike->computeOutputs();
      delete b;
    } else if (c == 'g') {
      bike->q2 = atof(optarg);
      bike->calcPitch();
      bike->eoms();
      bike->computeOutputs();
    } else if (c == 'o')
      strcpy(outfolder, optarg);
    else if (c == 's') {
      double * state = new double[10];
      readState(state, optarg);
      bike->setState(state);
      bike->eoms();
      bike->computeOutputs();
      delete [] state;
    } else if (c == 't')
      bike->tf = atof(optarg);
    else if (c == 'f')
      bike->fps = atof(optarg);
    else if (c == 'v')
      verbose_flag = true;
    else {
      cout << "Invalid option." << endl;
      abort();
    }
  } // while()

  if (verbose_flag) {
    bike->printParameters();
    bike->printState();
  }
} // processOptions()
