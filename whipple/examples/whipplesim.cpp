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

// Forward declarations
void processOptions(int argc,
                    char ** argv,
                    char * outfolder,
                    Whipple * bike);
std::ostream &operator<<(std::ostream &outfile, const Whipple * discs);
void writeSimRecord_dt(const char * filename);

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();
  char outfolder[512] = "";
  std::string filename;

  // Process command line options
  processOptions(argc, argv, outfolder, bb);

  // Write parameters
  filename = outfolder; filename += "simulation_parameters.txt";
  bb->writeParameters(filename.c_str());

  // Write data record file.
  filename = outfolder; filename += "sim_record.py";
  writeSimRecord_dt(filename.c_str());

  // Open data file
  filename = outfolder; filename += "simulation.dat";
  std::ofstream OutputFile(filename.c_str(), std::ios::binary);

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
  std::cout << "Simulation completed.  Simulation output written to "
       << outfolder << '\n';
  return 0;
} // main

/**
 * Process command line options.
*/
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
      std::cout <<
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
      std::cout << "Invalid option." << '\n';
      abort();
    }
  } // while()

  if (verbose_flag) {
    bike->printParameters();
    bike->printState();
  }
} // processOptions()

std::ostream &operator<<(std::ostream &file, const Whipple * discs)
{
  file.write((char *) &(discs->t), sizeof discs->t);
  file.write((char *) &discs->q0, sizeof discs->q0);
  file.write((char *) &discs->q1, sizeof discs->q1);
  file.write((char *) &discs->q2, sizeof discs->q2);
  file.write((char *) &discs->q3, sizeof discs->q3);
  file.write((char *) &discs->q4, sizeof discs->q4);
  file.write((char *) &discs->q5, sizeof discs->q5);
  file.write((char *) &discs->q6, sizeof discs->q6);
  file.write((char *) &discs->q7, sizeof discs->q7);
  file.write((char *) &discs->u0, sizeof discs->u0);
  file.write((char *) &discs->u1, sizeof discs->u1);
  file.write((char *) &discs->u2, sizeof discs->u2);
  file.write((char *) &discs->u3, sizeof discs->u3);
  file.write((char *) &discs->u4, sizeof discs->u4);
  file.write((char *) &discs->u5, sizeof discs->u5);
  file.write((char *) discs->no_fn, sizeof discs->no_fn);
  file.write((char *) &discs->Rx, sizeof discs->Rx);
  file.write((char *) &discs->Ry, sizeof discs->Ry);
  file.write((char *) &discs->Rz, sizeof discs->Rz);
  file.write((char *) &discs->Fx, sizeof discs->Fx);
  file.write((char *) &discs->Fy, sizeof discs->Fy);
  file.write((char *) &discs->Fz, sizeof discs->Fz);
  file.write((char *) &discs->ke, sizeof discs->ke);
  file.write((char *) &discs->pe, sizeof discs->pe);
  file.write((char *) &discs->fa_yaw, sizeof discs->fa_yaw);
  file.write((char *) &discs->fa_lean, sizeof discs->fa_lean);
  file.write((char *) &discs->fa_pitch, sizeof discs->fa_pitch);
  file.write((char *) discs->constraints, sizeof discs->constraints);
  return file;
} // ostream &operator<<()

/**
 * Write a Python file which defines a Numpy data type compatible with binary simulation data
 * data type to enable easy manipulation of simulation results in Python
 *
 * @pre filename should be a string with a valid filename
 * @post Python module written to filename
*/
void writeSimRecord_dt(const char * filename)
{
  std::ofstream fp(filename, std::ios::out);
  if (fp.is_open()) {
    fp << "import numpy as np\n";
    fp << "sim_dt = np.dtype([('t', np.float64),\n"
          "                   ('q0', np.float64),\n"
          "                   ('q1', np.float64),\n"
          "                   ('q2', np.float64),\n"
          "                   ('q3', np.float64),\n"
          "                   ('q4', np.float64),\n"
          "                   ('q5', np.float64),\n"
          "                   ('q6', np.float64),\n"
          "                   ('q7', np.float64),\n"
          "                   ('u0', np.float64),\n"
          "                   ('u1', np.float64),\n"
          "                   ('u2', np.float64),\n"
          "                   ('u3', np.float64),\n"
          "                   ('u4', np.float64),\n"
          "                   ('u5', np.float64),\n"
          "                   ('fnx', np.float64),\n"
          "                   ('fny', np.float64),\n"
          "                   ('fnz', np.float64),\n"
          "                   ('Rx', np.float64),\n"
          "                   ('Ry', np.float64),\n"
          "                   ('Rz', np.float64),\n"
          "                   ('Fx', np.float64),\n"
          "                   ('Fy', np.float64),\n"
          "                   ('Fz', np.float64),\n"
          "                   ('ke', np.float64),\n"
          "                   ('pe', np.float64),\n"
          "                   ('fa_yaw', np.float64),\n"
          "                   ('fa_lean', np.float64),\n"
          "                   ('fa_pitch', np.float64),\n"
          "                   ('nh1', np.float64),\n"
          "                   ('nh2', np.float64),\n"
          "                   ('nh3', np.float64)])";
    fp.close();
  } else {
    std::cerr << "Unable to open " << filename << "for writing.\n";
    std::cerr << "Aborting.\n";
    exit(0);
  }
} // writeSimRecord_dt()
