/* whipplesteady_origin.cpp
 *
 * Copyright (C) 2010-2011 Dale Lukas Peterson
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
#include <fstream>
#include "whipple.h"
#include "gslVecUtils.h"
#include "OBDConfig.h"

using namespace std;

struct record {
  double lean, steer, v, torque, ke, pe, te;
};
bool dominantreal(gsl_vector_complex * evals);

int main(int argc, char ** argv)
{
  Whipple * bike = new Whipple;

  bike->u1 = bike->u3 = 0.0;
  bike->u5 = -4.29238253634111/(bike->rf + bike->rft);
  bike->eoms();
  bike->computeOutputs();
  cout << "Weave total energy: " << bike->ke + bike->pe << "\n";
  bike->u5 = -6.02426201538837/(bike->rf + bike->rft);
  bike->eoms();
  bike->computeOutputs();
  cout << "Capsize total energy: " << bike->ke + bike->pe << "\n\n";


  if (argc != 7) {
    cout << "Must specify N_phi N_delta phi_min delta_min phi_max delta_max\n";
    return 1;
  }
  int N_phi, N_delta;
  double phi_min, delta_min, phi_max, delta_max;
  N_phi = atoi(argv[1]);
  N_delta = atoi(argv[2]);
  phi_min = atof(argv[3]);
  delta_min = atof(argv[4]);
  phi_max = atof(argv[5]);
  delta_max = atof(argv[6]);

  double dphi = (phi_max-phi_min) / double(N_phi-1);
  double ddelta = (delta_max-delta_min) / double(N_delta-1);
  int records = 0;
  struct record data;

  ofstream datafile("steadyturning.dat", ofstream::binary);
  ofstream stable_r("stablereal.dat", ofstream::binary);
  ofstream stable_i("stableimag.dat", ofstream::binary);

  bike->q1=bike->q3=0;
  bike->calcPitch();
  bike->steadyEqns();
  cout.precision(16);

  cout << "dphi = " << dphi << "\n";
  cout << "ddelta = " << ddelta << "\n";

  for (int j = 0; j < N_delta; ++j) {
    bike->q3 = delta_min + j*ddelta;  // Set steer
    for (int i = 0; i < N_phi; ++i) {
      bike->q1 = phi_min + i*dphi;      // Set lean
      bike->calcPitch();    // Calculate pitch
      bike->steadyEqns();   // Evaluate steady eoms

      // Write the turn conditions to record struct
      data.lean =  bike->q1*180.0/M_PI;      // lean
      data.steer = bike->q3*180.0/M_PI;      // steer
      data.v = pow(abs(bike->u5s_s), 0.5)*(bike->rf + bike->rft);   // fw rate squared
      ++records;

      // Set applied torque to steady turning torque
      bike->Ts = data.torque = bike->Ts_s;
      // Set fw rate to steady front wheel rate
      bike->u5 = -sqrt(bike->u5s_s);
      bike->eoms();             // evaluate motion equations
      bike->computeOutputs();   // compute output quantities (ke, pe, te);
      data.ke = bike->ke;
      data.pe = bike->pe;
      data.te = bike->ke + bike->pe;

      bike->calcEvals();
      if (stable(bike->evals)) {
        if (dominantreal(bike->evals)) {
          stable_r.write((char *) &bike->q1, sizeof(double));
          stable_r.write((char *) &bike->q3, sizeof(double));
        } else {
          stable_i.write((char *) &bike->q1, sizeof(double));
          stable_i.write((char *) &bike->q3, sizeof(double));
        }
      }      // Write one record to file
      datafile.write((char *)&data, sizeof(record));
    } // for i
  } // for j

  cout << "Points processed: " << records << "\n";

  bike->q1 = 5.0*M_PI/180.;
  for (int j = 0;; ++j) {
    bike->q3 = delta_min + j*ddelta;
    if (bike->q3 > M_PI)
      break;
    bike->calcPitch();
    bike->steadyEqns();
    if (bike->Ts_s > 0) {
      bike->Ts = bike->Ts_s;
      bike->u5 = -sqrt(bike->u5s_s);
      bike->calcEvals();
      if (stable(bike->evals)) {
        bike->printCfgCon();
        cout << "Steer torque = " << bike->Ts_s << "\n"
             << "Front wheel rate = " << bike->u5 << "\n"
             << "Velocity = " << -bike->u5*(bike->rf+bike->rft) << "\n";
        bike->printEvals();
        cout << "\n";
      }
    }
  }

  datafile.close();
  stable_r.close();
  stable_i.close();
  delete bike;
  return 0;
} // main

bool dominantreal(gsl_vector_complex * evals)
{
  double realpart = -DBL_MAX;
  int index = 0;
  // find the eigenvalue closest to the j-w axis
  for (int i = 0; i < 4; ++i)
    if (GSL_REAL(gsl_vector_complex_get(evals, i)) > realpart) {
      index = i;
      realpart = GSL_REAL(gsl_vector_complex_get(evals, i));
    }

  if (GSL_IMAG(gsl_vector_complex_get(evals, index)))
    return false;
  return true;
}
