#include <iostream>
#include <string>
#include <getopt.h>
#include "whippleutils.h"
#include "whipple.h"
#include "gslVecUtils.h"

#define N_MESH 21

int main(void)
{
  Whipple * bb = new Whipple();

  double lean_min = M_PI/9.0 + .09, lean_max = M_PI/6.0 - 0.02,
         steer_min = M_PI/9.0 + .09, steer_max = M_PI/6.0 - .02,
         delta_l = (lean_max - lean_min) / (N_MESH - 1),
         delta_s = (steer_max - steer_min) / (N_MESH - 1),
         lean, steer;
  double u5;

  bb->u1 = bb->u3 = bb->u5 = 0.0;

  std::cout.precision(16);
  for (int i = 0; i < N_MESH; ++i) {
    lean = lean_min + i*delta_l;
    for (int j = 0; j < N_MESH; ++j) {
      steer = steer_min + j*delta_s;
      bb->q1 = lean;
      bb->q3 = steer;
      bb->calcPitch();
      bb->steadyEqns(); // Deterimine equilibrium speed and steer torque

      // if u5^2 < 0.0, not a feasible turn
      if (bb->u5s_s < 0.0)
        continue;

      // Set front wheel rate
      u5 = sqrt(bb->u5s_s);
      // Set applied steer torque
      bb->Ts = bb->Ts_s;

      // set state to positive rotation rate, compute e-vals
      bb->u5 = u5;
      bb->calcEvals();
      if (stable(bb->evals)) {
        if (bb->Ts_s >= 0.0) {
          std::cout << "Stable steady turn:\n"
               << "lean = " << lean << "\n"
               << "pitch = " << bb->q2 << "\n"
               << "steer = " <<  steer << "\n"
               << "steer torque = " << bb->Ts_s << "\n"
               << "fw rate = " << bb->u5 << "\n"
               << "v = " << -bb->u5*bb->rf << "\n"
               << "fw contact height = " << hc_f(bb->q2, bb) << "\n";
          bb->printEvals();
          std::cout << "\n";
        }
      } //

      // set state to negative rotation rate, compute e-vals
      bb->u5 = -u5;
      bb->calcEvals();
      if (stable(bb->evals)) {
        if (bb->Ts_s >= 0.0) {
          std::cout << "Stable steady turn:\n"
               << "lean = " << lean << "\n"
               << "pitch = " << bb->q2 << "\n"
               << "steer = " <<  steer << "\n"
               << "steer torque = " << bb->Ts_s << "\n"
               << "fw rate = " << bb->u5 << "\n"
               << "v = " << -bb->u5 * bb->rf << "\n"
               << "fw contact height = " << hc_f(bb->q2, bb) << "\n";
          bb->printEvals();
          std::cout << "\n";
        }
      } //
    } // for j
  } // for i

  delete bb;
  return 0;
} // main
