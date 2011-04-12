#include <iostream>
#include <string>
#include <getopt.h>
#include "whippleutils.h"
#include "whipple.h"
#include "gslVecUtils.h"

#define N_MESH 101

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

  for (int i = 0; i < N_MESH; ++i) {
    lean = lean_min + i*delta_l;
    for (int j = 0; j < N_MESH; ++j) {
      steer = steer_min + j*delta_s;
      bb->q1 = lean;
      bb->q3 = steer;
      bb->calcPitch();
      bb->steadyEqns();

      if (bb->u5s_s < 0.0)
        continue;

      u5 = -sqrt(bb->u5s_s);

      // set state to positive rotation rate, compute e-vals
      bb->u5 = u5;
      bb->calcEvals();
      std::cout.precision(16);
      if (stable(bb->evals)) {
        if (bb->Ts_s >= 0.0) {
          std::cout << "Stable:\nlean = " << lean << "\n"
               << "steer = " <<  steer << "\n"
               << "steer torque = " << bb->Ts_s << "\n"
               << "speed = " << -u5*bb->rf << "\n";
          bb->printEvals();
        }
      } //
    } // for j
  } // for i


  delete bb;
  return 0;
} // main
