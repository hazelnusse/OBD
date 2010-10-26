#include <iostream>

#include "whipple.h"

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();
  ofstream OutputFile("simulation.dat", ios::binary);
  OutputFile << bb;
  
  double tj, state[10] = {bb->q0, bb->q1, bb->q3, bb->q4, bb->q5, bb->q6,
                          bb->q7, bb->u1, bb->u3, bb->u5};
  for (int j = 1; j < bb->fps*bb->tf + 1; ++j) {
    tj = ((double) j) / ((double) bb->fps);
    while (bb->t < tj)
      bb->evolve(tj, state);
    bb->computeOutputs(); // Compute outputs only at desired time points
    OutputFile << bb;     // Write quantities to file
  } // for j

  delete bb;
  OutputFile.close();
  cout << "Simulation completed." << endl;
  return 0;
} // main
