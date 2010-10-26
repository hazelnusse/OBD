#include "whipple.h"
#include "gslVecUtils.h"

// TODO:  add command line parser to parse for:
// speed range
// state       (must be a dynamic equilibrium)
// parameters

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();
  ofstream OutputFile("eigenvalues.dat", ios::binary);
  
  int N = 201;

  // Vector to store range of speeds to calculate eigenvalues
  gsl_vector * speed = linspaceN(0.0, 10.0, N);

  bb->u1 = 0.0;
  bb->u3 = 0.0;
  bb->rrt = bb->rft = 0.00;
  bb->evalConstants();

  for (int i = 0; i < N; ++i) {
    bb->u5 = -speed->data[i]/(bb->rf+bb->rft);
    bb->calcEvals();
    OutputFile.write((char *) &speed->data[i], sizeof(double));
    OutputFile.write((char *) bb->fourValues, 4*sizeof(double));
  } // for i

  OutputFile.close();
  delete bb;
  return 0;
} // main
