/* whipple.cpp
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

#include "whipple.h"
#include "whippleutils.h"
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE

int eomwrapper(double t, const double x[10], double f[10], void * params)
{
  Whipple * p = (Whipple *) params;
  // Assign the states of the Whipple object
  p->setState(x);
  // Evaluate the RHS of the ODE's representing the equations of motion
  p->eoms();
  // Assign the right hand sides of the ODE's to the array passed
  f[0] = p->q0p;
  f[1] = p->q1p;
  f[2] = p->q3p;
  f[3] = p->q4p;
  f[4] = p->q5p;
  f[5] = p->q6p;
  f[6] = p->q7p;
  f[7] = p->u1p;
  f[8] = p->u3p;
  f[9] = p->u5p;
  // Return the status
  return GSL_SUCCESS;
} // eomwrapper()

double hc_f(double q2, void * params)
{
  Whipple * p = (Whipple *) params;
  return
      p->rft + p->ls*cos(p->q1)*cos(q2) +
      p->lf*(sin(p->q1)*sin(p->q3)-sin(q2)*cos(p->q1)*cos(p->q3)) +
      p->rf*sqrt((1-pow((sin(p->q1)*cos(p->q3)+sin(q2)*sin(p->q3)*cos(p->q1)),2))) -
      p->rrt - p->rr*cos(p->q1) - p->lr*sin(q2)*cos(p->q1);
} // hc_f()

double hc_df(double q2, void * params)
{
  Whipple * p = (Whipple *) params;
  return
      -cos(p->q1)*(p->lr*cos(q2)+p->ls*sin(q2)+p->lf*cos(q2)*cos(p->q3)+p->rf*sin(p->q3)*cos(q2)*(
                  sin(p->q1)*cos(p->q3)+sin(q2)*sin(p->q3)*cos(p->q1))/sqrt((1-pow((sin(p->q1)*cos(p->q3)+sin(q2)*
                              sin(p->q3)*cos(p->q1)),2))));
} // hc_df()

void hc_fdf(double q2, void * params, double * f, double * df)
{
  *f = hc_f(q2, params);
  *df =  hc_df(q2, params);
} // hc_fdf()

void Whipple::printState(void) const
{
  std::cout.precision(16);
  std::cout << "q0 = " << q0 << " (ignorable)\n"
       << "q1 = " << q1 << '\n'
       << "q2 = " << q2 << " (dependent)\n"
       << "q3 = " << q3 << '\n'
       << "q4 = " << q4 << " (ignorable)\n"
       << "q5 = " << q5 << " (ignorable)\n"
       << "q6 = " << q6 << " (ignorable)\n"
       << "q7 = " << q7 << " (ignorable)\n"
       << "u0 = " << u0 << " (dependent)\n"
       << "u1 = " << u1 << '\n'
       << "u2 = " << u2 << " (dependent)\n"
       << "u3 = " << u3 << '\n'
       << "u4 = " << u4 << " (dependent)\n"
       << "u5 = " << u5 << '\n';
} // printState()

void Whipple::printParameters(void) const
{
  std::cout.precision(16);
  std::cout << "rr   = " << rr
       << "\nrrt  = " << rrt
       << "\nrf   = " << rf
       << "\nrft  = " << rft
       << "\nlr   = " << lr
       << "\nls   = " << ls
       << "\nlf   = " << lf
       << "\nmr   = " << mr
       << "\nmf   = " << mf
       << "\nICyy = " << ICyy
       << "\nIDxx = " << IDxx
       << "\nIDyy = " << IDyy
       << "\nIDzz = " << IDzz
       << "\nIDxz = " << IDxz
       << "\nIExx = " << IExx
       << "\nIEyy = " << IEyy
       << "\nIEzz = " << IEzz
       << "\nIExz = " << IExz
       << "\nIFyy = " << IFyy
       << "\nlrx  = " << lrx
       << "\nlrz  = " << lrz
       << "\nlfx  = " << lfx
       << "\nlfz  = " << lfz
       << "\ng    = " << g << '\n';
} // printParameters()

void Whipple::writeState(const char * filename) const
{
  std::ofstream fp(filename, std::ios::out);
  fp.precision(16);
  if (fp.is_open()) {
  fp << "q0 = " << q0 << " (ignorable coordinate)" << '\n'
     << "q1 = " << q1 << '\n'
     << "q2 = " << q2 << " (dependent coordinate)" << '\n'
     << "q3 = " << q3 << '\n'
     << "q4 = " << q4 << " (ignorable coordinate)" << '\n'
     << "q5 = " << q5 << " (ignorable coordinate)" << '\n'
     << "q6 = " << q6 << " (ignorable coordinate)" << '\n'
     << "q7 = " << q7 << " (ignorable coordinate)" << '\n'
     << "u0 = " << u0 << " (dependent speed)" << '\n'
     << "u1 = " << u1 << '\n'
     << "u2 = " << u2 << " (dependent speed)" << '\n'
     << "u3 = " << u3 << '\n'
     << "u4 = " << u4 << " (dependent speed)" << '\n'
     << "u5 = " << u5 << '\n';
  } else {
    std::cerr << "Unable to open " << filename << "for writing." << '\n';
    std::cerr << "Aborting." << '\n';
    exit(0);
  }
} // writeState()

void Whipple::writeParameters(const char * filename) const
{
  std::ofstream fp(filename, std::ios::out);
  fp.precision(16);
  if (fp.is_open()) {
  fp << "rr   = " << rr << '\n'
     << "rrt  = " << rrt << '\n'
     << "rf   = " << rf << '\n'
     << "rft  = " << rft << '\n'
     << "lr   = " << lr << '\n'
     << "ls   = " << ls << '\n'
     << "lf   = " << lf << '\n'
     << "mr   = " << mr << '\n'
     << "mf   = " << mf << '\n'
     << "ICyy = " << ICyy << '\n'
     << "IDxx = " << IDxx << '\n'
     << "IDyy = " << IDyy << '\n'
     << "IDzz = " << IDzz << '\n'
     << "IDxz = " << IDxz << '\n'
     << "IExx = " << IExx << '\n'
     << "IEyy = " << IEyy << '\n'
     << "IEzz = " << IEzz << '\n'
     << "IExz = " << IExz << '\n'
     << "IFyy = " << IFyy << '\n'
     << "lrx  = " << lrx << '\n'
     << "lrz  = " << lrz << '\n'
     << "lfx  = " << lfx << '\n'
     << "lfz  = " << lfz << '\n'
     << "g    = " << g << '\n';
  } else {
    std::cerr << "Unable to open " << filename << "for writing." << '\n';
    std::cerr << "Aborting." << '\n';
    exit(0);
  }
} // writeParameters()

Whipple::Whipple()
{
  // Setup the root finder and the numerical integrator
  initRootFinder();
  initODESolver();

  // Camera settings
  theta = M_PI / 4.0;
  phi = 0.0;
  d = 1.0;
  ctx = .35;
  cty = .35;
  ctz = 0.0;

  // Default inputs, parameters, and initial state
  Trw = Tfw = Ts = 0.0;
  setBenchmarkParameters();
  evalConstants();          // evaluate the constant z's
  setBenchmarkState();      // set the state
  eoms();                   // compute du_i/dt, and intermediate z's
  computeOutputs();         // compute remaining z's and all outputs

  // Initialize the eigenvalues and eigenvectors
  evals = gsl_vector_complex_alloc(4);
  evecs = gsl_matrix_complex_alloc(4, 4);
  m = gsl_matrix_alloc(4, 4);
  gsl_matrix_set_zero(m);
  gsl_matrix_set(m, 0, 2, 1.0);
  gsl_matrix_set(m, 1, 3, 1.0);
  w = gsl_eigen_nonsymmv_alloc(4);
} // constructor

Whipple::~Whipple()
{
  gsl_root_fdfsolver_free(fdf_s);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  gsl_vector_complex_free(evals);
  gsl_matrix_complex_free(evecs);
  gsl_matrix_free(m);
  gsl_eigen_nonsymmv_free(w);
} // destructor

bool Whipple::setParameters(WhippleParams * p)
{
  bool validparameters = true;

  // Check that distances are non-negative
  if (p->rr < 0.0) {
    std::cerr << "rr must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->rf < 0.0) {
    std::cerr << "rf must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->rrt < 0.0) {
    std::cerr << "rrt must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->rft < 0.0) {
    std::cerr << "rft must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->lr < 0.0) {
    std::cerr << "lr must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->ls < 0.0) {
    std::cerr << "ls must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->lf < 0.0) {
    std::cerr << "lf must be greater than or equal to zero.\n";
    validparameters = false;
  }

  // Check axle to axle offset is enough to avoid wheel to wheel overlap
  if (validparameters) {
    double minimumAxleOffset = p->rr + p->rf + p->rrt + p->rft,
           d_zero = sqrt(pow(p->lr + p->lf, 2.0) + pow(p->ls, 2.0)),
           d_pi   = sqrt(pow(p->lr - p->lf, 2.0) + pow(p->ls, 2.0));

    if (d_zero < minimumAxleOffset) {
      overlap[0] = true;
      std::cerr << "Tire overlap will occur when steer = 0\n";
    } else {
      overlap[0] = false;
    }
    if (d_pi < minimumAxleOffset) {
      overlap[1] = true;
      std::cerr << "Tire overlap will occur when steer = pi\n";
    } else {
      overlap[1] = false;
    }
    if (overlap[0] && overlap[1]) {
      validparameters = false;
      std::cerr << "Tire overlap occurs for steer = 0 and steer = pi,"
              "non realistic model.\n";
    }
  }

  // Check that masses are non-negative
  if (p->mr < 0.0) {
    std::cerr << "mr must be greater than or equal to zero.\n";
    validparameters = false;
  }
  if (p->mf < 0.0) {
    std::cerr << "mf must be greater than or equal to zero.\n";
    validparameters = false;
  }

  // Check that principal moments of inertia are all positive and that the
  // triangle inequalities are satisfied.

  // Rear wheel
  if (p->ICyy < 0.0) {
    std::cerr << "Rear wheel spin inertia must be non-negative.\n";
    validparameters = false;
  }
  // Front wheel
  if (p->IFyy < 0.0) {
    std::cerr << "Front wheel spin inertia must be non-negative.\n";
    validparameters = false;
  }
  // Rear assembly
  if (!validInertia(p->IDxx, p->IDyy, p->IDzz, p->IDxz)) {
    std::cerr << "Rear assembly inertia invalid.\n";
    validparameters = false;
  }
  // Front assembly
  if (!validInertia(p->IDxx, p->IDyy, p->IDzz, p->IDxz)) {
    std::cerr << "Rear assembly inertia invalid.\n";
    validparameters = false;
  }

  ICyy = p->ICyy;
  IDxx = p->IDxx;
  IDxz = p->IDxz;
  IDyy = p->IDyy;
  IDzz = p->IDzz;
  IExx = p->IExx;
  IExz = p->IExz;
  IEyy = p->IEyy;
  IEzz = p->IEzz;
  IFyy = p->IFyy;
  g = p->g;
  lf = p->lf;
  lfx = p->lfx;
  lfz = p->lfz;
  lr = p->lr;
  lrx = p->lrx;
  lrz = p->lrz;
  ls = p->ls;
  mr = p->mr;
  mf = p->mf;
  rf = p->rf;
  rft = p->rft;
  rr = p->rr;
  rrt = p->rrt;

  return validparameters;
} // setParameters()

void Whipple::setState(const double state[10])
{
  q0 = state[0];  // Yaw
  q1 = state[1];  // Lean
  q3 = state[2];  // Steer
  calcPitch();    // Solve holonomic constraint for pitch angle
  q4 = state[3];
  q5 = state[4];
  q6 = state[5];
  q7 = state[6];
  u1 = state[7];
  u3 = state[8];
  u5 = state[9];
} // setState()

// TODO: fix inertias and cm distances
void Whipple::setBenchmarkParameters(void)
{
  MJWhippleParams * bin = new MJWhippleParams;
  WhippleParams * bout = new WhippleParams;
  // From Meijaard2007
  bin->rr = 0.3;
  bin->rrt = 0.0;
  bin->rf = 0.35;
  bin->rft = 0.0;
  bin->w = 1.02;
  bin->c = 0.08;
  bin->lambda = M_PI/10.0;
  bin->mr = 2.0;
  bin->mb = 85.0;
  bin->mh = 4.0;
  bin->mf = 3.0;
  bin->IRxx = 0.0603;
  bin->IRyy = 0.12;
  bin->IBxx = 9.2;
  bin->IByy = 11.0;
  bin->IBzz = 2.8;
  bin->IBxz = 2.4;
  bin->IHxx = 0.05892;
  bin->IHyy = 0.06;
  bin->IHzz = 0.00708;
  bin->IHxz = -0.00756;
  bin->IFxx = 0.1405;
  bin->IFyy = 0.28;
  bin->xb = 0.3;
  bin->zb = -0.9;
  bin->xh = 0.9;
  bin->zh = -0.7;
  bin->g = 9.81;
  // Convert to gyrostat parameters
  convertParameters(bout, bin);
  setParameters(bout);

  delete bin;
  delete bout;
} // setBenchmarkParameters()

void Whipple::setBenchmarkState(void)
{
  q0 = q1 = q3 = q4 = q5 = q6 = q7 = 0.0;
  q2 = M_PI/10.0;
  u1 = 0.5;
  u3 = 0.0;
  u5 = -4.6/(rf+rft);
} // setBenchmarkState

void Whipple::initRootFinder(void)
{
  // Newton-Raphson setup
  FDF.f = &hc_f;
  FDF.df = &hc_df;
  FDF.fdf = &hc_fdf;
  FDF.params = this;
  fdf_T = gsl_root_fdfsolver_newton;
  fdf_s = gsl_root_fdfsolver_alloc(fdf_T);
} // initRootFinder()

void Whipple::initODESolver(void)
{
  // Set integration settings
  t = 0.0;
  tf = 5.0;
  h = 0.001;
  fps = 100;
  T = gsl_odeiv_step_rk8pd;
  s = gsl_odeiv_step_alloc(T, 10);
  c = gsl_odeiv_control_y_new(1e-6, 1e-9);
  e = gsl_odeiv_evolve_alloc(10);
  sys.function = eomwrapper;
  sys.jacobian = NULL;
  sys.dimension = 10;
  sys.params = this;
} // initODESolver()

void Whipple::calcPitch(void)
{
  iter = 0;
  // Newton's method
  gsl_root_fdfsolver_set(fdf_s, &FDF, q2);
  do
  {
    gsl_root_fdfsolver_iterate(fdf_s);
    status = gsl_root_test_residual(
              hc_f(gsl_root_fdfsolver_root(fdf_s), this), 1e-16);
  } while(status == GSL_CONTINUE && ++iter < 100);
  q2 = gsl_root_fdfsolver_root(fdf_s);
} // calcPitch

void Whipple::calcEvals(void)
{
  gsl_matrix_set(m, 0, 0, 0.0);
  gsl_matrix_set(m, 0, 1, 0.0);
  gsl_matrix_set(m, 0, 3, 0.0);
  gsl_matrix_set(m, 1, 0, 0.0);
  gsl_matrix_set(m, 1, 1, 0.0);
  gsl_matrix_set(m, 1, 2, 0.0);
  gsl_matrix_set(m, 0, 2, 1.0);
  gsl_matrix_set(m, 1, 3, 1.0);
  // Evaluate the EOMS
  eoms();
  // Evaluate the 10x10 A matrix
  computeOutputs();

  // Get the 4x4 sub matrix of the 10x10 A matrix
  gsl_matrix_set(m, 2, 0, A[71]);
  gsl_matrix_set(m, 2, 1, A[72]);
  gsl_matrix_set(m, 2, 2, A[77]);
  gsl_matrix_set(m, 2, 3, A[78]);
  gsl_matrix_set(m, 3, 0, A[81]);
  gsl_matrix_set(m, 3, 1, A[82]);
  gsl_matrix_set(m, 3, 2, A[87]);
  gsl_matrix_set(m, 3, 3, A[88]);

  // Get the eigenvalues
  gsl_eigen_nonsymmv(m, evals, evecs, w);
  gsl_eigen_nonsymmv_sort(evals, evecs, GSL_EIGEN_SORT_ABS_ASC);
  getFourValues();
} // calcEvals()

void Whipple::getFourValues(void)
{
  int i;
  double valr, vali;
  if (evalCase() == 0) {
    for (i = 0; i < 4; ++i)
      fourValues[i] = GSL_REAL(gsl_vector_complex_get(evals, i));
    return;
  } else
    for (i = 0; i < 4; ++i) {
      vali = GSL_IMAG(gsl_vector_complex_get(evals, i));
      valr = GSL_REAL(gsl_vector_complex_get(evals, i));
      if (vali == 0.0) {
        fourValues[i] = valr;
        continue;
      }
      fourValues[i] = valr;
      fourValues[++i] = vali;
    }
} // getFourValues()

int Whipple::evalCase(void) const
{
  int i, real = 0;
  for (i = 0; i < 4; ++i)
   if (GSL_IMAG(gsl_vector_complex_get(evals, i)) == 0.0)
     ++real;

  if (real == 4)
    return 0;
  if (real == 2)
    return 1;
  // Final case, real == 0
  return 0;
} // evalCase()

void Whipple::evolve(double tj, double * state)
{
  gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tj, &h, state);
} // evolve()

void Whipple::printEvals(void) const
{
  std::cout.precision(16);
  std::cout << "Eigenvalues:\n";
  for (int i = 0; i < 4; ++i) {
    double real = GSL_REAL(gsl_vector_complex_get(evals, i));
    double imag = GSL_REAL(gsl_vector_complex_get(evals, i));
    std::cout << "lambda_" << i << " = "
              << real;
    if (imag > 0.0)
      std::cout << " + " << imag << "j\n";
    else
      std::cout << " - " << -imag << "j\n";
  }
}

bool Whipple::validInertia(double Ixx, double Iyy, double Izz, double Ixz) const
{
  bool isvalid = true;
  double sqrt_term = sqrt(pow(Ixx - Izz, 2.0) + 4.0*Ixz*Ixz);
  double I[3];

  I[0] = Iyy;
  I[1] = (Ixx + Izz - sqrt_term)/2.0;
  I[2] = (Ixx + Izz + sqrt_term)/2.0;

  insertionSort(3, I);

  for (int i = 0; i < 3; ++i) {
    if (I[i] < 0.0) {
      std::cerr << "I" << i+1 << " is negative.\n";
      isvalid = false;
    }
  } // for i

  if (I[0] > I[1] + I[2]) {
    std::cerr << "Inertia inequality not satsified: I1 > I2 + I3\n";
    isvalid = false;
  }
  if (I[1] > I[0] + I[2]) {
    std::cerr << "Inertia inequality not satsified: I2 > I1 + I3\n";
    isvalid = false;
  }
  if (I[2] > I[0] + I[1]) {
    std::cerr << "Inertia inequality not satsified: I3 > I1 + I2\n";
    isvalid = false;
  }

  return isvalid;
}

void Whipple::insertionSort(int N, double ar[]) const
{
  int j;
  double temp;
  for (int i = 1; i < N; ++i)
  {
    temp = ar[i];
    for (j = i - 1; j >= 0 && ar[j] > temp; j--)
      ar[j + 1] = ar[j];
    ar[j + 1] = temp;
  } // for i
}
