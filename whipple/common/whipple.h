/* whipple.h
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
#ifndef WHIPPLE_H
#define WHIPPLE_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>

#define Z_MAX 987
#define ZS_MAX 652

using namespace std;

extern "C" {
  inline int eomwrapper(double t, const double x[10], double f[10], void * params);
}

// structure use to store command line options passed to steady turning code
typedef struct {
  char outfolder[512];// output folder name
  size_t N;           // # of points to mesh steer in [0, pi]
  bool all;           // controls whether to mesh whole feasible region
  gsl_vector * iso_v, * iso_t, * iso_mew;
} steadyOpts_t;

// Constant parameters used in the model derivation
typedef struct {
  double ICyy,IDxx,IDxz,IDyy,IDzz,IExx,IExz,IEyy,IEzz,IFyy;
  double g,lf,lfx,lfz,lr,lrx,lrz,ls,mr,mf,rf,rft,rr,rrt;
} WhippleParams;

// Constant parameters used in the Meijaard derivation
typedef struct {
  double rr, rrt, rf, rft, w, c, lambda;  // Wheel and frame geometric parameters
  double mr, mb, mh, mf;          // R. whl, rear frame, fork, F. whl mass
  double IRxx, IRyy;              // Rear wheel inertia scalars, IRxx == IRzz
  double IBxx, IByy, IBzz, IBxz;  // Rear frame and rider inertia scalars
  double IHxx, IHyy, IHzz, IHxz;  // Front fork and handlebar inertia scalars
  double IFxx, IFyy;              // Front wheel inertia scalars, IFxx == IFzz
  double xb, zb, xh, zh;          // COM locations relative to rear contact
  double g;
} MJWhippleParams;

class Whipple {
  public:
    double t, tf, h;
    // Constant parameters
    double ICyy,IDxx,IDxz,IDyy,IDzz,IExx,IExz,IEyy,IEzz,IFyy;
    double g,lf,lfx,lfz,lr,lrx,lrz,ls,mr,mf,rf,rft,rr,rrt;
    // State variables, and their derivatives
    double q0,q1,q2,q3,q4,q5,q6,q7,u1,u3,u5;
    double q0p,q1p,q2p,q3p,q4p,q5p,q6p,q7p,u0p,u1p,u2p,u3p,u4p,u5p;
    // Inputs torques
    double Tfw,Trw,Ts;

    // Output quantities
    double z[Z_MAX], no_fn[3], cn_cm[3], h2_cl[3], constraints[3];
    double fa_yaw, fa_lean, fa_pitch;
    double A[100], B[30],C[50];
    //double T_CN[16], T_CL[16], T_B[16], T_C[16], T_D[16], T_E[16], T_F[16];
    //double T_G[16], T_FN[16], T_CGL[16], T_FGL[16];
    double Fx,Fy,Fz,Rx,Ry,Rz,u0,u2,u4,ke,pe;
    // Camera variables
    double theta, phi, d, ctx, cty, ctz;

    // Numerical integrator variables
    gsl_odeiv_evolve * e;
    gsl_odeiv_control * c;
    gsl_odeiv_step * s;
    gsl_odeiv_system sys;
    const gsl_odeiv_step_type * T;
    int fps; // frames per second, controls output data time interval

    // Variables used for Newton-Raphson root finding algorithm
    gsl_root_fdfsolver * fdf_s;
    const gsl_root_fdfsolver_type * fdf_T;
    gsl_function_fdf FDF;
    int status, iter;

    // Variables for calculating eigenvalues and eigenvectors
    gsl_vector_complex * evals;
    gsl_matrix_complex * evecs;
    gsl_matrix * m;
    gsl_eigen_nonsymmv_workspace * w;
    double fourValues[4];

    // Variables associated with steady turning
    double z_s[ZS_MAX];
    // Coefficient of friction required for front and rear wheels
    double mew_f, mew_r;
    // These attributes are set when RigidRiderSteady is called
    double u5s_s, Ts_s, Ry_s, Rz_s, Fy_s, Fz_s;
    double F[11], dF[36];

    // Member functions
    Whipple();
    ~Whipple();

    // Mutators
    void calcEvals(void);
    void calcPitch(void);
    void computeOutputs(void);
    void eoms(void);
    void evalConstants(void);
    void evolve(double tj, double * state);
    void getFourValues(void);
    void initRootFinder(void);
    void initODESolver(void);
    void setBenchmarkParameters(void);
    void setBenchmarkState(void);
    bool setParameters(WhippleParams * p);
    void setState(const double state[10]);

    // Steady turning related functions
    void steadyEqns(void);
    void steadyCalcs(steadyOpts_t * options);

    // Accessors
    void writeEvalRecord_dt(const char * filename) const;
    void writeSimRecord_dt(const char * filename) const;
    void writeParameters(const char * filename) const;
    void writeState(const char * filename) const;
    void printState(void) const;
    void printParameters(void) const;
    void printEvals (void) const;
    int evalCase(void) const;

    // Wrapper functions to interface with GSL required calling conventions
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
    friend double hc_f(double q2, void * params);
    friend double hc_df(double q2, void * params);
    friend void hc_fdf(double q2, void * params, double * f, double * df);
    friend ostream &operator<<(ostream &file, const Whipple * discs);

  private:
    bool validInertia(double Ixx, double Iyy, double Izz, double Ixz) const;
    void insertionSort(int N, double ar[]) const;

};
#endif
