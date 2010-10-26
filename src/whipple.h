#ifndef WHIPPLE_H
#define WHIPPLE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>

#define Z_MAX 1127

using namespace std;

extern "C" {
  inline int eomwrapper(double t, const double x[10], double f[10], void * params);
}

// Constant parameters used in the model derivation
typedef struct {
  double ICxx,ICyy,IDxx,IDxz,IDyy,IDzz,IExx,IExz,IEyy,IEzz,IFxx,IFyy;
  double g,lf,lfx,lfz,lr,lrx,lrz,ls,mc,md,me,mf,rf,rft,rr,rrt;
} WhippleParams;

class Whipple {
  private:

  public:
    double t, tf, h;
    // Constant parameters
    double ICxx,ICyy,IDxx,IDxz,IDyy,IDzz,IExx,IExz,IEyy,IEzz,IFxx,IFyy;
    double g,lf,lfx,lfz,lr,lrx,lrz,ls,mc,md,me,mf,rf,rft,rr,rrt;
    // State variables, and their derivatives
    double q0,q1,q2,q3,q4,q5,q6,q7,u1,u3,u5;
    double q0p,q1p,q2p,q3p,q4p,q5p,q6p,q7p,u1p,u3p,u5p;
    // Inputs
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
    gsl_vector_complex *evals;
    gsl_matrix_complex *evecs;
    gsl_matrix * m;
    gsl_eigen_nonsymmv_workspace * w;
    double fourValues[4];

    // Member functions
    Whipple();
    ~Whipple();

    // Mutators
    void initRootFinder(void);
    void initODESolver(void);
    void setBenchmarkParameters(void);
    void setBenchmarkState(void);
    void calcPitch(void);
    void calcEvals(void);
    void setState(const double state[10]);
    void setParameters(WhippleParams * p);
    void evalConstants(void);
    void eoms(void);
    void computeOutputs(void);
    void getFourValues(void);
    int evalCase(void);

    // Accessors
    void writeRecord_dt(void) const;
    void printState(void) const;
    void printParameters(void) const;

    // Wrapper functions to interface with GSL required calling conventions
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
    friend double hc_f(double q2, void * params);
    friend double hc_df(double q2, void * params);
    friend void hc_fdf(double q2, void * params, double * f, double * df);
    friend ostream &operator<<(ostream &file, const Whipple * discs);
}; 
#endif
