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

/**  \brief Whipple bicycle class.
 *
 * This class represents the Whipple bicycle model and provides means to:
 *   - Simulate nonlinear equations of motion
 *   - Calculate steady turning equilibrium conditions
 *   - Calculate eigenvalues of equilibrium
*/

// Forward declarations
struct GyrostatParams;
struct MeijaardParams;

class Whipple {
  public:
    /**
     * Default constructor
     */
    Whipple();

    /**
     * Constructor
     */
    Whipple(const struct GyrostatParams & params);

    /**
     * Constructor
     */
    Whipple(const struct MeijaardParams & params);


    /**
     * Destructor
     */
    ~Whipple();

  private:
    /**
     * Rear Gyrostat Constants
     */
    double Jr, Irxx, Iryy, Irzz, Irxz, mr;
    double lr, xr, zr, Rr, rr;

    /**
     * Front Gyrostat Constants
     */
    double Jf, Ifxx, Ifyy, Ifzz, Ifxz, mf;
    double lf, xf, zf, Rf, rf;

    /**
     * Steer axis offset and gravity
     */
    double ls, g;

    /**
     * \brief Complete State of bicycle
     *
     * This array includes both the independent and dependent state variables,
     * as well as variables which are ignorable with respect to the dynamic
     * equations.
     *
     * The ordering of this state vector is:
     *   - Yaw
     *   - Lean
     *   - Pitch
     *   - Steer
     *   - Rear wheel angle
     *   - Front wheel angle
     *   - Rear wheel ground contact in n_x direction
     *   - Rear wheel ground contact in n_y direction
     *   - omegax
     *   - omegay
     *   - omegaz
     *   - omegas
     *   - omegar
     *   - omegas
     */
    double x_complete[14];

    /**
     * \brief Minimal State of bicycle
     *
     * This array includes only the independent state variables.
     * The ordering of this state vector is:
     *   - \f$q_1\f$
     *   - \f$q_2\f$
     *   - \f$u_1\f$
     *   - \f$u_2\f$
     *   - \f$u_3\f$
     *   - \f$\psi\f$ Yaw
     *   - \f$\theta_r\f$Rear wheel angle
     *   - \f$\theta_f\f$Front wheel angle
     *   - \f$x\f$ Rear wheel ground contact in \f$\mathbf{n}_x\f$ direction
     *   - \f$y\f$ Rear wheel ground contact in \f$\mathbf{n}_y\f$ direction
     *
     * q1, q2 represent the two independent configuration variables, they will
     * be either lean and pitch, lean and steer, or pitch and steer.
     *
     * u1, u2, u3 represent the two independent configuration variables, they
     * will be selected from the set wx, wy, wz, ws, wr, wf.  There are 20
     * possible choices of independent generalized speeds (6!/(6-3)!/3!)
     */
    double x[10];

    /**
     * \brief Indices of independent generalized coordinates
     *
     * The holonomic constraint imposes constrains the lean, pitch and steer
     * angles.  The most numerically reliable way to evalue which of these
     * three should be dependent is to compute the unit normal to the surface
     * of constraint, and take the dependent variable to be the one whose
     * component in the unit normal is largest.
     *
     * For motion near upright, pitch should be chosen as dependent.  The
     * farther the bicycle leans, the less this becomes true.
     *
     */
    int qIndependent[2];

    /**
     * \brief Indices of independent generalized speeds
     */
    int uIndependent[3];



    /*
  public:
    double t, tf, h;
    // Constant parameters
    // State variables, and their derivatives
    double q0,q1,q2,q3,q4,q5,q6,q7,u1,u3,u5;
    double q0p,q1p,q2p,q3p,q4p,q5p,q6p,q7p,u0p,u1p,u2p,u3p,u4p,u5p;
    // Inputs torques
    double Tfw,Trw,Ts;

    // Boolean array to detect when wheel overlap will occur in either
    // handlebar reversed or handlebar forward configurations.
    bool overlap[2];

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


    // Mutators
    void calcEvals(void);
    void calcPitch(void);
    void computeOutputs(void);
    void eoms(void);

    /
     * Method for evaluating z's which are constant with respect to the bicycle
     * state.
     *
     * @pre Physical parameters of bicycle model must be set.
     *
     * @post z's which are constant with respect to the bicycle state are evaluated.
     *
     * @return None
     * /
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
    //void writeEvalRecord_dt(const char * filename) const;
    void writeParameters(const char * filename) const;
    void writeState(const char * filename) const;
    void printState(void) const;
    void printCfgCon(void) const;
    void printParameters(void) const;
    void printEvals (void) const;
    int evalCase(void) const;

    // Wrapper functions to interface with GSL required calling conventions
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
    friend double hc_f(double q2, void * params);
    friend double hc_df(double q2, void * params);
    friend void hc_fdf(double q2, void * params, double * f, double * df);
    friend std::ostream &operator<<(std::ostream &file, const Whipple * discs);

  private:
    bool validInertia(double Ixx, double Iyy, double Izz, double Ixz) const;
    void insertionSort(int N, double ar[]) const;
  */
};
#endif
