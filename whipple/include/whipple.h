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

#include "whippleparams.h"

/**  \brief Whipple bicycle class.
 *
 * This class represents the Whipple bicycle model and provides means to:
 *   - Simulate nonlinear equations of motion
 *   - Calculate steady turning equilibrium conditions
 *   - Calculate eigenvalues of equilibrium
*/

class Whipple {
  public:
    /**
     * Default constructor
     *
     * Create an instance of the Whipple model using benchmark parameter
     * numerical values.  Initial condition is upright, zero steer, zero speed.
     *
     */
    Whipple();

    /**
     * Constructor
     *
     * Create an instance of the Whipple model using gyrostat parameter
     * numerical values.  Initial condition is upright, zero steer, zero speed.
     *
     */
    Whipple(const GyrostatParams & params);

    /**
     * Constructor
     *
     * Create an instance of the Whipple model using Meijaard parameter
     * numerical values.  Initial condition is upright, zero steer, zero speed.
     *
     */
    Whipple(const MeijaardParams & params);


    /**
     * Destructor
     */
    ~Whipple();

    /**
     * Externally applied inputs: steer torque, rear wheel torque, front wheel
     * torque.
     *  \param u An array of doubles to store \f$T_\delta\f$,
     *  \f$T_{\theta_r}\f$, \f$T_{\theta_f}\f$
     *  \param t Time
     *  \param x Complete state of bicycle.
     *
     *  By default, this function sets all input quantities to zero.  Overload
     *  this function to implement your own feedback controller or applied
     *  inputs to the bicycle model.
     */
    void inputs(double u[3], double t, const double x[14]) const { u[0] = u[1] = u[2] = 0.0;} ;

  private:
    GyrostatParams RearGyrostat,   /**< Rear gyrostat */
                   FrontGyrostat;  /**< Front gyrostat */

    double ls, /**< \f$l_s\f$, Steer axis offset */
           g;  /**< \f$g\f$, gravitational constant*/

    /**
     * \brief Complete state of bicycle.
     *
     * This vector includes both the independent and dependent state variables.
     *
     * The ordering of this state vector is:
     *   - \f$\psi\f$ Yaw
     *   - \f$\phi\f$ Lean
     *   - \f$\theta\f$ Pitch
     *   - \f$\delta\f$ Steer
     *   - \f$\theta_r\f$ Rear wheel angle
     *   - \f$\theta_f\f$ Front wheel angle
     *   - \f$x\f$ Rear wheel ground contact in n_x direction
     *   - \f$y\f$ Rear wheel ground contact in n_y direction
     *   - \f$\omega_x\f$
     *   - \f$\omega_y\f$
     *   - \f$\omega_z\f$
     *   - \f$\omega_s\f$
     *   - \f$\omega_r\f$
     *   - \f$\omega_f\f$
     */
    double x_complete[14];

    /**
     * \brief Minimal state of bicycle.
     *
     * This vector includes only the independent state variables.  There are 7
     * independent generalized coordinates and 3 indepedent generalized
     * speeds, for a total of a 10 dimensional space.  Of the 7 generalized
     * coordinates, 5 are ignorable with respect to the dynamic equations.
     *
     * This vector is used during numerical integration (simulation) as a means
     * to ensure that the state of the bicycle is evolved in a way that is
     * consistent with the constraint equations.  This approach is sometimes
     * referred to as "constraint embedding" in the literature, and it ensures
     * that the internal steps taken by the numerical ODE integrator are
     * consistent with the constraints.
     *
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
     * \f$q_1, q_2\f$ are the two independent configuration variables,
     * either lean and pitch, lean and steer, or pitch and steer.  The
     * dependent configuration variable is selected by computing the gradient
     * of the holonomic constraint equation and selecting the configuration
     * variable whose component in the gradient is largest.
     *
     * \f$u_1, u_2, u_3\f$ are the three independent generalized speeds, selected
     * from the set \f$\{\omega_x, \omega_y, \omega_z, \omega_s, \omega_r,
     * \omega_f\}\f$.  By the binomial coefficient, there are 20 possible
     * choices of independent generalized speeds (6!/(6-3)!/3!).  The choice of
     * independent generalized speeds is automatically selected using a
     * singular value decomposition of the nonholonomic constraint coefficient
     * matrix.
     */
    double x[10];

    /**
     * \brief Indices of independent generalized coordinates
     *
     * The holonomic constraint relates the lean, pitch and steer
     * angles.  The most numerically reliable way to evaluate which of these
     * three should be dependent is to compute the unit normal to the surface
     * of constraint, and take the dependent variable to be the one whose
     * component in the unit normal is largest.
     */
    int qIndependent[2];

    /**
     * \brief Indices of independent generalized speeds
     *
     * The nonholonomic constraint equations can be written as
     * \f[
     *   B \omega = 0
     * \f]
     * Where \f$B \in \mathbf{R}^{3 \times 6}, \omega \in \mathbf{R}^{6 \times
     * 1} \f$.  Once the independent speeds are selected, this can be rewritten
     * as
     * \f[
     *   B_i \omega_i + B_d \omega_d = 0
     * \f]
     */
    int uIndependent[3];
};
#endif
