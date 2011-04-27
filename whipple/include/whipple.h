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

class GyroStat;

class Whipple {

  public:
    /**
     * Default constructor
     */
    Whipple();

    /**
     * Default destructor
     */
    ~Whipple();

    /**
     * Evaluate the holonomic constraint and it's gradient
     */
    void holonomicConstraint(double lps[3], double * f, double df[3]) const;

    /**
     * Print model parameters
     */
    void printParameters(void) const;

    /**
     * Set model parameters to values equivalent to those in Meijaard2007
     */
    void setBenchmarkParameters(void);


  private:
    /**
     * Front and Rear gyrostat
     */
    GyroStat *R, *F;

    /**
     * Gyrostat steer axis offset
     */
    double ls;

    /**
     * Gravitational constant
     */
    double g;

    /**
     * Configuration variables and their time derivatives
     */
    double psi, phi, theta, delta, theta_r, theta_f, x, y;
    double psip, phip, thetap, deltap, theta_rp, theta_fp, xp, yp;
    double psipp, phipp, thetapp, deltapp, theta_rpp, theta_fpp, xpp, ypp;

    /**
     * Generalized speeds and their time derivatives
     */
    double wx, wy, wz, ws, wr, wf;
    double wxp, wyp, wzp, wsp, wrp, wfp;

    /**
     * Applied torques:  Rear wheel, Front wheel, Steer
     */
    double Trw, Tfw, Ts;

    // Output quantities
    double no_fn[3], cn_cm[3], h2_cl[3], constraints[3];

    /**
     * Rear wheel ground reaction forces
     */
    double Rx, Ry, Rz;

    /**
     * Front wheel ground reaction forces
     */
    double Fx, Fy, Fz;

    /**
     * Kinetic and potential energy
     */
    double ke, pe;

    /**
     * Boolean array to indicate wheel overlap will occur in either
     * handlebar reversed or handlebar forward configurations.
     */
    bool overlap[2];
};
#endif
