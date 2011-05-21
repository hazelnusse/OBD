/* whippleparams.h
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
#ifndef WHIPPLEPARAMS_H
#define WHIPPLEPARAMS_H

/** \file whippleparams.h
 * \brief Definitions of structures for model parameterizations.
 */

/** \struct GyrostatParams
 *
 * \brief Parameters used for the front and rear assemblies
 *
 * The front and rear assemblies of the bicycle may each be viewed as a
 * cylindrical gyrostat with a torus representing each wheel and a laterally
 * symmetric rigid body representing the frame and/or fork.  Eleven parameters
 * are need to describe each gyrostat.
 */
struct GyrostatParams {
  double J,    /**< Wheel spin moment of inertia  */
         Ixx,  /**< Gyrostat xx moment of inertia */
         Iyy,  /**< Gyrostat yy moment of inertia */
         Izz,  /**< Gyrostat zz moment of inertia */
         Ixz,  /**< Gyrostat xz moment of inertia */
         m,    /**< Gyrostat total mass */
         x,    /**< Distance from wheel center to mass center in body
                  fixed x direction */
         z,    /**< Distance from wheel center to mass center in body
                  fixed z direction */
         l,    /**< Distance from wheel center to steer axis */
         R,    /**< Wheel major radius */
         r;    /**< Wheel minor radius */
};

/** \struct MeijaardParams
 *
 * \brief Parameters used in the Meijaard derivation
 *
 * In the reference configuration, the bicycle is upright, with zero steer,
 * and is pointed in the positive \f$x\f$ direction relative to the inertial
 * frame.  The wheelbase \f$w\f$, trail \f$c\f$, and steer axis tilt
 * \f$ \lambda \f$, and the mass center locations \f$ x_B, z_B, x_H, z_H \f$, are
 * all defined relative to an inertial frame with origin at the rear wheel
 * contact, when the bicycle is in the reference configuration.  Similarly, the
 * inertia scalars are all defined with respect to body fixed frames which are
 * aligned with the inertial frame when the bicycle is in the reference
 * configuration.
 */
struct MeijaardParams {
  double w,      /**< \f$w\f$, Wheel base. */
         c,      /**< \f$c\f$, Trail. */
         lambda, /**< \f$\lambda\f$, Steer axis tilt (\f$ \pi/2 \f$ - head
                   angle). */
         g,       /**< \f$g\f$, Gravitational constant. */
         rR,     /**< \f$ r_{\textrm{R}} \f$, Rear wheel radius. */
         mR,     /**< \f$m_{\textrm{R}}\f$, Rear wheel mass. */
         IRxx,   /**< \f$I_{\textrm{R}xx}\f$, central moment of inertia of
                   \f$\textrm{R}\f$ about any line in plane of symmetry. */
         IRyy,   /**< \f$I_{\textrm{R}xx}\f$, central moment of inertia of
                   \f$\textrm{R}\f$ for \f$\mathbf{n}_y\f$. */
         xB,     /**< \f$x_{\textrm{B}}\f$, rear frame mass center location in
                   \f$\mathbf{n}_x\f$
                   direction. */
         zB,     /**< \f$z_{\textrm{B}}\f$, rear frame mass center location in
                   \f$\mathbf{n}_z\f$
                   direction. */
         mB,     /**< \f$m_{\textrm{B}}\f$, rear frame mass. */
         IBxx,   /**< \f$I_{\textrm{B}xx}\f$ central moment of inertia of
                   \f$\textrm{B}\f$ for \f$\mathbf{n}_x\f$. */
         IByy,   /**< \f$I_{\textrm{B}yy}\f$ central moment of inertia of
                   \f$\textrm{B}\f$ for \f$\mathbf{n}_y\f$. */
         IBzz,   /**< \f$I_{\textrm{B}zz}\f$ central moment of inertia of
                   \f$\textrm{B}\f$ for \f$\mathbf{n}_z\f$. */
         IBxz,   /**< \f$I_{\textrm{B}xz}\f$ central product of inertia of
                   \f$\textrm{B}\f$ for \f$\mathbf{n}_x\f$ and
                   \f$\mathbf{n}_z\f$. */
         xH,     /**< \f$x_{\textrm{H}}\f$, front frame mass center location in
                   \f$\mathbf{n}_x\f$. */
         zH,     /**< \f$z_{\textrm{H}}\f$, front frame mass center location in
                   \f$\mathbf{n}_z\f$. */
         mH,     /**< \f$m_{\textrm{H}}\f$, front frame mass. */
         IHxx,   /**< \f$I_{\textrm{H}xx}\f$ central moment of inertia of
                   \f$\textrm{H}\f$ for \f$\mathbf{n}_x\f$. */
         IHyy,   /**< \f$I_{\textrm{H}yy}\f$ central moment of inertia of
                   \f$\textrm{H}\f$ for \f$\mathbf{n}_y\f$. */
         IHzz,   /**< \f$I_{\textrm{H}zz}\f$ central moment of inertia of
                   \f$\textrm{H}\f$ for \f$\mathbf{n}_z\f$. */
         IHxz,   /**< \f$I_{\textrm{H}xz}\f$ central product of inertia of
                   \f$\textrm{H}\f$ for \f$\mathbf{n}_x\f$ and
                   \f$\mathbf{n}_z\f$ */
         rF,     /**< \f$r_\textrm{F}\f$, Front wheel radius. */
         mF,     /**< \f$m_\textrm{F}\f$, Front wheel mass. */
         IFxx,   /**< \f$I_{\textrm{F}xx}\f$, central moment of inertia of
                   \f$\textrm{R}\f$ about any line in plane of symmetry. */
         IFyy;   /**< \f$I_{\textrm{F}xx}\f$, central moment of inertia of
                   \f$\textrm{R}\f$ for \f$\mathbf{n}_y\f$. */
};
#endif
