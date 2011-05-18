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
/** Parameters used in the gyrostat derivation
 */
struct GyrostatParams {
  double Jr,    /**< Rear wheel spin moment of inertia  */
         Irxx,  /**< Rear gyrostat xx moment of inertia */
         Iryy,  /**< Rear gyrostat yy moment of inertia */
         Irzz,  /**< Rear gyrostat zz moment of inertia */
         Irxz,  /**< Rear gyrostat xz moment of inertia */
         mr,    /**< Rear gyrostat total mass */
         xr,    /**< Distance from rear wheel center to mass center in body
                  fixed x direction */
         zr,    /**< Distance from rear wheel center to mass center in body
                  fixed z direction */
         lr,    /**< Distance from rear wheel center to steer axis */
         Rr,    /**< Rear wheel major radius */
         rr,    /**< Rear wheel minor radius */

         Jf,    /**< Front wheel spin moment of inertia  */
         Ifxx,  /**< Front gyrostat xx moment of inertia */
         Ifyy,  /**< Front gyrostat yy moment of inertia */
         Ifzz,  /**< Front gyrostat zz moment of inertia */
         Ifxz,  /**< Front gyrostat xz moment of inertia */
         mf,    /**< Front gyrostat total mass */
         xf,    /**< Distance from front wheel center to mass center in body
                  fixed x direction */
         zf,    /**< Distance from front wheel center to mass center in body
                  fixed z direction */
         lf,    /**< Distance from front wheel center to steer axis */
         Rf,    /**< Front wheel major radius */
         rf,    /**< Front wheel minor radius */
         ls,    /**< Steer axis offset */
         g;     /**< Gravitational constant */
};

/** \struct MeijaardParams
 *
 * Parameters used in the Meijaard derivation
 *
 */
struct MeijaardParams {
  double w,         /**< Wheel base */
         c,         /**< Trail */
         lambda,    /**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         g,         /**< Gravitational constant */
         rR,        /**< Rear wheel radius */
         mR,        /**< Rear wheel mass */
         IRxx,      /**< Rear wheel transverse mass moment of inertia */
         IRyy,      /**< Rear wheel spin moment of inertia */
         xB,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         zB,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         mB,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IBxx,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IByy,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IBzz,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IBxz,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         xH,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         zH,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         mH,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IHxx,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IHyy,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IHzz,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IHxz,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         rF,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         mF,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IFxx,/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
         IFyy;/**< Steer axis tilt (\f$ \pi/2 \f$ - head angle) */
};
#endif
