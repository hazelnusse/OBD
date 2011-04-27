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

#include <iostream>
#include "whipple.h"
#include "GyroStat.h"

Whipple::Whipple()
{
  // Allocate rear and front gyrostats
  R = new GyroStat;
  F = new GyroStat;
  // Default inputs, parameters, and initial state
  Trw = Tfw = Ts = 0.0;
  setBenchmarkParameters();
  // setBenchmarkState();

} // constructor

Whipple::~Whipple()
{
} // destructor

void Whipple::holonomicConstraint(double lps[3], double *f, double df[3]) const
{
} // holonomicConstraint

void Whipple::printParameters(void) const
{
  std::cout << "Rear gyrostat parameters:\n";
  R->printParameters();
  std::cout << "Front gyrostat parameters:\n";
  F->printParameters();
  std::cout << "Steer axis offset and gravity:\n"
               "ls   = " << ls << "\n"
               "g    = " << g  << "\n";
} // printParameters()

void Whipple::setBenchmarkParameters(void)
{
  /*
  rr   = 0.3;
  rrt  = 0.0;
  lr   = 0.9534570696121847;
  xr  = 0.4599058376856175;
  zr  = -0.4669419422355365;
  mr   = 87.0;
  Irs  = 0.12;
  Irxx = 7.684799791449107;
  Iryy = 11.87931034482759;
  Irzz = 5.315110553378478;
  Irxz = 4.262158094617231;

  rf   = 0.35;
  rft  = 0.0;
  lf   = 0.0320714267276193;
  mf   = 7.0;
  xf  = -0.003411905099535338;
  zf  = -0.2114010400161699;
  IExx = 0.4335379755311008;
  IEyy = 0.2946857142857143;
  IEzz = 0.1481477387546136;
  IExz = 0.005332503757935522;
  IFyy = 0.28;

  ls   = 0.2676445084476887;
  g    = 9.81;
  */
} // setBenchmarkParameters()
