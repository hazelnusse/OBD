#include <iostream>
#include "GyroStat.h"

GyroStat::GyroStat()
{

}

GyroStat::~GyroStat()
{

}

bool GyroStat::setTorusRadii(double majorRadius, double minorRadius)
{
  if (majorRadius < 0 || minorRadius < 0)
    return false;
  params[r] = minorRadius;
  params[R] = majorRadius;
  return true;
} // setTorusRadii

bool GyroStat::setHingeOffset(double offset)
{
  if (offset < 0.0)
    return false;
  params[l] = offset;
  return true;
} // setHingeOffset

void GyroStat::setRotorCenter(double dx, double dz)
{
  params[x] = dx;
  params[z] = dz;
} // setRotorCenter

bool GyroStat::setMass(double mass)
{
  if (mass < 0.0)
    return false;
  params[m] = mass;
  return true;
} // setMass

bool GyroStat::setSpinInertia(double SpinInertia)
{
  if (SpinInertia < 0.0)
    return false;
  params[J] = SpinInertia;
  return true;
} // setSpinInertia

bool GyroStat::setCarrierInertia(double InertiaScalars[4])
{
  if (InertiaScalars[1] < 0.0)
    return false;
  // TODO Determine conditions for inertial feasibility given a cylindrical
  // gyrostat
  params[Ixx] = InertiaScalars[0];
  params[Iyy] = InertiaScalars[1];
  params[Izz] = InertiaScalars[2];
  params[Ixz] = InertiaScalars[3];
  return true;
}

bool GyroStat::setParameters(double (&parameters)[11])
{
  params[r] = parameters[r];
  params[R] = parameters[R];
  params[l] = parameters[l];
  params[x] = parameters[x];
  params[z] = parameters[z];
  params[m] = parameters[m];
  params[J] = parameters[J];
  params[Ixx] = parameters[Ixx];
  params[Iyy] = parameters[Iyy];
  params[Izz] = parameters[Izz];
  params[Ixz] = parameters[Ixz];
  return true;
}

void GyroStat::printParameters(void) const
{
  std::cout <<
       "r   = " << params[r] << "\n"
       "R  = " << params[R] << "\n"
       "l   = " << params[l] << "\n"
       "x   = " << params[x] << "\n"
       "z   = " << params[z] << "\n"
       "m   = " << params[m] << "\n"
       "J  = " << params[J] << "\n"
       "Ixx = " << params[Ixx] << "\n"
       "Iyy = " << params[Iyy] << "\n"
       "Izz = " << params[Izz] << "\n"
       "Ixz = " << params[Ixz] << "\n";
} // printParameters
