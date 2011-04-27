#ifndef GYROSTAT_H
#define GYROSTAT_H
class GyroStat {
  public:
    GyroStat();
    ~GyroStat();

    //! enum for GyroStat parameters
    enum parameters {
      r,  /*!< Rotor minor radius. */
      R,  /*!< Rotor major radius. */
      l,  /*!< Non-negative distance between mass center and hinge axis */
      x,  /*!< Distance from gyrostat CM to rotor center in x direction */
      z,  /*!< Distance from gyrostat CM to rotor center in z direction */
      m,  /*!< Gyrostat total mass (carrier and rotor). */
      J,  /*!< Rotor moment of inertia about spin axis. */
      Ixx,  /*!< Carrier xx moment of inertia. */
      Iyy,  /*!< Carrier yy moment of inertia. */
      Izz,  /*!< Carrier zz moment of inertia. */
      Ixz  /*!< Carrier xz product of inertia. */
    };
    // Mutators
    bool setTorusRadii(double majorRadius, double minorRadius);
    bool setHingeOffset(double offset);
    void setRotorCenter(double dx, double dz);
    bool setMass(double mass);
    bool setSpinInertia(double SpinInertia);
    bool setCarrierInertia(double CarrierInertiaScalars[4]);
    bool setParameters(double (&parameters)[11]);

    // Accessors
    void printParameters(void) const;

    // Dynamic methods


  private:
    double params[11];
    double omega_a_n[3];
    double omega_b_n[3];

};
#endif
