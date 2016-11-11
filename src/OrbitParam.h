#ifndef ORBITPARAM_H
#define ORBITPARAM_H

namespace IdealRocketEquation
{

// parameters used to define the initial orbit of the spacecraft
class OrbitParam
{
    public:
        OrbitParam();
        virtual ~OrbitParam();

        double ecc;     // the orbital eccentricity of the spacecraft's current orbit       (n/a)
        double mu;      // the 'GM' parameter for the central gravitating body              (m^3 s^-2)
        double smaj;    // the semi-major axis of the spacecraft's current orbit            (m)

    protected:

    private:
};
}
#endif // ORBITPARAM_H
