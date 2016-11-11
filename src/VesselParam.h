#ifndef VESSELPARAM_H
#define VESSELPARAM_H

namespace IdealRocketEquation
{
// parameters used to define the initial state of the spacecraft
class VesselParam
{
    public:
        VesselParam();
        virtual ~VesselParam();

        double m0;      // the spacecraft's initial mass                                    (kg)
        double gamma;   // the fuel mass flow rate                                           (kg/s)
        double ve;      // the effective velocity of engine exhaust gases                   (m/s)

    protected:

    private:
};
}
#endif // VESSELPARAM_H
