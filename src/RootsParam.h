#ifndef ROOTSPARAM_H
#define ROOTSPARAM_H

namespace IdealRocketEquation
{
class RootsParam
{
    public:
        RootsParam();
        virtual ~RootsParam();

        double r;              // the orbital radius at the start of the orbit insertion burn
        double tm;             // the duration of the orbit insertion burn (Initial guess)
        double dr;             // the small distance displacement for calculating derivatives of the final eccentricity vector with resepct to r;

    protected:

    private:
};
}
#endif // ROOTSPARAM_H
