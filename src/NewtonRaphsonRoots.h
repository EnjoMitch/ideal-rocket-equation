#ifndef NEWTONRAPHSONROOTS_H
#define NEWTONRAPHSONROOTS_H

#include <vector>

namespace IdealRocketEquation
{
class OrbitParam;
class VesselParam;
class RootsParam;

class NewtonRaphsonRoots
{
    public:
        NewtonRaphsonRoots();
        virtual ~NewtonRaphsonRoots();

        struct Result
        {
            Result()
            {
                tm = r = ecc = r_c = tbp = 0;
            }
            double tm;
            double r;
            double ecc;
            double r_c;
            double tbp;
        };

        Result FindRoots(const RootsParam & rpar, const OrbitParam & opar, const VesselParam & vpar) const;

    protected:

    private:
};
}
#endif // NEWTONRAPHSONROOTS_H
