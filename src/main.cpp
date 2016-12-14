//  main.c

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "NewtonRaphsonRoots.h"
#include "OrbitParam.h"
#include "VesselParam.h"
#include "RootsParam.h"

using namespace IdealRocketEquation;

int main(int argc, const char * argv[])
{
    OrbitParam opar;
    opar.ecc = 0.960122;
    opar.mu = 4.9048695e12;
    opar.smaj = 84781532.99;

    VesselParam vpar;
    vpar.m0 = 87200.0;
    vpar.gamma = 64.507;
    vpar.ve = 33000.0;

    RootsParam rpar;
    rpar.tm = 19.586;
    rpar.dr = 1.0;

    // A routine to calculate an initial 'rpar.r' from BTC's estimate of the time to periapsis, 'tbtc'.
    // To be used on the first pass-through.  Thereafter, the previous estimate calculated by the routine should be used.
    double M, dE, E, tbtc;
    tbtc = 10.0;

    if (opar.smaj >0) {
        M = -tbtc * sqrt( opar.mu / opar.smaj / opar.smaj / opar.smaj);
        E = M;
        dE = 1.0e10;
        while( fabs(dE) > 1.0e-10) {
            dE = - (E - opar.ecc * sin(E) - M)/(1 - opar.ecc * cos(E));
            E += dE;
        }
        rpar.r = opar.smaj * (1 - opar.ecc * cos(E));
    } else {
        M = -tbtc * sqrt( - opar.mu / opar.smaj / opar.smaj / opar.smaj);
        E = M;
        dE = 1.0e10;
        while( fabs(dE) > 1.0e-10) {
            dE = - (M + E - opar.ecc * sinh(E))/(1 - opar.ecc * cosh(E));
            E += dE;
        }
        rpar.r = opar.smaj * (1 - opar.ecc * cosh(E));//
    }

    NewtonRaphsonRoots nrr;
    const NewtonRaphsonRoots::Result & res = nrr.FindRoots(rpar, opar, vpar);

    // print the solution
    // 'r'   is the orbital radius of the spacecraft at the start of the orbit insertion burn (this can be converted to a time measure.)
    // 'tm'  is the duration (seconds) of the orbit insertion burn
    // 'ecc' is the final orbital eccentricity at the end of the orbit insertion burn.  This should be zero for a circular orbit
    // 'r_c' is the orbital radius at the end of the burn (metres).  For entry into a circular orbit, it also the radius of the final circular orbit.
    // 'tbp' is the time before (current) periapsis when the orbit insertion burn should start.

    printf("r       %2.10f\n", res.r  );
    printf("tm      %2.10f\n", res.tm );
    printf("\n");
    printf("ecc     %2.10f\n", res.ecc );
    printf("r_c     %2.10f\n", res.r_c );
    printf("tbp     %2.10f\n", res.tbp );
    printf("\n");

    return 0;
}
