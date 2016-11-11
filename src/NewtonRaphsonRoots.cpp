#include "NewtonRaphsonRoots.h"
#include "TaylorIntegration.h"
#include "OrbitParam.h"
#include "RootsParam.h"
#include <cmath>
#include <cstdio>

using namespace IdealRocketEquation;

NewtonRaphsonRoots::NewtonRaphsonRoots(){}
NewtonRaphsonRoots::~NewtonRaphsonRoots(){}

NewtonRaphsonRoots::Result NewtonRaphsonRoots::FindRoots(const RootsParam & rpar, const OrbitParam & opar, const VesselParam & vpar) const
{
    const double SMAJ = opar.smaj;
    const double ECC = opar.ecc;
    const double MU = opar.mu;
    // INPUTS
    double r = rpar.r;
    const double dr = rpar.dr;
    double tm = rpar.tm;

    // OUTPUTS
    double   e[4];                     // the radial and transverse components of the eccentricity vector at the end of the orbital insertion burn
                                       // e[0] holds the radial component of the eccentricity vector
                                       // e[1] holds the transverse component of the eccentricity vector
                                       // e[2] holds the derivative of the radial component with respect to time
                                       // e[3] holds the derivative of the transverse component with respect to time

    double   g[3];                     // a vector to hold other useful information to be used once the Newton-Raphson iteration has been completed
                                       // g[0] holds the the angle variable after orbit insertion
                                       // g[1] holds the spacecraft acceleration due to thrust immediatel prior to the end of the insertion burn
                                       // g[2] holds the spacecraft's specific angular momentum

    double f0[2];  f0[0] = f0[1] = 0;  // a 2-vector declaration to hold the values of the orbital eccentricity vector at the current N-R iteration
    double jac[2][2];                  // a 2 x 2 array declaration to hold the Jacobian
    double det;                        // the determinant of the Jacobian
    int    flag = 1;


    // perform the Newton-Raphson (N-R) root-finding algorithm
    int     status;
    TaylorIntegration tin;


    while ((sqrt(f0[0] * f0[0] + f0[1] * f0[1]) > 1.e-8) || flag == 1) {

        flag      = 0;

        status    = tin.fun(opar, vpar, r + dr, tm, e, g);
        jac[0][0] = e[0];
        jac[1][0] = e[1];

        status    = tin.fun(opar, vpar, r - dr, tm, e, g);
        jac[0][0]-= e[0]; jac[0][0] = jac[0][0] / dr / 2;
        jac[1][0]-= e[1]; jac[1][0] = jac[1][0] / dr / 2;

        status    = tin.fun(opar, vpar, r, tm, e, g);
        jac[0][1] = e[2];
        jac[1][1] = e[3];

        f0[0]     = e[0];
        f0[1]     = e[1];

        det       = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        r         = r  - (+ jac[1][1] * f0[0] - jac[0][1] * f0[1]) / det;
        if (r <= SMAJ * (1 - ECC)) {
            r = SMAJ * (1 - ECC) + 100;
        }
        if (r >= SMAJ * (1 + ECC) && SMAJ > 0) {
            r = SMAJ * (1 + ECC) - 100;
        }
        tm        = tm - (- jac[1][0] * f0[0] + jac[0][0] * f0[1]) / det;
    }

    status    = tin.fun(opar, vpar, r, tm, e, g);
    f0[0]     = e[0];
    f0[1]     = e[1];

    Result res;
    res.r = r;
    res.tm = tm;
    res.ecc = sqrt(f0[0] * f0[0] + f0[1] * f0[1]);
    res.r_c = g[2] * g[2] / MU / (1 + f0[0]);
    if (SMAJ > 0) {
        res.tbp = +sqrt(+SMAJ/MU)*(sqrt(-SMAJ*SMAJ*(1 - ECC*ECC) + 2*SMAJ*r - r*r) - 2*SMAJ* atan(sqrt((SMAJ*(-1 + ECC) + r)/(SMAJ + SMAJ*ECC - r))));
    } else {
        res.tbp = -sqrt(-SMAJ/MU)*(sqrt(+r*r - 2*r*SMAJ + (1 - ECC*ECC)*SMAJ*SMAJ) + 2*SMAJ*atanh(sqrt((-r + (1 - ECC)*SMAJ)/(-r + (1 + ECC)*SMAJ))));
    }
    return res;
}
