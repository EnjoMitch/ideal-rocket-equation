#include "TaylorIntegration.h"

#include <cmath>
#include "OrbitParam.h"
#include "VesselParam.h"

//  A procedure to use 'automatic differentiation' to allow fast computation of the derivtaives
//  to arbitrarily high order - first used in celestial mechanics problems
//  For further details see, for example, http://www.maia.ub.edu/~angel/taylor/taylor.pdf


// parameters used to define numerical integration order and step-size
#define DELTAT      20.0                        // the integration step-size                                        (s)
#define N           20                          // the order of the Taylor series expansion of the integration

using namespace IdealRocketEquation;

TaylorIntegration::TaylorIntegration(){}
TaylorIntegration::~TaylorIntegration(){}

int taylorIntegation( const OrbitParam & opar, const VesselParam & vpar,
                      double er[N+1],           // time derivatives of the radial component of the orbital eccentricity vector
                      double et[N+1],           // time derivatives of the transbsrese component of the orbital eccentricity vector
                      double th[N+1],           // time derivatives of the angular component of the spacecraft's position
                      double h [N+1],           // time derivatives of the spacecrafts specific angular momentum
                      double f [N+1]            // time derivatives of the spacecraft's acceleration
                     ) {
    const double MU = opar.mu;
    const double VE = vpar.ve;

    // declare the working arrays
    double u01[N], u02[N], u05[N], u06[N], u07[N], u08[N], u09[N], u10[N], u11[N], u12[N], u13[N], u14[N], u15[N];
    double v01[N], v02[N];
    double w01[N], w02[N], w04[N], w05[N], w06[N], w07[N], w08[N], w09[N], w10[N];

    // declare miscellaeous
    double alpha;

    // the integration step - building Taylor series expansions up to 20th order.
    for ( int n = 0; n < N; n++) {

        u01[n] = er[n];
        u02[n] = et[n];
        v01[n] = f [n];
        w01[n] = h [n];


        // calculate f'(t) - (the derivatives of the spacecraft's thrust function)
        v02[n] = 0;
        for (int j = 0; j <= n; j++) {          // f(t)^2
            v02[n] += v01[j] * v01[n-j];
        }

        v02[n] = v02[n] / VE;                   // f(t)^2 / ve

        f[n+1] = v02[n] / (n+1);                // f'(t) = f(t)^2 / ve


        // calculate h'(t) and higher order derivatives
        w02[n] = 0;
        for (int j = 0; j <= n; j++) {          // h(t)^2
            w02[n] += w01[j] * w01[n-j];
        }

        w04[n] = 0;
        for (int j = 0; j <= n; j++) {          // f(t) * h(t)^2
            w04[n] += w02[j] * v01[n-j];
        }

        w05[n] = u01[n];
        if (n == 0) {                           // 1 + e_r(t)
            w05[n] += 1;
        }

        w06[n] = 0;
        for (int j = 0; j <= n; j++) {          // (1 + e_r(t))^2
            w06[n] += w05[j] * w05[n-j];
        }

        w07[n] = 0;
        for (int j = 0; j <= n; j++) {          // e_t(t)^2
            w07[n] += u02[j] * u02[n-j];
        }

        w08[n] = w06[n] + w07[n];               // A = (1 + e_r(t))^2 + e_t(t)^2

        alpha = -0.5;
        if (n == 0) {
            w09[n] = pow( w08[0], alpha );      // A^-1/2
        } else {
            w09[n] = 0;
            for (int j = 0; j <= n-1; j++) {
                w09[n] += (n * alpha - j * (1 + alpha)) * w08[n-j] * w09[j];
            }
            w09[n] = w09[n] / n / w08[0];
        }

        w10[n] = 0;
        for (int j = 0; j <= n; j++) {
            w10[n] += w09[j] * w04[n-j];
        }
        w10[n] = - w10[n] / MU;                 // - mu^-1 * f(t) * h(t)^2 * A^-1/2

        h[n+1] = w10[n] / (n+1);


        // calculate theta'(t) and higher order derivatives
        alpha = -3.0;
        if (n == 0) {
            u05[n] = pow( w01[0], alpha );      // h(t)^-3
        } else {
            u05[n] = 0;
            for (int j = 0; j <= n-1; j++) {
                u05[n] += (n * alpha - j * (1 + alpha)) * w01[n-j] * u05[j];
            }
            u05[n] = u05[n] / n / w01[0];
        }

        u06[n] = 0;
        for (int j = 0; j <= n; j++) {
            u06[n] += u05[j] * w06[n-j];
        }
        u06[n] = u06[n] * MU * MU;              // mu^2 * (1 + e_r(t))^2 * h(t)^-3

        th[n+1] = u06[n] / (n+1);               // th'(t) = u06 / (n+1)


        // calculate e_r'(t) and higher order derivatives
        u07[n] = 0;
        for (int j = 0; j <= n; j++) {
            u07[n] += u06[j] * u02[n-j];        // mu^2 * (1 + e_r(t))^2 * h(t)^-3 * e_t(t)
        }

        u08[n] = 0;
        for (int j = 0; j <= n; j++) {          // - mu^-1 * f(t) * h(t)^2 * A^-1/2 * (1 + e_r(t))
            u08[n] += w10[j] * w05[n-j];
        }

        alpha = -1.0;
        if (n == 0) {
            u09[n] = pow( w01[0], alpha );      // h(t)^-1
        } else {
            u09[n] = 0;
            for (int j = 0; j <= n-1; j++) {
                u09[n] += (n * alpha - j * (1 + alpha)) * w01[n-j] * u09[j];
            }
            u09[n] = u09[n] / n / w01[0];
        }

        u10[n] = 0;
        for (int j = 0; j <= n; j++) {          // - mu^-1 * f(t) * h(t) * A^-1/2 * (1 + e_r(t))
            u10[n] += u09[j] * u08[n-j];
        }

        u11[n] = + u07[n] + 2 * u10[n];         // mu^2*(1+e_r(t))^2*h(t)^-3*e_t(t) - 2*mu^-1*f(t)*h(t)*A^-1/2*(1 + e_r(t))

        er[n+1] = u11[n] / (n+1);               // e_r'(t) = u11 / (n+1)


        //calculate e_t'(t) and higher order derivatives
        u12[n] = 0;
        for (int j = 0; j <= n; j++) {          // mu^2 * (1 + e_r(t))^2 * h(t)^-3 * e_r(t)
            u12[n] += u06[j] * u01[n-j];
        }

        u13[n] = 0;
        for (int j = 0; j <= n; j++) {
            u13[n] += w10[j] * u02[n-j];        // - mu^-1 * f(t) * h(t)^2 * A^-1/2 * e_t(t)
        }

        u14[n] = 0;
        for (int j = 0; j <= n; j++) {          // - mu^-1 * f(t) * h(t) * A^-1/2 * e_t(t)
            u14[n] += u13[j] * u09[n-j];
        }

        u15[n] = - u12[n] + 2 * u14[n];         // - mu^2*(1 + e_r(t))^2*h(t)^-3*e_r(t) - 2*mu^-1*f(t)*h(t)*A^-1/2*e_t(t)

        et[n+1] = u15[n] / (n+1);               // e_t'(t) = u15 / (n+1)
    }

    return 0;
}



int updateInitialValues( double er[N+1],         // time derivatives of the radial component of the orbital eccentricity vector
                         double et[N+1],         // time derivatives of the transbsrese component of the orbital eccentricity vector
                         double th[N+1],         // time derivatives of the angular component of the spacecraft's position
                         double h [N+1],         // time derivatives of the spacecrafts specific angular momentum
                         double f [N+1],         // time derivatives of the spacecraft's acceleration
                         double      dt          // the time step
) {
    double temp;

    temp = 0; for (int j = N; j >=0; j--) {temp = (temp * dt + er[j]);}; er[0] = temp;
    temp = 0; for (int j = N; j >=0; j--) {temp = (temp * dt + et[j]);}; et[0] = temp;
    temp = 0; for (int j = N; j >=0; j--) {temp = (temp * dt + th[j]);}; th[0] = temp;
    temp = 0; for (int j = N; j >=0; j--) {temp = (temp * dt + f [j]);}; f [0] = temp;
    temp = 0; for (int j = N; j >=0; j--) {temp = (temp * dt + h [j]);}; h [0] = temp;
    return 0;
}


int TaylorIntegration::fun ( const OrbitParam & opar, const VesselParam & vpar,
         double  R,                      // the orbital radius at the start of the orbit insertion burn
        double tm,                      // the duration of the orbit insertion burn
        double e[4],                    // the radial and transverse components of the eccentricity vector at the end of the orbital inserion burn
        double g[3]
       ) {

    const double SMAJ = opar.smaj;
    const double ECC = opar.ecc;
    const double MU = opar.mu;

    const double M0 = vpar.m0;
    const double GAMMA = vpar.gamma;
    const double VE = vpar.ve;
    // declare result arrays - normalised Taylor series derivatives up to and including 10th order
    double er[N+1];                     // the component of the eccentricity vector in the radial direction
    double et[N+1];                     // the component of the eccentricity vector in the 'theta' direction
    double th[N+1];                     // the spatial angular coordinate of the spacecraft
    double h [N+1];                     // the specific angular momentum of the spacecraft
    double f [N+1];                     // the spacecraft's thrust function g * ve / (m0 - g * t)

    // initalise the values to start the integration processteration process
    er[0] = SMAJ * (1 - ECC * ECC) / R - 1;
    et[0] = sqrt( (1 - ECC * ECC) * (SMAJ * ECC + SMAJ - R) * (SMAJ * ECC - SMAJ + R)) / R;
    th[0] = 0.0;
    h [0] = sqrt(MU * SMAJ * (1 - ECC * ECC));
    f [0] = GAMMA * VE / M0;

    // declare miscellanous quantities
    int    status;
    double t, temp;

    t      = tm;
    // carry out the integration step nine times - with each integration step spanning 1000 seconds
    while (t >= DELTAT) {
        status = taylorIntegation    (opar, vpar, er, et, th, h, f);
        status = updateInitialValues             (er, et, th, h, f, DELTAT);
        t     -= DELTAT;
    }
    status = taylorIntegation    (opar, vpar, er, et, th, h, f);
    status = updateInitialValues             (er, et, th, h, f, t);

    // calculate the components of the eccentricity vector at the end of the orbit insertion burn
    e[0]  = er[0];
    e[1]  = et[0];

    // calculate the firs derivatives with respect to time of the eccentricity vector at the end of the orbit insertion burn
    temp = 0; for (int j = N; j >=1; j--) {temp = (temp * t + j * er[j]);}; e[2] = temp;
    temp = 0; for (int j = N; j >=1; j--) {temp = (temp * t + j * et[j]);}; e[3] = temp;

    g[0] = th[0];
    g[1] = f [0];
    g[2] = h [0];

    return 0;
}
