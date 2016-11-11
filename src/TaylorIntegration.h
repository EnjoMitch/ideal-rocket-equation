#ifndef TAYLOR_INTEGRATION_H
#define TAYLOR_INTEGRATION_H

namespace IdealRocketEquation
{
    class OrbitParam;
    class VesselParam;

    class TaylorIntegration
    {
        public:
            TaylorIntegration();
            virtual ~TaylorIntegration();

            int fun ( const OrbitParam & opar, const VesselParam & vpar, double  R, double tm,  double e[4], double g[3]);

        protected:

        private:
    };
}
#endif // TAYLOR_INTEGRATION_H
