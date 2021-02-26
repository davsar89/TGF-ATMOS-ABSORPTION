#pragma once

#include <cmath>

// A bunch of static functions to convert coordinates

// x y z are in meters, lat lon are in degrees, alt in m
namespace geod_conv {
    static double kSmA = 6378137;
    static double kSemiminorAxis = 6356752.3142;
    static double kFES = 6.69437999014 * 0.001;
    static double kSES = 6.73949674228 * 0.001;

    static double Esq = kSmA * kSmA - kSemiminorAxis * kSemiminorAxis;
    static double kFESxE = kFES * Esq;
    static double kFES_x_kFES = kFES * kFES;
    static double kSmaAxkSmA = kSmA * kSmA;
    static double kSmA2 = kSemiminorAxis * kSemiminorAxis;
    static double kSmAxkAokSmA = kSmA2 / kSmA;

    static float kSemimajorAxis_f = 6378137.0f;
    static float kSemiminorAxis_f = 6356752.3142f;
    static float kFirstEccentricitySquared_f = 6.69437999014f * 0.001f;
    static float kSecondEccentricitySquared_f = 6.73949674228f * 0.001f;
    static float Esq_f = kSemimajorAxis_f * kSemimajorAxis_f - kSemiminorAxis_f * kSemiminorAxis_f;
    //static double kFlattening                = 1 / 298.257223563;

    class GeodeticConverter {
    public:

        GeodeticConverter() {
        }

        ~GeodeticConverter() {
        }

        // Default copy constructor and assignment operator are OK.

        static void
        geodetic2ecef(const double &latitude, const double &longitude, const double &altitude, double &x, double &y,
                      double &z)
        // x y z are in meters, lat lon are in degrees, alt in m
        {
            // Convert geodetic coordinates to ECEF.
            // http://code.google.com/p/pysatel/source/browse/trunk/coord.py?r=22
            double lat_rad = deg2Rad(latitude);
            double lon_rad = deg2Rad(longitude);
            double coslat = 0, sinlat = 0, sinlon = 0, coslon = 0;

#if defined(__linux__) && defined(__GNUG__) // if linux and GCC
            sincos(lat_rad, &sinlat, &coslat);
            sincos(lon_rad, &sinlon, &coslon);
#else
            sinlat = sin(lat_rad);
            coslat = cos(lat_rad);
            sinlon = sin(lon_rad);
            coslon = cos(lon_rad);
#endif

            double xi = sqrt(1 - kFES * sinlat * sinlat);
            x = (kSmA / xi + altitude) * coslat * coslon;
            y = (kSmA / xi + altitude) * coslat * sinlon;
            z = (kSmA / xi * (1 - kFES) + altitude) * sinlat;
        }

        static void
        ecef2Geodetic(const double &x, const double &y, const double &z, double &latitude, double &longitude,
                      double &altitude) {
            // Convert ECEF coordinates to geodetic coordinates.
            // J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates
            // to geodetic coordinates," IEEE Transactions on Aerospace and
            // Electronic Systems, vol. 30, pp. 957-961, 1994.

            double z2 = z * z;

            double r2 = x * x + y * y;

            double r = sqrt(r2);

            double F = 54.0 * kSmA2 * z2;

            double G = r * r + (1.0 - kFES) * z2 - kFESxE;

            double C = (kFES_x_kFES * F * r2) / pow(G, 3.0);

            double S = cbrt(1.0 + C + sqrt(C * C + 2.0 * C));

            double P = F / (3.0 * pow((S + 1.0 / S + 1.0), 2.0) * G * G);

            double Q = sqrt(1.0 + 2.0 * kFES_x_kFES * P);

            double r_0 = -(P * kFES * r) / (1.0 + Q)
                         + sqrt(
                    0.5 * kSmaAxkSmA * (1.0 + 1.0 / Q) -
                    P * (1.0 - kFES) * z2 / (Q * (1.0 + Q)) - 0.5 * P * r2);

            double U = sqrt(pow((r - kFES * r_0), 2.0) + z2);

            double V = sqrt(pow((r - kFES * r_0), 2.0) + (1.0 - kFES) * z2);

            double Z_0 = kSmAxkAokSmA * z / V;

            altitude = U * (1.0 - kSmA2 / (kSmA * V));

            latitude = rad2Deg(atan((z + kSES * Z_0) / r));

            longitude = rad2Deg(atan2(y, x));
        }

        static void
        ecef2Geodetic_float(const float &x, const float &y, const float &z, float &latitude, float &longitude,
                            float &altitude) {
            // Convert ECEF coordinates to geodetic coordinates.
            // J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates
            // to geodetic coordinates," IEEE Transactions on Aerospace and
            // Electronic Systems, vol. 30, pp. 957-961, 1994.
            float r = sqrt(x * x + y * y);
            float F = 54.0 * kSemiminorAxis_f * kSemiminorAxis_f * z * z;
            float G = r * r + (1.0 - kFirstEccentricitySquared_f) * z * z - kFirstEccentricitySquared_f * Esq;
            float C = (kFirstEccentricitySquared_f * kFirstEccentricitySquared_f * F * r * r) / pow(G, 3);
            float S = cbrt(1.0 + C + sqrt(C * C + 2.0 * C));
            float P = F / (3.0 * pow((S + 1.0 / S + 1.0), 2.0) * G * G);
            float Q = sqrt(1.0 + 2.0 * kFirstEccentricitySquared_f * kFirstEccentricitySquared_f * P);
            float r_0 = -(P * kFirstEccentricitySquared_f * r) / (1.0 + Q) + sqrt(
                    0.5 * kSemimajorAxis_f * kSemimajorAxis_f * (1.0 + 1.0 / Q) -
                    P * (1.0 - kFirstEccentricitySquared_f) * z * z / (Q * (1.0 + Q)) - 0.5 * P * r * r);
            float U = sqrt(pow((r - kFirstEccentricitySquared_f * r_0), 2.0) + z * z);
            float V = sqrt(pow((r - kFirstEccentricitySquared_f * r_0), 2.0) + (1.0 - kFirstEccentricitySquared_f) * z * z);
            float Z_0 = kSemiminorAxis_f * kSemiminorAxis_f * z / (kSemimajorAxis_f * V);
            altitude = U * (1.0 - kSemiminorAxis_f * kSemiminorAxis_f / (kSemimajorAxis_f * V));
            latitude = rad2Deg(atan((z + kSecondEccentricitySquared_f * Z_0) / r));
            longitude = rad2Deg(atan2(y, x));
        }

    private:

        static inline double rad2Deg(const double &radians) {
            return (radians / 3.14159265359) * 180.0;
        }

        static inline double deg2Rad(const double &degrees) {
            return (degrees / 180.0) * 3.14159265359;
        }
    };

    // class GeodeticConverter
} // namespace geodetic_conv
