#pragma once

#include <Eigen/Geometry>

#define PI double(EIGEN_PI)

using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

/* WGS84 Ellipsoidal model parameter */
const double WGS84_WIE = 7.2921151467E-5;       /* Earth Rotation Rate*/
const double WGS84_F   = 0.0033528106647474805; /* flattening */
const double WGS84_RA  = 6378137.0000000000;    /* major semi-axis a */
const double WGS84_RB  = 6356752.3142451793;    /* minor semi-axis b */
const double WGS84_GM0 = 398600441800000.00;    /* Gravitational Constant of the earth */
const double WGS84_E1  = 0.0066943799901413156; /* 1st Eccentricity square */
const double WGS84_E2  = 0.0067394967422764341; /* 2st Eccentricity square */


class Earth {

public:
    /* compute normal Gravity */
    static double gravity(const Vector3d& blh) {

        double sin2 = sin(blh[0]);
        sin2 *= sin2;

        return 9.7803267715 * (1 + 0.0052790414 * sin2 + 0.0000232718 * sin2 * sin2) +
            blh[2] * (0.0000000043977311 * sin2 - 0.0000030876910891) + 0.0000000000007211 * blh[2] * blh[2];
    }

    static Eigen::Vector2d meridianPrimeVerticalRadius(double lat) {
        double tmp, sqrttmp;

        tmp = sin(lat);
        tmp *= tmp;
        tmp = 1 - WGS84_E1 * tmp;
        sqrttmp = sqrt(tmp);

        return { WGS84_RA * (1 - WGS84_E1) / (sqrttmp * tmp), WGS84_RA / sqrttmp }; // [Rm, Rn]
    }

    static double RN(double lat) {
        double sinlat = sin(lat);
        return WGS84_RA / sqrt(1.0 - WGS84_E1 * sinlat * sinlat);   // only Rn
    }

    /* DCM Cne */
    static Matrix3d cne(const Vector3d& blh) { // C_downn_upe
        double coslon, sinlon, coslat, sinlat;

        sinlat = sin(blh[0]);
        sinlon = sin(blh[1]);
        coslat = cos(blh[0]);
        coslon = cos(blh[1]);

        Matrix3d dcm;
        dcm(0, 0) = -sinlat * coslon;
        dcm(0, 1) = -sinlon;
        dcm(0, 2) = -coslat * coslon;

        dcm(1, 0) = -sinlat * sinlon;
        dcm(1, 1) = coslon;
        dcm(1, 2) = -coslat * sinlon;

        dcm(2, 0) = coslat;
        dcm(2, 1) = 0;
        dcm(2, 2) = -sinlat;

        return dcm;
    }

    /* qne from n frame to e frame */
    static Quaterniond qne(const Vector3d& blh) {
        Quaterniond quat;

        double coslon, sinlon, coslat, sinlat;

        coslon = cos(blh[1] * 0.5);
        sinlon = sin(blh[1] * 0.5);
        coslat = cos(-PI * 0.25 - blh[0] * 0.5);
        sinlat = sin(-PI * 0.25 - blh[0] * 0.5);

        quat.w() = coslat * coslon;
        quat.x() = -sinlat * sinlon;
        quat.y() = sinlat * coslon;
        quat.z() = coslat * sinlon;

        return quat;
    }

    /* get BLH from qne and height */
    static Vector3d blh(const Quaterniond& qne, double height) {
        return { -2 * atan(qne.y() / qne.w()) - PI * 0.5, 2 * atan2(qne.z(), qne.w()), height };
    }

    
    static Matrix3d DRi(const Vector3d& blh) {  // Gongmin Yan 4.1.57
        Matrix3d dri = Matrix3d::Zero();

        Eigen::Vector2d rmn = meridianPrimeVerticalRadius(blh[0]);

        dri(0, 0) = 1.0 / (rmn[0] + blh[2]);
        dri(1, 1) = 1.0 / ((rmn[1] + blh[2]) * cos(blh[0]));
        dri(2, 2) = -1;
        return dri;
    }
};