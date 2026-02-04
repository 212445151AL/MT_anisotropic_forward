// Simple namespace and EM tools for EM modeling

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

//  ************************** Changelog **************************
// 2021/07/22
// (1) Add function "double tand(double deg)”

#ifndef _EM_H
#define _EM_H

#include <complex>
#include <cmath>

namespace EM 
{
   typedef std::complex<double>    Dcomplex;
   static const Dcomplex   II    = Dcomplex(0.,1.);
   static const double     pi    = 3.1415926535897932384626433832795;
   static const double     mu0   = 4.0*pi*1e-7; 
}
using namespace EM;



inline double cosd(double deg) 
{ 
   return cos(deg/180.0*pi); 
}



inline double sind(double deg) 
{
   return sin(deg/180.0*pi); 
}



inline double tand(double deg) 
{
   return tan(deg/180.0*pi); 
}



// acosd returns [0,180], which is the colatitude range [0 180]
inline double acosd(double rad)
{
   return acos(rad)/pi*180.0;
}



// asind returns degrees [-90,90]
inline double asind(double rad)
{
   return asin(rad)/pi*180.0;
}



// atand returns degrees [-90,90]
inline double atand(double rad)
{
   return atan(rad)/pi*180.0;
}



// atan2d(y,x) computes atan(y/x), returns degrees [-180,180]
inline double atan2d(double y, double x)
{
   return atan2(y,x)/pi*180.0;
}



/**
 * @brief Convert cartesian (x,y,z) to geocentric spherical (r,theta,phi) coordinates
 * 
 * @param[in] sph[0] = r[m], sph[1] = theta[degree], sph[2] = phi[degree]
 * @param[out] car[0] = x[m], car[1] = y[m], car[2] = z[m]
 */
inline void sph2car(const double sph[], double car[])
{
   double r = sph[0];
   double theta = sph[1];
   double phi = sph[2];
   car[0] = r*sind(theta)*cosd(phi);
   car[1] = r*sind(theta)*sind(phi);
   car[2] = r*cosd(theta);
}



/**
 * @brief Convert geocentric spherical (r,theta,phi) to cartesian (x,y,z) coordinates
 * 
 * @param[in] car[0] = x[m], car[1] = y[m], car[2] = z[m]
 * @param[out] sph[0] = r[m], sph[1] = theta[degree], sph[2] = phi[degree]
 */
inline void car2sph(const double car[], double sph[])
{
   double x = car[0];
   double y = car[1];
   double z = car[2];
   sph[0] = sqrt(x*x+y*y+z*z);
   sph[1] = acosd(z/sph[0]);
   sph[2] = atan2d(y,x); // returns degrees [-180,180]
   if (sph[2]<0) sph[2] = sph[2]+360; // convert to [0 360]
}

#endif // _EM_H
