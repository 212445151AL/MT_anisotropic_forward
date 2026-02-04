// Class for configuration of external libraries

// Author:     Hongbo Yao
// Institute:  School of Geosciences and Info-Physics,
//             Central South University (CSU)
// Email:      yaohongbo@csu.edu.cn
// Date:       2021/12/09

// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _CONFIG_H
#define _CONFIG_H

// https://github.com/Nek5000/gslib
// Robust interpolation, I use GSLIB-FindPoints to find the elements 
// that contain the electromagnetic site(station) point
// If you want to use the default FindPoints function in MFEM library, 
// please comment out this line
#define USE_GSLIB

// Todo:
// Add mumps direct solver

#endif // _CONFIG_H