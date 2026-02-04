// Class for configuration of external libraries

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

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
