// Error estimator that drives the adaptive FEM solver

// Error estimator is based on the face-jumps of the normal component of the 
// current density, total field formulation (J=sigma*E) is proposed by Zhengyong 
// Ren et al. Here we adapt it for secondary field formulation 
// J = sigma*E_s + (sigma-sigma_0)*E_0

// Reference: 
// Zhengyong Ren, Thomas Kalscheuer, Stewart Greenhalgh and Hansruedi Maurer 
// (2013). A goal-oriented adaptive finite-element approach for plane wave 3-D 
// electromagnetic modelling. Geophysical Journal International, 194, 700-718.

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#ifndef _ERROR_ESTIMATORS_H
#define _ERROR_ESTIMATORS_H


#include "post.h"
#include "param_handler.h"
#include "mfem.hpp"
using namespace mfem;

class Estimator_nJ
{
public:
   Estimator_nJ();
   ~Estimator_nJ();

public:
   // Compute estimated error of solution U (E field or it's dual field)
   Vector compute_estimated_error(ParamHandler& param_handler, 
                                  ParComplexGridFunction& U);

   // Get non-goal-oriented estimated error
   Vector get_estimated_error(ParamHandler& param_handler,  
                              ParComplexGridFunction& U);

   // Get goal-oriented estimated error
   Vector get_goal_estimated_error(ParamHandler& param_handler,  
                                   ParComplexGridFunction& U, 
                                   ParComplexGridFunction& W);
};


/* Tool functions */
// Compute primary field 
//##########delete by lrj ################



// Interpolate field on interface
inline void compute_field_on_surface(ParComplexGridFunction& U,
                                     ParFiniteElementSpace&  pfes,
                                     ElementTransformation&  tet_trans,
                                     const FiniteElement&    tet_fe,
                                     int                     tet_id,
                                     double                  U_re[], 
                                     double                  U_im[],
                                     bool                    MPI_neighbor=false)
{
   int ndof = tet_fe.GetDof();
   int dim = tet_fe.GetDim();
   DenseMatrix shape(ndof, dim);
   tet_fe.CalcPhysVShape(tet_trans, shape);
   Array<int> vdofs;
   Vector dofs_re, dofs_im;
   if (!MPI_neighbor)
   {
      pfes.GetElementVDofs(tet_id, vdofs);
      U.real().GetSubVector(vdofs, dofs_re);
      U.imag().GetSubVector(vdofs, dofs_im);
   }
   else
   {
      pfes.GetFaceNbrElementVDofs(tet_id, vdofs);
      U.real().FaceNbrData().GetSubVector(vdofs, dofs_re);
      U.imag().FaceNbrData().GetSubVector(vdofs, dofs_im);
   }
   
   // U = \sum_{j=1}^{ndof} N_{j}*U_{j}
   for (int j=0; j<ndof; j++)
   {
      U_re[0] += shape(j,0)*dofs_re[j];
      U_re[1] += shape(j,1)*dofs_re[j];
      U_re[2] += shape(j,2)*dofs_re[j];
      U_im[0] += shape(j,0)*dofs_im[j];
      U_im[1] += shape(j,1)*dofs_im[j];
      U_im[2] += shape(j,2)*dofs_im[j];
   }
}

#endif // _ERROR_ESTIMATORS_H
