// Parallel post-processing. Compute magnetic induction B from 
// Nedelec electric field E

// Author:     Hongbo Yao
// Institute:  School of Geosciences and Info-Physics,
//             Central South University (CSU)
// Email:      yaohongbo@csu.edu.cn
// Date:       2021/08/03

// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _POST_H
#define _POST_H

#include <fstream>
#include "em.h"
#include "param_handler.h"
#include "mfem.hpp"
using namespace mfem;

enum SpaceRelation{ IN, OUT, ONSURFACE};  //点与四面体关系
class Post
{
public:
   // Parameter handler
   ParamHandler*                 param_handler;

   // Parallel mesh
   ParMesh*                      pmesh;

   // Parallel H(curl) finite-element space
   ParFiniteElementSpace*        pfes;

   // Nedelec solutions
   ParGridFunction               U_re, U_im;

   // Angular frequency
   double omega;
   
   // Local sites Vertex
   Array<double>                 local_sites_r, local_sites_theta, local_sites_phi; // in degrees

   //----------- add by lrj ------------------
   int n_sites;
   Array<double> global_sites_r, global_sites_theta, global_sites_phi;
   Array<double> global_sites_x, global_sites_y, global_sites_z;





   //
   // Bradius, Btheta, Bphi, Re, Im
   Array<double>                 Br_re, Br_im, Bt_re, Bt_im, Bp_re, Bp_im,Er_re, Er_im, Et_re, Et_im, Ep_re, Ep_im;   

public:
   Post(ParamHandler&            param_handler_,
        ParMesh&                 pmesh_,
        ParFiniteElementSpace&   pfes_,
        ParGridFunction&         U_re_,
        ParGridFunction&         U_im_,
        double                   omega_);

   ~Post();

public:
   // Post-processing
   void post_process();
   
   SpaceRelation TestPointInTet(Vector point_, int &tet_id);
   // Operator of points
   Vector CrossProduct(Vector point_1, Vector point_2);
   Vector MinusVector(Vector point_1, Vector point_2);
   double DotProduct(Vector point_1, Vector point_2);
   double Length(Vector point_);
   
   // Save results of local sites
   void save(std::ofstream& out);

   // Save results of global sites in rank 0 (root rank)
   void save_as_one(std::ofstream& out);
};

#endif // _POST_H
