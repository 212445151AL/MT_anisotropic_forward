// class SolverParm： Solver's parameters handler. Including MFEM's FGMRES, 
// Hypre's PCG, and Hypre's auxiliary-space Maxwell solver (AMS)

// class FGMRES: Flexible GMRES iterative solver. This class is almost an exact 
// copy of MFEM's FGMRESSolver in linalg/solvers.cpp, see below for details

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#ifndef _SOLVER_H
#define _SOLVER_H

#include <string>
#include <mpi.h>
#include "mfem.hpp"
using namespace mfem;

// class SolverParm： Solver's parameters handler. Including MFEM's FGMRES, 
// Hypre's PCG, and Hypre's auxiliary-space Maxwell solver (AMS)
class FGMRES;
class SolverParm
{
public:
   //------------------------------ FGMRES options -----------------------------
   int fgmres_maxit;
   double fgmres_primal_tol;  // relative tolerance for primal problem, e.g. 1e-6
   double fgmres_dual_tol;    // relative tolerance for dual problem, can be relaxed, e.g. 1e-3
   int fgmres_restart;
   int fgmres_print_level;

   //------------------------------ preconditioner options ---------------------
   int prec_type;             // 0-AMS, 1-PCG-AMS

   //------------------------------ PCG options --------------------------------
   int pcg_maxit;             // maximum iterations, e.g. 20
   double pcg_tol;            // relative tolerance, e.g. 1e-3
   int pcg_print_level;

   //--------------------- auxiliary-space Maxwell solver options --------------
   // Solver options
   int ams_cycle_type;        // AMS circle type
   int ams_maxit;             // maximum AMS circles ,1 means use as a preconditioner
   double ams_tol;            // the convergence tolerance
   int ams_print_level;       // The default is 1 (print residual norm at each step)

   // Smoothing options for A
   int A_relax_type;          // 1 = l1-scaled Jacobi, 2 = l1-scaled block symmetric 
                              // Gauss-Seidel/SSOR, 3 = Kaczmarz, 4 = truncated 
                              // version of l_1-scaled block symmetric 
                              // Gauss-Seidel/SSOR 16 - Chebyshev 
   int A_relax_sweeps;        // relaxation sweeps (relax_times) on each level
   double A_relax_weight;     // damping parameter
   double A_relax_omega;      // SSOR coefficient

   // AMG options for B_Pi, B_G, we set the same options 
   // 1. AMG coarsening options
   int amg_coarsen_type;      // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP
   int amg_agg_levels;        // number of aggressive coarsening levels
   double amg_theta;          // strength threshold: 0.25, 0.5, 0.8
   // 2. AMG interpolation options
   int amg_interp_type;       // 6 = extended+i, 0 = classical
   int amg_Pmax;              // max number of elements per row in P
   // 3. AMG relaxation options
   int amg_relax_type;        // 8 = l1-GS, 6 = symm. GS, 3 = GS, 18 = l1-Jacobi

public:
   SolverParm(std::string& solver_parm_file);
   ~SolverParm();

public:
   // Load solver's parameters from given file
   void load_solver_parm(std::string& solver_parm_file);

   // Output to file
   void output_to_file(std::string out="input_solver_parm.log");

   // Set PCG solver, called when prec_type=1
   void set_pcg_parm(HyprePCG& pcg);

   // Set AMS solver, called when prec_type=0
   void set_ams_parm(HypreAMS& ams);
};




// class FGMRES: Flexible GMRES iterative solver. This class is almost an exact 
// copy of MFEM's FGMRESSolver in linalg/solvers.cpp. I do this for extending 
// the FGMRESSolver::Mult function, so I can compute and output the relative 
// residual during the iteration procedure.

// The following three inline functions are the direct copies of MFEM's FGMRESSolver 
// in linalg/solvers.cpp. For the FGMRESSolver::Mult function, I just did a minor
// modifications to compute and output the relative residual.

// The following are a copy of the header's introduction in MFEM's linalg/solvers.cpp:

// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.

// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.

// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

class FGMRES: public FGMRESSolver
{
public:
   FGMRES();
   FGMRES(MPI_Comm comm_);
   
   void Mult(const Vector &b, Vector &x) const;
};



// ---------------------------------------------------------------------------------
// The following three inline functions are the direct copies of MFEM's FGMRESSolver 
// in linalg/solvers.cpp.
inline void GeneratePlaneRotation(double &dx, double &dy,
                                  double &cs, double &sn)
{
   if (dy == 0.0)
   {
      cs = 1.0;
      sn = 0.0;
   }
   else if (fabs(dy) > fabs(dx))
   {
      double temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
   }
   else
   {
      double temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
   }
}



inline void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
   double temp = cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}



inline void Update(Vector &x, int k, DenseMatrix &h, Vector &s,
                   Array<Vector*> &v)
{
   Vector y(s);

   // Backsolve:
   for (int i = k; i >= 0; i--)
   {
      y(i) /= h(i,i);
      for (int j = i - 1; j >= 0; j--)
      {
         y(j) -= h(j,i) * y(i);
      }
   }

   for (int j = 0; j <= k; j++)
   {
      x.Add(y(j), *v[j]);
   }
}

#endif // _SOLVER_H
