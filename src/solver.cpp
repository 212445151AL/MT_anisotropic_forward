// class SolverParm： Solver's parameters handler. Including MFEM's FGMRES, 
// Hypre's PCG, and Hypre's auxiliary-space Maxwell solver (AMS)

// class FGMRES: Flexible GMRES iterative solver. This class is almost an exact 
// copy of MFEM's FGMRESSolver in linalg/solvers.cpp, see below for details

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#include "solver.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// class SolverParm： Solver's parameters handler. Including MFEM's FGMRES, 
// Hypre's PCG, and Hypre's auxiliary-space Maxwell solver (AMS)
SolverParm::SolverParm(std::string& solver_parm_file)
{
   this->load_solver_parm(solver_parm_file);
   // this->output_to_file();
}



SolverParm::~SolverParm()
{

}



void SolverParm::load_solver_parm(std::string& solver_parm_file)
{
   std::ifstream in_stream(solver_parm_file);
   std::string tmp;

   //------------------------------ FGMRES options -----------------------------
   in_stream >> tmp >> fgmres_maxit;
   in_stream >> tmp >> fgmres_primal_tol;
   in_stream >> tmp >> fgmres_dual_tol;
   in_stream >> tmp >> fgmres_restart;
   in_stream >> tmp >> fgmres_print_level;

   //------------------------------ preconditioner options ---------------------
   in_stream >> tmp >> prec_type;
   if (prec_type!=0 && prec_type!=1 && prec_type!=2)
      throw std::invalid_argument("In SolverParm::load_solver_parm(): Wrong preconditioner type!\n");

   //------------------------------ PCG options --------------------------------
   in_stream >> tmp >> pcg_maxit;
   in_stream >> tmp >> pcg_tol;
   in_stream >> tmp >> pcg_print_level;

   //--------------------- auxiliary-space Maxwell solver options --------------
   // Solver options
   in_stream >> tmp >> ams_cycle_type;
   in_stream >> tmp >> ams_maxit;
   in_stream >> tmp >> ams_tol;
   in_stream >> tmp >> ams_print_level;
   
   // Smoothing options for A
   in_stream >> tmp >> A_relax_type;
   in_stream >> tmp >> A_relax_sweeps;
   in_stream >> tmp >> A_relax_weight;
   in_stream >> tmp >> A_relax_omega;

   // AMG options for B_Pi, B_G, we set the same options 
   in_stream >> tmp >> amg_coarsen_type;
   in_stream >> tmp >> amg_agg_levels;
   in_stream >> tmp >> amg_theta;
   in_stream >> tmp >> amg_interp_type;
   in_stream >> tmp >> amg_Pmax;
   in_stream >> tmp >> amg_relax_type;

   in_stream.close();
}

void SolverParm::output_to_file(std::string out)
{
   std::ofstream out_stream(out);

   //------------------------------ FGMRES options -----------------------------
   out_stream << fgmres_maxit << "\n";
   out_stream << fgmres_primal_tol << "\n";
   out_stream << fgmres_dual_tol << "\n";
   out_stream << fgmres_restart << "\n";
   out_stream << fgmres_print_level << "\n\n";

   //------------------------------ preconditioner options ---------------------
   out_stream << prec_type << "\n\n";

   //------------------------------ PCG options --------------------------------
   out_stream << pcg_maxit << "\n";
   out_stream << pcg_tol << "\n";
   out_stream << pcg_print_level << "\n\n";

   //--------------------- auxiliary-space Maxwell solver options --------------
   // Solver options
   out_stream << ams_cycle_type << "\n";
   out_stream << ams_maxit << "\n";
   out_stream << ams_tol << "\n";
   out_stream << ams_print_level << "\n\n";
   
   // Smoothing options for A
   out_stream << A_relax_type << "\n";
   out_stream << A_relax_sweeps << "\n";
   out_stream << A_relax_weight << "\n";
   out_stream << A_relax_omega << "\n\n";

   // AMG options for B_Pi, B_G, we set the same options 
   out_stream << amg_coarsen_type << "\n";
   out_stream << amg_agg_levels << "\n";
   out_stream << amg_theta << "\n";
   out_stream << amg_interp_type << "\n";
   out_stream << amg_Pmax << "\n";
   out_stream << amg_relax_type << "\n";

   out_stream.close();
}

void SolverParm::set_pcg_parm(HyprePCG& pcg)
{
   pcg.SetMaxIter(pcg_maxit);
   pcg.SetTol(pcg_tol);
   pcg.SetPrintLevel(pcg_print_level);
}

/*https://hypre.readthedocs.io/en/latest/solvers-ams.html*/
void SolverParm::set_ams_parm(HypreAMS& ams)
{
   HYPRE_AMSSetCycleType(ams, ams_cycle_type);
   HYPRE_AMSSetMaxIter(ams, ams_maxit);
   HYPRE_AMSSetTol(ams, ams_tol);
   HYPRE_AMSSetPrintLevel(ams, ams_print_level);
   // set additional AMS options
   HYPRE_AMSSetSmoothingOptions(ams, 
                                A_relax_type, 
                                A_relax_sweeps, 
                                A_relax_weight, 
                                A_relax_omega);
   HYPRE_AMSSetAlphaAMGOptions(ams, 
                               amg_coarsen_type, 
                               amg_agg_levels, 
                               amg_relax_type,
                               amg_theta, 
                               amg_interp_type, 
                               amg_Pmax);
   HYPRE_AMSSetBetaAMGOptions(ams, 
                              amg_coarsen_type, 
                              amg_agg_levels, 
                              amg_relax_type,
                              amg_theta, 
                              amg_interp_type, 
                              amg_Pmax);

   HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams, amg_relax_type);
   HYPRE_AMSSetBetaAMGCoarseRelaxType(ams, amg_relax_type);

}

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

FGMRES::FGMRES():FGMRESSolver()
{

}



FGMRES::FGMRES(MPI_Comm comm_):FGMRESSolver(comm_)
{

}


void FGMRES::Mult(const Vector &b, Vector &x) const
{
   DenseMatrix H(m+1,m);
   Vector s(m+1), cs(m+1), sn(m+1);
   Vector r(b.Size());

   int i, j, k;


   if (iterative_mode)
   {
      oper->Mult(x, r);
      subtract(b,r,r);
   }
   else
   {
      x = 0.;
      r = b;
   }
   double beta = Norm(r);  // beta = ||r||
   MFEM_ASSERT(IsFinite(beta), "beta = " << beta);

   final_norm = std::max(rel_tol*beta, abs_tol);

   if (beta <= final_norm)
   {
      final_norm = beta;
      final_iter = 0;
      converged = 1;
      return;
   }

   /*
   if (print_level == 1)
   {
      mfem::out << "   Pass: " << std::setw(2) << 1
                << "   Iteration: " << std::setw(3) << 0
                << "  || r || = " << beta << std::endl;
   }
   */
   
   // Hongbo Yao, 2020/09/06
   double norm_b = Norm(b);
   if (print_level == 1)
   {
      mfem::out << "   Pass: " << std::setw(2) << 1
                << "   Iteration: " << std::setw(3) << 0
                << "   ||r|| = " << beta << "\t"
                << "||r||/||b|| = " << beta/norm_b << std::endl;
   }

   Monitor(0, beta, r, x);

   Array<Vector*> v(m+1);
   Array<Vector*> z(m+1);
   for (i= 0; i<=m; i++)
   {
      v[i] = NULL;
      z[i] = NULL;
   }

   j = 1;
   while (j <= max_iter)
   {
      if (v[0] == NULL) { v[0] = new Vector(b.Size()); }
      (*v[0]) = 0.0;
      v[0] -> Add (1.0/beta, r);   // v[0] = r / ||r||
      s = 0.0; s(0) = beta;

      for (i = 0; i < m && j <= max_iter; i++, j++)
      {

         if (z[i] == NULL) { z[i] = new Vector(b.Size()); }
         (*z[i]) = 0.0;

         if (prec)
         {
            prec->Mult(*v[i], *z[i]);
         }
         else
         {
            (*z[i]) = (*v[i]);
         }
         oper->Mult(*z[i], r);

         for (k = 0; k <= i; k++)
         {
            H(k,i) = Dot( r, *v[k]); // H(k,i) = r * v[k]
            r.Add(-H(k,i), (*v[k])); // r -= H(k,i) * v[k]
         }

         H(i+1,i)  = Norm(r);       // H(i+1,i) = ||r||
         if (v[i+1] == NULL) { v[i+1] = new Vector(b.Size()); }
         (*v[i+1]) = 0.0;
         v[i+1] -> Add (1.0/H(i+1,i), r); // v[i+1] = r / H(i+1,i)

         for (k = 0; k < i; k++)
         {
            ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
         }

         GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

         double resid = fabs(s(i+1));
         MFEM_ASSERT(IsFinite(resid), "resid = " << resid);
         /*
         if (print_level == 1)
         {
            mfem::out << "   Pass: " << std::setw(2) << (j-1)/m+1
                      << "   Iteration: " << std::setw(3) << j
                      << "  || r || = " << resid << std::endl;
         }
         */
         // Hongbo Yao, 2020/09/06
         if (print_level == 1)
         {
            mfem::out << "   Pass: " << std::setw(2) << (j-1)/m+1
                      << "   Iteration: " << std::setw(3) << j
                      << "   ||r|| = " << resid << "\t"
                      << "||r||/||b|| = " << resid/norm_b << std::endl;
         }
         Monitor(j, resid, r, x, resid <= final_norm);

         if (resid <= final_norm)
         {
            Update(x, i, H, s, z);
            final_norm = resid;
            final_iter = j;
            converged = 1;

            if (print_level == 2)
            {
               mfem::out << "Number of FGMRES iterations: " << final_iter << std::endl;
            }

            for (i= 0; i<=m; i++)
            {
               if (v[i]) { delete v[i]; }
               if (z[i]) { delete z[i]; }
            }
            return;
         }
      }

      if (print_level == 1)
      {
         mfem::out << "Restarting..." << std::endl;
      }

      Update(x, i-1, H, s, z);

      oper->Mult(x, r);
      subtract(b,r,r);
      beta = Norm(r);
      MFEM_ASSERT(IsFinite(beta), "beta = " << beta);
      if (beta <= final_norm)
      {
         final_norm = beta;
         final_iter = j;
         converged = 1;

         if (print_level == 2)
         {
            mfem::out << "Number of FGMRES iterations: " << final_iter << std::endl;
         }

         for (i= 0; i<=m; i++)
         {
            if (v[i]) { delete v[i]; }
            if (z[i]) { delete z[i]; }
         }
         return;
      }
   }

   for (i = 0; i <= m; i++)
   {
      if (v[i]) { delete v[i]; }
      if (z[i]) { delete z[i]; }
   }
   converged = 0;

   if (print_level >= 0)
   {
      mfem::out << "FGMRES: No convergence!" << std::endl;
   }

   return;
}
