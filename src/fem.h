// Global motional induction finite-element forward modeling

// Author:     Hongbo Yao
// Institute:  School of Geosciences and Info-Physics,
//             Central South University (CSU)
// Email:      yaohongbo@csu.edu.cn
// Date:       2021/08/05

// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _FEM_H
#define _FEM_H



#include "config.h"
#include <mpi.h>
#include "param_handler.h"
#include "motional_current.h"
#include "solver.h"
#include "mfem.hpp"
#include "PwMatrixCoefficient.h"
//----------- add by lrj ------------

#include <fstream>
#include "error_estimators.h"
#include <vector>
using namespace mfem;

class FEM
{
public:
   // MPI variables
   int myid;
   int n_procs;
   MPI_Comm comm;

   // Class ParamHandler for handlering global EM input parameters
   ParamHandler* param_handler;
   
   // Motional induction current source
   MotionalCurrent* motional_current;

   // Parallel mesh
   ParMesh* pmesh;

   // Parallel H(curl) finite-element space
   ParFiniteElementSpace* pfes;

   // Ccefficients of governing curl-curl equation
   Coefficient* stiffness_coef;  // curl-curl term
   Coefficient* mass_coef;       // mass term
   PwMatrixCoefficient* mass_coef_tensor;
   Coefficient* neg_stiffness_coef;
   Coefficient* neg_mass_coef;

   // Weak formulation B(U,V) = D(V)
   ParSesquilinearForm* B;
   ParComplexLinearForm* D;
   ParComplexGridFunction* Ugf;

   // ------------------ add by lrj -------------------
   // Weak formulation B(W,V) = L(V) for dual problem 
   ParComplexLinearForm* L;
   ParComplexGridFunction* Wgf;


   // Linear system KU = F
   OperatorHandle K;
   Vector F, U;

   Vector C, W;


   // Auxiliary-space block-diagonal preconditioner
   ParBilinearForm* precond;  // preconditioning bilinearform
   HypreParMatrix A;          // preconditioning matrix

   // Dirichlet boundary condition on the outer boundary
   Array<int> dirichlet_bdr;        // bool marker for Dirichlet bc
   Array<int> dirichlet_tdof_list;  // dofs id for Dirichlet bc

   //------------------add by lrj------------------------------
   bool solve_dual_problem;
   Vector estimated_errors;
   Vector estimated_solution_errors;





public:
   FEM(ParamHandler&    para_handler_,
       MotionalCurrent& motional_current_,
       ParMesh&         pmesh_);
   virtual ~FEM();

   void initialize(int p_order=1);

public:
   // Setup ccefficients of governing curl-curl equation
   void setup_coefficients();

   // Print time in seconds(sec), minutes(min), and hours(hr)
   void print_time(double time_in_seconds, std::string time_str="MPI_Wtime: ");

   // Problem size
   HYPRE_Int get_problem_dofs();
   void print_problem_size();

   // Assemble linear form (source term)
   void assemble_linearform();

   // Assemble all bilinearforms and linearforms and form linear system
   void form_linear_system();

   // Setup auxiliary-space block-diagonal preconditioner
   void form_preconditioner(SolverParm& solver_parm);

   // Solve the linear system
   void solve();

   //----------add by lrj --------------------------
   void assemble_system_matrix();
   void assemble_rhs();
   // no need void assemble_primal_rhs();  =====      void assemble_linearform();
   void assemble_dual_rhs(); // L(V) = 1/Omega_s*\iiint_Omega_s (V dot I), I=(1,1,1) when ip in \Omega_s
   // Local sites and tets that contains local sites
   Array<int> local_sites_tets;
   //Array<int> local_sites_id;
   Array<double> local_sites_r, local_sites_theta, local_sites_phi;

   // Find local sites tets
   void find_sites_tets();
   void find_sites_tets_with_gslib();
   void forward_modeling(int source_id);
   void error_estimating();
   void update();
   // Print VTK mesh and estimated error
   void print_vtk_as_one(std::ofstream& out);









};

#endif // _FEM_H