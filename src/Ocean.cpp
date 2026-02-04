// Global motional induction forward modeling using high-order 
// tetrahedral Nedelec elements

// Author:     Hongbo Yao
// Institute:  School of Geosciences and Info-Physics,
//             Central South University (CSU)
// Email:      yaohongbo@csu.edu.cn
// Date:       2021/08/05

// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#include <mpi.h>
#include <sys/stat.h> // int mkdir(const char*, __mode_t)
#include "em.h"
#include "param_handler.h"
#include "motional_current.h"
#include "fem.h"
#include "post.h"
#include "mfem.hpp"
using namespace mfem;

int main(int argc, char *argv[])
{
   // MPI handler
   int n_procs, myid;
   MPI_Init(&argc, &argv);// MPI progress initial
   MPI_Comm_size(MPI_COMM_WORLD, &n_procs);//MPI progress numbers
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);//MPI progress identification, uniqueness

   double main_tic, main_toc;
   main_tic = MPI_Wtime();//Watch Timer ,here is starting

   // Check program usage
   if (argc < 2) 
   {
      if (myid==0) 
      {
         std::cout << "Usage: " << argv[0] << " input_model_filename\n"; 
      }
      MPI_Finalize();//MPI progress finalize
      return 1;    
   }

   /**
    * Step 1: Input model parameters
    * 
    */
   ParamHandler param_handler(argv[1]);//parameter initial,read file and so on,set conductivity

   MotionalCurrent motional_current1(param_handler.source_file1);//read the id and amplitude of tide source in global element
   MotionalCurrent motional_current2(param_handler.source_file2);//read the id and amplitude of tide source in global element
   MotionalCurrent motional_current3(param_handler.source_file3);//read the id and amplitude of tide source in global element

   int p_order = param_handler.p_order;//finite element's order
   if (myid==0)
   {
      std::cout << "\n*************** Motional induction forward modeling ***************\n";
      std::cout << motional_current1.current_name << "\t" << motional_current1.period_in_hours << " h\n";
      std::cout << motional_current2.current_name << "\t" << motional_current2.period_in_hours << " h\n";
      std::cout << motional_current3.current_name << "\t" << motional_current3.period_in_hours << " h\n";
      std::cout << "p_order: " << p_order << "\n\n";
      if (myid==0) int is_create = mkdir("Forward_solutions",S_IRWXU);
   } 
   /**
    * Step 2: Read mesh
    * 
    */
   double mesh_tic = MPI_Wtime();
   if (myid==0) std::cout << "Reading mesh...\n";
   Mesh* mesh = new Mesh(param_handler.mesh_file.c_str(), 1, 1);

   if (mesh->GetNodes() == NULL)
   {
      mesh->EnsureNodes();
      // mesh->SetCurvature (1);
      
   }   
   // Set attributes by ourself, so we can handle realistic conductivity models
   if (param_handler.marker_type=="element_id")
   {
      for (int i=0; i<mesh->GetNE(); i++)
      {
         mesh->SetAttribute(i,param_handler.marker_vec[i]);
      }
      mesh->SetAttributes();
   }

   // Generate parallel mesh by partitioning a serial mesh using METIS
   ParMesh* pmesh1 = new ParMesh(MPI_COMM_WORLD, *mesh); 
   ParMesh* pmesh2 = new ParMesh(MPI_COMM_WORLD, *mesh); 
   ParMesh* pmesh3 = new ParMesh(MPI_COMM_WORLD, *mesh); 

   delete mesh;
   // Save mesh statistics
   std::ostringstream stats_name;
   stats_name << "Forward_solutions/FEM_mesh" << ".stats";
   std::ofstream stats_ofs(stats_name.str().c_str());
   stats_ofs.precision(8);

   pmesh1->PrintInfo(stats_ofs);
   pmesh2->PrintInfo(stats_ofs);
   pmesh3->PrintInfo(stats_ofs);

   double mesh_toc = MPI_Wtime();
   if (myid==0) std::cout << "MPI_Wtime: " << mesh_toc-mesh_tic << " (s)\n\n";

   /**
    * Step 3: Call parallel finite element solver
    * 
    */
   if (p_order>1) pmesh1->ReorientTetMesh(); //Note: Refinement does not work after a call to this method!
   if (p_order>1) pmesh2->ReorientTetMesh();
   if (p_order>1) pmesh3->ReorientTetMesh();
   // FEM
   //Source1
   FEM fem1(param_handler, motional_current1, *pmesh1);
   fem1.initialize(p_order);
   // fem.print_problem_size();
   // fem.form_linear_system();
   fem1.forward_modeling(1);
   //Source2
   FEM fem2(param_handler, motional_current2, *pmesh2);
   fem2.initialize(p_order);
   // fem.print_problem_size();
   // fem.form_linear_system();
   fem2.forward_modeling(2);
   //Source3
   FEM fem3(param_handler, motional_current3, *pmesh3);
   fem3.initialize(p_order);
   // fem.print_problem_size();
   // fem.form_linear_system();
   fem3.forward_modeling(3);





   // Post processing
   // if (myid==0) std::cout << "\nPost processing...\n";
   // double post_tic = MPI_Wtime();
   
   // std::ostringstream sol_name;
   // sol_name << "Forward_solutions/" << motional_current.current_name 
   //          << "_fem_solution_p" << p_order << ".sol";
   // std::ofstream sol_ofs(sol_name.str().c_str());
   // sol_ofs.precision(8);

   // Post post(param_handler, *pmesh, *(fem.pfes), fem.Ugf->real(), fem.Ugf->imag(), motional_current.omega);
   // post.post_process();
   // post.save_as_one(sol_ofs);

   // double post_toc = MPI_Wtime();
   // if (myid==0) std::cout << "MPI_Wtime: " << post_toc-post_tic << " (s)\n\n";

   // main_toc = MPI_Wtime();
   // if (myid==0) std::cout << "\nTotal MPI_Wtime: " << main_toc-main_tic << " (s)\n";

   delete pmesh1;
   delete pmesh2;
   delete pmesh3;

   MPI_Finalize();
   return 0;
}