// Global motional induction finite-element forward modeling

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#include "fem.h"
#include "em.h"
#include <stdexcept>
#include <algorithm> 
#include "post.h"
#include "PwMatrixCoefficient.h"


FEM::FEM(ParamHandler&  para_handler_,
         MotionalCurrent& motional_current_,
         ParMesh&      pmesh_):
         myid(0),
         n_procs(1),
         param_handler(&para_handler_),
         motional_current(&motional_current_),
         pmesh(&pmesh_),
         pfes(NULL),
         stiffness_coef(NULL),
         mass_coef(NULL),
         mass_coef_tensor(NULL),
         B(NULL),
         D(NULL),
         Ugf(NULL),                                  //  L Wgf add by lrj
         L(NULL),
         Wgf(NULL),
         precond(NULL)
{
   if (param_handler->amr_method == 0) // goal-oriented AMR   //  add by lrj 
      solve_dual_problem = true;
   else // non-goal-oriented and global AMR
      solve_dual_problem = false;
}

FEM::~FEM()
{
   delete pfes;
   delete B;
   delete D;
   delete Ugf;
   delete L;                                         // add by lrj :delete L Wgf
   delete Wgf;
   delete precond;
   delete stiffness_coef;
   delete mass_coef;
   delete mass_coef_tensor;
}

void FEM::setup_coefficients()
{
   // Stiffness coefficient of curl-curl term, 1.0
   stiffness_coef = new ConstantCoefficient(1.0);
   neg_stiffness_coef = new ConstantCoefficient(-1.0);

   // Mass coefficient of massterm, omega*mu0*sigma
   int n_attris = pmesh->attributes.Size();
   int dim = pmesh->SpaceDimension();
   
   // 创建分片常数矩阵系数
   mass_coef_tensor = new PwMatrixCoefficient(dim, dim);
   
   // 为不同的网格属性设置不同的矩阵系数
   for (int i = 0; i < n_attris; i++) 
   {
      int attr = pmesh->attributes[i];
      double omega = motional_current->omega;
      
      // 获取该属性的电导率张量
      DenseMatrix sigma_tensor(dim, dim);
      param_handler->get_elem_conductivity_tensor(attr, sigma_tensor);
      
      // 缩放：ω * μ₀ * σ_tensor
      sigma_tensor *= (omega * EM::mu0);
      
      // 创建矩阵系数并添加到分片系数中
      MatrixConstantCoefficient* coef = new MatrixConstantCoefficient(sigma_tensor);
      mass_coef_tensor->AddCoefficient(attr, coef);
   } 
}

void FEM::initialize(int p_order)
{
   // MPI variables
   comm = pmesh->GetComm();
   MPI_Comm_size(comm, &n_procs);
   MPI_Comm_rank(comm, &myid);
   
   // Parallel H(curl) finite-element space
   int dim = pmesh->SpaceDimension();
   pfes = new ParFiniteElementSpace(pmesh, new ND_FECollection(p_order, dim));
   
   // Setup coefficients
   this->setup_coefficients();



   // Bilinear form
   B = new ParSesquilinearForm(pfes, ComplexOperator::BLOCK_SYMMETRIC);
   B->AddDomainIntegrator(new CurlCurlIntegrator(*stiffness_coef), NULL);
   B->AddDomainIntegrator(NULL, new VectorFEMassIntegrator(*mass_coef_tensor));
   
   // Linear form
   D = new ParComplexLinearForm(pfes, ComplexOperator::BLOCK_SYMMETRIC);
   D->Vector::operator=(0.0);

   // Gridfunction of primary problem
   Ugf = new ParComplexGridFunction(pfes);  *Ugf = 0.0;

   // ------------- add by lrj --------------------------
   // Linear form of dual problem
   L = new ParComplexLinearForm(pfes, ComplexOperator::BLOCK_SYMMETRIC);
   L->Vector::operator=(0.0);

   // Gridfunction of dual problem
   Wgf = new ParComplexGridFunction(pfes);  *Wgf = 0.0;



   // Preconditioner blinearform
   precond = new ParBilinearForm(pfes);
   precond->AddDomainIntegrator(new CurlCurlIntegrator(*stiffness_coef));
   precond->AddDomainIntegrator(new VectorFEMassIntegrator(*mass_coef_tensor));
   // precond->AddDomainIntegrator(new VectorFEMassIntegrator(*neg_mass_coef));

   // Boundary elements attributes, extracted from Gmsh mesh generator
   int n_bdr_attributes = pmesh->bdr_attributes.Size();
   if (n_bdr_attributes<1 && myid==0)
      throw std::invalid_argument("In RegionalFEM::initialize(): The meshes should have at least one boundary attribute!\n");
   // Dirichlet boundary condition on the outer boundary
   int dirichlet_bdr_location = -1;
   for (int i=0; i<n_bdr_attributes; i++)
   {
      if (pmesh->bdr_attributes[i] == param_handler->bdr_marker)
      {
         dirichlet_bdr_location = i;
         break;
      }
   }
   if (dirichlet_bdr_location==-1 && myid==0)
      throw std::invalid_argument("In FEM::initialize(): Can not find dirichlet boundary attribute!\n");
   dirichlet_bdr.SetSize(n_bdr_attributes);
   dirichlet_bdr = 0;  
   dirichlet_bdr[dirichlet_bdr_location] = 1;
   pfes->GetEssentialTrueDofs(dirichlet_bdr, dirichlet_tdof_list);
}

void FEM::print_time(double time_in_seconds, std::string time_str)
{
   // In hour-min-sec
   // int seconds = (int)time_in_seconds;
   // int hr = seconds/3600;
   // seconds %= 3600;
   // int min = seconds/60;
   // double sec = time_in_seconds-hr*3600-min*60;
   // std::cout << time_str << hr << " hr " << min << " min " << sec <<  " sec\n";

   // In seconds
   std::cout << time_str << time_in_seconds << " (s)\n";
}

HYPRE_Int FEM::get_problem_dofs()
{
   return 2*pfes->GlobalTrueVSize();
}

void FEM::print_problem_size()
{
   HYPRE_Int size = pfes->GlobalTrueVSize();
   long n_global_elems = pmesh->GetGlobalNE();
   if (myid == 0)
   {
      std::cout << "Number of tetrahedral elements: " << n_global_elems << "\n";
      std::cout << "Number of complex-valued unknowns: " << size << "\n";
      std::cout << "Degrees of Freedom (DOFs): " << 2*size << "\n";
   }
}

//------------  add by lrj -----------
void FEM::assemble_system_matrix()
{
   B->Assemble();
   B->Finalize();
}


//------------  add by lrj -----------
void FEM::assemble_rhs()
{
   D->Vector::operator=(0.0); // VERY IMPORTANT!
   assemble_linearform();
   if (solve_dual_problem)
   {
      L->Vector::operator=(0.0); // VERY IMPORTANT!
      assemble_dual_rhs();
   }
}

void FEM::assemble_linearform()
{
   int n_current_tets = motional_current->n_current_tets;
   double omega = motional_current->omega;
   for (int tet_id=0;  tet_id<pmesh->GetNE();  tet_id++) {
       int global_elem_id = pmesh->GetAttribute(tet_id);
       auto it = motional_current->current_tet_id_to_index.find(global_elem_id);
       int idx;
       if (it != motional_current->current_tet_id_to_index.end()) {
            idx = it->second;
       }
       else {
           continue;
       }
//      int idx;
//      std::vector<int>::iterator it;
//      it = std::find(motional_current->current_tets_id.begin(),motional_current->current_tets_id.end(),global_elem_id);
//      if (it != motional_current->current_tets_id.end()){
//         idx = it - motional_current->current_tets_id.begin();
//        }
//    else continue;

      double J_re[3], J_im[3];
      J_re[0] = motional_current->Jx_re[idx];
      J_re[1] = motional_current->Jy_re[idx];
      J_re[2] = motional_current->Jz_re[idx];
      J_im[0] = motional_current->Jx_im[idx];
      J_im[1] = motional_current->Jy_im[idx];
      J_im[2] = motional_current->Jz_im[idx];

      // Do volume integration
      const FiniteElement* fe = pfes->GetFE(tet_id);
      ElementTransformation* eltrans = pfes->GetElementTransformation(tet_id);
      const IntegrationRule* ir = &IntRules.Get(fe->GetGeomType(), 3*fe->GetOrder() + eltrans->OrderW());
//      double elem_volume = pmesh->GetElementVolume(tet_id);
      // Retrieve the number of basis functions
      const int ndof = fe->GetDof();
      const int dim = fe->GetDim();
      // Allocate a vector to hold the values of each basis function
      DenseMatrix shape(ndof, dim);
      Vector elemvect_re(ndof), elemvect_im(ndof);
      elemvect_re = 0.0;
      elemvect_im = 0.0;
      // Loop over each quadrature point in the reference element
      for (int q=0; q<ir->GetNPoints(); ++q)
      {
         // Extract the current quadrature point from the integration rule
         const IntegrationPoint& ip = ir->IntPoint(q);
         // Prepare to evaluate the coordinate transformation at the current
         // quadrature point
         eltrans->SetIntPoint(&ip);
         // Compute the Jacobian determinant at the current integration point
         double WJ = eltrans->Weight()*ip.weight;
         // Evaluate the basis functions at the point ip
         fe->CalcPhysVShape(*eltrans, shape);
         for (int l=0; l<ndof; l++) // loop all dofs
         {
            // // e^iwt
            elemvect_re[l] += WJ*omega*EM::mu0*(shape(l,0)*J_im[0]+shape(l,1)*J_im[1]+shape(l,2)*J_im[2]);
            elemvect_im[l] += -WJ*omega*EM::mu0*(shape(l,0)*J_re[0]+shape(l,1)*J_re[1]+shape(l,2)*J_re[2]);

            // e^-iwt
            // elemvect_re[l] += -WJ*omega*EM::mu0*(shape(l,0)*J_im[0]+shape(l,1)*J_im[1]+shape(l,2)*J_im[2]);
            // elemvect_im[l] += WJ*omega*EM::mu0*(shape(l,0)*J_re[0]+shape(l,1)*J_re[1]+shape(l,2)*J_re[2]);
         }
      }
      Array<int> vdofs;
      pfes->GetElementVDofs(tet_id, vdofs);
      D->real().AddElementVector(vdofs, elemvect_re);
      D->imag().AddElementVector(vdofs, elemvect_im);
   }
}



//-------------- add by lrj -----------------
void FEM::assemble_dual_rhs()
{
   for (int i=0; i<local_sites_tets.Size(); i++)
   {
      int tet_id = local_sites_tets[i];
      const FiniteElement* fe = pfes->GetFE(tet_id);
      const IntegrationRule* ir = &IntRules.Get(fe->GetGeomType(), 2*fe->GetOrder());
      ElementTransformation* eltrans = pfes->GetElementTransformation(tet_id);

      double elem_volume = pmesh->GetElementVolume(tet_id);
      double inv_volume = 1.0/elem_volume;
      const int ndof = fe->GetDof();
      const int dim = fe->GetDim();
      DenseMatrix shape(ndof, dim);
      Array<int> vdofs;
      pfes->GetElementVDofs(tet_id, vdofs);
      Vector elemvect(ndof);
      elemvect = 0.0;
      for (int q=0; q<ir->GetNPoints(); ++q) 
      {
         const IntegrationPoint& ip = ir->IntPoint(q);
         eltrans->SetIntPoint(&ip);
         fe->CalcPhysVShape(*eltrans, shape);
         double WJ = eltrans->Weight()*ip.weight;
         for (int l=0; l<ndof; l++) // loop all dofs
         {
            elemvect[l] += WJ*inv_volume*(shape(l,0)*1.0+shape(l,1)*1.0+shape(l,2)*1.0); 
         }
      }
      L->real().AddElementVector(vdofs, elemvect);
   }
}

//-------------- add by lrj -----------------
void FEM::find_sites_tets()
{
   local_sites_tets.DeleteAll();
   //local_sites_id.DeleteAll();
   local_sites_r.DeleteAll();
   local_sites_theta.DeleteAll();
   local_sites_phi.DeleteAll();

   //------------- add by lrj ---------------
   // Load global sites coordinates
   std::ifstream sites_stream(param_handler->sites_file);
   assert(sites_stream.good());
   int n_sites;
   sites_stream >> n_sites;
   Array<double> global_sites_r(n_sites), global_sites_theta(n_sites), global_sites_phi(n_sites);
   Array<double> global_sites_x(n_sites), global_sites_y(n_sites), global_sites_z(n_sites);
   for (int i=0; i<n_sites; i++) 
   {
      double r, theta, phi;
      sites_stream >> r >> theta >> phi;
      assert(r>=0);                    // r:       [0, +inf]
      assert(theta>=0 && theta<=180);  // theta:   [0, 180] => [0, pi]
      assert(phi>=0 && phi<=360);      // phi      [0, 360] => [0, 2pi]
      global_sites_r[i] = r;
      global_sites_theta[i] = theta;
      global_sites_phi[i] = phi;
      theta = theta/180.0*EM::pi;
      phi = phi/180.0*EM::pi;
      global_sites_x[i] = r*sin(theta)*cos(phi);
      global_sites_y[i] = r*sin(theta)*sin(phi);
      global_sites_z[i] = r*cos(theta);
   }
   // Find local sites and tets which contain the local sites
   DenseMatrix point_mat(3,n_sites);
   Array<int> sites_tets;
   Array<IntegrationPoint> ips;
   for (int i=0; i<n_sites; i++) 
   {
      point_mat(0,i) = global_sites_x[i];
      point_mat(1,i) = global_sites_y[i];
      point_mat(2,i) = global_sites_z[i];
   }
   pmesh->FindPoints(point_mat, sites_tets, ips);
   for (int i=0; i<sites_tets.Size(); i++) 
   {
      int tet_id = sites_tets[i];
      if (tet_id == -1 && myid==0)
         throw std::invalid_argument("FEM::find_sites_tets(): some sites were not found! Please check your mesh.\n");
      if (tet_id != -2) // -2 means find at other processor
      {
         local_sites_tets.Append(tet_id);
         //local_sites_id.Append(global_sitpost->es_id[i]);
         local_sites_r.Append(global_sites_r[i]);
         local_sites_theta.Append(global_sites_theta[i]);
         local_sites_phi.Append(global_sites_phi[i]);
      }
   }
}

//-------------- add by lrj -----------------
void FEM::find_sites_tets_with_gslib()
{
   #ifdef USE_GSLIB
   local_sites_tets.DeleteAll();
   //local_sites_id.DeleteAll();
   local_sites_r.DeleteAll();
   local_sites_theta.DeleteAll();
   local_sites_phi.DeleteAll();

   // Set sites coordinates
   //int n_sites = param_handler->n_sites;

   std::ifstream sites_stream(param_handler->sites_file);
   assert(sites_stream.good());
   int n_sites;
   sites_stream >> n_sites;
   std::vector<double> global_sites_r(n_sites), global_sites_theta(n_sites), global_sites_phi(n_sites);
   std::vector<double> global_sites_x(n_sites), global_sites_y(n_sites), global_sites_z(n_sites);
//   Array<double> global_sites_r(n_sites), global_sites_theta(n_sites), global_sites_phi(n_sites);
//   Array<double> global_sites_x(n_sites), global_sites_y(n_sites), global_sites_z(n_sites);
   for (int i=0; i<n_sites; i++) 
   {
      double r, theta, phi;
      sites_stream >> r >> theta >> phi;
      assert(r>=0);                    // r:       [0, +inf]
      assert(theta>=0 && theta<=180);  // theta:   [0, 180] => [0, pi]
      assert(phi>=0 && phi<=360);      // phi      [0, 360] => [0, 2pi]
      global_sites_r[i] = r;
      global_sites_theta[i] = theta;
      global_sites_phi[i] = phi;
      theta = theta/180.0*EM::pi;
      phi = phi/180.0*EM::pi;
      global_sites_x[i] = r*sin(theta)*cos(phi);
      global_sites_y[i] = r*sin(theta)*sin(phi);
      global_sites_z[i] = r*cos(theta);
   }


   Vector vxyz(n_sites*3);
   vxyz.SetSize(3 * n_sites);

   for (int i=0; i<n_sites; i++) 
   {
      vxyz(i) = global_sites_x[i];
      vxyz(n_sites+i) = global_sites_y[i];
      vxyz(2*n_sites+i) = global_sites_z[i];
   }

   // Find tets that contain xyz coordinates
   FindPointsGSLIB finder(pmesh->GetComm());
   finder.Setup(*pmesh);
   finder.FindPoints(vxyz);
   Array<unsigned int> code = finder.GetCode();
   Array<unsigned int> rank = finder.GetProc();
   Array<unsigned int> elem = finder.GetElem();

   // Check if there are some points not found
   int flag = code.Find(2); // code=2 means the point is not found
   if (flag != -1) // != -1 means code=2 is found, indicating that some points were not found
   {
      throw std::invalid_argument("FEM::find_sites_tets_with_gslib()：some sites were not found! Please check your mesh.\n");
   }

   for (int i=0; i<n_sites; i++) 
   {
      if (rank[i]==myid)
      {
         local_sites_tets.Append(elem[i]);
         //local_sites_id.Append(world_sites_id[i]);
         local_sites_r.Append(global_sites_r[i]);
         local_sites_theta.Append(global_sites_theta[i]);
         local_sites_phi.Append(global_sites_phi[i]);
      }
   }

   finder.FreeData();
   #endif
}


















void FEM::form_linear_system()
{
   assemble_system_matrix();
   assemble_rhs();
   
   double tic, toc;
   if (myid == 0) 
   {
      std::cout << "\nStep 1: Forming linear system...\n"; 
   }
   // Bilinear form B(U,V)
   // tic = MPI_Wtime();
   // B->Assemble();
   // B->Finalize();
   // toc = MPI_Wtime();
   // if (myid == 0)
   // {
   //    print_time(toc-tic, "System matrix assembly: ");
   // }
   // Assemble linear form D(V)
   // assemble_linearform();
   // Impose dirichlet boundary condition
   Vector zero_vec(pmesh->SpaceDimension()); 
   zero_vec = 0.0;
   VectorConstantCoefficient zero_vec_coef(zero_vec);
   Ugf->ProjectBdrCoefficientTangent(zero_vec_coef, zero_vec_coef, dirichlet_bdr);

   // Linear system of primal problem KU=F
   tic = MPI_Wtime();
   B->FormLinearSystem(dirichlet_tdof_list, *Ugf, *D, K, U, F);

   //------------ add by lrj -----------------
   if (solve_dual_problem) Wgf->ProjectBdrCoefficientTangent(zero_vec_coef, zero_vec_coef, dirichlet_bdr);
   if (solve_dual_problem)
   {
      B->FormLinearSystem(dirichlet_tdof_list, *Wgf, *L, K, W, C);
   }
   //------------ add by lrj -----------------

   toc = MPI_Wtime();
   if (myid == 0) 
   {
      print_time(toc-tic, "FormLinearSystem: ");
   }
}

void FEM::form_preconditioner(SolverParm& solver_parm)
{
   double tic, toc;
   if (myid == 0) 
   {
      std::cout << "\nStep 2: Setup auxiliary-space block-diagonal preconditioner...\n"; 
      if (solver_parm.prec_type == 0)
         std::cout << "PrecondType: AMS\n";
      else if (solver_parm.prec_type == 1)
         std::cout << "PrecondType: PCG-AMS\n";
      else
         std::cout << "PrecondType: Multigrid\n";

   }
   
   // Assemble bilinearform of BlockDiagonalPreconditioner
   tic = MPI_Wtime();
   precond->Assemble();
   precond->Finalize();
   toc = MPI_Wtime();
   if (myid == 0) 
   {
      print_time(toc-tic, "Preconditioning matrix assembly: ");
   }

   tic = MPI_Wtime();
   precond->FormSystemMatrix(dirichlet_tdof_list, A);
   toc = MPI_Wtime();
   if (myid == 0) 
   {
      print_time(toc-tic, "FormSystemMatrix: ");
   }
}

void FEM::solve()
{
   // // // ######### Solve using a direct #######
   // {
   //    HypreParMatrix *A1 = K.As<ComplexHypreParMatrix>()->GetSystemMatrix();
   //    MUMPSSolver mumps;
   //    mumps.SetPrintLevel(0);
   //    mumps.SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_INDEFINITE );
   //    mumps.SetOperator(*A1);
   //    mumps.Mult(F,U);
   //    delete A1;
   // }

   // //######### Solve using an iterative solver #######
   {
      SolverParm solver_parm(param_handler->linear_opts_file);
      this->form_preconditioner(solver_parm);
      // Construct BlockDiagonalPreconditioner
      Array<int> blockTrueOffsets;
      blockTrueOffsets.SetSize(3);
      blockTrueOffsets[0] = 0;
      blockTrueOffsets[1] = A.Height();
      blockTrueOffsets[2] = A.Height();
      blockTrueOffsets.PartialSum();
      BlockDiagonalPreconditioner block_prec(blockTrueOffsets);

      if (solver_parm.prec_type == 0) // Plain AMS
      {
         HypreAMS *pc_r = new HypreAMS(A, pfes);
         solver_parm.set_ams_parm(*pc_r);
         HypreAMS *pc_i = new HypreAMS(A, pfes);
         solver_parm.set_ams_parm(*pc_i);
         block_prec.SetDiagonalBlock(0, pc_r);
         block_prec.SetDiagonalBlock(1, pc_i);
      }
      else if (solver_parm.prec_type == 1) // PCG-AMS
      {
         HypreAMS *ams = new HypreAMS(A, pfes);
         HyprePCG *pc_r = new HyprePCG(A);
         pc_r->SetPreconditioner(*ams);
         solver_parm.set_pcg_parm(*pc_r);
         HyprePCG *pc_i = new HyprePCG(A);
         pc_i->SetPreconditioner(*ams);
         solver_parm.set_pcg_parm(*pc_i);
         block_prec.SetDiagonalBlock(0, pc_r);
         block_prec.SetDiagonalBlock(1, pc_i);
      }
      else if (solver_parm.prec_type == 2) // Multigrid
      {
         Operator *pc_r = nullptr;
         Operator *pc_i = nullptr;
         OperatorHandle PCOp;
         precond->FormSystemMatrix(dirichlet_tdof_list, PCOp);
         pc_r = new HypreAMS(*PCOp.As<HypreParMatrix>(), pfes);
         pc_i = new ScaledOperator(pc_r,
                                   (ComplexOperator::BLOCK_SYMMETRIC == ComplexOperator::HERMITIAN) ? -1.0 : 1.0);
         block_prec.SetDiagonalBlock(0, pc_r);
         block_prec.SetDiagonalBlock(1, pc_i);
      }
      block_prec.owns_blocks = 1;

      // // Setup GMRES solver
      // GMRESSolver gmres(MPI_COMM_WORLD);
      // gmres.SetPrintLevel(0);
      // gmres.SetKDim(50);
      // gmres.SetMaxIter(2000);
      // gmres.SetRelTol(1e-5);
      // gmres.SetAbsTol(0.0);
      // gmres.SetOperator(*K.Ptr());
      // gmres.SetPreconditioner(block_prec);
      // // Solve linear system
      // double tic, toc;
      // if (myid == 0)
      // {
      //    std::cout << "\nStep 3: Solving linear system...\n";
      // }
      // tic = MPI_Wtime();
      // gmres.Mult(F, U);
      // toc = MPI_Wtime();
      // if (myid == 0)
      // {
      //    print_time(toc-tic);
      // }

      // Setup FGMRES solver
      FGMRES fgmres(MPI_COMM_WORLD);
      fgmres.SetOperator(*K.Ptr());
      fgmres.SetPreconditioner(block_prec);
      fgmres.SetMaxIter(solver_parm.fgmres_maxit);
      fgmres.SetRelTol(solver_parm.fgmres_primal_tol);
      fgmres.SetKDim(solver_parm.fgmres_restart); // restart
      // 1 - print norm of each iterations, 2 - print number of iterations
      fgmres.SetPrintLevel(solver_parm.fgmres_print_level);
      // Solve linear system
      double tic, toc;
      if (myid == 0)
      {
         std::cout << "\nStep 3: Solving linear system...\n";
      }
      tic = MPI_Wtime();
      fgmres.Mult(F, U);
      toc = MPI_Wtime();
      if (myid == 0)
      {
         print_time(toc - tic);
      }
      // Recover primary solution
      B->RecoverFEMSolution(U, *D, *Ugf);

      if (solve_dual_problem)
      {
         fgmres.SetRelTol(solver_parm.fgmres_dual_tol);
         if (myid == 0)
         {
            std::cout << " Solve dual problem...\n";
         }
            
         tic = MPI_Wtime();
         fgmres.Mult(C, W);
         toc = MPI_Wtime();
         if (myid == 0)
         {
            std::cout << " MPI_Wtime: " << toc - tic << " (s)\n";
         }
            
         B->RecoverFEMSolution(W, *L, *Wgf);
      }
   }
}

// ----------- add by lrj ---------------------
void FEM::error_estimating()
{
   Estimator_nJ estimator;
   // estimated_errors = estimator.get_goal_estimated_error(*param, *gem1d, source_origin, *Ugf, *Wgf);

   // Compute goal-oriented estimated error
   estimated_solution_errors = estimator.get_estimated_error(*param_handler, *Ugf);
   int Nt = estimated_solution_errors.Size();
   estimated_errors.SetSize(Nt);

   if (param_handler->amr_method==0) // goal-oriented AMR
   {
      Vector estimated_influence_errors = estimator.get_estimated_error(*param_handler, *Wgf);
      for (int t=0; t<Nt; t++)
      {
         estimated_errors[t] = estimated_solution_errors[t]*estimated_influence_errors[t];
      }
   }
   else // non-goal-oriented and global AMR
   {
      for (int t=0; t<Nt; t++)
      {
         estimated_errors[t] = estimated_solution_errors[t];
      }
   }
}

//------------ add by lrj ---------------------
void FEM::update()
{
   pfes->Update();

   Ugf->Update(); *Ugf = 0.0;
   Wgf->Update(); *Wgf = 0.0;

   B->Update();
   D->Update();
   L->Update();
   precond->Update();
}

//------------ add by lrj ---------------------
void FEM::forward_modeling(int source_id)
{
   // For printing cumulative modeling time in adaptive FEM
   double modeling_tic, modeling_toc;
   modeling_tic = MPI_Wtime();

   int iter = 0;
   double tic, toc;
   bool reach_max_dofs = false;

   // else if (param_handler->amr_method==1) // non-goal-oriented AMR
   // {
   //    print_info = "Non-goal-oriented hp-adaptive refinement loop #";
   // }
   // else if (param_handler->amr_method==2) // global AMR
   // {
   //    print_info = "Global hp-adaptive refinement loop #";
   // }

   // Step 1: h-refinement
   initialize();
   for (iter = 0; iter < param_handler->h_iter; iter++)
   {
      if (myid == 0) // goal-oriented AMR
         std::cout << "-------Goal-oriented hp-adaptive refinement loop #" << iter + 1 << "\n";
      int mesh_squence = pmesh->GetSequence();
      if (myid == 0)
      {
         std::cout << "p_order: 1\n";
      }

      print_problem_size();

      if (myid == 0)
         std::cout << "Find station elements...\n";
      tic = MPI_Wtime();
     #ifdef USE_GSLIB
     find_sites_tets_with_gslib();
     #else
      find_sites_tets();
     #endif
      toc = MPI_Wtime();
      if (myid == 0)
         std::cout << "MPI_Wtime: " << toc - tic << " (s)\n";

      if (myid == 0)
         std::cout << "Form linear system and preconditioner...\n";
      tic = MPI_Wtime();
      form_linear_system();
      toc = MPI_Wtime();
      if (myid == 0)
         std::cout << "MPI_Wtime: " << toc - tic << " (s)\n";

      solve();

      if (myid == 0)
         std::cout << "\n###############Post processing...###############\n";
      double post_tic = MPI_Wtime();
      int p_order = param_handler->p_order;
      std::ostringstream sol_name;
      sol_name << "Forward_solutions/" << "MT" << "_source" << source_id << "_iter" << iter
               << ".sol";
      std::ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);

      Post post(*param_handler, *pmesh, *pfes, Ugf->real(), Ugf->imag(), motional_current->omega);
      post.post_process();
      post.save_as_one(sol_ofs);

      // if (myid==0) *log_stream << " Post-processing...\n";
      // tic = MPI_Wtime();
      // std::ostringstream sol_name;
      // sol_name << "Forward_solutions/period" << period_id << "_amr" << mesh_squence << "_p1.sol";
      // std::ofstream sol_ofs(sol_name.str().c_str());
      // Post post(*param, *pfes, Ugf->real(), Ugf->imag(), local_sites_tets,
      //          local_sites_id, local_sites_r, local_sites_theta, local_sites_phi);
      // post.post_process(*gem1d, source_origin);
      // post.save_results(period_id, sol_ofs);
      // toc = MPI_Wtime();
      // if (myid==0) *log_stream << " MPI_Wtime: " << toc-tic << " (s)\n";

      // Goal-oriented error estimation
      error_estimating();
// Save vtk mesh for showing
      if (param_handler->print_vtk)
      {
         std::ostringstream vtk_name;
         vtk_name << "Forward_solutions/period" << "_amr" << mesh_squence << "_p1.vtk";
         std::ofstream vtk_ofs(vtk_name.str().c_str());
         print_vtk_as_one(vtk_ofs);
      }


      if (param_handler->amr_method == 0 || param_handler->amr_method == 1) // goal-oriented or non goal-oriented AMR
      {
         // Obtain global maximum estimated error
         double local_max_err = estimated_errors.Max();
         double global_max_err;
         MPI_Allreduce(&local_max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, pmesh->GetComm());

         // Refine the mesh based on error estimator
         double threshold = param_handler->beta * global_max_err;
         pmesh->RefineByError(estimated_errors, threshold);
      }

      update();

      if (get_problem_dofs() > param_handler->max_dofs)
      {
         reach_max_dofs = true;
         if (myid == 0)
         {
            std::cout << "\n Stop due to reached max number of dofs\n";
            std::cout << " Dofs: " << get_problem_dofs() << "\n";
            std::cout << " Dofs: " << param_handler->max_dofs << "\n";
            modeling_toc = MPI_Wtime();
            std::cout << " Refinement cumulative MPI_Wtime: " << modeling_toc - modeling_tic << " (s)\n";
         }
         break;
      }

      modeling_toc = MPI_Wtime();
      if (myid == 0)
         std::cout << "Refinement cumulative MPI_Wtime: " << modeling_toc - modeling_tic << " (s)\n";
      if (myid == 0)
         std::cout << "\n ---------------------Next loop--------------------\n";
   }
}

void FEM::print_vtk_as_one(std::ofstream& out)
{
   // Note: works for tetrahedral meshes
   if (myid==0)
   {
      out << "# vtk DataFile Version 3.0\n"
             "Generated by Hongbo Yao\n"
             "ASCII\n"
             "DATASET UNSTRUCTURED_GRID\n";
   }

   MPI_Status status;
   MPI_Comm comm = pmesh->GetComm();
   int comm_size;
   MPI_Comm_size(comm, &comm_size);
   int space_dim = pmesh->SpaceDimension();

   // Points
   int nv = pmesh->GetNV();
   int nv_glob;
   MPI_Reduce(&nv, &nv_glob, 1, MPI_INT, MPI_SUM, 0, comm);
   double* vert; // nv*space_dim
   if (myid==0)
   {
      out << "POINTS " << nv_glob << " double\n";
      for (int i=0; i<pmesh->GetNV(); i++)
      {
         double* point = pmesh->GetVertex(i);
         for (int j=0; j<space_dim; j++)
         {
            out << point[j] << " ";
         }
         out << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&nv, 1, MPI_INT, p, 101, comm, &status);
         vert = new double[nv*space_dim];
         MPI_Recv(vert, nv*space_dim, MPI_DOUBLE, p, 102, comm, &status);
         for (int i=0; i<nv; i++)
         {
            for (int j=0; j<space_dim; j++)
            {
               out << vert[i*space_dim+j] << " ";
            }
            out << "\n";
         }
         delete vert;
      }
   }
   else
   {
      Vector vert_coord(nv*space_dim);
      int k = 0;
      for (int i=0; i<nv; i++)
      {
         double* elem_vert = pmesh->GetVertex(i);
         for (int j=0; j<space_dim; j++)
         {
            vert_coord[k] = elem_vert[j];
            k++;
         }
      }
      vert = vert_coord.GetData();
      MPI_Send(&nv, 1, MPI_INT, 0, 101, comm);
      MPI_Send(vert, nv*space_dim, MPI_DOUBLE, 0, 102, comm);
   }
   
   // Cells
   int ne = pmesh->GetNE();
   int ne_glob;
   MPI_Reduce(&ne, &ne_glob, 1, MPI_INT, MPI_SUM, 0, comm);
   const int nv_per_element = pmesh->GetElement(0)->GetNVertices();
   int* vert_id; // ne*nv_per_element
   if (myid==0)
   {
      out << "CELLS " << ne_glob << " " << ne_glob*(nv_per_element+1) << "\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         const int *v = pmesh->GetElement(i)->GetVertices();
         out << nv_per_element;
         for (int j=0; j<nv_per_element; j++)
         {
            out << " " << v[j];
         }
         out << "\n";
      }
      int nv_sum = pmesh->GetNV();
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 103, comm, &status);
         vert_id = new int[ne*nv_per_element];
         MPI_Recv(vert_id, ne*nv_per_element, MPI_INT, p, 104, comm, &status);
         MPI_Recv(&nv, 1, MPI_INT, p, 105, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << nv_per_element;
            for (int j=0; j<nv_per_element; j++)
            {
               out << " " << nv_sum+vert_id[i*nv_per_element+j];
            }
            out << "\n";
         }
         nv_sum += nv;
         delete vert_id;
      }
   }
   else
   {
      Array<int> vert_id_;
      for (int i=0; i<ne; i++)
      {
         Array<int> elem_vert_id;
         pmesh->GetElementVertices(i,elem_vert_id);
         vert_id_.Append(elem_vert_id);
      }
      vert_id = vert_id_.GetData();
      MPI_Send(&ne, 1, MPI_INT, 0, 103, comm);
      MPI_Send(vert_id, ne*nv_per_element, MPI_INT, 0, 104, comm);
      MPI_Send(&nv, 1, MPI_INT, 0, 105, comm);
   }

   // Cell types
   if (myid==0)
   {
      out << "CELL_TYPES " << ne_glob << "\n";
      for (int i=0; i<ne_glob; i++)
      {
         int vtk_cell_type = 10; // tet
         out << vtk_cell_type << "\n";
      }
   }

   // Element error indicator
   double local_max_err = estimated_solution_errors.Max();
   double global_max_err;
   MPI_Allreduce(&local_max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, pmesh->GetComm());
   double* element_error_indicator;
   if (myid==0)
   {
      out << "CELL_DATA " << ne_glob << "\n";
      out << "SCALARS element_relative_error_indicator double 1\n"
         << "LOOKUP_TABLE default\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         out << estimated_solution_errors[i]/global_max_err << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 106, comm, &status);
         element_error_indicator = new double[ne];
         MPI_Recv(element_error_indicator, ne, MPI_DOUBLE, p, 107, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << element_error_indicator[i]/global_max_err << "\n";
         }
      }
      delete element_error_indicator;
   }
   else
   {
      MPI_Send(&ne, 1, MPI_INT, 0, 106, comm);
      element_error_indicator = estimated_solution_errors.GetData();
      MPI_Send(element_error_indicator, ne, MPI_DOUBLE, 0, 107, comm);
   }

   // Element marker
   int* element_marker;
   if (myid==0)
   {
      out << "SCALARS element_marker int 1\n"
          << "LOOKUP_TABLE default\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         out << pmesh->GetAttribute(i) << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 108, comm, &status);
         element_marker = new int[ne];
         MPI_Recv(element_marker, ne, MPI_INT, p, 109, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << element_marker[i] << "\n";
         }
      }
      delete element_marker;
   }
   else
   {
      Array<int> element_marker_(ne);
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         element_marker_[i] = pmesh->GetAttribute(i);
      }
      element_marker = element_marker_.GetData();
      MPI_Send(&ne, 1, MPI_INT, 0, 108, comm);
      MPI_Send(element_marker, ne, MPI_INT, 0, 109, comm);
   }

   // Element process id
   int element_process_id;
   if (myid==0)
   {
      out << "SCALARS element_process_id int 1\n"
          << "LOOKUP_TABLE default\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         out << pmesh->GetMyRank() << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 110, comm, &status);
         MPI_Recv(&element_process_id, 1, MPI_INT, p, 111, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << element_process_id << "\n";
         }
      }
   }
   else
   {
      element_process_id = pmesh->GetMyRank();
      MPI_Send(&ne, 1, MPI_INT, 0, 110, comm);
      MPI_Send(&element_process_id, 1, MPI_INT, 0, 111, comm);
   }

   // Element volume
   double* element_volume;
   if (myid==0)
   {
      out << "SCALARS element_volume double 1\n"
          << "LOOKUP_TABLE default\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         out << pmesh->GetElementVolume(i) << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 112, comm, &status);
         element_volume = new double[ne];
         MPI_Recv(element_volume, ne, MPI_DOUBLE, p, 113, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << element_volume[i] << "\n";
         }
      }
      delete element_volume;
   }
   else
   {
      Array<double> element_volume_(ne);
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         element_volume_[i] = pmesh->GetElementVolume(i);
      }
      element_volume = element_volume_.GetData();
      MPI_Send(&ne, 1, MPI_INT, 0, 112, comm);
      MPI_Send(element_volume, ne, MPI_DOUBLE, 0, 113, comm);
   }

    // Element conductivity tensor - 修改为输出张量形式
   double* element_conductivity_tensor;
   if (myid==0)
   {
      out << "TENSORS element_conductivity_tensor double\n";
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         DenseMatrix sigma_tensor(3, 3);
         param_handler->get_elem_conductivity_tensor(pmesh->GetAttribute(i), sigma_tensor);
         out << sigma_tensor(0,0) << " " << sigma_tensor(0,1) << " " << sigma_tensor(0,2) << "\n";
         out << sigma_tensor(1,0) << " " << sigma_tensor(1,1) << " " << sigma_tensor(1,2) << "\n";
         out << sigma_tensor(2,0) << " " << sigma_tensor(2,1) << " " << sigma_tensor(2,2) << "\n";
      }
      for (int p=1; p<comm_size; p++)
      {
         MPI_Recv(&ne, 1, MPI_INT, p, 114, comm, &status);
         element_conductivity_tensor = new double[ne * 9]; // 每个单元9个分量
         MPI_Recv(element_conductivity_tensor, ne * 9, MPI_DOUBLE, p, 115, comm, &status);
         for (int i=0; i<ne; i++)
         {
            out << element_conductivity_tensor[i*9 + 0] << " " 
                << element_conductivity_tensor[i*9 + 1] << " " 
                << element_conductivity_tensor[i*9 + 2] << "\n";
            out << element_conductivity_tensor[i*9 + 3] << " " 
                << element_conductivity_tensor[i*9 + 4] << " " 
                << element_conductivity_tensor[i*9 + 5] << "\n";
            out << element_conductivity_tensor[i*9 + 6] << " " 
                << element_conductivity_tensor[i*9 + 7] << " " 
                << element_conductivity_tensor[i*9 + 8] << "\n";
         }
         delete [] element_conductivity_tensor;
      }
   }
   else
   {
      Array<double> element_conductivity_tensor_(ne * 9); // 每个单元9个分量
      for (int i=0; i<pmesh->GetNE(); i++)
      {
         DenseMatrix sigma_tensor(3, 3);
         param_handler->get_elem_conductivity_tensor(pmesh->GetAttribute(i), sigma_tensor);
         // 按行优先顺序存储张量
         element_conductivity_tensor_[i*9 + 0] = sigma_tensor(0,0);
         element_conductivity_tensor_[i*9 + 1] = sigma_tensor(0,1);
         element_conductivity_tensor_[i*9 + 2] = sigma_tensor(0,2);
         element_conductivity_tensor_[i*9 + 3] = sigma_tensor(1,0);
         element_conductivity_tensor_[i*9 + 4] = sigma_tensor(1,1);
         element_conductivity_tensor_[i*9 + 5] = sigma_tensor(1,2);
         element_conductivity_tensor_[i*9 + 6] = sigma_tensor(2,0);
         element_conductivity_tensor_[i*9 + 7] = sigma_tensor(2,1);
         element_conductivity_tensor_[i*9 + 8] = sigma_tensor(2,2);
      }
      element_conductivity_tensor = element_conductivity_tensor_.GetData();
      MPI_Send(&ne, 1, MPI_INT, 0, 114, comm);
      MPI_Send(element_conductivity_tensor, ne * 9, MPI_DOUBLE, 0, 115, comm);
   }

   if (myid==0)
   {
      out.close();
   }
}
