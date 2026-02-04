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

#include "error_estimators.h"

Estimator_nJ::Estimator_nJ()
{

}



Estimator_nJ::~Estimator_nJ()
{

}



Vector Estimator_nJ::compute_estimated_error(ParamHandler& param_handler, 
                                             ParComplexGridFunction& U)
{
   ParFiniteElementSpace* pfes = U.ParFESpace();
   ParMesh* pmesh = pfes->GetParMesh();
   int dim = pmesh->SpaceDimension();
   int Nt = pmesh->GetNE();
   Vector estimated_error(Nt);
   estimated_error = 0.0;


   // Loop interior faces
   int n_faces = pmesh->GetNFaces();
   for (int f=0; f<n_faces; f++)
   {
      FaceElementTransformations* face_elem_trans = pmesh->GetInteriorFaceTransformations(f);
      if (face_elem_trans == NULL)
      {
         continue;
      }
      DenseMatrix sigma_tensor1(dim, dim);
      DenseMatrix sigma_tensor2(dim, dim);
      // 创建电导率张量
      for (int i = 0; i < 3; i++) {
         for (int j = 0; j < 3; j++) {
            sigma_tensor1(i, j) = 0.0;
            sigma_tensor2(i, j) = 0.0;
         }
      }

      ElementTransformation* face_trans = face_elem_trans->Face;
      const IntegrationRule* face_ir = &IntRules.Get(face_trans->GetGeometryType(), 2*face_trans->Order()-1);
      
      // Tetrahedron 1
      int tet1_id = face_elem_trans->Elem1No;
      const FiniteElement* tet1_fe = pfes->GetFE(tet1_id);
      ElementTransformation* tet1_trans = face_elem_trans->Elem1;
      
      // 修改这里：获取张量电导率
      param_handler.get_elem_conductivity_tensor(tet1_trans->Attribute, sigma_tensor1);
      
      // Tetrahedron 2
      int tet2_id = face_elem_trans->Elem2No;
      const FiniteElement* tet2_fe = pfes->GetFE(tet2_id);
      ElementTransformation* tet2_trans = face_elem_trans->Elem2;
      
      // 修改这里：获取张量电导率
      param_handler.get_elem_conductivity_tensor(tet2_trans->Attribute, sigma_tensor2);

      // Face integration
      for (int q=0; q<face_ir->GetNPoints(); q++) 
      {
         const IntegrationPoint& ip = face_ir->IntPoint(q);
         face_trans->SetIntPoint(&ip);
         
         Vector normalJ(dim);
         CalcOrtho(face_trans->Jacobian(), normalJ);
         Vector normal(dim);
         normal.Set(1.0/normalJ.Norml2(), normalJ);
         double WJ = face_trans->Weight()*ip.weight;

         // Compute secondary fields on face's tet1
         IntegrationPoint tet1_ip;
         face_elem_trans->Loc1.Transform(ip, tet1_ip);
         tet1_trans->SetIntPoint(&tet1_ip);
         double U1_re[3] = {0.}, U1_im[3] = {0.};
         compute_field_on_surface(U, *pfes, *tet1_trans, *tet1_fe, tet1_id, U1_re, U1_im);

         // Compute secondary fields on face's tet2
         IntegrationPoint tet2_ip;
         face_elem_trans->Loc2.Transform(ip, tet2_ip);
         tet2_trans->SetIntPoint(&tet2_ip);
         double U2_re[3] = {0.}, U2_im[3] = {0.};
         compute_field_on_surface(U, *pfes, *tet2_trans, *tet2_fe, tet2_id, U2_re, U2_im);

         // Complex value
         Dcomplex U1[3] = {0.}, U2[3] = {0.};
         for (int c=0; c<3; c++)
         {
            U1[c] = Dcomplex(U1_re[c], U1_im[c]);
            U2[c] = Dcomplex(U2_re[c], U2_im[c]);
         }

         // 修改这里：计算张量形式的电流密度 J = σ · E
         Dcomplex Js1[3] = {0.}, Js2[3] = {0.};
         for (int i=0; i<3; i++)
         {
            for (int j=0; j<3; j++)
            {
                Js1[i] += sigma_tensor1(i, j) * U1[j];
                Js2[i] += sigma_tensor2(i, j) * U2[j];
            }
         }

         // Compute n dot (J1-J2)
         Dcomplex nJ = 0.0;
         for (int c=0; c<3; c++)
         {
            nJ += normal[c]*(Js1[c]-Js2[c]);
         }

         // Integration
         double integral_value = 0.5*WJ*std::pow(std::abs(nJ),2);
         estimated_error[tet1_id] += integral_value;
         estimated_error[tet2_id] += integral_value;
      }
   }

   // 同样修改 MPI-shared faces 部分
   U.real().ExchangeFaceNbrData();
   U.imag().ExchangeFaceNbrData();
   for (int f=0; f<pmesh->GetNSharedFaces(); f++)
   {
      FaceElementTransformations* face_elem_trans = pmesh->GetSharedFaceTransformations(f);
      ElementTransformation* face_trans = face_elem_trans->Face;
      const IntegrationRule* face_ir = &IntRules.Get(face_trans->GetGeometryType(), 2*face_trans->Order()-1);
      // 在MPI共享面循环中也定义张量
      DenseMatrix sigma_tensor1(dim, dim);
      DenseMatrix sigma_tensor2(dim, dim);
      for (int i = 0; i < 3; i++) {
         for (int j = 0; j < 3; j++) {
            sigma_tensor1(i, j) = 0.0;
            sigma_tensor2(i, j) = 0.0;
         }
      }

      // Tetrahedron 1
      int tet1_id = face_elem_trans->Elem1No;
      const FiniteElement* tet1_fe = pfes->GetFE(tet1_id);
      ElementTransformation* tet1_trans = face_elem_trans->Elem1;
      
      // 修改这里
      param_handler.get_elem_conductivity_tensor(tet1_trans->Attribute, sigma_tensor1);

      // Tetrahedron 2
      int tet2_id = face_elem_trans->Elem2No-pmesh->GetNE(); 
      const FiniteElement* tet2_fe = pfes->GetFaceNbrFE(tet2_id);
      ElementTransformation* tet2_trans = face_elem_trans->Elem2;
      
      // 修改这里
      param_handler.get_elem_conductivity_tensor(tet2_trans->Attribute, sigma_tensor2);

      // Face integration
      for (int q=0; q<face_ir->GetNPoints(); ++q) 
      {
         const IntegrationPoint& ip = face_ir->IntPoint(q);
         face_trans->SetIntPoint(&ip);
         
         Vector normalJ(dim);
         CalcOrtho(face_trans->Jacobian(), normalJ);
         Vector normal(dim);
         normal.Set(1.0/normalJ.Norml2(), normalJ);
         double WJ = face_trans->Weight()*ip.weight;

         // Compute fields on face's tet1
         IntegrationPoint tet1_ip;
         face_elem_trans->Loc1.Transform(ip, tet1_ip);
         tet1_trans->SetIntPoint(&tet1_ip);
         double U1_re[3] = {0.}, U1_im[3] = {0.};
         compute_field_on_surface(U, *pfes, *tet1_trans, *tet1_fe, tet1_id, U1_re, U1_im);

         // Compute fields on face's tet2
         IntegrationPoint tet2_ip;
         face_elem_trans->Loc2.Transform(ip, tet2_ip);
         tet2_trans->SetIntPoint(&tet2_ip);
         double U2_re[3] = {0.}, U2_im[3] = {0.};
         compute_field_on_surface(U, *pfes, *tet2_trans, *tet2_fe, tet2_id, U2_re, U2_im, true);

         // Complex value
         Dcomplex U1[3] = {0.}, U2[3] = {0.};
         for (int c=0; c<3; c++)
         {
            U1[c] = Dcomplex(U1_re[c], U1_im[c]);
            U2[c] = Dcomplex(U2_re[c], U2_im[c]);
         }

         // 修改这里：计算张量形式的电流密度
         Dcomplex Js1[3] = {0.}, Js2[3] = {0.};
         for (int i=0; i<3; i++)
         {
            for (int j=0; j<3; j++)
            {
                Js1[i] += sigma_tensor1(i, j) * U1[j];
                Js2[i] += sigma_tensor2(i, j) * U2[j];
            }
         }

         // Compute n dot (J1-J2)
         Dcomplex nJ = 0.0;
         for (int c=0; c<3; c++)
         {
            nJ += normal[c]*(Js1[c]-Js2[c]);
         }

         // Integration
         double integral_value = 0.5*WJ*std::pow(std::abs(nJ),2);
         estimated_error[tet1_id] += integral_value;
      }
   }

   Vector final_estimated_error(Nt);
   final_estimated_error = 0.0;
   for (int t=0; t<Nt; t++)
   {
      final_estimated_error[t] = std::sqrt(estimated_error[t]);
   }

   return final_estimated_error;
}



Vector Estimator_nJ::get_estimated_error(ParamHandler& param_handler,  
                                         ParComplexGridFunction& U)
{
   return compute_estimated_error(param_handler, U);
}



Vector Estimator_nJ::get_goal_estimated_error(ParamHandler& param_handler,  
                                              ParComplexGridFunction& U, 
                                              ParComplexGridFunction& W)
{
   Vector estimated_error_U = compute_estimated_error(param_handler, U);
   Vector estimated_error_W = compute_estimated_error(param_handler, W);
   int Nt = estimated_error_U.Size();
   Vector goal_estimated_error(Nt);
   for (int t=0; t<Nt; t++)
   {
      goal_estimated_error[t] = estimated_error_U[t]*estimated_error_W[t];
   }
   return goal_estimated_error;
}
