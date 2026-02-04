// Piecewise Constant Matrix Coefficient Class, returns different matrix coefficients based on grid cell properties

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#ifndef PW_MATRIX_COEFFICIENT_H
#define PW_MATRIX_COEFFICIENT_H

#include "mfem.hpp"
#include <map>

using namespace mfem;

/**
 * @brief 
 */
class PwMatrixCoefficient : public MatrixCoefficient
{
private:
   std::map<int, MatrixCoefficient*> attr_to_coef;  // Mapping of Attributes to Matrix Coefficients
   DenseMatrix default_mat;                       

public:
   /// @brief 
   /// @param height 
   /// @param width 
   PwMatrixCoefficient(int height, int width)
      : MatrixCoefficient(height, width), default_mat(height, width)
   {
      // Set the default matrix to the identity matrix
      default_mat = 0.0;
      for (int i = 0; i < height && i < width; i++) {
         default_mat(i, i) = 1.0;
      }
   }

   /// @brief 
   /// @param attr 
   /// @param coef 
   void AddCoefficient(int attr, MatrixCoefficient* coef)
   {
      attr_to_coef[attr] = coef;
   }

   /// @brief 
   /// @param mat 
   void SetDefaultMatrix(const DenseMatrix& mat)
   {
      default_mat = mat;
   }

   /// @brief
   virtual void Eval(DenseMatrix &M, ElementTransformation &T,
                    const IntegrationPoint &ip)
   {
      int attr = T.Attribute;
      
      auto it = attr_to_coef.find(attr);
      if (it != attr_to_coef.end()) {
         it->second->Eval(M, T, ip);
      } else {
         M = default_mat;
      }
   }

   virtual ~PwMatrixCoefficient()
   {
      for (auto& pair : attr_to_coef) {
         delete pair.second;
      }
   }
};

#endif // PW_MATRIX_COEFFICIENT_H
