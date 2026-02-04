#ifndef PW_MATRIX_COEFFICIENT_H
#define PW_MATRIX_COEFFICIENT_H

#include "mfem.hpp"
#include <map>

using namespace mfem;

/**
 * @brief 分片常数矩阵系数类，根据网格单元属性返回不同的矩阵系数
 */
class PwMatrixCoefficient : public MatrixCoefficient
{
private:
   std::map<int, MatrixCoefficient*> attr_to_coef;  // 属性到矩阵系数的映射
   DenseMatrix default_mat;                         // 默认矩阵

public:
   /// @brief 构造函数
   /// @param height 矩阵高度
   /// @param width 矩阵宽度
   PwMatrixCoefficient(int height, int width)
      : MatrixCoefficient(height, width), default_mat(height, width)
   {
      // 设置默认矩阵为单位矩阵
      default_mat = 0.0;
      for (int i = 0; i < height && i < width; i++) {
         default_mat(i, i) = 1.0;
      }
   }

   /// @brief 添加属性对应的矩阵系数
   /// @param attr 网格属性
   /// @param coef 矩阵系数
   void AddCoefficient(int attr, MatrixCoefficient* coef)
   {
      attr_to_coef[attr] = coef;
   }

   /// @brief 设置默认矩阵
   /// @param mat 默认矩阵
   void SetDefaultMatrix(const DenseMatrix& mat)
   {
      default_mat = mat;
   }

   /// @brief 重写 Eval 方法
   virtual void Eval(DenseMatrix &M, ElementTransformation &T,
                    const IntegrationPoint &ip)
   {
      int attr = T.Attribute;
      
      // 查找对应的矩阵系数
      auto it = attr_to_coef.find(attr);
      if (it != attr_to_coef.end()) {
         it->second->Eval(M, T, ip);
      } else {
         // 使用默认矩阵
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