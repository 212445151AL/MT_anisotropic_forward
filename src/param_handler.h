// A class for handering input parameters

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#ifndef _PARAM_HANDLER_H
#define _PARAM_HANDLER_H

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "mfem.hpp"

namespace mfem {
class DenseMatrix;
}

class ParamHandler
{
public:
   /* Input parameters in *.config file */
   // Inducing source file
   std::string             source_file1;
   std::string             source_file2;
   std::string             source_file3;
   // order of polynomial of last amr iteration
   int                     p_order;          
   // FGMRES and AMS options
   std::string             linear_opts_file;
   // Conductivity model file
   std::string             conductivity_model_file;
   // Mesh file
   std::string             mesh_file; 
   // Geomagnetic sites file
   std::string             sites_file;
   // Upper and lower boundary marker nxE=0
   int                     bdr_marker;   

   //----------------------add by lrj ---------------------------
   // AMR method: 0-goal-oriented AMR, 1-non-goal-oriented AMR, 2-global mesh refinement, default = 0
   int amr_method;
   // Iterations of adaptive h-refinement (mesh refinement) 
   int h_iter;
   // Max_dofs
   int max_dofs;
   // Marking parameter that controls refinement elements
   double beta;
   // Print VTK mesh and estimated error for each element
   int print_vtk;


   /* Detailed parameters */
   // Model parameters
   std::string             marker_type;   // "region_marker" or "element_id"
   int                     n_regions;     // number of markers (regions or elements)
   std::vector<int>        marker_vec;    // the same as map's key, but more easier to visit
   std::map<int, std::vector<double>>   region_conductivity_tensor;   
   std::map<int, int>   region_tag;//in order to find measure point (Yangcong)

public:
   ParamHandler(char* parameter_file);
   ~ParamHandler();

public:
   // Skip comments, modified from the original line_process() function provided by 
   // Yiyuan Zhong, Central South University.
   void skip_comments(std::istream& in, std::vector<std::string>& para_str_vec, 
                      std::string comment_str="#");
   
   // Read model information from file
   void read_parameter_file(char* parameter_file);

   // Get element tensor_conductivity according to marker (attribute)
   void get_elem_conductivity_tensor(int marker, mfem::DenseMatrix& sigma_tensor);

   // Get element tag according to marker (attribute)
   double get_elem_tag(int marker);

   

   // Output to file
   void output(std::string out="input_parameter_list.log");
};

#endif // _PARAM_HANDLER_H
