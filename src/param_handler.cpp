// A class for managing model input information

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#include "param_handler.h"
#include "em.h"
#include <cassert>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <sstream>

ParamHandler::ParamHandler(char* parameter_file)
{
   this->read_parameter_file(parameter_file);
}



ParamHandler::~ParamHandler()
{

}



void ParamHandler::skip_comments(std::istream& in, std::vector<std::string>& para_str_vec, 
                                std::string comment_str)
{
   std::string line;
   while (std::getline(in, line))
   {
      for (char &c : line) // loop string in C++11 grammar
      {
         if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
         {
            c = ' ';
         }
      }
      line.erase(0, line.find_first_not_of(" ")); // delete the space at the beginning of the line 
      line.erase(line.find_last_not_of(" ") + 1); // delete the space at the ending of the line   
      int n_comment_start = line.find_first_of(comment_str); 
      if (n_comment_start != std::string::npos)
      {
         line.erase(n_comment_start); // delete comments
      }
      if (line.empty()) continue;
      para_str_vec.push_back(line);
   }
}



void ParamHandler::read_parameter_file(char* parameter_file)
{
   std::ifstream in_stream(parameter_file);
   assert(in_stream.good());
   
   int para_lines = 14; // number of parameters' lines
   std::vector<std::string> para_str_vec;
   skip_comments(in_stream, para_str_vec);
   assert(para_str_vec.size()==para_lines);
   in_stream.close();
   std::stringstream ss;
   
   // Read parameters
   ss << para_str_vec[0]; ss >> source_file1; ss.clear();
   ss << para_str_vec[1]; ss >> source_file2; ss.clear();
   ss << para_str_vec[2]; ss >> source_file3; ss.clear();
   ss << para_str_vec[3]; ss >> p_order; ss.clear();
   ss << para_str_vec[4]; ss >> linear_opts_file; ss.clear();
   ss << para_str_vec[5]; ss >> conductivity_model_file; ss.clear();
   ss << para_str_vec[6]; ss >> mesh_file; ss.clear();
   ss << para_str_vec[7]; ss >> sites_file; ss.clear();
   ss << para_str_vec[8]; ss >> bdr_marker; ss.clear();
   ss << para_str_vec[9]; ss >> amr_method; ss.clear();
   ss << para_str_vec[10]; ss >> h_iter; ss.clear();
   ss << para_str_vec[11]; ss >> max_dofs; ss.clear();
   ss << para_str_vec[12]; ss >> beta; ss.clear();
   ss << para_str_vec[13]; ss >> print_vtk; ss.clear();

   // Output to file for checking
   output();

   
   // Read conductivity model
   std::ifstream cond_in_stream(conductivity_model_file.c_str());
   assert(cond_in_stream.good());
   cond_in_stream >> marker_type;
   std::transform(marker_type.begin(), marker_type.end(), marker_type.begin(), tolower);
   if (marker_type!="region_marker" && marker_type!="element_id")
   {
      std::cout << "parahandler.cpp function ParamHandler::read_parameter_file(): wrong marker_type!\n";
      std::abort();
   }
   cond_in_stream >> n_regions;
   for (int i=0; i<n_regions; i++)
      {
         int marker;
         std::vector<double> cond_tensor(9);
         double tag;
         cond_in_stream >> marker ;
         for (int i=0;i<9;i++)
         {
            cond_in_stream >> cond_tensor[i] ;
         }
         cond_in_stream >> tag;
         marker_vec.push_back(marker);
         region_conductivity_tensor[marker] = cond_tensor;
         region_tag[marker] = tag;
         }
   cond_in_stream.close();

   // Check file stream
   std::ifstream source_in_stream(source_file1.c_str());
   assert(source_in_stream.good());
   source_in_stream.close();
   

   std::ifstream msh_stream(mesh_file);
   assert(msh_stream.good());
   msh_stream.close();

   std::ifstream linear_opts_stream(linear_opts_file);
   assert(linear_opts_stream.good());
   linear_opts_stream.close();
   
   std::ifstream sites_stream(sites_file);
   assert(sites_stream.good());
   sites_stream.close();
}


void  ParamHandler::get_elem_conductivity_tensor(int marker, mfem::DenseMatrix& sigma_tensor)
{  
   assert(!region_conductivity_tensor.empty());
   std::map<int, std::vector<double>>::iterator it = region_conductivity_tensor.find(marker);
   assert(it!=region_conductivity_tensor.end());
   std::vector<double>& vec = it->second;
   sigma_tensor(0,0) = vec[0];
   sigma_tensor(0,1) = vec[1];
   sigma_tensor(0,2) = vec[2];
   sigma_tensor(1,0) = vec[3];
   sigma_tensor(1,1) = vec[4];
   sigma_tensor(1,2) = vec[5];
   sigma_tensor(2,0) = vec[6];
   sigma_tensor(2,1) = vec[7];
   sigma_tensor(2,2) = vec[8];
}

double ParamHandler::get_elem_tag(int marker)
{
   assert(!region_tag.empty());
   std::map<int, int>::iterator it = region_tag.find(marker);
   assert(it!=region_tag.end());
   return (*it).second;
}



void ParamHandler::output(std::string out)
{
   std::ofstream out_stream(out);
   out_stream << source_file1 << "\n";
   out_stream << source_file2 << "\n";
   out_stream << source_file3 << "\n";
   out_stream << p_order << "\n";
   out_stream << linear_opts_file << "\n";
   out_stream << conductivity_model_file << "\n";
   out_stream << mesh_file << "\n";
   out_stream << sites_file << "\n";
   out_stream << bdr_marker << "\n";
   out_stream << amr_method << "\n";
   out_stream << h_iter << "\n";
   out_stream << max_dofs << "\n";
   out_stream << beta << "\n";
   out_stream << print_vtk << "\n";
}



