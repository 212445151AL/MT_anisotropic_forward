// Motional induction current

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#ifndef _MOTIONAL_CURRENT_H
#define _MOTIONAL_CURRENT_H

#include <string>
#include <vector>
#include <map>

class MotionalCurrent
{
public:
   std::string          current_file;
   std::string          current_name; // e.g. M2_tide
   double               period_in_hours; // in hours
   double               omega; // angular frequency
   int                  n_current_tets;
   std::vector<int>     current_tets_id;
   std::map<int, int>   current_tet_id_to_index;
//   std::set<int>        current_tets_id;
   std::vector<double>  Jx_re, Jx_im;
   std::vector<double>  Jy_re, Jy_im;
   std::vector<double>  Jz_re, Jz_im;

public:
   MotionalCurrent(std::string& current_file_);
   ~MotionalCurrent();

public:
   void load_motional_current();
};

#endif // _MOTIONAL_CURRENT_H
