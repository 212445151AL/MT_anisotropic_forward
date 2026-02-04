// Motional induction current

// Author:     Liangyu Xie
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#include "motional_current.h"
#include "em.h"
#include <fstream>
#include <cassert>

MotionalCurrent::MotionalCurrent(std::string& current_file_): current_file(current_file_)
{
   load_motional_current();
}



MotionalCurrent::~MotionalCurrent()
{

}



void MotionalCurrent::load_motional_current()
{
   std::ifstream current_stream(current_file);
   assert(current_stream.good());
   current_stream >> current_name >> period_in_hours;
   double T = period_in_hours*3600.0;
   omega = 2*EM::pi/T;
   current_stream >> n_current_tets;
   current_tets_id.resize(n_current_tets);
   Jx_re.resize(n_current_tets);
   Jx_im.resize(n_current_tets);
   Jy_re.resize(n_current_tets);
   Jy_im.resize(n_current_tets);
   Jz_re.resize(n_current_tets);
   Jz_im.resize(n_current_tets);
   for (int i=0; i<n_current_tets; i++)
   {
      current_stream >> current_tets_id[i]
                     >> Jx_re[i] >> Jx_im[i]
                     >> Jy_re[i] >> Jy_im[i]
                     >> Jz_re[i] >> Jz_im[i];
      this->current_tet_id_to_index.insert(std::make_pair(current_tets_id[i], i));
   }
   current_stream.close();
}
