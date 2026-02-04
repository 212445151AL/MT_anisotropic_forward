// Parallel post-processing. Compute magnetic induction B from 
// Nedelec electric field E

// Author:     Liangyu Xie,Renjie Li,Chaojian Chen
// Institute:  Central South University (CSU)            
// Email:      8211221219@csu.edu.cn
// Date:       2026/02/04

// GitHub Page: https://github.com/212445151AL

#include "post.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>

Post::Post(ParamHandler&            param_handler_,
           ParMesh&                 pmesh_,
           ParFiniteElementSpace&   pfes_,
           ParGridFunction&         U_re_,
           ParGridFunction&         U_im_,
           double                   omega_):
           param_handler            (&param_handler_),
           pmesh                    (&pmesh_),
           pfes                     (&pfes_),
           U_re                     (U_re_),
           U_im                     (U_im_),
           omega                    (omega_)
{

}

Post::~Post()
{

}

SpaceRelation Post::TestPointInTet(Vector po, int &tet_id)  
{ 
   //assert(point.Size()==3);
	// int global_elem_id = pmesh->GetNode(tet_id)

   Array<int> Vertices;
   pmesh->GetElementVertices(tet_id, Vertices);

   double volume;
   volume = pmesh->GetElementVolume(tet_id);
   Array<double *> Vertices_coord(4);
   for (int i = 0; i < Vertices.Size(); i++)
   {
      Vertices_coord[i] = pmesh->GetVertex(Vertices[i]);
      // std::cout<< "Vertices: "<<Vertices[i]<<"\t coordinate[x,y,z]:"<< Vertices_coord[i][0] <<"\t"<<
      //    Vertices_coord[i][1] <<"\t"<<Vertices_coord[i][2] <<"\n";
   }
   // std::cout<< Vertices_coord.Size()<<"\n";
   Vector pa(3),pb(3),pc(3),pd(3);
   pa = Vertices_coord[0];
   pb = Vertices_coord[1];
   pc = Vertices_coord[2];
   pd = Vertices_coord[3];

   Vector ab(3),ba(3),ac(3),ad(3),bc(3),cb(3),bd(3),db(3),cd(3);
   Vector oa(3),ob(3);
   double Residual = 1e-8;
   ab = MinusVector(pb,pa);
   ba = MinusVector(pa,pb);
   ac = MinusVector(pc,pa);
   ad = MinusVector(pd,pa);
   bc = MinusVector(pc,pb);
   cb = MinusVector(pb,pc);
   bd = MinusVector(pd,pb);
   db = MinusVector(pb,pd);
   cd = MinusVector(pd,pc);
   oa = MinusVector(pa,po);
   ob = MinusVector(pb,po);

   Vector nabc(3),nacd(3),nadb(3),ncbd(3);
   // calulate plane abc's outside normal vector
   nabc = CrossProduct(ab,ac);
   // get the value of the ad with nabc vector
   double valu;
   valu = DotProduct(nabc,ad);
   if(valu > 0.) nabc.Neg();

   // calulate plane acd's outside normal vector
   nacd = CrossProduct(ac,ad);
   // get the value of the ad with nacd vector
   valu = DotProduct(nacd,ab);
   if(valu > 0.) nacd.Neg();

   // calulate plane adb's outside normal vector
   nadb = CrossProduct(ab,ad);
   // get the value of the ad with nadb vector
   valu = DotProduct(nadb,ac);
   if(valu > 0.) nadb.Neg();

   // calulate plane cbd's outside normal vector
   ncbd = CrossProduct(bc,bd);
   // get the value of the ad with ncbd vector
   valu = DotProduct(ncbd,ba);
   if(valu > 0.) ncbd.Neg();

   double vabc,vacd,vadb,vcbd;
   vabc = DotProduct(nabc,oa);
   vacd = DotProduct(nacd,oa);
   vadb = DotProduct(nadb,oa);
   vcbd = DotProduct(ncbd,ob);
   
   // if abc,acd,adb,cbd's value lager 0,mean the point is located the inside of element.
   if(vabc > 0. && vacd > 0. && vadb > 0. && vcbd > 0.)
      return IN;  
   else if (vabc >= -Residual && vacd >= -Residual && vadb >= -Residual && 
        vcbd >= -Residual && std::abs(vabc*vacd*vadb*vcbd) <= Residual)
      return ONSURFACE;  
   else return OUT; 
}  

Vector Post::CrossProduct(Vector p1, Vector p2)
{
   assert(p1.Size()==3);
   assert(p2.Size()==3);
   Vector Point_(3);
   Point_[0] = p1[1]*p2[2] - p1[2]*p2[1];
   Point_[1] = -p1[0]*p2[2] + p1[2]*p2[0];
   Point_[2] = p1[0]*p2[1] - p1[1]*p2[0];
   return Point_;
}
Vector Post::MinusVector(Vector p1, Vector p2)
{
   assert(p1.Size()==3);
   assert(p2.Size()==3);
   Vector Point_(3);
   Point_[0] = p2[0] - p1[0];
   Point_[1] = p2[1] - p1[1];
   Point_[2] = p2[2] - p1[2];
   return Point_;
}

double Post::Length(Vector p1)
{
   assert(p1.Size()==3);
   double value = std::pow(p1[0],2)+
                  std::pow(p1[1],2)+
                  std::pow(p1[2],2);
   return std::sqrt(value);
}

double Post::DotProduct(Vector p1, Vector p2)
{
   assert(p1.Size()==3);
   assert(p2.Size()==3);
   return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

void Post::post_process()
{
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


   
   // // // yangcong ,2023.2.27
   // // // the run velocity is too slow
   // //in order to find sites
   // Array<int> local_sites_tets;
	// for (int p_sites = 0; p_sites < n_sites; p_sites++)
	// {
   //    Vector point(3);
   //    point[0] = global_sites_x[p_sites];
   //    point[1] = global_sites_y[p_sites];
   //    point[2] = global_sites_z[p_sites];
   //    // std::ofstream of("check.txt");
	// 	//loop all element
	// 	for (int tet_id = 0; tet_id < pmesh->GetNE() ; tet_id++)
	// 	{
   //       int tet_tag = param_handler->get_elem_tag(pmesh->attributes[tet_id]);
   //       // of<<"tet_id: "<<tet_id << " tag: "<<tet_tag<<"\n";
	// 		// /*sites point in air layer or surface*/
	// 		if( tet_tag == 1 || tet_tag ==2 ) // 1--air layer || 2-surface layer
	// 		{
   //          Vector center(3);
   //          pmesh->GetElementCenter(tet_id,center);
   //          double p2length;
   //          p2length = Length(MinusVector(point , center));
   //          Array< int > edges, cor;
   //          pmesh->GetElementEdges(tet_id,edges,cor);
   //          Vector coord_1(3),coord_2(3);
   //          Vector diameter_tp(6);
   //          double diameter;
   //          for (int i = 0; i < edges.Size(); i++)
   //          {
   //             Array< int > vert(2);
   //             pmesh->GetEdgeVertices(edges[i],vert);
   //             pmesh->GetNode(vert[0],coord_1);
   //             pmesh->GetNode(vert[1],coord_2);
   //             diameter_tp[i] = Length(MinusVector(coord_1 , coord_2));
   //          }
   //          diameter = diameter_tp.Max();
	// 			if(p2length < diameter)
	// 			{//in order to reduce calculate time cost
	// 				//check if the station point inside of the element or not;
	// 				if(SpaceRelation(TestPointInTet(point, tet_id))!= 1)
	// 				{
   //                // std::cout<<p_sites<<" is find:"<< tet_id<<"\n";
   //                local_sites_r.Append(global_sites_r[p_sites]);
   //                local_sites_theta.Append(global_sites_theta[p_sites]);
   //                local_sites_phi.Append(global_sites_phi[p_sites]);
   //                local_sites_tets.Append(tet_id);
	// 					break;
	// 				}	
	// 			}
	// 		}	
	// 	}
   //    // of.close();
   //    // return ;
	// }
   

   // // // Find local sites and tets which contain the local sites
   // //Note:: This method is not 100 percent reliable, i.e. it is not guaranteed to find a point,
   // // even if it lies inside a mesh element. 
   // DenseMatrix point_mat(3, n_sites);
   // Array<int> sites_tets;
   // Array<IntegrationPoint> ips;
   // for (int i=0; i<n_sites; i++) 
   // {
   //    point_mat(0,i) = global_sites_x[i];
   //    point_mat(1,i) = global_sites_y[i];
   //    point_mat(2,i) = global_sites_z[i];
   // }
   // pmesh->FindPoints(point_mat, sites_tets, ips);
   // Array<int> local_sites_tets;
   // for (int i=0; i<sites_tets.Size(); i++) 
   // {
   //    int tet_id = sites_tets[i];
   //    assert(tet_id != -1);
   //    if (tet_id != -2) // find at other processor
   //    {
   //       local_sites_r.Append(global_sites_r[i]);
   //       local_sites_theta.Append(global_sites_theta[i]);
   //       local_sites_phi.Append(global_sites_phi[i]);
   //       local_sites_tets.Append(tet_id);
   //    }
   // }


   /* Using gslib */
   // Hongbo Yao, 2021/12/27
   Vector vxyz(n_sites*3);
   for (int i=0; i<n_sites; i++) 
   {
      vxyz(i) = global_sites_x[i];
      vxyz(n_sites+i) = global_sites_y[i];
      vxyz(2*n_sites+i) = global_sites_z[i];
   }
   // Find tets that contain xyz coordinates
   int myrank;
   MPI_Comm_rank(pmesh->GetComm(), &myrank);
   FindPointsGSLIB finder(pmesh->GetComm());
   finder.Setup(*pmesh);
   finder.FindPoints(vxyz);
   Array<int> code = finder.GetCode();
   Array<int> rank = finder.GetProc();
   Array<int> elem = finder.GetElem();

   Array<int> local_sites_tets;
   // Check if there are some points not found
   int flag = code.Find(2); // code=2 means the point is not found, return -1 if not found
   if (flag != -1) // != -1 means code=2 is found, indicating that some points were not found
   {
      throw std::invalid_argument("GOAFEM::find_sites_tets_with_gslib(): some sites were not found! Please check your mesh.\n");
   }
   for (int i=0; i<n_sites; i++) 
   {
      if (rank[i]==myrank)
      {
         local_sites_tets.Append(elem[i]);
         local_sites_r.Append(global_sites_r[i]);
         local_sites_theta.Append(global_sites_theta[i]);
         local_sites_phi.Append(global_sites_phi[i]);
      }
   }
   finder.FreeData();

   // Compute magnetic induction B = 1.0/(-i*omega) Curl E
   int n_local_sites = local_sites_r.Size();
   Array<double> Bx_re(n_local_sites); Bx_re = 0.0;
   Array<double> Bx_im(n_local_sites); Bx_im = 0.0;
   Array<double> By_re(n_local_sites); By_re = 0.0;
   Array<double> By_im(n_local_sites); By_im = 0.0;
   Array<double> Bz_re(n_local_sites); Bz_re = 0.0;
   Array<double> Bz_im(n_local_sites); Bz_im = 0.0;
   Array<double> Ex_re(n_local_sites); Ex_re = 0.0;
   Array<double> Ex_im(n_local_sites); Ex_im = 0.0;
   Array<double> Ey_re(n_local_sites); Ey_re = 0.0;
   Array<double> Ey_im(n_local_sites); Ey_im = 0.0;
   Array<double> Ez_re(n_local_sites); Ez_re = 0.0;
   Array<double> Ez_im(n_local_sites); Ez_im = 0.0;

   for (int i=0; i<n_local_sites; i++) 
   {
      double r = local_sites_r[i];
      double theta = local_sites_theta[i]/180.0*EM::pi;
      double phi = local_sites_phi[i]/180.0*EM::pi;

      int tet_id = local_sites_tets[i];
      ElementTransformation* tr = pfes->GetElementTransformation(tet_id);
      const FiniteElement* fe = pfes->GetFE(tet_id);
      Vector pt;
      pt.SetSize(3);
      pt[0] = r*sin(theta)*cos(phi);
      pt[1] = r*sin(theta)*sin(phi);
      pt[2] = r*cos(theta);
      IntegrationPoint ip;
      tr->TransformBack(pt, ip);
      tr->SetIntPoint(&ip);

      
      // Compute Curl E by myself
      int ndof = fe->GetDof();
      int dim = fe->GetDim();
      DenseMatrix curl_shape(ndof, dim);
      fe->CalcPhysCurlShape(*tr, curl_shape);
      DenseMatrix shape(ndof, dim);
      fe->CalcPhysVShape(*tr, shape);
      Array<int> vdofs;
      pfes->GetElementVDofs(tet_id, vdofs);
      Vector U_re_elemvect, U_im_elemvect;
      U_re.GetSubVector(vdofs, U_re_elemvect);
      U_im.GetSubVector(vdofs, U_im_elemvect);

      Vector B_re, B_im;
      B_re.SetSize(3); B_re = 0.0;
      B_im.SetSize(3); B_im = 0.0;
      Vector E_re, E_im;
      E_re.SetSize(3); E_re = 0.0;
      E_im.SetSize(3); E_im = 0.0;

      // B_re = -1.0/omega*Curl E_imag
      // B_im =  1.0/omega*Curl E_real
      double inv_omega = 1.0/omega;
      for (int l=0; l<ndof; l++) 
      {
         Vector curl_shape_l;
         Vector shape_l;
         shape.GetRow(l, shape_l);
         curl_shape.GetRow(l, curl_shape_l);
         for (int n=0; n<3; n++) 
         {
            E_re[n] += shape_l[n]*U_re_elemvect[l];
            E_im[n] += shape_l[n]*U_im_elemvect[l];
            B_re[n] += -inv_omega*curl_shape_l[n]*U_im_elemvect[l];
            B_im[n] += inv_omega*curl_shape_l[n]*U_re_elemvect[l];
         }
      }

      Bx_re[i] = B_re[0];
      Bx_im[i] = B_im[0];
      By_re[i] = B_re[1];
      By_im[i] = B_im[1];
      Bz_re[i] = B_re[2];
      Bz_im[i] = B_im[2];
      
      Ex_re[i] = E_re[0];
      Ex_im[i] = E_im[0];
      Ey_re[i] = E_re[1];
      Ey_im[i] = E_im[1];
      Ez_re[i] = E_re[2];
      Ez_im[i] = E_im[2];
       /*
      // Compute Curl E using MFEM's GridFunction's GetCurl
      // Note that these two results are the same
      Vector curl_E_r, curl_E_i;
      U_re.GetCurl(*tr, curl_E_r);
      U_im.GetCurl(*tr, curl_E_i);
      // // e^iwt
      // // // B_re = -1.0/omega*Curl E_imag
      // // // B_im =  1.0/omega*Curl E_real
      double inv_omega = 1.0/omega;
      Bx_re[i] = -inv_omega*curl_E_i[0];
      Bx_im[i] =  inv_omega*curl_E_r[0];
      By_re[i] = -inv_omega*curl_E_i[1];
      By_im[i] =  inv_omega*curl_E_r[1];
      Bz_re[i] = -inv_omega*curl_E_i[2];
      Bz_im[i] =  inv_omega*curl_E_r[2];
      */
      

      // // e^-iwt
      // // B_re = 1.0/omega*Curl E_imag
      // // B_im =  -1.0/omega*Curl E_real
      // double inv_omega = 1.0/omega;
      // Bx_re[i] = inv_omega*curl_E_i[0];
      // Bx_im[i] =  -inv_omega*curl_E_r[0];
      // By_re[i] = inv_omega*curl_E_i[1];
      // By_im[i] =  -inv_omega*curl_E_r[1];
      // Bz_re[i] = inv_omega*curl_E_i[2];
      // Bz_im[i] =  -inv_omega*curl_E_r[2];
      // std::cout<<Bx_re[i]<<"\t";
      // std::getchar();
   }

   // Convert to geocentric spherical components
   Br_re.SetSize(n_local_sites); Br_re = 0.0;
   Br_im.SetSize(n_local_sites); Br_im = 0.0;
   Bt_re.SetSize(n_local_sites); Bt_re = 0.0;
   Bt_im.SetSize(n_local_sites); Bt_im = 0.0;
   Bp_re.SetSize(n_local_sites); Bp_re = 0.0;
   Bp_im.SetSize(n_local_sites); Bp_im = 0.0;
   Er_re.SetSize(n_local_sites); Er_re = 0.0;
   Er_im.SetSize(n_local_sites); Er_im = 0.0;
   Et_re.SetSize(n_local_sites); Et_re = 0.0;
   Et_im.SetSize(n_local_sites); Et_im = 0.0;
   Ep_re.SetSize(n_local_sites); Ep_re = 0.0;
   Ep_im.SetSize(n_local_sites); Ep_im = 0.0;
   for (int i=0; i<local_sites_r.Size(); i++) 
   {
      double r = local_sites_r[i];
      double theta = local_sites_theta[i];
      double phi = local_sites_phi[i];

      Br_re[i] = sind(theta)*cosd(phi)*Bx_re[i] + sind(theta)*sind(phi)*By_re[i] + cosd(theta)*Bz_re[i];
      Br_im[i] = sind(theta)*cosd(phi)*Bx_im[i] + sind(theta)*sind(phi)*By_im[i] + cosd(theta)*Bz_im[i];
      Bt_re[i] = cosd(theta)*cosd(phi)*Bx_re[i] + cosd(theta)*sind(phi)*By_re[i] - sind(theta)*Bz_re[i];
      Bt_im[i] = cosd(theta)*cosd(phi)*Bx_im[i] + cosd(theta)*sind(phi)*By_im[i] - sind(theta)*Bz_im[i];
      Bp_re[i] = -sind(phi)*Bx_re[i] + cosd(phi)*By_re[i];
      Bp_im[i] = -sind(phi)*Bx_im[i] + cosd(phi)*By_im[i];
   
      Er_re[i] = sind(theta)*cosd(phi)*Ex_re[i] + sind(theta)*sind(phi)*Ey_re[i] + cosd(theta)*Ez_re[i];
      Er_im[i] = sind(theta)*cosd(phi)*Ex_im[i] + sind(theta)*sind(phi)*Ey_im[i] + cosd(theta)*Ez_im[i];
      Et_re[i] = cosd(theta)*cosd(phi)*Ex_re[i] + cosd(theta)*sind(phi)*Ey_re[i] - sind(theta)*Ez_re[i];
      Et_im[i] = cosd(theta)*cosd(phi)*Ex_im[i] + cosd(theta)*sind(phi)*Ey_im[i] - sind(theta)*Ez_im[i];
      Ep_re[i] = -sind(phi)*Ex_re[i] + cosd(phi)*Ey_re[i];
      Ep_im[i] = -sind(phi)*Ex_im[i] + cosd(phi)*Ey_im[i];
      // T to nT
      Br_re[i] *= 1e9;
      Br_im[i] *= 1e9;
      Bt_re[i] *= 1e9;
      Bt_im[i] *= 1e9;
      Bp_re[i] *= 1e9;
      Bp_im[i] *= 1e9;
   }
}



void Post::save(std::ofstream& out)
{
   for (int i=0; i<local_sites_r.Size(); i++) 
   {
      out.setf(std::ios::scientific);
      out.precision(6);
      out 	<< std::setw(15) << local_sites_r[i]/1000.0 << "\t" 
            << std::setw(15) << local_sites_theta[i] << "\t" 
            << std::setw(15) << local_sites_phi[i] << "\t" 
            << std::setw(15) << Br_re[i] << "\t" 
            << std::setw(15) << Br_im[i] << "\t" 
            << std::setw(15) << Bt_re[i] << "\t" 
            << std::setw(15) << Bt_im[i] << "\t" 
            << std::setw(15) << Bp_re[i] << "\t" 
            << std::setw(15) << Bp_im[i] << "\t"
            << std::setw(15) << Er_re[i] << "\t" 
            << std::setw(15) << Er_im[i] << "\t"
            << std::setw(15) << Et_re[i] << "\t"
            << std::setw(15) << Et_im[i] << "\t"
            << std::setw(15) << Ep_re[i] << "\t"
            << std::setw(15) << Ep_im[i] << "\n";
   }
}

void Post::save_as_one(std::ofstream& out)
{
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Status status;
   MPI_Comm mycomm = pfes->GetComm();
   int myid, n_procs;
   MPI_Comm_size(mycomm, &n_procs);
   MPI_Comm_rank(mycomm, &myid);

   int count;
   double* msg_r;
   double* msg_theta;
   double* msg_phi;
   double* msg_Br_re;
   double* msg_Br_im;
   double* msg_Bt_re;
   double* msg_Bt_im;
   double* msg_Bp_re;
   double* msg_Bp_im;

   double* msg_Er_re;
   double* msg_Er_im;
   double* msg_Et_re;
   double* msg_Et_im;
   double* msg_Ep_re;
   double* msg_Ep_im;

   if (myid == 0) 
   {
      this->save(out);
      for (int p=1; p<n_procs; p++) 
      {
         MPI_Recv(&count, 1, MPI_INT, p, 101, mycomm, &status);
         msg_r = new double[count];
         msg_theta = new double[count];
         msg_phi = new double[count];
         msg_Br_re = new double[count];
         msg_Br_im = new double[count];
         msg_Bt_re = new double[count];
         msg_Bt_im = new double[count];
         msg_Bp_re = new double[count];
         msg_Bp_im = new double[count];

         msg_Er_re = new double[count];
         msg_Er_im = new double[count];
         msg_Et_re = new double[count];
         msg_Et_im = new double[count];
         msg_Ep_re = new double[count];
         msg_Ep_im = new double[count];
         
         MPI_Recv(msg_r, count, MPI_DOUBLE, p, 102, mycomm, &status);
         MPI_Recv(msg_theta, count, MPI_DOUBLE, p, 103, mycomm, &status);
         MPI_Recv(msg_phi, count, MPI_DOUBLE, p, 104, mycomm, &status);
         MPI_Recv(msg_Br_re, count, MPI_DOUBLE, p, 105, mycomm, &status);
         MPI_Recv(msg_Br_im, count, MPI_DOUBLE, p, 106, mycomm, &status);
         MPI_Recv(msg_Bt_re, count, MPI_DOUBLE, p, 107, mycomm, &status);
         MPI_Recv(msg_Bt_im, count, MPI_DOUBLE, p, 108, mycomm, &status);
         MPI_Recv(msg_Bp_re, count, MPI_DOUBLE, p, 109, mycomm, &status);
         MPI_Recv(msg_Bp_im, count, MPI_DOUBLE, p, 110, mycomm, &status);
         MPI_Recv(msg_Er_re, count, MPI_DOUBLE, p, 111, mycomm, &status);
         MPI_Recv(msg_Er_im, count, MPI_DOUBLE, p, 112, mycomm, &status);
         MPI_Recv(msg_Et_re, count, MPI_DOUBLE, p, 113, mycomm, &status);
         MPI_Recv(msg_Et_im, count, MPI_DOUBLE, p, 114, mycomm, &status);
         MPI_Recv(msg_Ep_re, count, MPI_DOUBLE, p, 115, mycomm, &status);
         MPI_Recv(msg_Ep_im, count, MPI_DOUBLE, p, 116, mycomm, &status);


         for (int i=0; i<count; i++) 
         {
            out.setf(std::ios::scientific);
            out.precision(6);
            out 	<< std::setw(15) << msg_r[i]/1000.0 << "\t" 
                  << std::setw(15) << msg_theta[i] << "\t" 
                  << std::setw(15) << msg_phi[i] << "\t" 
                  << std::setw(15) << msg_Br_re[i] << "\t" 
                  << std::setw(15) << msg_Br_im[i] << "\t" 
                  << std::setw(15) << msg_Bt_re[i] << "\t" 
                  << std::setw(15) << msg_Bt_im[i] << "\t" 
                  << std::setw(15) << msg_Bp_re[i] << "\t" 
                  << std::setw(15) << msg_Bp_im[i] << "\t"
                  << std::setw(15) << msg_Er_re[i] << "\t" 
                  << std::setw(15) << msg_Er_im[i] << "\t" 
                  << std::setw(15) << msg_Et_re[i] << "\t" 
                  << std::setw(15) << msg_Et_im[i] << "\t" 
                  << std::setw(15) << msg_Ep_re[i] << "\t" 
                  << std::setw(15) << msg_Ep_im[i] << "\n";                 
         }
         
         // delete 
         delete msg_r;
         delete msg_theta;
         delete msg_phi;
         delete msg_Br_re;
         delete msg_Br_im;
         delete msg_Bt_re;
         delete msg_Bt_im;
         delete msg_Bp_re;
         delete msg_Bp_im;
         delete msg_Er_re;
         delete msg_Er_im;
         delete msg_Et_re;
         delete msg_Et_im;
         delete msg_Ep_re;
         delete msg_Ep_im;
      }
   }
   else 
   {
      count = local_sites_r.Size();
      msg_r = local_sites_r.GetData();
      msg_theta = local_sites_theta.GetData();
      msg_phi = local_sites_phi.GetData();
      msg_Br_re = Br_re.GetData();
      msg_Br_im = Br_im.GetData();
      msg_Bt_re = Bt_re.GetData();
      msg_Bt_im = Bt_im.GetData();
      msg_Bp_re = Bp_re.GetData();
      msg_Bp_im = Bp_im.GetData();
      msg_Er_re = Er_re.GetData();
      msg_Er_im = Er_im.GetData();
      msg_Et_re = Et_re.GetData();
      msg_Et_im = Et_im.GetData();
      msg_Ep_re = Ep_re.GetData();
      msg_Ep_im = Ep_im.GetData();
         
      MPI_Send(&count, 1, MPI_INT, 0, 101, mycomm);
      MPI_Send(msg_r, count, MPI_DOUBLE, 0, 102, mycomm);
      MPI_Send(msg_theta, count, MPI_DOUBLE, 0, 103, mycomm);
      MPI_Send(msg_phi, count, MPI_DOUBLE, 0, 104, mycomm);
      MPI_Send(msg_Br_re, count, MPI_DOUBLE, 0, 105, mycomm);
      MPI_Send(msg_Br_im, count, MPI_DOUBLE, 0, 106, mycomm);
      MPI_Send(msg_Bt_re, count, MPI_DOUBLE, 0, 107, mycomm);
      MPI_Send(msg_Bt_im, count, MPI_DOUBLE, 0, 108, mycomm);
      MPI_Send(msg_Bp_re, count, MPI_DOUBLE, 0, 109, mycomm);
      MPI_Send(msg_Bp_im, count, MPI_DOUBLE, 0, 110, mycomm);
      MPI_Send(msg_Er_re, count, MPI_DOUBLE, 0, 111, mycomm);
      MPI_Send(msg_Er_im, count, MPI_DOUBLE, 0, 112, mycomm);
      MPI_Send(msg_Et_re, count, MPI_DOUBLE, 0, 113, mycomm);
      MPI_Send(msg_Et_im, count, MPI_DOUBLE, 0, 114, mycomm);
      MPI_Send(msg_Ep_re, count, MPI_DOUBLE, 0, 115, mycomm);
      MPI_Send(msg_Ep_im, count, MPI_DOUBLE, 0, 116, mycomm);
   }
}
