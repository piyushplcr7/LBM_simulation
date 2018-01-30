/**
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED
//#include "flagella.hpp"
#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include "elb.hpp"
//#include <iostream>
//#include <cmath>

namespace lb {

/**
 *  @brief Simulation class implementing LB
 *
 *  This class holds a lattice as member (see @ref simulation::l) and
 *  carries out the simulation steps on top of it. The main methods of
 *  this class are @ref simulation::advect() and
 *  @ref simulation::collide().
 */
class simulation
{
public: // ctor

	/**
	 *  @brief Construct from domain size and flow parameters
	 *  @param[in] nx    extent in x direction
	 *  @param[in] ny    extent in y direction
	 *  @param[in] _Re   Reynolds number
	 *  @param[in] _Vmax mean flow velocity
	 */
	simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax)
	: l(nx, ny),
	  shift(velocity_set().size),
	  Re(_Re),
	  Vmax(_Vmax),
	 	Cyl_radius(10),
	  visc(Vmax*Cyl_radius*2/Re),
	  beta(1./(6*visc+1)),
	  time(0),
	  file_output(false), // set to true if you want to write files
	  output_freq(100),
	  output_index(0),
		flag_moving_cyl(!true),
		using_entropic(!true),
		using_flagella(!true)
	{
		// define amount to shift populations for advection (according to the array model of domain)
		for (unsigned int i=0; i<velocity_set().size; ++i)
			shift[i]=velocity_set().c[0][i] + l.real_nx * velocity_set().c[1][i];
	}

	/**
	 *  @brief Initialize the flow field
	 *
	 *  Initialization includes defining initial density, velocity and
	 *  populations. You can use Taylor-Green vortex flow conditions.
	 */
	void initialize()
	{

		//Initialize Cyclinder
		Cyl_center_0[0] = l.nx/4+.5 ;//C yl_radius*10;
		Cyl_center_0[1] = l.ny/2 +.5;//Cyl_radius*10;
		Cyl_center[0] = Cyl_center_0[0];
		Cyl_center[1] = Cyl_center_0[1];
		Cyl_vel[0] = 0.0;
		Cyl_vel[1] = 0.0;
		unsigned int n_links = 1;

		partition = Cyl_center_0[0] + Cyl_radius - 10 ;

		if (!using_flagella)
			partition = l.nx;
		//force.open("Force.txt",std::ios::out);
		//force << std::setw(10) << "Fx: " << std::setw(10) << "Fy: " << "\n";
		std::cout << " Starting to mark solid nodes using add_wallCylinder" << std::endl;
		l.add_wallCylinder(Cyl_center, Cyl_vel, Cyl_radius, using_flagella, partition);

		//unsigned int I_1; unsigned int J_1;
		//I_1 = l.cylinder_fbn.front().coord.i;
		//J_1 = l.cylinder_fbn.front().coord.j;
		//std::cout << "size of cylinder_fbn: " << l.cylinder_fbn.size() << "& mps" << l.get_node(I_1,J_1).missing_populations.size() << std::endl;
		std::cout << "Cylinder nodes added" <<std::endl;


		//Init B.C.
		u_inlet = 0.05;
		runUptime = 0;
		rho_inlet = 1;
		bool flag_Complete_Init = true;
		//const float_type pi(std::acos(-1.0));

		//#pragma omp parallel for
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			double ux,uy,rho;
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				if( !(l.get_node(i,j).has_flag_property("solid")) )
				{
					//Initialize flow around Cylinder
					if(flag_Complete_Init) {ux = u_inlet;}
					else{ux = 0.0;}
					uy = 0.0;
					rho = 1;
					l.get_node(i,j).u()  = ux;
					l.get_node(i,j).v()  = uy;
					l.get_node(i,j).rho() = rho;
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).f(k)=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
					}
				}
				else
				{
					l.get_node(i,j).u()  = Cyl_vel[0];
					l.get_node(i,j).v()  = Cyl_vel[1];
					l.get_node(i,j).rho() = 1.0;
					//for inside the solid, initialize with 0 value for populations
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).f(k)=0.;
					}
				}
			}
		} //loop j ends

		
		//Reset vectors of boundary nodes again, because the old population are needed!!
		//resetBN(n_links);
		l.cylinder_fbn.clear();

		l.add_wallCylinder(Cyl_center, Cyl_vel, Cyl_radius,using_flagella,partition);
		//initialize flagella
		if (using_flagella)
		{
			coordinate<float> attach_point = {Cyl_center[0]+Cyl_radius,Cyl_center[1]};
			std::cout << "Attachment Point: " << attach_point.i << "   " << attach_point.j << std::endl;
			std::cout << "Partition Value: " << partition << std::endl;
			double length = 2*Cyl_radius, mass = 200, K=100, c=100;
			flg = new flagella(n_links, length, mass, K, c, attach_point.i, attach_point.j, 0, 0);
			std::cout << "Flagella created" << std::endl;
			//adding the flagella nodes and corresponding fluid boundary nodes to the lattice
			l.add_flagella_nodes(flg, Cyl_vel, Cyl_center, Cyl_radius, partition);
			std::cout << "Flagella nodes added" <<std::endl;
		}

		//Merge different fluid boundary nodes with different directions into one
		l.merge_into_fbn(using_flagella);

		std::cout << "Fluid Boundary Nodes merged " << std::endl;

		//std::cout << "&&&&size of cylinder_fbn: " << l.cylinder_fbn.size() << "& mps" << l.get_node(I_1,J_1).missing_populations.size() << std::endl;


	}

	/**
	 *  @brief advect the populations
	 *
	 *  Include periodic boundary conditions here also
	 */
	void advect() //advection to be done in parallel
	{
		for (unsigned int k=0; k<velocity_set().size; ++k)
		{

			if (shift[k] == 1 || shift[k] == -1)
			{
				#pragma omp parallel for schedule(static)
				for (int j=0; j<static_cast<int>(l.ny); ++j)
				{
					for (int i=0; i<static_cast<int>(l.nx); ++i)
					{
							//right top index or increasing index
							unsigned int rtindex= l.real_nx*(j+l.buffer_size) + i + l.buffer_size;
							//left bottom index or decreasing index
							unsigned int lbindex= l.real_nx*(l.ny-1-j+l.buffer_size) + l.nx-1 -i + l.buffer_size;
							//PUSH ADVECTION!
								if (shift[k]>0) //if shift is positive, using decreasing index
									l.f[k][lbindex+shift[k]]=l.f[k][lbindex];
								else
									l.f[k][rtindex+shift[k]]=l.f[k][rtindex];
					}
				}
			}

			else if (shift[k] > 1 || shift[k] < -1)
			{
				for (int j=0; j<static_cast<int>(l.ny); ++j)
				{
					#pragma omp parallel for schedule(static)
					for (int i=0; i<static_cast<int>(l.nx); ++i)
					{
							//right top index or increasing index
							unsigned int rtindex= l.real_nx*(j+l.buffer_size) + i + l.buffer_size;
							//left bottom index or decreasing index
							unsigned int lbindex= l.real_nx*(l.ny-1-j+l.buffer_size) + l.nx-1 -i + l.buffer_size;
							//PUSH ADVECTION!
								if (shift[k]>0) //if shift is positive, using decreasing index
									l.f[k][lbindex+shift[k]]=l.f[k][lbindex];
								else
									l.f[k][rtindex+shift[k]]=l.f[k][rtindex];
					}
				}
			}

		}
	}

		/**  @brief calculate fraction of fluid from node n to the solid boundary along the direction i */
	/*static*/ float_type get_qi(const node& n,const int& i) /*const*/{
		double x1 = (double)n.coord.i; double y1 = (double)n.coord.j;
		double x2 =  x1 + (double)velocity_set().c[0][i]; double y2 =  y1 + (double)velocity_set().c[1][i];
		assert(l.get_node((int)x2,(int)y2).has_flag_property("solid"));
		//calc distance c_i
		float_type dx = x2-x1;
		float_type dy = y2-y1;

		float_type c_i = sqrt((dx)*(dx)+(dy)*(dy));

		//std::cout << "dx " << dx << " dy " << dy << " c_i " << c_i << std::endl;

		//calc x_w,i
		//find root of quadratic equation
		//((x1+t*dx -l.center[0])*(x1+t*dx -l.center[0]) + (y1+t*dy -l.center[1])*(y1+t*dy -l.center[1])-l.radius*l.radius == 0)
		float_type a = (dx*dx + dy*dy);
		float_type b = 2*dx* (x1-Cyl_center[0]) + 2*dy * (y1 - Cyl_center[1]);
		float_type c = (x1-Cyl_center[0])*(x1-Cyl_center[0]) + (y1-Cyl_center[1])*(y1-Cyl_center[1]) - Cyl_radius*Cyl_radius;
		float_type t = solve_quadratic(a,b,c);
		float_type q_i = sqrt(t*dx *t*dx + t*dy *t*dy)/c_i;
		//std::cout << "Q_i: " << q_i << std::endl;

		return q_i;

	}


	/**  @brief calculate smaller root of quadratic function */
	inline float_type solve_quadratic(float_type a, float_type b, float_type c){
		return -b/2/a - sqrt(b*b-4*a*c)/2/a;
	}


	/**  @brief Calculate Eq Presure tensor */
	void calc_Peq(int i, int j, float_type (&Peq)[2][2]){
		double rho = l.get_node(i,j).rho(); double u = l.get_node(i,j).u();	double v = l.get_node(i,j).v();	double cs = velocity_set().cs;
		Peq[0][0] = rho*(cs*cs + u*u) ;
		Peq[1][0] = rho*(u*v);
		Peq[0][1] = Peq[1][0];
		Peq[1][1] = rho*(cs*cs + v*v);
	}


	/**  @brief Calculate NonEq for with direct influence of Boundary */
	void calc_Pneq(int i, int j/*, bool dir_solid[9]*/, float_type q_i[9], float_type (&Pneq)[2][2]){

		float_type du_x =0, du_y=0, dv_x =0, dv_y = 0;
		//calculate velocity gradients
		if( !(q_i[3]==0) and !(q_i[1] == 0) ){
			du_x = (l.get_node(i+1,j).u()-l.get_node(i-1,j).u())/(q_i[3] + q_i[1]);
			dv_x = (l.get_node(i+1,j).v()-l.get_node(i-1,j).v())/(q_i[3] + q_i[1]);
		}

		if( !(q_i[4] == 0) and !(q_i[2] == 0) ){
			du_y = (l.get_node(i,j+1).u()-l.get_node(i,j-1).u())/(q_i[4] + q_i[2]);
			dv_y = (l.get_node(i,j+1).v()-l.get_node(i,j-1).v())/(q_i[4] + q_i[2]);
		}
		double rho = l.get_node(i,j).rho(); double cs = velocity_set().cs;
		//Set non equilibrium pressure tensor
		Pneq[0][0] = - rho*cs*cs/beta*(du_x);
		Pneq[1][0] = -rho*cs*cs/beta/2*(du_y + dv_x);
		Pneq[0][1] = Pneq[1][0];
		Pneq[1][1] = - rho*cs*cs/beta*(dv_y);
	}

	/**  @brief Calculate NonEq for with direct influence of Boundary */
	void calc_Pneq(int i, int j, const std::vector<std::pair<double,double>>& uv, float_type q_i[9], float_type (&Pneq)[2][2]){

		float_type du_x =0, du_y=0, dv_x =0, dv_y = 0;
		//calculate velocity gradients
		if(1 or (!(q_i[3]==0) and !(q_i[1] == 0) ))
		{
			du_x = (uv[1].first - uv[3].first)/(q_i[3] + q_i[1]);
			dv_x = (uv[1].second - uv[3].second)/(q_i[3] + q_i[1]);
		}

		if(1 or (!(q_i[4] == 0) and !(q_i[2] == 0)) )
		{
			du_y = (uv[2].first - uv[4].first)/(q_i[2] + q_i[4]);
			dv_y = (uv[2].second - uv[4].second)/(q_i[2] + q_i[4]);
		}
		double rho = l.get_node(i,j).rho(); double cs = velocity_set().cs;
		//Set non equilibrium pressure tensor
		Pneq[0][0] = - rho*cs*cs/beta*(du_x);
		Pneq[1][0] = -rho*cs*cs/beta/2*(du_y + dv_x);
		Pneq[0][1] = Pneq[1][0];
		Pneq[1][1] = - rho*cs*cs/beta*(dv_y);
	}

	/**  @brief Calculate NonEq for without influence of Boundary */
	void calc_Pneq(int i, int j, float_type (&Pneq)[2][2]){

		//calculate velocity gradients with central difference
		float_type du_x = (l.get_node(i+1,j).u()-l.get_node(i-1,j).u())/2;
		float_type du_y = (l.get_node(i+1,j).u()-l.get_node(i-1,j).u())/2;
		float_type dv_x = (l.get_node(i+1,j).v()-l.get_node(i-1,j).v())/2;
		float_type dv_y = (l.get_node(i,j+1).v()-l.get_node(i,j-1).v())/2;


		//Set non equilibrium pressure tensor
		Pneq[0][0] = - l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/beta*(du_x);
		Pneq[1][0] = - l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/2/beta*(du_y + dv_x);
		Pneq[0][1] = Pneq[1][0];
		Pneq[1][1] = - l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/beta*(dv_y);

		return;
	}
	//function to invert the population index
	int inv_popl(const int& i) const{
		int inverse[9] = {0,3,4,1,2,7,8,5,6};
		return inverse[i];
	}

	unsigned int getIndex(unsigned int Node_index, std::vector<node>& bn){
		unsigned int i;
		for(auto it = bn.begin(); it != bn.end(); ++it){
			if(it->index == Node_index){
				i = std::distance(bn.begin(),it);
			}
		}
		return i;
	}

	struct predicate {
		predicate(unsigned int indx):ind_(indx){}
		bool operator()(const node& n) {
			return n.index==ind_;
		}
		unsigned int ind_;
	};


	/**  @brief apply wall boundary conditions */
	void curved_wall_bc()
	{
		std::ofstream scatter;
		scatter.open("Nodes_scatter.txt",std::ios::app);
		if(time==0) 
			scatter << "#xj " << std::setw(15) << " yj " << std::endl;
		//Curved Wall
		#pragma omp parallel for schedule(dynamic)
		for (unsigned int i=0; i<l.fluid_boundary_nodes.size(); ++i)
		{
			unsigned int x_i = l.fluid_boundary_nodes[i].coord.i;
			unsigned int y_j = l.fluid_boundary_nodes[i].coord.j;
			if(time ==0) 
				scatter << x_i << std::setw(15) << y_j << std::endl;
			/*unsigned int index = l.fluid_boundary_nodes[i].index;
			predicate p(index);
			auto M_size_cyl_it = std::find_if(l.cylinder_fbn.begin(),l.cylinder_fbn.end(),p);
			unsigned int M_size_cyl = 0; if (M_size_cyl_it != l.cylinder_fbn.end()) M_size_cyl = M_size_cyl_it->missing_populations.size();

			auto M_size_cyl_f_it= std::find_if(l.cylinder_fbn_f.begin(),l.cylinder_fbn_f.end(),p);	
			unsigned int M_size_cyl_f = 0; if (M_size_cyl_f_it != l.cylinder_fbn_f.end()) M_size_cyl_f = M_size_cyl_f_it->missing_populations.size();

			auto M_size_flag_it = std::find_if(l.flagella_nodes[0].begin(),l.flagella_nodes[0].end(),p);
			unsigned int M_size_flag = 0; if (M_size_flag_it != l.flagella_nodes[0].end() ) M_size_flag = M_size_flag_it->missing_populations.size();

			unsigned int M_size_node = l.get_node(x_i, y_j).missing_populations.size();
			std::cout << " Missing populations size in global node: " << M_size_node << "// in addition" << M_size_flag+ M_size_cyl + M_size_cyl_f << std::endl;
			assert(M_size_flag+ M_size_cyl + M_size_cyl_f == M_size_node);*/


			//std::cout << "Applying curved BC for "<<x_i << " " << y_j << std::endl;
			//std::cout << "(x_i,y_j) = (" << x_i<< " " << y_j << ") "  << std::endl;
			//std::cout << "l.get_node(x_i,y_j) = " << l.get_node(x_i,y_j).u() << "  "  << std::endl;//<< l.get_node(x_i,y_j).missing_populations.size() << std::endl;
			std::vector<int> i_d(l.get_node(x_i,y_j).return_missing_populations()); //contains the direction in which solid is encountered
			std::vector<float_type> q = l.get_node(x_i,y_j).q_i;												//contains the corresponding q values
			std::vector<std::pair<double,double> > uvw = l.get_node(x_i,y_j).uvw_i;	//contains the corresponding wall velocities
			float_type utgt[2] = {0,0};
			float_type rho_bb = 0;
			float_type rho_s = 0;
			float_type rho_0 = 1;
			float_type q_i[9] = {1,1,1,1,1,1,1,1,1};
			assert(i_d.size()==q.size() and i_d.size() == uvw.size());							//consistency check
			for (unsigned int j = 0; j < i_d.size(); ++j)
			{
				q_i[i_d[j]] = q[j];
			}
			std::vector<float_type> u_P(8);

			for (unsigned int j = 0; j < i_d.size(); ++j)
			{
				//Velocities away from the solid (for the missing populations D bar )
				int c_x = -velocity_set().c[0][i_d[j]];
				int c_y = -velocity_set().c[1][i_d[j]];
				//get q_i for the direction given by i_d[j]
				float_type u_fi[2];
				//get u_f,i for the fluid node (in the opposite direction from the one where solid is encountered)
				if ( q_i[inv_popl(i_d[j])] == 1. ){
					u_fi[0] = l.get_node(x_i+c_x,y_j+c_y).u(); u_fi[1]=l.get_node(x_i+c_x,y_j+c_y).v();
				}
				else
					{
						unsigned int tempit;
						for (tempit = 0; tempit<i_d.size() ; ++tempit)
							{
								if ( i_d[tempit] == inv_popl(i_d[j]) )
									break;
							}
						float_type utemp,vtemp;
						std::tie(utemp,vtemp) = uvw[tempit];
						u_fi[0] = utemp; u_fi[1] = vtemp;
					}
					//float_type u_fi[2] = {l.get_node(x_i+c_x,y_j+c_y).u(), l.get_node(x_i+c_x,y_j+c_y).v()};
				//std::cout << "q_i["<<i_d[j]<<"] = " << q_i[i_d[j]] << " Out of " << i_d.size() << " & u_fi = {" << u_fi[0] << "," << u_fi[1] << "}" ;
				//wall velocity

				float_type u_w_i,v_w_i;
				std::tie(u_w_i,v_w_i) = uvw[j];
				//calc utgt part
				//utgt[0] += (q_i[i_d[j]]* u_fi[0] + u_w_i)/(1+q_i[i_d[j]])/i_d.size();
				//utgt[1] += (q_i[i_d[j]]* u_fi[1] + v_w_i)/(1+q_i[i_d[j]])/i_d.size();
				if( !(q_i[inv_popl(i_d[j])] == 0) and !(q_i[i_d[j]]==0) ){ //ONly set to non zero if the distance is not zero for both direction
					utgt[0] += (q_i[i_d[j]]* u_fi[0] + q_i[inv_popl(i_d[j])]*u_w_i)/(q_i[inv_popl(i_d[j])]+q_i[i_d[j]])/i_d.size();
					utgt[1] += (q_i[i_d[j]]* u_fi[1] + q_i[inv_popl(i_d[j])]*v_w_i)/(q_i[inv_popl(i_d[j])]+q_i[i_d[j]])/i_d.size();
				}

				//get rho_s
				rho_s += 6*rho_0*velocity_set().W[i_d[j]]*(c_x*u_w_i + c_y*v_w_i);
				//adjusting the bounce-back populations (the D-bar populations, which are inverse of the missing_populations)
				l.get_node(x_i, y_j).f( inv_popl(i_d[j]) ) = l.get_node(x_i, y_j).f(i_d[j]);
				//std::cout << "Inversed Population:  " << l.get_node(x_i, y_j).f( inv_popl(i_d[j]) ) << std::endl;

			}

			for(unsigned int j = 0; j < velocity_set().size ; ++j){
				//rho_bb += l.get_node(x_i, y_j).f(j);
				rho_bb += l.get_node(x_i, y_j).f(j);
			}

			//
			std::vector<std::pair<double,double>> uv(1);
			for (unsigned int i=1 ; i<5 ; ++i)
			{
				if (q_i[i]<1)
					{
						auto found = std::find(std::begin(i_d),std::end(i_d),i);
						unsigned int req_ind = std::distance(i_d.begin(),found);
						uv.push_back(uvw[req_ind]);
					}
				else
					{
						double tempu = l.get_node(x_i+velocity_set().c[0][i] , y_j+velocity_set().c[1][i]).u();
						double tempv = l.get_node(x_i+velocity_set().c[0][i] , y_j+velocity_set().c[1][i]).v();
						uv.push_back(std::make_pair(tempu,tempv));
					}
			}


			//rho_tgt
			float_type rho_tgt = rho_bb + rho_s;

			l.get_node(x_i, y_j).rho() = rho_tgt;
			l.get_node(x_i, y_j).u() = utgt[0];
			l.get_node(x_i, y_j).v() = utgt[1];
			//Calculate missing populations
			float_type Peq[2][2];
			calc_Peq(x_i, y_j, Peq);
			float_type Pneq[2][2];
			calc_Pneq(x_i, y_j ,uv, q_i, Pneq);
			float_type f_new;

			//std::cout << "Peq: " << Peq[0][0] << " /" << Peq[0][1] << " /" << Peq[1][0] << " /" << Peq[1][1] << " /" << std::endl;

			for (unsigned int j = 0; j < i_d.size(); ++j)
			{
				//Reset population with first two parts of grads' approximation
				//f_new = velocity_set().W[i_d[j]]*(l.get_node(x_i, y_j).rho()*(1 + (velocity_set().c[0][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).u() + velocity_set().c[1][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).v())/velocity_set().cs/velocity_set().cs));
				f_new = (l.get_node(x_i, y_j).rho()*(1 + (velocity_set().c[0][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).u() + velocity_set().c[1][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).v())/velocity_set().cs/velocity_set().cs));
				//std::cout << "First part " << f_new << "// rho: " <<l.get_node(x_i, y_j).rho() << " // u: " << l.get_node(x_i, y_j).u() << " // v: " << l.get_node(x_i, y_j).v() << std::endl;
				//Update with Tensor part and summation over alpha and Beta
				for (int alpha = 0; alpha < 2; ++alpha){
					for( int beta = 0; beta < 2; ++beta){
						if( alpha == beta)
							f_new +=  1/2/pow(velocity_set().cs,4)*((Pneq[alpha][beta]+Peq[alpha][beta]-l.get_node(x_i, y_j).rho()*velocity_set().cs*velocity_set().cs)*(velocity_set().c[alpha][inv_popl(i_d[j])]*velocity_set().c[beta][inv_popl(i_d[j])] - velocity_set().cs*velocity_set().cs ));
						else
							f_new +=  1/2/pow(velocity_set().cs,4)*((Pneq[alpha][beta]+Peq[alpha][beta])*(velocity_set().c[alpha][inv_popl(i_d[j])]*velocity_set().c[beta][inv_popl(i_d[j])] ));
					}
				}
				f_new *= velocity_set().W[i_d[j]];
				l.get_node(x_i, y_j).f( inv_popl(i_d[j]) ) = f_new;
				//std::cout << "Fnew:  " << l.get_node(x_i, y_j).f( inv_popl(i_d[j]) ) << std::endl;

				//if (!std::isnan(f_new))
					//std::cout << "F_ new_ missing: " << f_new << std::endl;
			}

		}
		scatter.close();
		
	}




	void periodic_bc(){
	//using the buffers to implement periodic boundary conditions
		#pragma omp parallel for
		for (int i=0; i<static_cast<int>(l.nx); ++i)
		//iteration over the top and bottom buffers (filled due to up and down avection)
		{
			//south direction
			l.get_node(i,l.ny-1).f(4)=l.get_node(i,-1).f(4);
			//north direction
			l.get_node(i,0).f(2)=l.get_node(i,l.ny).f(2);
			//south east
			l.get_node(i,l.ny-1).f(8)=l.get_node(i,-1).f(8);
			//south west
			l.get_node(i,l.ny-1).f(7)=l.get_node(i,-1).f(7);
			//north east
			l.get_node(i,0).f(5)=l.get_node(i,l.ny).f(5);
			//north west
			l.get_node(i,0).f(6)=l.get_node(i,l.ny).f(6);
		}
		#pragma omp parallel for
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		//iteration over the left and right buffer (filled due to left right avection)
		{
			//east direction
			l.get_node(0,j).f(1)=l.get_node(l.nx,j).f(1);
			//west direction
			l.get_node(l.nx-1,j).f(3)=l.get_node(-1,j).f(3);
			//south east
			l.get_node(0,j).f(8)=l.get_node(l.nx,j).f(8);
			//south west
			l.get_node(l.nx-1,j).f(7)=l.get_node(-1,j).f(7);
			//north east
			l.get_node(0,j).f(5)=l.get_node(l.nx,j).f(5);
			//north west
			l.get_node(l.nx-1,j).f(6)=l.get_node(-1,j).f(6);
		}

		//corner adjustment for periodic boundary conditions

		//southeast corner(buffer) mapped to northwest corner of the domain
		l.get_node(0,l.ny-1).f(8)=l.get_node(l.nx,-1).f(8);

		//southwest corner(buffer) mapped to northeast corner of the domain
		l.get_node(l.nx-1,l.ny-1).f(7)=l.get_node(-1,-1).f(7);

		//northeast corner(buffer) mapped to southwest corner of the domain
		l.get_node(0,0).f(5)=l.get_node(l.nx,l.ny).f(5);

		//northwest corner(buffer) mapped to southeast corner of the domain
		l.get_node(l.nx-1,0).f(6)=l.get_node(-1,l.ny).f(6);
	}

	/** @brief Apply Boundary Conditions on top wall  */
	void top_wall_bc(){

		int j=l.ny - 1;
		bool flag_top_wall = false;

		//Wall, starts at i=1!
		if(flag_top_wall){
			#pragma omp parallel for
			for (int i=1; i<static_cast<int>(l.nx)-1; ++i)
			//iteration over the top buffers (filled due to up and down avection)
			{
				//south direction
				l.get_node(i,j).f(4)=l.get_node(i,j+1).f(2);
				//south east
				l.get_node(i,j).f(8)=l.get_node(i+1,j+1).f(5);
				//south west
				l.get_node(i,j).f(7)=l.get_node(i-1,j+1).f(6);
			}
		}
		//Free Slip
		else{
			#pragma omp parallel for schedule(static)
			for (int i=1; i<static_cast<int>(l.nx)-1; ++i)
			//iteration over the top buffers (filled due to up and down avection)
			{
				//IMMEDIATE BOUNCE BACK CONDITION
				//south direction
				l.get_node(i,j).f(4)=l.get_node(i,j).f(2);
				//south east
				l.get_node(i,j).f(8)=l.get_node(i,j).f(5);
				//south west
				l.get_node(i,j).f(7)=l.get_node(i,j).f(6);
			}
		}

	}

	/** @brief Apply Boundary Conditions on bottom wall  */
	void bottom_wall_bc(){

		int j=0;
		bool flag_bottom_wall = false;
		//Wall, starts at i = 1!
		if(flag_bottom_wall){
			#pragma omp parallel for
			for (int i=1; i<static_cast<int>(l.nx); ++i)
			//iteration over bottom buffers (filled due to up and down avection)
			{
				//south direction
				l.get_node(i,j).f(2)=l.get_node(i,j-1).f(4);
				//south east
				l.get_node(i,j).f(5)=l.get_node(i+1,j-1).f(8);
				//south west
				l.get_node(i,j).f(6)=l.get_node(i-1,j-1).f(7);
			}
		}
		//FreeSlip
		else{
			//IMMEDIATE BOUNCE BACK CONDITION
			#pragma omp parallel for schedule(static)
			for (int i=1; i<static_cast<int>(l.nx)-1; ++i)
			//iteration over bottom buffers (filled due to up and down advection)
			{
				//south direction
				l.get_node(i,j).f(2)=l.get_node(i,j).f(4);
				//south east
				l.get_node(i,j).f(5)=l.get_node(i,j).f(8);
				//south west
				l.get_node(i,j).f(6)=l.get_node(i,j).f(7);
			}
		}


	}

	/** @brief Apply Boundary Conditions on left wall  */
	void left_wall_bc(){
		//Inlet Conditions
		float_type u_x;

		if(time < runUptime)
			u_x = time * u_inlet / runUptime;
		else
			u_x = u_inlet;

		int i = 0;
		float_type uy = 0;
		#pragma omp parallel for schedule(dynamic)
		for(int j=0; j <static_cast<int>(l.ny) ; ++j){
			l.get_node(i,j).u()  = u_inlet;
			l.get_node(i,j).v()  = uy;
			l.get_node(i,j).rho() = rho_inlet;
			//initialize populations, rho initially 1 for doubly periodic shear layer
			for (unsigned int k=0; k<velocity_set().size; ++k)
			{
				l.get_node(i,j).f(k)=rho_inlet*velocity_set().W[k]*(2.-sqrt(1.+3.*u_x*u_x))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*u_x+sqrt(1.+3.*u_x*u_x))/(1.-u_x) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
			}
		}
	}

	/** @brief Apply Boundary Conditions on right wall  */
	void right_wall(){
		//No Boundary Condition
		const unsigned int i = l.nx -1;
		#pragma omp parallel for schedule(dynamic)
		for(int j = 1; j <static_cast<int>(l.ny)-1 ; ++j)
		{
			l.get_node(i,j).f(3) = l.get_node(i-1,j).f(3);
			l.get_node(i,j).f(6) = l.get_node(i-1,j/*+1*/).f(6);
			l.get_node(i,j).f(7) = l.get_node(i-1,j/*-1*/).f(7);
		}
	}

	/** @brief collide the populations */
	void collide()
	{
		//calculation rho,ux,uy at each lattice point then eqbm populations (for each element of velocity set)
		float_type ave_rho = 0;
		const int n_solid = l.solid_nodes.size();
		#pragma omp parallel for schedule(static)
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{	//calculation of rho, ux and uy for the node
				if( !(l.get_node(i,j).has_flag_property("solid")) )
				{
					double ux=0.,uy=0.,rho=0.,feq[9], alpha=2.;
					for (unsigned int temp=0; temp<velocity_set().size; ++temp)
						{
						rho+=l.get_node(i,j).f(temp);
						ux+=l.get_node(i,j).f(temp) * velocity_set().c[0][temp];
						uy+=l.get_node(i,j).f(temp) * velocity_set().c[1][temp];
						}
					ux=ux/rho;
					uy=uy/rho;
					l.get_node(i,j).rho() = rho;
					l.get_node(i,j).u()   = ux;
					l.get_node(i,j).v()   = uy;

					if (flag_moving_cyl)
					{
						#pragma omp atomic
						ave_rho += rho;
					}
					//#pragma omp parallel for
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						feq[k]=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
					}

					if (using_entropic)
					alpha = get_alpha(l.get_node(i,j), feq);
					//#pragma omp critical
					//std::cout << alpha << std::endl;

					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).f(k) += beta*alpha*(feq[k]-l.get_node(i,j).f(k));
					}

				}

			}
		}
		if (flag_moving_cyl)
		l.s_a_rho = ave_rho/(l.nx*l.ny-n_solid);
	}

	/** @brief Adaption of the Cylinder position  */
	void Adapt_Cyl()
	{
		//Move as as sinus oszillation
		float_type f = 0.0008;
		float_type y_move = Cyl_radius/2; // A: in Paper
		float_type omega = f*2*M_PI;
		Cyl_center[0] = Cyl_center_0[0];
		Cyl_center[1] = Cyl_center_0[1] + y_move * sin(omega*time);
		Cyl_vel[0] = 0.0;
		Cyl_vel[1] = y_move * omega * cos(omega*time);

		//std::cout << "Cylx: " << Cyl_center[0] << "Cyly: " << Cyl_center[1] << "Cyl v: " << Cyl_vel[1]<< std::endl;


		//Calculate solid equilibrium population with new properties
		//sim->l.f[0][sim->l.index(2,0)] = 3;
		float_type rho = l.s_a_rho,ux=Cyl_vel[0],uy=Cyl_vel[1];
		float_type f_solid[9];
		for (unsigned int k=0; k<velocity_set().size; ++k)
		{
			f_solid[k]=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
		}


		// Update with new Properties of all nodes(old position) inside the solid
		for(unsigned int k = 0 ; k< l.solid_nodes.size() ; ++k) {
			int i = l.solid_nodes[k].coord.i;
			int j = l.solid_nodes[k].coord.j;
			l.get_node(i,j).u()  = ux;
			l.get_node(i,j).v()  = uy;
			l.get_node(i,j).rho() = rho;
			for (unsigned int h=0; h<velocity_set().size; ++h)
				l.get_node(i,j).f(h) = f_solid[h];
		}


		//Delete the fluid boundary nodes from the step before
		l.delete_fluid_boundary_nodes();
		l.delete_solids();


		// Add new nodes of new position
		l.add_wallCylinder(Cyl_center, Cyl_vel, Cyl_radius, using_flagella, partition);


		// Update Properties of all nodes of new position inside the solid
		for(unsigned int k = 0 ; k< l.solid_nodes.size() ; ++k) {
			int i = l.solid_nodes[k].coord.i;
			int j = l.solid_nodes[k].coord.j;
			l.get_node(i,j).u()  = ux;
			l.get_node(i,j).v()  = uy;
			l.get_node(i,j).rho() = rho;
			for (unsigned int h=0; h<velocity_set().size; ++h)
				l.get_node(i,j).f(h) = f_solid[h];
		}

		//Steps to implement the boundary in every timestep in case of a moving solid

		//Calculate Force

		//Calculate Acceleration

		//Calculate and update Velocity

		//Update position of solid


		/*
		//Update Boundaries and node Properties
		l.delete_solids()
		l.delete_walls()

		l.add_wallCylinder(Cyl_center, Cyl_radius)
		*/
	}

	/** @brief Calculate force on the solid object in the flow */

	void eval_F_total (double &Fx_, double &Fy_) const {
		double Fx=0.,Fy=0.;
		#pragma omp parallel for schedule(dynamic)
		for (unsigned int it=0; it<l.fluid_boundary_nodes.size(); ++it) //summing for all fluid_boundary_nodes
		{
			unsigned int xb = l.fluid_boundary_nodes[it].coord.i; //coordinates of fluid boundary node
			unsigned int yb = l.fluid_boundary_nodes[it].coord.j;
			/*if(xb == 200){
				std::cout << "x = 200/ y =" << yb << " /// " << l.fluid_boundary_nodes.size() <<std::endl;
			}*/
			std::vector<int> D(l.get_node(xb,yb).missing_populations); //vector of intersecting populations, inverse of D-bar
			//std::vector<int> D(l.get_node(xb,yb).return_missing_populations());
			for (unsigned int i = 0; i < D.size(); ++i)
			{
				//unsigned int xs = xb + velocity_set().c[0][D[i]]; //coordinates of solid node
				//unsigned int ys = yb + velocity_set().c[1][D[i]];
				double common =  (l.get_node(xb,yb).f( inv_popl(D[i]) )
									//+l.get_node(xb,yb).s_di_populations[i] );
													+l.get_node(xb,yb).s_di_populations[i] );
				#pragma omp critical
				{
					Fx += velocity_set().c[0][D[i]] * common;
					Fy += velocity_set().c[1][D[i]] * common;
				}
			}
		}
		Fx_ = Fx; Fy_ =Fy;
	}
	/** @brief Calculate force on the solid object in the flow */

	void eval_F_general (const std::vector<node>& boundary_nodes, double &Fx_, double &Fy_) const {
		double Fx=0.,Fy=0.;
		#pragma omp parallel for schedule(dynamic)
		for (unsigned int it=0; it<boundary_nodes.size(); ++it) //summing for all fluid_boundary_nodes
		{
			unsigned int xb = boundary_nodes[it].coord.i; //coordinates of fluid boundary node
			unsigned int yb = boundary_nodes[it].coord.j;
			/*if(xb == 200){
				std::cout << "x = 200/ y =" << yb << " /// " << boundary_nodes.size() <<std::endl;
			}*/
			std::vector<int> D(boundary_nodes[it].missing_populations); //vector of intersecting populations, inverse of D-bar
			//std::vector<int> D(l.get_node(xb,yb).return_missing_populations());
			for (unsigned int i = 0; i < D.size(); ++i)
			{
				//unsigned int xs = xb + velocity_set().c[0][D[i]]; //coordinates of solid node
				//unsigned int ys = yb + velocity_set().c[1][D[i]];
				double common =  (l.get_node(xb,yb).f( inv_popl(D[i]) )
									//+l.get_node(xb,yb).s_di_populations[i] );
													+boundary_nodes[it].s_di_populations[i] );
				#pragma omp critical
				{
					Fx += velocity_set().c[0][D[i]] * common;
					Fy += velocity_set().c[1][D[i]] * common;
				}
			}
		}
		Fx_ = Fx; Fy_ =Fy;
	}

	std::vector<float_type> eval_M_flagella()
	{
		assert(using_flagella);
		unsigned int n_elements = flg->n;
		flg->updx(); //Update all x values inside the flagella class
		std::vector<float_type> Moments/*(n_elements)*/;


		for (unsigned int link_no = 0 ; link_no < n_elements ; ++link_no)
		{
			double moment = 0.;
			for (unsigned int it=0; it<l.flagella_nodes[link_no].size(); ++it)
			{
				double Fx=0.,Fy=0.;
				unsigned int xb = l.flagella_nodes[link_no][it].coord.i; //coordinates of fluid boundary node
				unsigned int yb = l.flagella_nodes[link_no][it].coord.j;
				std::vector<int> D(l.flagella_nodes[link_no][it].return_missing_populations()); //vector of intersecting populations, inverse of D-bar, only for flagella
				for (unsigned int i = 0; i < D.size(); ++i)
				{
					//unsigned int xs = xb + velocity_set().c[0][D[i]]; //coordinates of solid node
					//unsigned int ys = yb + velocity_set().c[1][D[i]];
					double common =  ( l.get_node(xb,yb).f( inv_popl(D[i]) )    //population of a fluid node from same time step
														+l.flagella_nodes[link_no][it].s_di_populations[i] );	//advected population from previous time step
														//+l.get_node(xs,ys).f( D[i] ) );
					//std::cout << "Common Forces: "  << common << " Populations: "<< inv_popl(D[i]) << "// " << l.get_node(xb,yb).f( inv_popl(D[i]) ) <<  std::endl;
					#pragma omp critical
					{
						Fx += velocity_set().c[0][D[i]] * common;
						Fy += velocity_set().c[1][D[i]] * common;
					}
				}
				moment += flg->eval_M(link_no,Fx,Fy,xb,yb);
			}
			Moments.push_back(moment);
		}
		return Moments;
	}

	void save_populations() {
		std::ofstream populations;
		populations.open("Populations.txt",std::ios::out);
		/* Format [f0] [f1] [f2] .......; where [f0] is column of population 0*/
		for (unsigned int i=0 ; i<l.real_size ; ++i )
		{
			for (unsigned int j=0 ; j<velocity_set().size ; ++j)
				populations << std::setw(15) << l.f[j][i];
			populations << "\n";
		}
		populations.close();
	}

	void resume(std::string filename) {
		std::ifstream filestream;
		filestream.open(filename,std::ios::in);
		for (unsigned int i=0 ; i<l.real_size ; ++i )
		{
			for (unsigned int j=0 ; j<velocity_set().size ; ++j)
			{
				double inp;
				filestream >> inp;
				l.f[j][i] = inp;
			}
		}
		filestream.close();
	}

	void resetBN(unsigned int n_links){
		for (auto it = l.cylinder_fbn_f.begin() ; it!=l.cylinder_fbn_f.end() ; ++it)
		{
			unsigned int i,j; i = it->coord.i; j = it->coord.j;
			l.get_node(i,j).missing_populations.clear();
			l.get_node(i,j).q_i.clear();
			l.get_node(i,j).uvw_i.clear();
			l.get_node(i,j).s_di_populations.clear();
			l.get_node(i,j).unset_flag_property("Fluid_Boundary_Node");
		}
		l.cylinder_fbn_f.clear();

		//if(using_flagella){
			for (unsigned int l_no = 0 ; l_no <  n_links ; ++l_no)
			{
				for (auto it = l.flagella_nodes[l_no].begin() ; it!=l.flagella_nodes[l_no].end() ; ++it)
				{
					unsigned int i, j; i = it->coord.i; j = it->coord.j;

					assert( (l.get_node(i,j).missing_populations.size() == l.get_node(i,j).q_i.size()) and  ( l.get_node(i,j).q_i.size()== l.get_node(i,j).uvw_i.size() ) and ( l.get_node(i,j).q_i.size() ==  l.get_node(i,j).s_di_populations.size() ) );
					//if(l.get_node(i,j).missing_populations.size()!=0){
						l.get_node(i,j).missing_populations.clear();
						l.get_node(i,j).q_i.clear();
						l.get_node(i,j).uvw_i.clear();
						l.get_node(i,j).s_di_populations.clear();
					//}
					l.get_node(i,j).unset_flag_property("Fluid_Boundary_Node");
				}
				l.flagella_nodes[l_no].clear();
			}
			l.flagella_nodes.clear();
		//}


		//probably the next loop is not required because cylinder_fbn is exclusive ie does not intersect with anyone else
		//only updation of the s_di_populations is necessary
		for (auto it = l.cylinder_fbn.begin() ; it!=l.cylinder_fbn.end() ; ++it)
		{
			unsigned int i, j, reset_size; 
			i = it->coord.i; 
			j = it->coord.j; 
			reset_size = it->missing_populations.size();
			unsigned int global_size = l.get_node(i,j).missing_populations.size();
			if(global_size != reset_size){
				std::cout << "Size difference in the: " << i << ", " << j << " // reset:" <<reset_size << " global: " << global_size << std::endl; 
			}
			//l.get_node(i,j).missing_populations.resize(reset_size);
			//l.get_node(i,j).q_i.resize(reset_size);
			//l.get_node(i,j).uvw_i.resize(reset_size);
			//l.get_node(i,j).s_di_populations.resize(reset_size);
			for (unsigned int temp = 0 ; temp < reset_size ; ++temp){
				l.get_node(i,j).s_di_populations[temp] = l.get_node(i,j).f(it->missing_populations[temp]);
				it->s_di_populations[temp] = l.get_node(i,j).f(it->missing_populations[temp]);
			}

			/*if ( !l.get_node(i,j).has_flag_property("Fluid_Boundary_Node") )
				l.get_node(i,j).set_flag_property("Fluid_Boundary_Node");*/
		}
		//std::cout << "changed the Nodes" << std::endl;

		l.fluid_boundary_nodes.resize(l.cylinder_fbn.size());
	}

	/** @brief Apply all Boundary Conditions */
	void Adapt_flagella(const std::vector<float_type>& Moments)
	{
		unsigned int n_links = Moments.size();
		resetBN(n_links);


		float_type dt = 1.0;
		std::cout << "Moments0:  " << Moments[0] << std::endl;
		//flg->step(Moments,dt);
		flg->writeAlphas();

		//std::cout << "Stepped the flagella" << std::endl;
		//adding the flagella nodes and corresponding fluid boundary nodes to the lattice
		l.add_flagella_nodes(flg,Cyl_vel, Cyl_center, Cyl_radius,partition);
		l.merge_into_fbn(using_flagella);

	}

	/** @brief Apply all Boundary Conditions */
	void wall_bc()
	{
		curved_wall_bc();
		//flagella_bc();
		top_wall_bc();
		bottom_wall_bc();
		left_wall_bc();
		right_wall();
		//periodic_bc();
	}

	void CalcForces(){
		//Calculate Forces

		/*std::tie(Fx_,Fy_) = */eval_F_total( Fx_, Fy_);
		if(using_flagella){
			std::vector<double> Fx_flag, Fy_flag;
			//for (auto it = l.flagella_nodes.begin();it !=l.flagella_nodes.end(); ++it){
				double FX, FY;
				/*std::tie(FX, FY) = */eval_F_general(l.flagella_nodes[0],FX, FY);
				Fx_flag.push_back(FX);
				Fy_flag.push_back(FY);
			//}
			Fx_flag_ = 0; Fy_flag_ = 0;
			for(auto& n: Fx_flag)
				Fx_flag_ += n;
			for(auto& n: Fy_flag)
				Fy_flag_ += n;
		}
		else{
			Fx_flag_ = 0;
			Fy_flag_ = 0;
		}
		
		double Fxtemp, Fytemp;
		eval_F_general(l.cylinder_fbn, Fxtemp, Fytemp);

		eval_F_general(l.cylinder_fbn_f, Fx_cyl_, Fy_cyl_);
		//std::tie(Fxtemp, Fytemp) = eval_F_general(l.cylinder_fbn);
		//std::tie(Fx_cyl_, Fy_cyl_) = eval_F_general(l.cylinder_fbn_f);

		std::cout << "Size of Fluid:" << l.fluid_boundary_nodes.size();
		if(using_flagella)
			std::cout << " // flagella : " << l.flagella_nodes[0].size();
		std::cout << " //Cylinder Fixed" << l.cylinder_fbn.size() << " Cylinder variable: " << l.cylinder_fbn_f.size() << std::endl;

		//std::cout << Fxtemp << "  " << Fx_cyl_ << "  //" << Fytemp << " " << Fy_cyl_ << std::endl;
		Fx_cyl_ += Fxtemp;
		Fy_cyl_ += Fytemp;
		std::cout << "Forcex difference: " << Fx_ - Fx_cyl_ - Fx_flag_ << "  // y-Force: " << Fy_ - Fy_cyl_ - Fy_flag_ << std::endl;
	}


	/** @brief LB step */
	void step()
	{
		std::cout << "*--------------------------------------*" << std::endl;
		std::cout << "Time; " << time << std::endl;
		//std::cout << "before advection" <<std::endl;
		advect();
		//std::cout << "before wall bc" <<std::endl;
		wall_bc();
		
		//std::cout << "before collision" <<std::endl;
		collide();

		CalcForces();


		//force << std::setw(10) << Fx_ << std::setw(10) << Fy_ << "\n";
		if(flag_moving_cyl)
			Adapt_Cyl();
		if(using_flagella)
		{
			//std::cout << "before Calculation of Moments for flagella" <<std::endl;
			std::vector<float_type> Moments(flg->n);
			Moments = eval_M_flagella();
			//std::cout << "Before adaption of flagella" << std::endl;X
			Adapt_flagella(Moments);
		}
		else{
			for (auto it = l.cylinder_fbn.begin() ; it!=l.cylinder_fbn.end() ; ++it)
			{
			unsigned int i, j, reset_size; 
			i = it->coord.i; 
			j = it->coord.j; 
			reset_size = it->missing_populations.size();
			unsigned int global_size = l.get_node(i,j).missing_populations.size();
			if(global_size != reset_size){
				std::cout << "Size difference in the: " << i << ", " << j << " // reset:" <<reset_size << " global: " << global_size << std::endl; 
			}
			//l.get_node(i,j).missing_populations.resize(reset_size);
			//l.get_node(i,j).q_i.resize(reset_size);
			//l.get_node(i,j).uvw_i.resize(reset_size);
			//l.get_node(i,j).s_di_populations.resize(reset_size);
			std::cout << l.get_node(i,j).s_di_populations.size() << "   " << reset_size << std::endl;
			for (unsigned int temp = 0 ; temp < reset_size ; ++temp){
				l.get_node(i,j).s_di_populations[temp] = l.get_node(i,j).f(it->missing_populations[temp]);
				it->s_di_populations[temp] = l.get_node(i,j).f(it->missing_populations[temp]);
			}

			/*if ( !l.get_node(i,j).has_flag_property("Fluid_Boundary_Node") )
				l.get_node(i,j).set_flag_property("Fluid_Boundary_Node");*/
		}
		}
		//std::cout << "step completed" << std::endl;


		// file io
		if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) )
		{
			write_fields();
			++output_index;
		}

		++time;
	}

public: // write to file

	/** write macroscopic variables to ascii file */
	void write_fields()
	{
		std::stringstream fns;
		fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
		l.write_fields(fns.str());
	}

public: // print

	/** print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const simulation& sim)
	{
		os << "simulation parameters\n"
		   << "---------------------\n";
		os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
		os << "Re:     " << sim.Re << "\n";
		os << "Vmax:   " << sim.Vmax << "\n";
		os << "visc:   " << sim.visc << "\n";
		os << "beta:   " << sim.beta << "\n";
		return os;
	}

public: // members

	lattice l;                 ///< lattice
	std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
	const float_type Re;       ///< Reynolds number
	const float_type Vmax;     ///< mean flow velocity
	const float_type Cyl_radius;
	const float_type visc;     ///< viscosity
	const float_type beta;     ///< LB parameter beta
	unsigned int time;         ///< simulation time
	bool file_output;          ///< flag whether to write files
	unsigned int output_freq;  ///< file output frequency
	unsigned int output_index; ///< index for file naming
	float_type Cyl_center[2];
	float_type Cyl_center_0[2];
	//float_type Cyl_radius;
	float_type Cyl_vel[2];
	float_type rho_inlet;
	float_type u_inlet;
	double Fx_,Fy_ , Fy_cyl_, Fx_cyl_, Fx_flag_, Fy_flag_;
	unsigned int runUptime;
	const bool flag_moving_cyl;
	const bool using_entropic;
	const bool using_flagella;
	flagella* flg;
	unsigned int partition;
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
