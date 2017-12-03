/**
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>
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
	  visc( /*filled in your code here*/ /*0.0015*/ Vmax*nx/Re),
	  beta( /*filled in your code here*/ /*0.9911*/ 1./(6*visc+1)),
	  time(0),
	  file_output(false), // set to true if you want to write files
	  output_freq(100),
	  output_index(0)
	{
		// define amount to shift populations for advection (according to the array model of domain)
		for (unsigned int i=0; i<velocity_set().size; ++i)
		{
			shift[i]=velocity_set().c[0][i] + l.real_nx * velocity_set().c[1][i];
		}
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
		Cyl_center[0] = 75.0;
		Cyl_center[1] = 75.0;
		Cyl_radius = 30.0;
		Cyl_vel[0] = 0.0;
		Cyl_vel[1] = 0.0;
		l.add_wallCylinder(Cyl_center, Cyl_radius);

		//Init B.C.
		u_inlet = 0.01;
		rho_inlet = 1;
		const float_type pi(std::acos(-1.0));

		//#pragma omp parallel for
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			double ux,uy,rho;
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				/*if( !(l.get_node(i,j).has_flag_property("solid")) )
				{
					//Initialize flow around Cylinder
					ux = 0.0;
					uy = 0.0;
					rho = 1;
					l.get_node(i,j).u()  = ux;
					l.get_node(i,j).v()  = uy;
					l.get_node(i,j).rho() = rho;
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						//l.get_node(i,j).f(k)=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
						l.get_node(i,j).f(k)=velocity_set().W[k]; //for ux=0 uy=0 rho=1
					}
				}
				else
				{
					//for inside the solid, initialize with negative values for populations
					ux = 0.;
					uy = 0.;
					rho = 10;
					l.get_node(i,j).u()  = ux;
					l.get_node(i,j).v()  = uy;
					l.get_node(i,j).rho() = rho;
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						//l.get_node(i,j).f(k)=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
						l.get_node(i,j).f(k)=-velocity_set().W[k]; //for ux=0 uy=0 rho=1
					}
				}*/

				//Initialize flow around Cylinder
				ux = 0.0;
				uy = 0.0;
				rho = 1;
				l.get_node(i,j).u()  = ux;
				l.get_node(i,j).v()  = uy;
				l.get_node(i,j).rho() = rho;
				for (unsigned int k=0; k<velocity_set().size; ++k)
				{
					l.get_node(i,j).f(k)=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
					//l.get_node(i,j).f(k)=velocity_set().W[k]; //for ux=0 uy=0 rho=1
				}
				
			}
		} //loop j ends

	}

	/**
	 *  @brief advect the populations
	 *
	 *  Include periodic boundary conditions here also
	 */
	void advect() //advection to be done
	{
		unsigned int rtindex,lbindex;
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				//advect everything, even the solid nodes, if bcs are not implemented correctly, if will be seen
				//if( !(l.get_node(i,j).has_flag_property("solid")) ) //advect normally if not a solid node (causes filling of the populations in the solid node)
				//{
					//right top index or increasing index
					rtindex= l.real_nx*(j+l.buffer_size) + i + l.buffer_size;
					//left bottom index or decreasing index
					lbindex= l.real_nx*(l.ny-1-j+l.buffer_size) + l.nx-1 -i + l.buffer_size;
					//PUSH ADVECTION!
					for (unsigned int k=0; k<velocity_set().size; ++k)
						{
						if (shift[k]>0) //if shift is positive, using decreasing index
							l.f[k][lbindex+shift[k]]=l.f[k][lbindex];
						else
							l.f[k][rtindex+shift[k]]=l.f[k][rtindex];
						}
				//}
			}
		}


	}

	/**  @brief calculate smaller root of quadratic function */
	inline float_type solve_quadratic(float_type a, float_type b, float_type c){
		return -b/2/a - sqrt(b*b-4*a*c)/2/a;
	}

	/**  @brief calculate fraction of fluid from node n to the solid boundary along the direction i */
	float_type get_qi(const node& n,const int& i){
		double x1 = (double)n.coord.i; double y1 = (double)n.coord.j;
		double x2 =  x1 + (double)velocity_set().c[0][i]; double y2 =  y1 + (double)velocity_set().c[1][i];
		//calc distance c_i
		float_type c_i = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

		float_type dx = x2-x1;
		float_type dy = y2-y1;

		//calc x_w,i
		//find root of quadratic equation
		//((x1+t*dx -l.center[0])*(x1+t*dx -l.center[0]) + (y1+t*dy -l.center[1])*(y1+t*dy -l.center[1])-l.radius*l.radius == 0)
		float_type a = (dx*dx + dy*dy);
		float_type b = dx* (x1-Cyl_center[0]) + dy * (y1 - Cyl_center[1]);
		float_type c = (x1-Cyl_center[0])*(x1-Cyl_center[0]) + (y1-Cyl_center[1])*(y1-Cyl_center[1]) - Cyl_radius*Cyl_radius;
		float_type t = solve_quadratic(a,b,c);


		float_type q_i = sqrt(t*dx *t*dx + t*dy *t*dy)/c_i;
		std::cout << "Q_i: " << q_i << std::endl; 

		return q_i;

	}
	/**  @brief Calculate Eq Presure tensor */
	void calc_Peq(int i, int j, float_type (&Peq)[2][2]){

		//Set non equilibrium pressure tensor
		Peq[0][0] = l.get_node(i,j).rho()*(velocity_set().cs*velocity_set().cs + l.get_node(i,j).u()* l.get_node(i,j).u()) ;
		Peq[1][0] = l.get_node(i,j).rho()*( l.get_node(i,j).u()* l.get_node(i,j).v());
		Peq[0][1] = Peq[1][0];
		Peq[1][1] = l.get_node(i,j).rho()*(velocity_set().cs*velocity_set().cs + l.get_node(i,j).v()* l.get_node(i,j).v());
		return;

	}


	/**  @brief Calculate NonEq for with direct influence of Boundary */
	void calc_Pneq(int i, int j/*, bool dir_solid[9]*/, float_type q_i[9], float_type (&Pneq)[2][2]){

		//calculate velocity gradients
		float_type du_x = (l.get_node(i+1,j).u()-l.get_node(i-1,j).u())/(q_i[3] + q_i[1]);
		float_type dv_x = (l.get_node(i+1,j).v()-l.get_node(i-1,j).v())/(q_i[3] + q_i[1]);

		float_type du_y = (l.get_node(i+1,j).u()-l.get_node(i-1,j).u())/(q_i[4] + q_i[2]);
		float_type dv_y = (l.get_node(i,j+1).v()-l.get_node(i,j-1).v())/(q_i[4] + q_i[2]);


		//Set non equilibrium pressure tensor
		Pneq[0][0] = - l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/beta*(du_x);
		Pneq[1][0] = -l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/beta/2*(du_y + dv_x);
		Pneq[0][1] = Pneq[1][0];
		Pneq[1][1] = - l.get_node(i,j).rho()*velocity_set().cs*velocity_set().cs/beta*(dv_y);

		return;
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
	inline int inv_popl(const int& i) {
		int inverse[9] = {0,3,4,1,2,7,8,5,6};
		return inverse[i];
	}

	/**  @brief apply wall boundary conditions */
	void curved_wall_bc()
	{
		//Curved Wall
		//#pragma omp parallel for
		std::cout << "Fluid boundary node: " <<  l.fluid_boundary_nodes.size() << std::endl;
		for (unsigned int i=0; i<l.fluid_boundary_nodes.size(); ++i)
		{
			//bool dir_solid[9] = {false, false, false, false, false, false, false, false, false};
			unsigned int x_i = l.fluid_boundary_nodes[i].coord.i;
			std::cout << x_i << std::endl;
			unsigned int y_j = l.fluid_boundary_nodes[i].coord.j;
			//the directions are already stored within the node
			/*/Get direction to solid according to the c-directions
			for( int j = 1; j < sizeof(dir_solid); ++j){
				if (l.get_node(x_i+velocity_set().c[0][j],y_j+velocity_set().c[1][j]).has_flag_property("solid") )
					dir_solid[j] = true;
			} */

			//Get fraction of fluid between current node and the adjacent node (1 if adjacent node is on or outside the solid)
			/*std::vector<int> i_d;  //the information recorded by this i_d is contained in missing_populations vector
			for (int j=0; j < sizeof(dir_solid); ++j)
			{
				if (dir_solid[j]){
					i_d.push_back(j);
				}
				// NO ADVECTION DURING BC!!!
				//Advect the direction away from Solid
				else{
					l.get_node(x_i+velocity_set().c[0][j], y_j + velocity_set().c[1][j]).f(j) = l.get_node(x_i, y_j).f(j);
				}
			}*/
			std::cout << l.get_node(x_i,y_j).missing_populations.size() << std::endl;
			std::vector<int> i_d(l.fluid_boundary_nodes[i].missing_populations); //contains the direction in which solid is encountered
			//utgt
			float_type utgt[2] = {0,0};
			float_type rho_bb = 0;
			float_type rho_s = 0;
			float_type rho_0 = 1;
			float_type q_i[9] = {1,1,1,1,1,1,1,1,1};

			std::cout << i_d.size() << std::endl;


			for (int j = 0; j < i_d.size(); ++j)
			{
				//Velocities away from the solid (for the missing populations D bar )
				int c_x = -velocity_set().c[0][i_d[j]];
				int c_y = -velocity_set().c[1][i_d[j]];

				//get q_i for the direction given by i_d[j]
				//q_i[i_d[j]] = get_qi(x_i, y_j, x_i - c_x, y_j - c_y);
				q_i[i_d[j]] = get_qi( l.fluid_boundary_nodes[i], i_d[j] );
				std::cout << j << " "<< std::endl;

				//get u_f,i
				float_type u_fi[2] = { l.get_node(x_i+c_x,y_j+c_y).u(), l.get_node(x_i+c_x,y_j+c_y).v()};

				//calc utgt part
				utgt[0] += (q_i[i_d[j]]* u_fi[0] + Cyl_vel[0])/(1+q_i[i_d[j]])/i_d.size();
				utgt[1] += (q_i[i_d[j]]* u_fi[1] + Cyl_vel[1])/(1+q_i[i_d[j]])/i_d.size();

				float_type u_w_i = Cyl_vel[0];
				float_type v_w_i = Cyl_vel[1];

				//get rho_s
				rho_s += 6*rho_0*velocity_set().W[i_d[j]]*(c_x*u_w_i + c_y*v_w_i);
				//adjusting the bounce-back populations (the D-bar populations, which are inverse of the missing_populations)
				l.fluid_boundary_nodes[i].f( inv_popl(i_d[j]) ) = l.get_node(x_i, y_j).f(i_d[j]);

			}

			for(int j = 0; j < velocity_set().size ; ++j){
				//rho_bb += l.get_node(x_i, y_j).f(j);
				rho_bb += l.fluid_boundary_nodes[i].f(j);
			}

			//rho_tgt
			float_type rho_tgt = rho_bb + rho_s;

			//Calculate missing populations
			float_type Peq[2][2];
			calc_Peq(x_i, y_j, Peq);
			float_type Pneq[2][2];
			calc_Pneq(x_i, y_j, /*dir_solid,*/ q_i, Pneq);
			float_type f_new;

			std::cout << "Peq: " << Peq[0][0] << " /" << Peq[0][1] << " /" << Peq[1][0] << " /" << Peq[1][1] << " /" << std::endl;
			std::cout << "Pneq: " << Pneq[0][0] << " /" << Pneq[0][1] << " /" << Pneq[1][0] << " /" << Pneq[1][1] << " /" << std::endl;
			std::cout << "q_i: " << q_i[0] << " /" << q_i[1] << " /" << q_i[2] << " /" << q_i[3] << " /" << q_i[4] << " /" << q_i[5] << " /" << q_i[6] << " /" << q_i[7] << " /" << q_i[8] << std::endl;

			for (int j = 0; j < i_d.size(); ++j)
			{
				//Reset population with first two parts of grads' approximation
				//f_new = velocity_set().W[i_d[j]]*(l.get_node(x_i, y_j).rho()*(1 + (velocity_set().c[0][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).u() + velocity_set().c[1][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).v())/velocity_set().cs/velocity_set().cs));
				f_new = (l.get_node(x_i, y_j).rho()*(1 + (velocity_set().c[0][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).u() + velocity_set().c[1][ inv_popl(i_d[j]) ]*l.get_node(x_i, y_j).v())/velocity_set().cs/velocity_set().cs));

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
				if (!std::isnan(f_new))
					std::cout << "F_ new_ missing: " << f_new << std::endl;
			}
		}
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
			//#pragma omp parallel for
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
			//#pragma omp parallel for
			for (int i=1; i<static_cast<int>(l.nx); ++i)
			//iteration over bottom buffers (filled due to up and down avection)
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
		unsigned int runUptime = 100;
		float_type u_x;

		if(time < runUptime)
			u_x = time * u_inlet / runUptime;
		else
			u_x = u_inlet;

		int i = 0;
		float_type uy = 0;
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
		unsigned int i = l.nx -1;

		for(int j = 1; j <static_cast<int>(l.ny)-1 ; ++j){
			l.get_node(i,j).f(3) = l.get_node(i-1,j).f(3);
			l.get_node(i,j).f(6) = l.get_node(i-1,j+1).f(6);
			l.get_node(i,j).f(7) = l.get_node(i-1,j-1).f(7);

		}
	}

	/** @brief collide the populations */
	void collide()
	{
		// **************************
		// * fill in your code here *
		// **************************

		//calculation rho,ux,uy at each lattice point then eqbm populations (for each element of velocity set)
		double ux,uy,rho,feq;

		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{	//calculation of rho, ux and uy for the node
				if( !(l.get_node(i,j).has_flag_property("solid")) )
				{
					rho=0.,ux=0.,uy=0.;
					for (unsigned int temp=0; temp<velocity_set().size; ++temp)
						{
						rho+=l.get_node(i,j).f(temp);
						ux+=l.get_node(i,j).f(temp) * velocity_set().c[0][temp];
						uy+=l.get_node(i,j).f(temp) * velocity_set().c[1][temp];
						}
					ux=ux/rho;
					uy=uy/rho;
					l.get_node(i,j).rho()=rho;
					l.get_node(i,j).u()   = ux;
					l.get_node(i,j).v()   = uy;

					//collide populations
					#pragma omp parallel for /*num_threads(8)*/
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						feq=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
						l.get_node(i,j).f(k)+=2.*beta*(feq-l.get_node(i,j).f(k));
					}
				}
			}
		}
	}

	/** @brief Adaption of the Cylinder position  */
	void Adapt_Cyl()
	{
		// **************************
		// * fill in your code here *
		// **************************

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
	/** @brief Apply all Boundary Conditions */
	void wall_bc()
	{
		curved_wall_bc();
		top_wall_bc();
		bottom_wall_bc();
		left_wall_bc();
		right_wall();

		//periodic_bc();
	}

	/** @brief LB step */
	void step()
	{
		advect();
		wall_bc();
		collide();
		Adapt_Cyl();

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
	const float_type visc;     ///< viscosity
	const float_type beta;     ///< LB parameter beta
	unsigned int time;         ///< simulation time
	bool file_output;          ///< flag whether to write files
	unsigned int output_freq;  ///< file output frequency
	unsigned int output_index; ///< index for file naming
	float_type Cyl_center[2];
	float_type Cyl_radius;
	float_type Cyl_vel[2];
	float_type rho_inlet;
	float_type u_inlet;
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
