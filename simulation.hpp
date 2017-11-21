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
			// **************************
			// * filled in your code here *
			// **************************
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
		// **************************
		// * filled in your code here *
		// **************************
		/*const float_type pi(std::acos(-1.0));
		double ux,uy,rho;
//		double k=80.,delta=0.05; //for doubly periodic shear layer
		double Kx=2*pi/l.nx,Ky=2*pi/l.ny,K=sqrt(pow(Kx,2)+pow(Ky,2)),Ma=Vmax*sqrt(3);*/

		//Initialize Cyclinder
		Cyl_center[0] = 30.0;
		Cyl_center[1] = 30.0;
		Cyl_radius = 5.0;
		Cyl_vel[0] = 0.0;
		Cyl_vel[1] = 0.0;
		l.add_wallCylinder(Cyl_center, Cyl_radius);

		#pragma omp parallel for
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			const float_type pi(std::acos(-1.0));
			double ux,uy,rho;
//			double k=80.,delta=0.05; //for doubly periodic shear layer
			double Kx=2*pi/l.nx,Ky=2*pi/l.ny,K=sqrt(pow(Kx,2)+pow(Ky,2)),Ma=Vmax*sqrt(3);

			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{	//initializing velocities and density for doubly periodic shear layer	
				if(!l.get_node(i,j).has_flag_property("solid"))
				{
					/*ux=Vmax*( std::tanh(k*( (float)j/(float)l.ny -.25 )) );
					uy=Vmax*delta*( std::sin(2.*pi*( (float)i/(float)l.nx + .25 ) ) );
					l.get_node(i,j).u()  = ux;
					l.get_node(i,j).v()  = uy;
					l.get_node(i,j).rho() = 1.0;
					if (j>static_cast<int>(l.ny/2))
					{
						ux=Vmax*( std::tanh(-k*( (float)j/(float)l.ny-.75 )) );
						l.get_node(i,j).u()  = ux;
					}*/

					//initialization for 2d taylor green vortex flow

					ux=-Vmax*Ky/K*sin(Ky*j)*cos(Kx*i);
					uy=Vmax*Kx/K*sin(Kx*i)*cos(Ky*j);
					rho=1-pow(Ma/K,2)/2*(pow(Ky,2)*cos(2*Kx*i)+pow(Kx,2)*cos(2*Ky*j));
					l.get_node(i,j).u()  = ux;
					l.get_node(i,j).v()  = uy;
					l.get_node(i,j).rho() = rho;
					//initialize populations, rho initially 1 for doubly periodic shear layer
					for (unsigned int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).f(k)=rho*velocity_set().W[k]*(2.-sqrt(1.+3.*ux*ux))*(2.-sqrt(1.+3.*uy*uy))*pow((2.*ux+sqrt(1.+3.*ux*ux))/(1.-ux) ,velocity_set().c[0][k])*pow((2.*uy+sqrt(1.+3.*uy*uy))/(1.-uy) ,velocity_set().c[1][k]);
					}
				}
			}
		}
	}
	
	/** 
	 *  @brief advect the populations
	 *  
	 *  Include periodic boundary conditions here also
	 */
	void advect()
	{
		// **************************
		// * filled in your code here *
		// **************************
		unsigned int rtindex,lbindex;
		for (int j=0; j<static_cast<int>(l.ny); ++j)  
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)   
			{	
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
			}
		}

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
	
	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{
		#pragma omp parallel for
		for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
		{
			// **************************
			// * fill in your code here *
			// **************************

			bool dir_solid[9] = {false, false, false, false, false, false, false, false, false};


			unsigned int x_i = l.wall_nodes[i].coord.i;
			unsigned int y_j = l.wall_nodes[i].coord.j;

			//Get direction to solid according to the c-directions
			if (l.get_node(x_i+1,y_j).has_flag_property("solid") )
				dir_solid[1] = true;
			if (l.get_node(x_i,y_j+1).has_flag_property("solid") )
				dir_solid[2] = true;
			if ( l.get_node(x_i-1,y_j).has_flag_property("solid") )
				dir_solid[3] = true;
			if ( l.get_node(x_i,y_j-1).has_flag_property("solid") )
				dir_solid[4] = true;
			if ( l.get_node(x_i+1,y_j+1).has_flag_property("solid") )
				dir_solid[5] = true;
			if ( l.get_node(x_i-1,y_j+1).has_flag_property("solid") )
				dir_solid[6] = true;
			if ( l.get_node(x_i-1,y_j-1).has_flag_property("solid") )
				dir_solid[7] = true;
			if ( l.get_node(x_i+1,y_j-1).has_flag_property("solid") )
				dir_solid[8] = true;

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
				if(l.get_node(i,j).has_flag_property("fluid"))
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
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
