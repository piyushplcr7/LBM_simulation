//#include "lattice.hpp"
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>
#include <sys/stat.h>



int main(int argc, char *argv[])
{
	//omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));
	//omp_set_num_threads(1);
	unsigned int nx = 400, ny = 400, Time = 50000;
	float Vmax = 0.05, Re = 100;
	lb::simulation* sim = new lb::simulation(nx,ny,Re,Vmax);
	sim->initialize();
	std::cout << *sim << std::endl;
	std::cout <<"Reynoldnumber: " << Re;
	//sim->resume("Populations.txt");

	#ifdef USE_OPENGL_VISUALIZATION

		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();

	#else

		// Here are some hints for getting aquainted with the lattice class
		// ================================================================

		// how to print the lattice:
		// -------------------------

		//std::cout << sim->l << std::endl;

		// how to access the lattice:
		// --------------------------

		// 1) access via node proxy
		//sim->l.get_node(1,0).f(0) = 2;

		// 2) access data directly (make sure you know what you're doing)
		//sim->l.f[0][sim->l.index(2,0)] = 3;

		// 3) using iterators to nodes
		//(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;



		//std::cout << sim->l << std::endl;
		std::string file_name_basic = "Result/Fields/out";
		mkdir("Result2", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		mkdir("Result2/Fields", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


		std::ofstream force,v1,v2,v3,v4;
		
		force.open("Force.txt",std::ios::out);
		v1.open("v_400_200.txt",std::ios::out);
		v2.open("v_400_250.txt",std::ios::out);
		v3.open("v_600_200.txt",std::ios::out);
		v4.open("v_600_250.txt",std::ios::out);
		force << std::setw(15) << "#Fx:" << std::setw(15) << "Fy:" << "\n";
		v1 << std::setw(15) << "#ux:" << std::setw(15) << "uy:" << "\n";
		v2 << std::setw(15) << "#ux:" << std::setw(15) << "uy:" << "\n";
		v3 << std::setw(15) << "#ux:" << std::setw(15) << "uy:" << "\n";
		v4 << std::setw(15) << "#ux:" << std::setw(15) << "uy:" << "\n";
		// use a loop like this to run the simulation

		for (unsigned int i=0; i<Time; ++i)
		{
			
			//if(i==1) {sim->l.write_fields(file_name);}			
			sim->step();
			force << std::setw(15) << sim->Fx_ << std::setw(15) << sim->Fy_ << "\n";
			v1 << std::setw(15) << sim->l.get_node(400,200).u() << std::setw(15) << sim->l.get_node(400,200).v() << "\n";
			v2 << std::setw(15) << sim->l.get_node(400,250).u() << std::setw(15) << sim->l.get_node(400,250).v() << "\n";
			v3 << std::setw(15) << sim->l.get_node(600,200).u() << std::setw(15) << sim->l.get_node(600,200).v() << "\n";
			v4 << std::setw(15) << sim->l.get_node(600,250).u() << std::setw(15) << sim->l.get_node(600,250).v() << "\n";
			//Write to files
			if(i > 30000 and i < 35000){
				std::string file_name = file_name_basic + "_" + std::to_string(i) + ".txt";			
				sim->l.write_fields(file_name);
			}

			//if (i == 14000)
				//sim->save_populations();
		}
		force.close();
		v1.close();
		v2.close();
		v3.close();
		v4.close();
	#endif

	return 0;
}
