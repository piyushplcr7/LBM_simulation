
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>



int main(int argc, char *argv[])
{
	//omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));
	omp_set_num_threads(1);
	unsigned int nx = 600, ny = 600, Time = 100;
	float Vmax = 0.05, Re = 500;
	lb::simulation* sim = new lb::simulation(nx,ny,Re,Vmax);
	sim->initialize();
	std::cout << *sim << std::endl;

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

		std::ofstream force;
		force.open("Force.txt",std::ios::out);
		force << std::setw(15) << "#Fx:" << std::setw(15) << "Fy:" << "\n";
		// use a loop like this to run the simulation

		for (unsigned int i=0; i<Time; ++i)
		{
			sim->step();
			force << std::setw(15) << sim->Fx_ << std::setw(15) << sim->Fy_ << "\n";
		}
		force.close();
	#endif

	return 0;
}
