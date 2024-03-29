/**
 *  @file
 *  @author Fabian Bösch
 *  @brief lattice and node
 */

#ifndef LB_LATTICE_HPP_INCLUDED
#define LB_LATTICE_HPP_INCLUDED

#include "velocity_set.hpp"
//#include "distance.hpp"
#include "property_array.hpp"
#include <vector>
#include <fstream>
#include<cassert>

namespace lb {

class lattice; // forward declaration

/**
 *  @brief Node representing one lattice site.
 *
 *  Easy access to bundled quantities and properties (works as proxy to
 *  the lattice class).
 */
struct node
{
public: // ctors

	/** @brief Default constructor */
	node() {}

	/**
	 *  @brief Construct from lattice and position
	 *  @param[in] lat Pointer to the lattice
	 *  @param[in] i   x coordinate
	 *  @param[in] j   y coordinate
	 *  @pre coordinates are in domain
	 */
	node(lattice* lat, int i, int j);

	node(const node&) = default;

public: // init

	/**
	 *  @brief Set lattice and position.
	 *  @param[in] lat Pointer to lattice
	 *  @param[in] i   x coordinate
	 *  @param[in] j   y coordinate
	 *  @pre coordinates are in domain
	 */
	void set(lattice* lat, int i, int j);

public: // access populations and macroscopic quantities

	/**
	 *  @brief Get population.
	 *  @param i Population index
	 *  @return Value of distribution function
	 *  @pre population index exists
	 */
	inline float_type f(unsigned int i) const;

	/**
	 *  @brief Get/set population.
	 *  @param i Population index
	 *  @return Reference to value of distribution function
	 *  @pre population index exists
	 */
	inline float_type& f(unsigned int i);

	/**
	 *  @brief Get density.
	 *  @return Local density
	 */
	inline float_type rho() const;

	/**
	 *  @brief Get/set density.
	 *  @return Reference to local density
	 */
	inline float_type& rho();

	/**
	 *  @brief Get x-velocity.
	 *  @return Local flow velocity in x direction
	 */
	inline float_type u() const;

	/**
	 *  @brief Get/set x-velocity.
	 *  @return Reference to local flow velocity in x direction
	 */
	inline float_type& u();

	/**
	 *  @brief Get y-velocity.
	 *  @return Local flow velocity in y direction
	 */
	inline float_type v() const;

	/**
	 *  @brief Get/set y-velocity.
	 *  @return Reference to local flow velocity in y direction
	 */
	inline float_type& v();

public: // query and access properties

	/**
	 *  @brief Query for flag property.
	 *  Query whether a flag is set for the node.
	 *  @param[in] name Flag name
	 *  @return True if flag is set, otherwise false
	 */
	inline bool has_flag_property(std::string name) const;

	/**
	 *  @brief Set a flag
	 *  Set the flag "name" to true
	 *  @param[in] name Flag name
	 *  @return True if flag exists, otherwise false
	 */
	inline bool set_flag_property(std::string name);

	/**
	 *  @brief Unset a flag
	 *  Set the flag "name" to false
	 *  @param[in] name Flag name
	 *  @return True if flag exists, otherwise false
	 */
	inline bool unset_flag_property(std::string name);

	/**
	 *  @brief Query for data property.
	 *  Query whether data property (object) is stored for the node.
	 *  @param[in] name Data property name
	 *  @return True if thre is such a data property, otherwise false
	 */
	inline bool has_data_property(std::string name) const;

	/**
	 *  @brief Store a data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @param[in] property Data property object
	 *  @return True if data property exists, otherwise false
	 */
	template<typename T>
	inline bool set_data_property(std::string name, const T& property);

	/**
	 *  @brief Delete a data property
	 *  @param[in] name Data property name
	 *  @return True if data property exists, otherwise false
	 */
	bool unset_data_property(std::string name);

	/**
	 *  @brief Get data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @return Reference to data property object
	 */
	template <typename T>
	T& get_data_property(std::string name);

	/**
	 *  @brief Get data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @return Reference to data property object
	 */
	template <typename T>
	const T& get_data_property(std::string name) const;
	/**
	*		@Additionaly implemented function
	*		@Adds to the missing_populations vector (only for the Fluid Boundary Nodes)
	*/
	void add_missing_populations(const int&);
	/**
	*		@Additionaly implemented function
	*		@Returns copy of the vector of missing populations
	*/
	std::vector<int> return_missing_populations() const;
	/**
	*		@Additionaly implemented function
	*		@Clears the missing_populations vector
	*/
	void clear_missing_populations();
public: // members

	lattice* l;            ///< Pointer to a lattice object
	unsigned int index;    ///< Index for looking up data in the lattice
	coordinate<int> coord; ///< Coordinate of node's position
	std::vector<int> missing_populations; ///<Contains the indices of missing populations (for Fluid Boundary Nodes, empty otherwise)
	std::vector<float_type> q_i; ///<Contains the indices of missing populations (for Fluid Boundary Nodes, empty otherwise)
};

/**
 *  @brief Lattice containing the populations.
 *
 *  The lattice is constructed using the function @ref velocity_set()
 *  which returns a velocity set object. Hence, the number of
 *  populations is defined through that function. Data structures are
 *  set up accordingly.
 *
 *  The basic data structure for the population and the macroscopic
 *  qunatities are one dimensional arrays (vectors) interpreted as two
 *  dimensional planes. The x (i) dimension varies first and the y (j)
 *  dimension last.
 *
 *  This class does provide access to the data through node iterators or
 *  through direct access of the public members. The node iterators
 *  return a @ref node object that provides easy access to all local
 *  quantities according to the 2d lattice coordinate.
 *
 *  There are buffer regions (extent is one in all directions) around
 *  the data to make the advection procedure easier.
 *
 *  The data is indexed in the range [0, nx-1][0, ny-1]; including
 *  buffers indices span the range [-1, nx][-1, ny], repectively.
 */
class lattice
{
public: // typedefs

	/** @brief Iterator type */
	typedef typename std::vector<node>::iterator node_iterator;
	/** @brief Const iterator type */
	typedef typename std::vector<node>::const_iterator const_node_iterator;
	/** @brief Reverse iterator type */
	typedef typename std::vector<node>::reverse_iterator reverse_node_iterator;
	/** @brief Const reverse iterator type */
	typedef typename std::vector<node>::const_reverse_iterator const_reverse_node_iterator;

public: // ctor

	/**
	 *  @brief Construct the lattice with given extent
	 *  @param[in] _nx Number of nodes in x direction
	 *  @param[in] _ny Number of nodes in y direction
	 */
	lattice(unsigned int _nx, unsigned int _ny);

public: // coordinates to index conversion

	/**
	 *  @brief Convert a coordinate to a unique index
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return unique index
	 *  @pre Coordinates are in the domain
	 */
	inline unsigned int index(int i, int j) const;

public: // node access

	/** @brief Iterator pointing to the beginning  @return iterator */
	node_iterator begin();
	/** @brief Const iterator pointing to the beginning @return const iterator */
	const_node_iterator begin() const;
	/** @brief Iterator pointing to the end @return iterator */
	node_iterator end();
	/** @brief Const iterator pointing to the end @return const iterator */
	const_node_iterator end() const;
	/** @brief Reverse iterator pointing to the end @return reverse iterator */
	reverse_node_iterator rbegin();
	/** @brief Const reverse iterator pointing to the end @return const reverse iterator */
	const_reverse_node_iterator rbegin() const;
	/** @brief Reverse iterator pointing to the beginning @return reverse iterator */
	reverse_node_iterator rend();
	/** @brief Const reverse iterator pointing to the beginning @return const reverse iterator */
	const_reverse_node_iterator rend() const;

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return reference to node at coordinate (i,j)
	 *  @pre coordinates are in domain
	 */
	inline node& get_node(int i, int j);

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return const reference to node at coordinate (i,j)
	 *  @pre coordinates are in domain
	 */
	inline const node& get_node(int i, int j) const;

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] idx unique node index
	 *  @return reference to node at coordinate (i,j)
	 *  @pre idx is between [0, @ref lattice::real_size )
	 */
	inline node& get_node(unsigned int idx);

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] idx unique node index
	 *  @return const reference to node at coordinate (i,j)
	 *  @pre idx is between [0, @ref lattice::real_size )
	 */
	inline const node& get_node(unsigned int idx) const;

public: // walls

	/**
	 *  @brief Add a solid wall
	 *
	 *  Creates wall flags in the coordinate rectangle defined by
	 *  min_coord and max_coord. The corresponding nodes get the flag
	 *  "wall" and they are also stored in the vector
	 *  @ref lattice:wall_nodes for convienience.
	 *
	 *  @param[in] min_coord minimum bounding rectangle corner
	 *  @param[in] max_coord maximum bounding rectangle corner
	 *  @pre (min_coord, max_coord) define a rectangle
	 *  @pre Both min_coord and max_coord are in the domain
	 */
	void add_wall(coordinate<int> min_coord, coordinate<int> max_coord);

	void add_wallCylinder(float_type center[2], float_type radius);

	/** @brief Delete all existing walls */
	void delete_walls();

	/** @brief Delete all existing solids */
	void delete_solids();
	/** @brief Delete all Fluid boundary nodes*/
	void delete_fluid_boundary_nodes();

public: // file dump

	/**
	 *  @brief Write fields to file
	 *
	 *  Write macroscopic variables to simple ascii file.
	 *
	 *  @param[in] file_name file name
	 */
	void write_fields(std::string file_name);

public: // print

	/** @brief print to output stream, useful for debugging only */
	friend std::ostream& operator<<(std::ostream& os, const lattice& l);

public: // members

	const unsigned int nx;                    ///< extent in x direction (excluding buffers)
	const unsigned int ny;                    ///< extent in y direction (excluding buffers)
	const unsigned int size;                  ///< total number of nodes (excluding buffers)
	const unsigned int buffer_size;           ///< buffer width (equal to one)
	const unsigned int real_nx;               ///< extent in x direction including buffers
	const unsigned int real_ny;               ///< extent in y direction including buffers
	const unsigned int real_size;             ///< total number of nodes including buffers
	const unsigned int n_populations;         ///< number of populations
	float_type s_a_rho;
	std::vector<std::vector<float_type> > f;  ///< population data
	std::vector<float_type> rho;              ///< density data
	std::vector<float_type> u;                ///< flow x-velocity data
	std::vector<float_type> v;                ///< flow y-velocity data
	std::vector<node> nodes;                  ///< array holding all node objects
	std::vector<node> wall_nodes;             ///< array holding node objects belonging to a solid wall
	std::vector<node> solid_nodes;
	std::vector<node> fluid_boundary_nodes;
	property_array properties;                ///< properties datastructure (can hold many different properties per node)
	const bool periodic_x;                    ///< flag whether to use periodicity in x direction
	const bool periodic_y;                    ///< flag whether to use periodicity in y direction
	std::vector<std::vector<int>> miss_pop;   ///< Vector of missing populations for each node
};



// implementation
// --------------


// node

node::node(lattice* lat, int i, int j)
: l(lat), index(l->real_nx*(j+l->buffer_size) + i + l->buffer_size), coord(i,j) { }

void node::set(lattice* lat, int i, int j)
{
	l = lat;
	index = l->real_nx*(j+l->buffer_size) + i + l->buffer_size;
	coord.i = i;
	coord.j = j;
}

inline float_type node::f(unsigned int i) const { return l->f[i][index]; }
inline float_type& node::f(unsigned int i) { return l->f[i][index]; }
inline float_type node::rho() const { return l->rho[index]; }
inline float_type& node::rho() { return l->rho[index]; }
inline float_type node::u() const { return l->u[index]; }
inline float_type& node::u() { return l->u[index]; }
inline float_type node::v() const { return l->v[index]; }
inline float_type& node::v() { return l->v[index]; }
//inline std::vector<int> node::missing_pop(const unsigned int index_fluid_node) const { return l->miss_pop[index_fluid_node]; }
//inline std::vector<int>& node::missing_pop(const unsigned int index_fluid_node){return l->miss_pop[index_fluid_node]; }


inline bool node::has_flag_property(std::string name) const { return l->properties.has_flag_property(name, index); }
inline bool node::set_flag_property(std::string name) { return l->properties.set_flag_property(name, index); }
inline bool node::unset_flag_property(std::string name) { return l->properties.unset_flag_property(name, index); }
inline bool node::has_data_property(std::string name) const { return l->properties.has_data_property(name,index); }
template<typename T>
inline bool node::set_data_property(std::string name, const T& property) { return l->properties.set_data_property(name, index, property); }
bool node::unset_data_property(std::string name) { return l->properties.unset_data_property(name, index); }
template <typename T>
T& node::get_data_property(std::string name) { return l->properties.get_data_property<T>(name, index); }
template <typename T>
const T& node::get_data_property(std::string name) const { return l->properties.get_data_property<T>(name, index); }
/**
*		@Additionaly implemented function
*		@Adds to the missing_populations vector (only for the Fluid Boundary Nodes)
*/
void node::add_missing_populations(const int& i) {
	if (has_flag_property("Fluid_Boundary_Node"))
		missing_populations.push_back(i);
	else
		std::cout << "No Missing Populations! Not a fluid boundary node!" << std::endl;
}
/**
*		@Additionaly implemented function
*		@Returns copy of the vector of missing populations
*/
std::vector<int> node::return_missing_populations() const {return missing_populations;}
/**
*		@Additionaly implemented function
*		@Clears the missing_populations vector
*/
void node::clear_missing_populations() {missing_populations.clear();}



// lattice

lattice::lattice(unsigned int _nx, unsigned int _ny)
: nx(_nx), ny(_ny), size(nx*ny), buffer_size(1), real_nx(nx+2*buffer_size), real_ny(ny+2*buffer_size),
  real_size(real_nx*real_ny), n_populations(velocity_set().size),
  f( n_populations, std::vector<float_type>(real_size, 0) ),
  rho(real_size, 0), u(real_size, 0), v(real_size, 0), nodes(real_size),
  properties(real_size), periodic_x(true), periodic_y(true)
{
	// register some properties
	properties.register_flag_property("fluid");
	properties.register_flag_property("buffer");
	properties.register_flag_property("wall");
	properties.register_flag_property("solid"); //used for the cylinder
	properties.register_flag_property("Fluid_Boundary_Node"); //used for the nodes adjacent to the cylinder


	// set up nodes and properties
	unsigned int k(0);
	for (unsigned int j=0; j<real_ny; ++j)
	{
		for (unsigned int i=0; i<real_nx; ++i)
		{
			nodes[k].set(this, static_cast<int>(i)-buffer_size, static_cast<int>(j)-buffer_size);
			if (i<buffer_size || i>=real_nx-buffer_size || j<buffer_size || j>=real_ny-buffer_size)
				properties.set_flag_property("buffer",nodes[k].index);
			else properties.set_flag_property("fluid",nodes[k].index);
			++k;
		}
	}
}

lattice::node_iterator lattice::begin() { return nodes.begin(); }
lattice::const_node_iterator lattice::begin() const { return nodes.begin(); }
lattice::node_iterator lattice::end() { return nodes.end(); }
lattice::const_node_iterator lattice::end() const { return nodes.end(); }
lattice::reverse_node_iterator lattice::rbegin() { return nodes.rbegin(); }
lattice::const_reverse_node_iterator lattice::rbegin() const { return nodes.rbegin(); }
lattice::reverse_node_iterator lattice::rend() { return nodes.rend(); }
lattice::const_reverse_node_iterator lattice::rend() const { return nodes.rend(); }

inline unsigned int lattice::index(int i, int j) const { return real_nx*(j+buffer_size) + i + buffer_size; }

inline node& lattice::get_node(int i, int j) { return nodes[real_nx*(j+buffer_size) + i + buffer_size]; }

inline const node& lattice::get_node(int i, int j) const { return nodes[real_nx*(j+buffer_size) + i + buffer_size]; }

inline node& lattice::get_node(unsigned int idx) { return nodes[idx]; }

inline const node& lattice::get_node(unsigned int idx) const { return nodes[idx]; }


std::ostream& operator<<(std::ostream& os, const lattice& l)
{
	for (unsigned int p=0; p<l.n_populations; ++p)
	{
		os << " f" << std::setw(2) << p << ":";
		os << "\n" << "   y  " << std::setw((l.nx+2*l.buffer_size)*12+1) << std::setfill('-') << "" << "\n" << std::setfill(' ');
		for (int j=static_cast<int>(l.ny+l.buffer_size-1); j>-static_cast<int>(l.buffer_size+1); --j)
		{
			os << std::setw(4) << j << " |";
			for (int i=-static_cast<int>(l.buffer_size); i<static_cast<int>(l.nx+l.buffer_size); ++i)
			{
				const unsigned int index = (j+l.buffer_size)*l.real_nx + i + l.buffer_size;
				if (i>=0 && i<static_cast<int>(l.nx) && j>=0 && j<static_cast<int>(l.ny))
				os << std::setw(12) << std::setprecision(5) /*<< std::scientific*/ << l.f[p][index];
				else
				os << std::setw(12) << "*";
			}
			os << " |" <<  "\n";
		}
		os << std::setw(6) << "" << std::setw((l.nx+2*l.buffer_size)*12+1) << std::setfill('-') << "" << "\n" << std::setfill(' ') << std::setw(6) << "";
		for (int i=-static_cast<int>(l.buffer_size); i<static_cast<int>(l.nx+l.buffer_size); ++i) os << std::setw(12) << i;
		os << " x\n";
	}
	os << l.properties;
	return os;
}

void lattice::add_wall(coordinate<int> min_coord, coordinate<int> max_coord)
{
	for (int j=min_coord.j; j<=max_coord.j; ++j)
	{
		for (int i=min_coord.i; i<=max_coord.i; ++i)
		{
			// check if node not yet labelled as wall
			if (!get_node(i,j).has_flag_property("wall"))
			{
				// set wall property
				get_node(i,j).set_flag_property("wall");
				wall_nodes.push_back( get_node(i,j) );
			}
		}
	}
}

void lattice::add_wallCylinder(float_type center[2], float_type radius) //function to mark the solid nodes (In and on the cylinder) & Fluid Boundary Nodes(some of whose populations come from solid)
{
	std::vector<node> not_solid;
	int x1 = floor(center[0] - radius)-1;  int y1 = floor(center[1] - radius)-1;
	int x2 = ceil(center[0] + radius)+1;  int y2 = ceil(center[1] + radius)+1;
	coordinate<int> min_coord = {x1,y1}; //bottom left corner of the bounding box
	coordinate<int> max_coord = {x2,y2}; //top right corner of the bounding box

	//differentiate between the solid nodes and the rest
	for (int j = min_coord.j; j<=max_coord.j; ++j)
	{
		for (int i=min_coord.i; i<=max_coord.i; ++i)
		{
			if ( (i-center[0])*(i-center[0]) + (j-center[1])*(j-center[1]) - radius*radius <=0 ) //if on or inside the circle
			{
				if ( !(get_node(i,j).has_flag_property("solid")) ) { // if on or inside the circle and and not a solid node, mark as a solid node
					get_node(i,j).set_flag_property("solid");
					solid_nodes.push_back(get_node(i,j));
				}

			}
			else  //if the point is completely outside the circle
				not_solid.push_back(get_node(i,j));
		}
	}
	//differentiate between the rest nodes and Fluid boundary nodes, updating the missing populations vector for the fluid boundary nodes
	for ( std::vector<node>::iterator it = not_solid.begin() ; it!=not_solid.end() ; ++it ) //iterating over the nodes inside bounding box which are not solid
	{
		unsigned int x = it->coord.i;	unsigned int y = it->coord.j;
		//checking the populations which intersect with the solid nodes
		for (unsigned int i = 0 ; i<velocity_set().size ; ++i)
		{
			if ( get_node(x + velocity_set().c[0][i] ,y + velocity_set().c[1][i] ).has_flag_property("solid") ) // if the adjacent node (according to Ci is a solid node)
			{
				if ( !(get_node(x,y).has_flag_property("Fluid_Boundary_Node") ) )
				{
					get_node(x,y).set_flag_property("Fluid_Boundary_Node");
					fluid_boundary_nodes.push_back(*it);
				}
				get_node(x,y).add_missing_populations(i);
			}
		}
	}
	not_solid.clear();
}

void lattice::delete_fluid_boundary_nodes()
{
	for (node n : fluid_boundary_nodes)
	{
		unsigned int x = n.coord.i;	unsigned int y = n.coord.j;
		get_node(x,y).unset_flag_property("Fluid_Boundary_Node");
		get_node(x,y).clear_missing_populations();
	}
	fluid_boundary_nodes.clear();

}

void lattice::delete_walls()
{
	for (node n : wall_nodes)
	{
		n.unset_flag_property("wall");
	}
	wall_nodes.clear();

}

void lattice::delete_solids()
{
	for (node n : solid_nodes)
	{
		n.unset_flag_property("solid");
	}
	solid_nodes.clear();
}

void lattice::write_fields(std::string file_name)
{
	std::ofstream ofs(file_name.c_str());
	if (ofs.is_open())
	{
		// write header (comment that part if necessary)
		ofs << "x y rho u v\n";
		// write body
		for (unsigned int j=0; j<ny; ++j)
		{
			for (unsigned int i=0; i<nx; ++i)
			{
				ofs << i << " " << j << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].rho() << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].u() << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].v() << "\n";
			}
		}
	}
	else throw std::runtime_error("could not write to file");
}

} // lb

#endif //LB_LATTICE_HPP_INCLUDED
