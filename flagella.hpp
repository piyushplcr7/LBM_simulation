// Calculation of te movement properties of the flagella



#ifndef FLAG_HPP
#define FLAG_HPP

#include "global.hpp"
#include <cmath>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>

typedef std::vector< double > state_type;
using namespace std;



class flagella{
	
	//Flagella
	using float_type = float;
	const int n;
	const std::vector<float_type> l, m, k_spring;
	const std::vector<float_type> B;
	const float_type x_0, y_0;
	const float_type y_move;
	const float_type omega;
	std::vector<float_type> theta;

	//variable dependent of time
	state_type alpha;
	std::vector<double> dd_alpha;

	float_type x, y;
	float_type x_d, y_d;
	float_type x_dd, y_dd;
	std::vector<float_type> Q;
	std::vector<std::vector<float_type>> x_vec;

public:
	//Vector Constructor
	flagella(int _n, std::vector<float_type> l_fla, std::vector<float_type> mass_fla, std::vector<float_type> k_fla, std::vector<float_type> B_fla, float_type x, float_type y, float_type _omega, float_type _ymove):
		n(_n),
		l(l_fla),
		m(mass_fla),
		k_spring(k_fla),
		B(B_fla),
		x_0(x),
		y_0(y),
		omega(_omega),
		y_move(_ymove),
		alpha (2*n,0),
		dd_alpha(n,0),
		Q(n,0),
		x_vec(n,{0,0})
		{
			x_d = .0;
			x_dd = .0;
			y_d = y_move*omega* cos(.0);
			y_dd = .0;
			for( int i = 0; i<n ; ++i){ theta.push_back(1/12*m[i]*l[i]*l[i]) ;}

			for( int i= 0; i<n; ++i){std::cout << std::setw(15) << "#alpha_" << std::to_string(i) << ": ";} std::cout <<std::endl;
		}

	//Uniform Constructor
	flagella(int _n, float_type l_fla, float_type mass_fla, float_type k_fla, float_type B_fla, float_type x, float_type y, float_type _omega, float_type _ymove):
		n(_n),
		l(n,l_fla),
		m(n,mass_fla),
		k_spring(n,k_fla),
		B(n,B_fla),
		x_0(x),
		y_0(y),
		omega(_omega),
		y_move(_ymove),
		alpha (2*n,0),
		dd_alpha(n,0),
		Q (n,0),
		x_vec(n,{0,0})
		{
			x_d = .0;
			x_dd = .0;
			y_d = y_move*omega* cos(.0);
			y_dd = .0;
			for( int i = 0; i<n ; ++i){ theta.push_back(1/12*m[i]*l[i]*l[i]) ;}
			for( int i= 0; i<n; ++i){std::cout << std::setw(15) << "#alpha_" << std::to_string(i) << ": ";} std::cout <<std::endl;
		}

	void InitAngle(int i, float_type j){alpha[2*i]= j;}
	void step(float_type dt);
	void step(std::vector<float_type> _Q, float_type dt){Q = _Q; step(dt);}
	void step(std::vector<float_type> Q){step(Q, 1.0);}

	void GetRHS(const state_type &x, state_type &dxdt, const double time);

	//Coordinates
	void updx();
	//Change speed of cylinder
	void updateX0(int time);

	//Access functions
	float_type getX(int i){updx(); return x_vec[i][0];};
	float_type getY(int j){updx(); return x_vec[0][j];};
	std::vector<float_type> getCoord(int i){updx(); return x_vec[i];}
	float_type sumMass(int first, int last);
};

//dt alpha_1, ddt alpha_1,
void flagella::GetRHS(const state_type &x, state_type &dxdt, const double time){
	
	int k = 0;
	//i iterates over the whole vector q_n = (alpha_1, alp_1_d, alp_1_dd, alp_2, alp_2_d, ...); [0, 2*n-1]
	//k iterates only over the elements of the flagella; [0, n-1]
	for(int i = 0; i < 2*n; ++i){
		float_type part1 = 0;
		float_type part2 = 0;
		float_type part3 = 0;
		float_type part4 = 0;
		float_type partk = 0;

		//alpha_k = x[2*k]
		//alpha_d_k = x[2*k+1]
		//alpha_dd_k = x[2*k + 2]

		// if i == alpha and dxdt = alpha_d
		if(i % 2 == 0){
			dxdt[i] = x[i+1];
		}

		//if i  == alpha_d and dxdt = alpha_dd
		if(i % 2 == 1){
			if(i==0){ // First Row
				dxdt[i] = Q[k] -k_spring[k]*(x[i-1])+k_spring[k+1]*(x[i+2]-x[i-1])- B[k]* x[i]- l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				part1 = 0; part2 = 0;
				//Part 2 and 4
				for(int j = k+1; j < n-1; ++j){
					part2 = part2 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
					part4 = part4 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
				}
			}
			if(i==2*n-1){ //Last Row
				dxdt[i] = Q[k] -k_spring[k]*(x[i-1]-x[i-4])- B[k]* x[i]- l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				part2 = 0; part4 = 0;
				//Part 1 and 2
				for(int j = 0; j < k-1; ++j){
					part1 = part1 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*cos(x[2*k]-x[2*j]);
					part3 = part3 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*sin(x[2*k]-x[2*j]);
				}
			}
			else{ // Normal Row
				dxdt[i] = Q[k] -k_spring[k]*(x[i-1]-x[i-4])+k_spring[k+1]*(x[i+2]-x[i-1])- B[k]* x[i] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				//Part 1 and 3
				for(int j = 0; j < k-1; ++j){
					part1 = part1 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*cos(x[2*k]-x[2*j]);
					part3 = part3 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*sin(x[2*k]-x[2*j]);
				}
				//Part 2 and 4
				for(int j = k+1; j < n-1; ++j){
					part2 = part2 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
					part4 = part4 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
				}
			}

			partk = l[k]*l[k]*(0.25*m[k] + sumMass(k+1, n))+ theta[k];
			dxdt[i] = (dxdt[i] -part1 - part2 -part2- part4)/partk;
			//std::cout << "dd_alpha" << dd_alpha[0] <<  "//" << dd_alpha[1] <<  "//" <<dd_alpha[2] << std::endl;
			//std::cout << part1 << "//" << part2 <<  "//" << part3 << "//" << part4 << "//" << partk << std::endl;
			dd_alpha[k] = dxdt[i];
			k=k+1;
		}
		//if i == alpha_dd and dxdt = alpha_ddd
		//if(i % 2 == 2){
		//	dxdt[i] = 0;
		//}
	}

	//Save second derivatives of the angle
	//for(int k = 0; k <n; k++){
	//	dd_alpha[k] = dxdt[2*k+1];
	//}
}


//calculates one step forward of the body movement equation
void flagella::step(float_type delta_t){
	int split = 100;

	boost::numeric::odeint::integrate(std::bind(&flagella::GetRHS, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
		
	//Update second derivatives of the angle
	//for(int k = 0; k <n; k++){
	//	alpha[2*k+1] = dd_alpha[k];
	//}
	for( int i = 0; i<n; ++i)
		std::cout << std::setw(15) << alpha[2*i];
	std::cout << std::endl;
}

//Update Coordinates
void flagella::updx(){
	x_vec[0][0] = x_0 + l[0] * cos(alpha[0]);
	x_vec[0][1] = y_0 + l[0] * sin(alpha[0]);

	for( int i = 1; i < n; ++i){
		x_vec[i][0] = x_vec[i-1][0] + l[1] * cos(alpha[1]);
		x_vec[i][1] = x_vec[i-1][1] + l[1] * sin(alpha[1]);
	}
}

flagella::float_type flagella::sumMass(int first, int last){
	float_type sum = 0;
	for(int i = first; i < last; ++i){
		sum = sum + m[i];
	}
	return sum;
}

//Updates the velocity of the attached part to the cylinder
void flagella::updateX0(int time){
	y = y_0 + y_move* sin(omega*time);
	x = x_0;
	y_d = y_move*omega* cos(omega*time);
	y_dd = -y_move*omega*omega* sin(omega*time);
}


#endif // FlAG_HPP
