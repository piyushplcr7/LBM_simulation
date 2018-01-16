// Calculation of te movement properties of the flagella



#ifndef FLAG_HPP
#define FLAG_HPP

#include "global.hpp"
#include <cmath>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>

#include<fstream>

typedef std::vector< double > state_type;
using namespace std;

class flagella{
	
	using float_type = double;

	//Flagella
	const int n;
	const std::vector<float_type> l, m, k_spring;
	const std::vector<float_type> B;
	const float_type x_0, y_0;
	const float_type y_move;
	const float_type omega;
	std::vector<float_type> theta;
	std::ofstream txt_a;

	//variable dependent of time
	float_type t;
	float_type alpha_1_0;
	state_type alpha;
	std::vector<double> dd_alpha;
	std::vector<double> dd_alpha_old;

	float_type x, y;
	float_type x_d, y_d;
	float_type x_dd, y_dd;
	std::vector<float_type> Q;
	std::vector<std::vector<float_type>> x_vec;
	bool flag_out;

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
		dd_alpha_old(n,0),
		Q(n,0),
		x_vec(n,{0,0}),
		flag_out(true),
		t(0.0)
		{
			x_d = .0;
			x_dd = .0;
			y_d = y_move*omega* cos(.0);
			y_dd = .0;
			for( int i = 0; i<n ; ++i){ theta.push_back(1/12*m[i]*l[i]*l[i]) ;}

			txt_a.open("Alpha.txt",std::ios::out);
			for( int i= 0; i<n; ++i)
				txt_a << std::setw(15) << "#alpha_" << std::to_string(i) << ": "; 
			txt_a <<std::endl;
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
		dd_alpha_old(n,0),
		Q (n,0),
		x_vec(n,{0,0}),
		flag_out(true),
		t(0.0)
		{
			x_d = .0;
			x_dd = .0;
			y_d = y_move*omega* cos(.0);
			y_dd = .0;
			for( int i = 0; i<n ; ++i){ theta.push_back(1/12*m[i]*l[i]*l[i]) ;}

			txt_a.open("Alpha.txt",std::ios::out);
			for( int i= 0; i<n; ++i)
				txt_a << std::setw(15) << "#alpha_" << std::to_string(i) << ": "; 
			txt_a <<std::endl;
			}

	void InitAngle(int i, float_type j){alpha[2*i]= j; if( i == 0) alpha_1_0 = j;}
	void step(float_type dt);
	void step(std::vector<float_type> _Q, float_type dt){Q = _Q; step(dt);}
	void step(std::vector<float_type> Q){step(Q, 1.0);}

	void GetRHS(const state_type &x, state_type &dxdt, const double time);
	void getRHS2D(const state_type &x, state_type &dxdt, const double time);

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

void flagella::getRHS2D(const state_type &x, state_type &dxdt, const double time){
	int k = 0;
	float_type part1, part2, part3, part4, partk, RHS;

	dxdt[0] = x[1];
	//k=1
	RHS = Q[0] - k_spring[0]*(x[0]) + k_spring[1]*(x[2]-x[0]) - B[0] * x[1] + dd_alpha[1] * l[0]* l[1] *( 0.5 * m[1]) *cos(x[0]-x[2]) - x[3]*x[3] * l[2]* l[1] *( 0.5 * m[1])*sin(x[0]-x[2]) ;
	partk = l[0]*l[0]*(0.25*m[0] + theta[0])+l[1]*l[1]*m[1];
	dxdt[1] = RHS/partk;
	dd_alpha[0] = dxdt[1];

	dxdt[2] = x[3];
	//k = 2
	RHS = Q[1] -k_spring[1]*(x[2]-x[0]) - B[1] * x[3] - dd_alpha[0] * l[0]* l[1] *( 0.5 * m[1]) *cos(x[0]-x[2]) + x[1]*x[1] * l[2]* l[1] *( 0.5 * m[1])*sin(x[0]-x[2]) ;
	partk = l[1]*l[1]*(0.25*m[1] + theta[1]);
	dxdt[3] = RHS/partk;
	dd_alpha[1] = dxdt[3];
 
};	

//dt alpha_1, ddt alpha_1,
void flagella::GetRHS(const state_type &x, state_type &dxdt, const double time){
	
	int k = 0;
	float_type part1, part2, part3, part4, partk, RHS;

	//i iterates over the whole vector q_n = (alpha_1, alp_1_d,, alp_2, alp_2_d, ...); [0, 2*n-1]
	//k iterates only over the elements of the flagella; [0, n-1]
	for(int i = 0; i < 2*n; ++i){
		part1 = 0; part2 = 0; part3 = 0; part4 = 0; partk = 0;

		//alpha_k = x[2*k]
		//alpha_d_k = x[2*k+1]

		// if i == alpha and dxdt = alpha_d
		if(i % 2 == 0){
			dxdt[i] = x[i+1];
		}

		//if i  == alpha_d and dxdt = alpha_dd
		else if(i % 2 == 1){
			if(i==1){ // First Row
				RHS = Q[k] - k_spring[k]*(x[i-1])+k_spring[k+1]*(x[i+1]-x[i-1]) - B[k] * x[i] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				part1 = 0; part3 = 0;
				//Part 2 and 4
				for(int j = k+1; j < n-1/*-1*/; ++j){
					part2 = part2 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
					part4 = part4 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
				}
			}
			else if(i==2*n-1){ //Last Row
				RHS = Q[k] -k_spring[k]*(x[i-1]-x[i-3])- B[k]* x[i]- l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				part2 = 0; part4 = 0;
				//Part 1 and 2
				for(int j = 0; j < k-1; ++j){
					part1 = part1 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*cos(x[2*k]-x[2*j]);
					part3 = part3 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*sin(x[2*k]-x[2*j]);
				}
			}
			else{ // Normal Row
				RHS = Q[k] -k_spring[k]*(x[i-1]-x[i-3])+k_spring[k+1]*(x[i+1]-x[i-1])- B[k]* x[i] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
				//Part 1 and 3
				for(int j = 0; j < k-1; ++j){
					part1 = part1 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*cos(x[2*k]-x[2*j]);
					part3 = part3 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[k]+ sumMass(k+1, n))*sin(x[2*k]-x[2*j]);
				}
				//Part 2 and 4
				for(int j = k+1; j < n-1/*-1*/; ++j){
					part2 = part2 + dd_alpha[j] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
					part4 = part4 + x[j*2 + 1]*x[j*2 + 1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
				}
			}
			
			partk = l[k]*l[k]*(0.25*m[k] + sumMass(k+1, n))+ theta[k];
			dxdt[i] = (RHS -part1 - part2 -part2- part4)/partk;
			dd_alpha[k] = dxdt[i];
			if(!(abs(part1/partk) < 1. and abs(part2/partk) < 1. and abs(part3/partk) < 1. and abs(part4/partk) < 1.))
				flag_out = false;
			if(flag_out){
				std::cout << RHS << " // " << part1 << " // " << part2 <<  " // " << part3 << " // " << part4 << " // " << partk << std::endl;
				std::cout << "dd_Alpha " << std::to_string(k) << ": " << dd_alpha[k] << " / d_Alpha " << std::to_string(k) << ": " << x[i] << " / Alpha " << std::to_string(k) << ": " << x[i-1] << std::endl;
			}
			k=k+1;
		}
	}
	dd_alpha_old = dd_alpha;
}


//calculates one step forward of the body movement equation
void flagella::step(float_type delta_t){
	int split = 100;
	///*
	if(n==2){
		boost::numeric::odeint::integrate(std::bind(&flagella::getRHS2D, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
	}
	else{
		boost::numeric::odeint::integrate(std::bind(&flagella::GetRHS, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
	}
	//*/
	

	//boost::numeric::odeint::integrate(std::bind(&flagella::GetRHS, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
	//Update second derivatives of the angle
	//for(int k = 0; k <n; k++){
	//	alpha[2*k+1] = dd_alpha[k];
	//}
	t = t + delta_t;
	for( int i = 0; i<n; ++i){
		txt_a << std::setw(15) << alpha[2*i];
		if(n==1){
			float_type m_star = (l[0]*l[0]*0.25*m[0]+theta[0]);
			float_type delta = B[0]/m_star/2;
			float_type omega_0 = sqrt(k_spring[0]/m_star);
			float_type omega = sqrt(omega_0*omega_0-delta*delta);
			//std::cout << alpha_1_0 << "  // " << m_star << " // " << omega_0 << " // " << delta << std::endl;
			txt_a << std::setw(15) << alpha_1_0*exp(-delta*t)*cos(omega*t);
		}
	}
	txt_a << std::endl;
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
