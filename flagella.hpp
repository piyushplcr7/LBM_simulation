// Calculation of te movement properties of the flagella



#ifndef FLAG_HPP
#define FLAG_HPP

#include "global.hpp"
#include <cmath>

#include <boost/numeric/odeint/config.hpp>
#include "boost/numeric/odeint.hpp"

typedef std::vector< double > state_type;
using namespace std;



class flagella{
	
	//Flagella
	using float_type = float;
	const float_type l, m, theta_1, theta_2, k_1, k_2;
	const float_type x_0, y_0;
	const float_type y_move;
	const float_type omega;

	//variable dependent of time
	float_type alp_1, alp_2;
	float_type alp_d_1, alp_d_2;
	float_type alp_dd_1, alp_dd_2;
	float_type x, y;
	float_type x_d, y_d;
	float_type x_dd, y_dd;
	float_type Q_1, Q_2;
	const float_type B_1, B_2;

public:
	flagella(float_type l_fla, float_type mass_fla, float_type k_1_fla, float_type k_2_fla, float_type x, float_type y, float_type _omega, float_type _ymove)
	: l(l_fla),
	m(mass_fla),
	k_1(k_1_fla),
	k_2(k_2_fla),
	B_1(0.2*k_1),
	B_2(0.2*k_2),
/*
	B_1(0)
	B_2(0)
*/
	x_0(x),
	y_0(y),
	theta_1( 1/12 * m * l*l), theta_2(theta_1),
	omega(_omega),
	y_move(_ymove)
	{
		alp_1 = 3.14/4; alp_2 = 3.14/4;
		alp_d_1 = .0; alp_d_2 = .0;
		alp_dd_1 = .0; alp_dd_2 = .0;
		x_d = .0;
		x_dd = .0;
		y_d = y_move*omega* cos(.0);
		y_dd = .0;
		Q_1 = .0; Q_2 = .0;
		std::cout << std::setw(15) << "#alpha_1: " << std::setw(15) <<  " alpha_2: " << std::endl;
	}



	void updateX(int time);
	void step(float_type Q_1, float_type Q_2, float_type dt);
	void step(float_type Q_1, float_type Q_2);

	void GetRHS1(const state_type &x1, state_type &dxdt1, const double time);
	void GetRHS2(const state_type &x2 , state_type &dxdt2 , const double time);
	void ODE(float_type dt);


	//Access functions
	float_type getX0(){return x;};
	float_type getY0(){return y;};
	float_type getX1(){return x+cos(alp_1)*l;};
	float_type getY1(){return y+sin(alp_1)*l;};
	float_type getX2(){return x+cos(alp_1)*l+cos(alp_2)*l;};
	float_type getY2(){return y+ sin(alp_1) * l + sin(alp_2)*l;};
	void uptCoord(float_type *coord);

};

/*flagella::flagella(float_type l_fla, float_type mass_fla, float_type k_1_fla, float_type k_2_fla, float_type x, float_type y ){
	
	theta_1 = 1/12 * m * l*l; theta_2 = theta_1;

	//Time dependent variables
	alp_1 = 0; alp_2 = 0;
	alp_d_1 = 0; alp_d_2 = 0;
	alp_dd_1 = 0; alp_dd_2 = 0;
	x = x_0;
	x_d = 0;
	x_dd = 0;
	y = y_0;
	y_d = y_move*omega* cos(omega*time);
	y_dd = 0;

} */

//Updates the velocity of the attached part to the cylinder
void flagella::updateX(int time){
	y = y_0 + y_move* sin(omega*time);
	x = x_0;
	y_d = y_move*omega* cos(omega*time);
	y_dd = -y_move*omega*omega* sin(omega*time);
	/*Cyl_center[1] = Cyl_center_0[1] + y_move * sin(omega*time);
	Cyl_vel[1] = y_move * omega * cos(omega*time); */
}

//RHS side of dt alpha_2, ddt alpha_2 
void flagella::GetRHS2(const state_type &x2 , state_type &dxdt2 , const double time){
	//Set matrix coefficient
	float_type a11 = 1.25*m*l*l + theta_1;
	float_type a12 = -0.5*m *l*l* cos(alp_1+alp_2);
	float_type a21 = 0.5*m*l*l*cos(alp_2-alp_1);
	float_type a22 = 0.5*m*l*l + theta_2;

	//Set right hand side of equation 
	float_type RHS1 = -0.5 * m *l*l * alp_d_2*alp_d_2* sin(alp_1-alp_2) - 1.5*m *l*(-sin(alp_1)*x_dd+cos(alp_1)*y_dd)-k_1*alp_1-k_2*(alp_2-alp_1)+Q_1- alp_d_1*B_1;
	float_type RHS2 = -0.5 * m *l*l * alp_d_1*alp_d_1* sin(alp_2-alp_1) + 0.5*m *l*(sin(alp_2)*x_dd+cos(alp_2)*y_dd)-k_2*(alp_2-alp_1)+Q_2- alp_d_2*B_2;

	//std::cout << "a11: " << a11 << " / a12: " << a12 << " / a21: " << a21 << "/ a22: " << a22 << std::endl;
	//std::cout << "RHS1: " << RHS1 << " RHS2 " << RHS2 << "/// damping1: " << alp_d_1*B_1 << "  damping2: " << alp_d_2*B_2 <<   std::endl;
	
	dxdt2[0] = x2[1];
	dxdt2[1] = 1/a22*(RHS2 - a21/a11 * RHS1);
	alp_dd_2 = dxdt2[1];
	alp_d_2 = x2[1];
	alp_2 = x2[0];

	//std::cout << "alpha_2 : " << x2[0] << " alpha_2' " << x2[1] << std::endl;

}

//dt alpha_1, ddt alpha_1,
void flagella::GetRHS1(const state_type &x1, state_type &dxdt1, const double time){
	//Set matrix coefficient
	float_type a11 = 1.25*m*l*l + theta_1;
	float_type a12 = -0.5*m *l*l* cos(alp_1+alp_2);
	float_type a21 = 0.5*m*l*l*cos(alp_2-alp_1);
	float_type a22 = 0.5*m*l*l + theta_2;

	//Set right hand side of equation 
	float_type RHS1 = -0.5 * m *l*l * alp_d_2*alp_d_2* sin(alp_1-alp_2) - 1.5*m *l*(-sin(alp_1)*x_dd+cos(alp_1)*y_dd)-k_1*alp_1-k_2*(alp_2-alp_1)+Q_1 - alp_d_1*B_1;
	float_type RHS2 = -0.5 * m *l*l * alp_d_1*alp_d_1* sin(alp_2-alp_1) + 0.5*m *l*(sin(alp_2)*x_dd+cos(alp_2)*y_dd)-k_2*(alp_2-alp_1)+Q_2 - alp_d_2*B_2;		
	
	//std::cout << "a11: " << a11 << " / a12: " << a12 << " / a21: " << a21 << "/ a22: " << a22 << std::endl;
	//std::cout << "RHS1: " << RHS1 << " RHS2 " << RHS2 << "/// damping1: " << alp_d_1*B_1 << "  damping2: " << alp_d_2*B_2 <<   std::endl;	
	
	dxdt1[0] = x1[1];
	dxdt1[1] = 1/a11*(RHS1 - a21 * alp_dd_2);

	alp_dd_1 = dxdt1[1];
	alp_d_1 = dxdt1[0];
	alp_1 = x1[0];
	
	//std::cout << "alpha_1 : " << x1[0] << " alpha_1' " << x1[1]<<std::endl;

}

void flagella::ODE(float_type dt){
	state_type alpha_1(2);
	state_type alpha_2(2);
	
	alpha_1[0] = alp_1;
	alpha_1[1] = alp_d_1;
	alpha_2[0] = alp_2;
	alpha_2[1] = alp_d_2;
	
	
	boost::numeric::odeint::integrate(std::bind(&flagella::GetRHS2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha_2, 0.0, 1.0, 0.2);
	//std::cout << "alpha_2 : " << alpha_2[0] << " alpha_2' " << alpha_2[1]<<std::endl;

	boost::numeric::odeint::integrate(std::bind(&flagella::GetRHS1, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha_1, 0.0, 1.0, 0.2);
	//std::cout << "alpha_1 : " << alpha_1[0] << " alpha_1' " << alpha_1[1]<<std::endl;
}

//calculates one step forward of the body movement equation
void flagella::step(float_type _Q_1, float_type _Q_2, float_type dt){
	Q_1 = _Q_1;
	Q_2 = _Q_2;
 /*	
	//Set left hand side of equation 
	float_type LHS1 = -0.5 * m *l*l * alp_d_2*alp_d_2* sin(alp_1-alp_2) - 1.5*m *l*(-sin(alp_1)*x_dd+cos(alp_1)*y_dd)-k_1*alp_1-k_2*(alp_2-alp_1)+Q_1;
	float_type LHS2 = -0.5 * m *l*l * alp_d_1*alp_d_1* sin(alp_2-alp_1) + 0.5*m *l*(sin(alp_2)*x_dd+cos(alp_2)*y_dd)-k_2*(alp_2-alp_1)+Q_2;
	
	//Set matrix coefficient
	float_type a11= 1.25* m * l*l + theta_1;
	float_type a12 = -0.5*m *l*l *cos(alp_1+alp_2);
	float_type a21 = 0.5*m*l*l*cos(alp_2-alp_1);
	float_type a22 = 0.5*m*l*l+theta_2;

	//calculate angular accelaration
	alp_dd_2 = 1/a22*(LHS2 -a21/a11 * LHS1);
	alp_dd_1 = 1/a11*(LHS1 - a21 * alp_dd_2);

	//update angular velocities
	alp_d_1 = alp_d_1 + alp_dd_1 * dt;
	alp_d_2 = alp_d_2 + alp_dd_2 *dt;

	//update angular position
	alp_1 = alp_1 + alp_d_1 * dt;
	alp_2 = alp_2 + alp_d_2 * dt; */
	
	ODE(dt);

	std::cout << std::setw(15) << alp_1 << std::setw(15) << alp_2 << std::endl;
}

//dt == 1 
void flagella::step(float_type Q_1, float_type Q_2){
	step(Q_1, Q_2, 1);

}

//Updates the passed array with the new coordinates of the flagella
void flagella::uptCoord(float_type *coord){
	coord[0] = x;
	coord[1] = y;
	coord[2] = x+cos(alp_1)*l;
	coord[3] = y+sin(alp_1)*l;
	coord[4] = x+cos(alp_1)*l+cos(alp_2)*l;
	coord[5] = y+ sin(alp_1) * l + sin(alp_2)*l;
}

#endif // FlAG_HPP
