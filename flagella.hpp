// Calculation of te movement properties of the flagella



#ifndef FLAG_HPP
#define FLAG_HPP

#include "global.hpp"
#include <cmath>
#include <utility>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <lapacke.h>
#include <eigen3/Eigen/Dense>

#include<fstream>

//typedef std::vector< double > state_type;
typedef boost::numeric::ublas::vector<double> state_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;
using namespace std;

namespace lb{

class flagella{
	
	//using float_type = double;

	//Flagella
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
	state_type dd_alpha;

	float_type x, y;
	float_type x_d, y_d;
	float_type x_dd, y_dd;
	std::vector<float_type> Q;
	std::vector<std::vector<float_type>> x_vec;
	bool flag_out;

public:
	const int n;
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
		//dd_alpha_old(n,0),
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
		//dd_alpha_old(n,0),
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
	state_type stepRet(std::vector<float_type> _Q, float_type dt){step(_Q, dt); return alpha;}

	std::vector<float> stepRetvec(std::vector<float_type> _Q, float_type dt){step(_Q, dt); return alpha2vec();}
	std::vector<float> alpha2vec();

	void GetRHS(const state_type &x, state_type &dxdt, const double time);
	void getRHS1D(const state_type &x, state_type &dxdt, const double time);
	void getRHS2D(const state_type &x, state_type &dxdt, const double time);
	void getRHSMatrix(const state_type &x, state_type &dxdt, const double time);
	void getJacobiT(const state_type &x, matrix_type &J, const double time, state_type &dfdt);
	void getJacobi(const state_type &x, matrix_type &J, const double time);

	void writeOut(boost::numeric::ublas::matrix<double> M);
	coordinate<float_type> getX0(){coordinate<float_type> x0; x0.i = x; x0.j=y; return x0; }
	

	float_type eval_M(int link_no, double Fx, double Fy, unsigned int xb, unsigned int yb);
	std::pair<coordinate<int>, coordinate<int>> get_bbox();

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

void flagella::getRHS1D(const state_type &x, state_type &dxdt, const double time){
	int k = 0;
	float_type partk, RHS;

	dxdt[0] = x[1];
	//k=1
	RHS = Q[0] - k_spring[0]*(x[0]) - B[0] * x[1]- l[0]*0.5*m[0]* (-x_dd*sin(x[0])+ y_dd*cos(x[0]));
	partk = l[0]*l[0]*(0.25*m[0] + theta[0]);
	dxdt[1] = RHS/partk;
	dd_alpha[0] = dxdt[1];
 
};	

void flagella::getRHS2D(const state_type &x, state_type &dxdt, const double time){
	int k = 0;
	float_type partk, RHS;

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

void flagella::getRHSMatrix(const state_type &x, state_type &dxdt, const double time){
	boost::numeric::ublas::matrix<double>  Matrix (n,n);
	boost::numeric::ublas::vector<double>  RHS (n);


	int i = 0;
	for( int k= 0; k<n; ++k){
		i=2*k;
		Matrix(k,k) = l[k]*l[k]*(0.25*m[k] + sumMass(k+1, n))+ theta[k];
		if(k==0){
			RHS[k] = Q[k] - k_spring[k]*(x[i])+k_spring[k+1]*(x[i+2]-x[i]) - B[k] * x[i+1] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
	        for (int j =k+1; j < n-2; ++j){
	            Matrix(k,j) = l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
	            RHS(k) = RHS(k) - x[j*2+1]*x[j*2+1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
			}
		}
		else if(k == n-1){
			RHS(k) = Q[k] - k_spring[k]*(x[i]-x[i-2]) - B[k]*x[i+1] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
	        for (int j =0; j < k-1; ++j){
	            Matrix(k,j) = l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
	            RHS(k) = RHS(k) - x[j*2+1]*x[j*2+1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
			}	
		}
		else{
			RHS(k) = Q[k] - k_spring[k]*(x[i]-x[i-2]) +k_spring[k+1]*(x[i+2]-x[i])- B[k]*x[i+1] - l[k]* ( 0.5 *m[k]+ sumMass(k+1, n))*(-x_dd *sin(x[2*k])+y_dd*cos(x[2*k]));
	        for (int j =0; j < k-1; ++j){
	            Matrix(k,j) = l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
	            RHS(k) = RHS(k) - x[j*2+1]*x[j*2+1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
			}
			for (int j =k+1; j < n-2; ++j){
	            Matrix(k,j) = l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*cos(x[2*k]-x[2*j]);
	            RHS(k) = RHS(k) - x[j*2+1]*x[j*2+1] * l[j]* l[k] *( 0.5 * m[j]+ sumMass(j+1, n))*sin(x[2*k]-x[2*j]);
			}
		}
	}

	//writeOut(Matrix);
	//std::cout << "before" << std::endl;
	bool Eigen = true ;
	if(Eigen){
		Eigen::MatrixXd Mat(n,n);
		Eigen::VectorXd RHS_d(n);
		for(int i = 0; i<n; ++i){
			RHS_d(i) = RHS(i);
			for( int j = 0; j <n ; ++j){
				Mat(i,j) = Matrix(i,j);
			}
		}
		Eigen::VectorXd x_star;
		if (n<=4){
			x_star = Mat.inverse()*RHS_d;
		}
		else{
			x_star = Mat.fullPivLu().solve(RHS_d);
		}

		for(int i = 0; i<n; ++i){
			RHS(i) = RHS_d(i);
		}
		dd_alpha = RHS;
	}
	else{
		boost::numeric::ublas::permutation_matrix <std::size_t> piv(Matrix.size1());
		boost::numeric::ublas::lu_factorize(Matrix, piv);
		dd_alpha = RHS;
		boost::numeric::ublas::lu_substitute(Matrix, piv, dd_alpha);
	}

	//writeOut(Matrix);
	//std::cout << "afterwards" << std::endl;	

	for(int i = 0; i<2*n;++i){
		if(i%2==0)
			dxdt(i) = x(i+1);
		else if(i%2 == 1)
			dxdt(i) = dd_alpha((i-1)/2);
	}
	//std::cout << time << std::endl;
}

void flagella::writeOut(boost::numeric::ublas::matrix<double> M){
	int round = 5;
	std::cout << std::endl << std::setw(15);
	for(int i = 0; i < M.size1(); ++i){
		for(int j = 0; j < M.size2(); ++j)
			std::cout << std::round(pow(10,round)*M(i,j))*pow(10,-round) << std::setw(15);
	std::cout << std::endl;
	}
	std::cout << std::endl;
}

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
	//dd_alpha_old = dd_alpha;
}

void flagella::getJacobiT(const state_type &x, matrix_type &J, const double time, state_type &dfdt){

	getJacobi(x,J,time);
	for(int j=0; j <2*n; ++j)
		dfdt[j] = 0.;
}

void flagella::getJacobi(const state_type &x, matrix_type &J, const double time){

	int k = 0;
	for(int j = 0; j < 2*n; ++j){

		if(j%2 == 0){//dfdx of Angles
			J(j,j+1) = 1.;
			}
		else if(j%2 == 1){ // dfdx of Angular Velocities
			//From angle
			//i<k
			for(int i = 0; i<k-1; ++i){
				J(j,2*i) = - x[2*i+1]*x[2*i+1]*l[k]*l[i]*(0.5*m[k]+sumMass(k+1,n))*cos(x[2*k]- x[2*(i)]);
			}
			if(k>0){
				J(j,2*(k-1)) = - x[2*k+1]*x[2*k+1]*l[k]*l[k-1]*(0.5*m[k]+sumMass(k+1,n))*cos(x[2*k]- x[2*(k-1)]) - k_spring[k];
			}
	
			//i=k
			J(j,2*k) = k_spring[k]+k_spring[k+1] - l[k]*(0.5*m[k]+ sumMass(k+1,n))*(x_dd*cos(x[2*k])+y_dd*sin(x[2*k]));
			for(int i=0; i < k-1; ++i){
				J(j,2*k) = J(k,k) + x[2*i+1]*x[2*i+1]*l[k]*l[i]*(0.5*m[k]+sumMass(k+1,n))*cos(x[2*k]- x[2*i]);
			}
			for(int i = k+1; i < n; ++i){
				J(j,2*k) = J(k,k) - x[2*i+1]*x[2*i+1]*l[k]*l[i]*(0.5*m[k]+sumMass(i+1,n))*sin(x[2*k]- x[2*i]);
			}
	
			//i>k
			if(k<n-1){
				J(j,2*(k+1))= - x[2*k+1]*x[2*(k+1)+1]*l[k]*l[k+1]*(0.5*m[k]+sumMass(k+2,n))*cos(x[2*k]- x[2*(k+1)]) +k_spring[k+1];
			}
			for( int i = k+2; i < n; ++i){
				J(j,2*i) = - x[2*i+1]*x[2*i+1]*l[k]*l[i]*(0.5*m[k]+sumMass(i+1,n))*cos(x[2*k]- x[2*(i)]);
			}
		
			//From angular vel.
			//i<k
			for(int i = 0; i<k; ++i){
				J(j, 2*i+1) = 2*x[2*i+1]*l[i]*l[k]*(0.5*m[k]+ sumMass(k+1, n))*sin(x[2*k]-x[2*i]);
			}
			//i==k
			J(j,2*k+1) = B[k];
			//i<k
			for( int i = k+1; i<n; ++i){
				J(j,2*i+1) = 2*x[2*i+1]*l[i]*l[k]*(0.5*m[k]+ sumMass(i+1, n))*sin(x[2*k]-x[2*i]);
			}

			k=k+1;
		}
	}
}


//calculates one step forward of the body movement equation
void flagella::step(float_type delta_t){
	int split = 100;
	///*
	boost::numeric::odeint::runge_kutta_cash_karp54<state_type> stepper; //seems to be the best
	//boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
	//boost::numeric::odeint::bulirsch_stoer<state_type> stepper;
	
	if(n==1){
		boost::numeric::odeint::integrate(std::bind(&flagella::getRHS1D, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
	}
	else if(n==2){
		boost::numeric::odeint::integrate(std::bind(&flagella::getRHS2D, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
	}
	else{
		boost::numeric::odeint::integrate_adaptive(stepper, std::bind(&flagella::getRHSMatrix, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), alpha, float_type(0.0), delta_t, delta_t/split);
		/*boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_dense_output<boost::numeric::odeint::rosenbrock4<double>>(1.06e-6, 1.06e-6),
		 std::make_pair(std::bind(&flagella::getRHSMatrix, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3),std::bind(&flagella::getJacobiT, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)),
		 alpha, float_type(0.0), delta_t, delta_t/split);*/
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

float_type flagella::eval_M(int link_no, double Fx, double Fy, unsigned int xb, unsigned int yb){
	float_type dx, dy, x1, y1, moment; //x2, y2, t,
	assert(link_no <= n);
	if (link_no == 0){
		x1 = x_0; y1 = y_0;
		//x2 = x_vec[0][0]; y2 = x_vec[0][1];
	}
	else{
		x1 = x_vec[link_no-1][0]; y1 = x_vec[link_no-1][1];
		//x2= x_vec[link_no][0]; y2= x_vec[link_no][1];
	}
	/*dx = x2 -x1;
	dy = y2-y1;
	t = ( (dx)*(xb-x1)+(dy)*(yb-y1) )/( dx*dx + dy*dy); //find t of line equation for shortest distance
	moment = (xb-x1)*t*Fy -(yb-y1)*t*Fx; //calculate moment from forces on the link*/
	moment = (xb-x1)*Fy - (yb-y1)*Fx;
	return moment;
}

std::pair<coordinate<int>, coordinate<int>> flagella::get_bbox(){
	coordinate<int> min,max;
	float_type xmin, xmax, ymin, ymax;
	int safety = 2;

	x_vec[0][0] = x_0 + l[0] * cos(alpha[0]);
	x_vec[0][1] = y_0 + l[0] * sin(alpha[0]);
	xmin = x_vec[0][0]; xmax = x_vec[0][0];
	ymin = x_vec[0][1]; ymax = x_vec[0][1];

	for( int i = 1; i < n; ++i){
		x_vec[i][0] = x_vec[i-1][0] + l[1] * cos(alpha[1]);
		x_vec[i][1] = x_vec[i-1][1] + l[1] * sin(alpha[1]);
		if(xmin > x_vec[i][0])
			xmin = x_vec[i][0];
		else if(xmax < x_vec[i][0])
			xmax = x_vec[i][0];
		if(ymin > x_vec[i][1])
			ymin = x_vec[i][1];
		else if(ymax < x_vec[i][1])
			ymax = x_vec[i][1];
	}

	min.i = (int) xmin - safety;
	min.j = (int) ymin - safety;
	max.i = (int) xmax + safety;
	max.j = (int) ymax + safety;

	return std::make_pair(min, max);
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

std::vector<float> flagella::alpha2vec(){
	std::vector<float> a_vec(2*n);
	for(int i = 0; i<2*n;++i)
		a_vec[i] = alpha[i];
	return a_vec;
}


float_type flagella::sumMass(int first, int last){
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

}

#endif // FlAG_HPP

