#include<iostream>
#include<array>
#include"flagella.hpp"
#include<cmath>



int main() {
	int time = 10000;
	int ele = 3;
	double l = 5;
	double m = 0.1;
	double k = 0.05;
	double b = 0.2;

	/*std::vector<double> l_fla, m_fla, k_fla, B_fla;

	for(int i = 0; i<ele;++i){
		l_fla.push_back(l);
		m_fla.push_back( m *(1.5-float(i)/ele));
		k_fla.push_back( k *(1.5-float(i)/ele));
		B_fla.push_back( b * (0.5+float(i)/ele));
	}

				// n,  l , m  ,		k ,	B,   x, y, ymove, omega
	flagella test(ele, l_fla, m_fla, k_fla, B_fla, 0, 0, 0, 0);
	*/

	flagella test(ele, l, m, k, b, 0,0,0,0);

	std::vector<double> Q = {0,0};
	for(int i = 0; i<ele;++i){
		test.InitAngle(i,0.3);
		//test.InitAngle(i, sin(float(i)/ele*2*M_PI)*0.1);
		//if(i%2 == 0) test.InitAngle(i, 0.1);
		//if(i%2 == 1) test.InitAngle(i,-0.001);
	}

	for(int i = 0; i< time;++i){
		//Q[0] = 1*sin(i*100/time);
		test.step(Q, 1);
	}
}
