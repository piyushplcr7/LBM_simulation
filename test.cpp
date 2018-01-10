#include<iostream>
#include<array>
#include"flagella.hpp"

int main() {
	int time = 10000;
	int ele = 2;
	flagella test(ele, 5, 10, 0.05, 0.05, 0, 0, 0, 0);

	std::vector<float> Q = {0,0};
	for(int i = 0; i<ele;++i)
		test.InitAngle(i, 3.14/2);
	
	for(int i = 0; i< time;++i){
		//Q[0] = 1*sin(i*100/time);
		test.step(Q, 1);
	}
}
