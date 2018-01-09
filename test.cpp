#include<iostream>
#include<array>
#include"flagella.hpp"

int main() {
	flagella test(5, 10, 0.05, 0.05, 0.005, 0.005, 0, 0, 0, 0);

	for(int i = 0; i< 10000;++i){
		test.step(0,0, 1.0);
	}
}
