#include<iostream>
#include<array>
#include"flagella.hpp"

int main() {
	flagella test(1, 1, 0.05, 0.05, 0, 0, 0, 0);

	for(int i = 0; i< 1000;++i){
		test.step(0,0, 0.05);
	}
}
