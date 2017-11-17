#include<iostream>
#include<array>

int main() {
const std::array<std::array<int, 9>, 2> c = {{{{       0,       1,       0,      -1,       0,       1,      -1,      -1,       1}}, 
	                                              {{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}}};

for (int i=0; i<9;++i)
	std::cout << c[0][i] << " " ;
std::cout<<std::endl;
for (int i=0; i<9;++i)
	std::cout << c[1][i] << " " ;
return 0;
}
