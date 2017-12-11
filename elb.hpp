// implementation of entropic lattice boltzmann
//function for calculating the entropic equilibrium using multivariate newton rhapson
//function for calculating the alpha for entropic relaxation

#ifndef ELB_HPP
#define ELB_HPP

#include"velocity_set.hpp"
#include "lattice.hpp"
//#include<Eigen/Dense>
#include<iostream>
#include<cmath>
using namespace lb;
//function to get the entropic lattice boltzmann coefficients A Bx and By
//and store them in the array whose pointer has been passed as an argument

/*void get_elbeq_coeffs(double* x, const double& p, const double& pux, const double& puy) {
    double Eabs=1e-14, Erel=1e-10;
    Eigen::VectorXd X(3); X.setZero();
    Eigen::VectorXd delX(3);
    Eigen::VectorXd b(3);
    Eigen::VectorXd f(9);
    Eigen::MatrixXd LDLT(3,3);

    do {
    //populating the populations vector
    for (int i=0 ; i<9 ; ++i)
      f(i) = velocity_set().W[i] * std::exp(X(0) + X(1)*velocity_set().c[0][i] + X(2)*velocity_set().c[1][i]);

    //creating the LDLT matrix
    double f00 = 0., f01=0. , f02= 0.,f11=0.,f12=0.,f22=0.;
    for (int i=0 ; i<9 ; ++i) {
      f00+=f(i);
      f01+=f(i)*velocity_set().c[0][i];
      f02+=f(i)*velocity_set().c[1][i];
      f11+=f(i)*velocity_set().c[0][i]*velocity_set().c[0][i];
      f12+=f(i)*velocity_set().c[0][i]*velocity_set().c[1][i];
      f22+=f(i)*velocity_set().c[1][i]*velocity_set().c[1][i];
    }
    LDLT(0,0)= f00; LDLT(0,1)= f01; LDLT(0,2)=f02;
    LDLT(1,0) = f01 ; LDLT(1,1)= f11; LDLT(1,2)=f12;
    LDLT(2,0)=f02; LDLT(2,1)= f12; LDLT(2,2)= f22;

    //constructing the b vector
    b(0) = f00-p ; b(1) = f01-pux ; b(2) = f02-puy;

    delX = LDLT.lu().solve(b);
    X = X - delX;
    } while (delX.norm() > Eabs + X.norm() * Erel);

    //convergence loop exited
    for (int i=0 ; i<3 ; ++i)
      x[i] = X(i);
}
*/

//function for finding alpha
double get_alpha(const node& n, double* feq) {
  //std::cout << "Test 1!" <<std::endl;
  double Erel=1e-3;
  int stupidcount =0;
  //double* X=new double[3];
  double alphamax = 10;
  //get_elbeq_coeffs(X,n.rho(),n.rho()*n.u(),n.rho()*n.v());
  /*Eigen::VectorXd feq(9);*/ std::vector<double> delta(9);
  for (int i=0 ; i<9 ; ++i) {
  //  feq(i) = velocity_set().W[i] * std::exp(X[0] + X[1]*velocity_set().c[0][i] + X[2]*velocity_set().c[1][i]);
    delta[i] = feq[i] - n.f(i);
    if (delta[i]<0)
    {
      double alphaimax = fabs(n.f(i) / delta[i]);
      if (alphaimax < alphamax)
        alphamax = alphaimax;
    }
    //if(fabs(delta[i]/n.f(i)) < Erel)
    //if(fabs(delta[i]/feq[i]) < Erel)
    if(fabs(delta[i]) < 1e-3) //condition for being cose to equilibrium 1e-3 working
      ++stupidcount;
  }
  if (stupidcount == 9)
    return 2.;
  if (alphamax <= 2 )
    return alphamax;

//const auto alphamax2 = alphamax;
    double alpha=2;
  //else get the value by newton rhapson
  double Hf=0.;
  for (int i=0 ; i<9 ; ++i){
    Hf+= n.f(i)*std::log(n.f(i)/velocity_set().W[i]);
  }
  bool flag = false; double troublesaver = 2;
  int steps = 1;
  //double s=1.0;
  do {
    ++steps;
    alphamax = alpha;
    double Haf=0.,Hpr=0.;
    for (int i=0 ; i<9 ; ++i){
      Haf+= (n.f(i) + alpha*delta[i])*std::log( (n.f(i) + alpha*delta[i])/velocity_set().W[i]);
      Hpr+= delta[i]*std::log( (n.f(i) + alpha*delta[i])/velocity_set().W[i]);
    }
  //  std::cout << "H prime value : " << Hpr << std::endl;
    if (Hpr < 1e-8) { //condition for smallness of h prime
      //if (flag == true )
        //troublesaver = ;
      flag = true;
      std::cout << "H too low! " << std::endl;
      //alpha=troublesaver + std::pow(-1/2,steps);
      break;
      //continue;
    }
    if(steps > 12){
      flag = true;
      break;
    }
    alpha = alpha - /*s*/(Haf-Hf)/Hpr ;

    //if (alpha > alphamax2)
  //  {
  //    alpha = 0.98*alphamax2;
  //    s*=0.99;
  //  }
    //std::cout << "Alpha: " << alpha << "Error: " << fabs((alpha-alphamax)/alphamax) <<std::endl;
  //  std::cout << "In the newton rhapson loop!" << std::endl;
  } while( fabs((alpha-alphamax)/alphamax) > Erel);
  if (flag)
    return 2.;
//std::cout << "Alpha calculated by newton rhapson" << std::endl;
//  if (alpha < 10 && alpha > 0)
  return alpha;
//  else
//  return 10.;
}



#endif //for the header
