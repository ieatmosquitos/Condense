#include <vector>
#include <math.h>
#include "Eigen/Core"
#include "Eigen/SVD"

#ifndef _STABILITY_H
#define _STABILITY_H


// every function is a function in the form ax + by = c 
class linSystem{
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  
 public:
  void addConstraint(double x, double y, double theta);
  void addConstraint(double x, double y, double robot_theta, double bearing);
  
  bool checkStability();
};

// we here take advantage of the correspondance with the form    y = mx + q    , where m is the angular coefficient of the line, and q is a constant
void linSystem::addConstraint(double x, double y, double theta){
  double m = tan(theta);
  double q = y - m*x;
  
  a.push_back(-m);
  b.push_back(1);
  c.push_back(q);
}


void linSystem::addConstraint(double x, double y, double robot_theta, double bearing){
  double theta = robot_theta + bearing;
  if(theta > M_PI){
    theta = -M_PI + (theta - M_PI);
  }
  else if(theta <= -M_PI){
    theta = M_PI - (theta + M_PI);
  }
  
  this->addConstraint(x, y, theta);
}


bool linSystem::checkStability(){
  Eigen::Matrix<double, Eigen::Dynamic, 2> A;
  A.resize(a.size(), 2);
  
  for(unsigned int i = 0; i<a.size(); i++){
    A(i,0) = a[i];
    A(i,1) = b[i];
  }
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
  Eigen::MatrixXd eigenvalues = svd.singularValues();
  
  if(eigenvalues(1,0) < 10e-9){
    return false;
  }
  
  if(eigenvalues(0,0) / eigenvalues(1,0) > 7){
    return false;
  }
  
  return true;
}



#endif // _STABILITY_H
