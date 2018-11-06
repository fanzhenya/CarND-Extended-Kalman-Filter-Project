#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd ret(4);
  ret << 0, 0, 0, 0;
  if (estimations.size() == 0 || ground_truth.size() != estimations.size()) return ret;
  
  for (size_t i = 0; i < estimations.size(); ++i) {
    VectorXd t = estimations[i] - ground_truth[i];
    t = t.array() * t.array();
    ret = ret + t;
  }
  ret /= estimations.size();
  return ret.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0), py = x_state(1), vx = x_state(2), vy = x_state(3);
  float t1 = px * px + py * py, t2 = sqrt(t1), t3 = t1 * t2;

  //check division by zero
  if(fabs(t1) < 0.0001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/t2), (py/t2), 0, 0,
       -(py/t1), (px/t1), 0, 0,
       py*(vx * py - vy * px) / t3, px * (px * vy - py * vx) / t3, px / t2, py / t2;
  return Hj;
}
