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
  VectorXd sum_squared = VectorXd(4);
  sum_squared << 0, 0, 0, 0;
  
  for (int i=0; i<estimations.size(); ++i) {
	  sum_squared = sum_squared.array() + (estimations[i] - ground_truth[i]).array().square();
  }

  VectorXd rmse = sum_squared / estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}