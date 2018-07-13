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
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if (estimations.size() != ground_truth.size()){
        cout << "Error - Estimation and ground truth are not the same size" << endl;
        return rmse;
    }
    if (estimations.size() == 0){
        cout << "Error - estimation vector is empty" << endl;
        return rmse;
    }
    if (estimations[0].size() != ground_truth[0].size()) {
    	cout << "Error - estimation and ground_truth are not the same size" << endl;
    	return rmse;
    }
    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd resid = estimations[i] - ground_truth[i];
		resid = resid.array() * resid.array();
        rmse += resid;
	}

	//calculate the mean
	// ... your code here
    rmse /= estimations.size();
    
	//calculate the squared root
	// ... your code here
    rmse  = rmse.array().sqrt();

    return rmse;
    
}