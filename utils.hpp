//
//  utils.hpp
//  ExtendedKalmanFilter
//
//  Created by Karan on 4/6/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include "Eigen/Dense"

Eigen::MatrixXd calculate_joacobian(const Eigen::VectorXd& v, const double dt);

#endif /* utils_hpp */
