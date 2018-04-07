//
//  tools.hpp
//  ExtendedKalmanFilter
//
//  Created by Karan on 4/6/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef tools_hpp
#define tools_hpp

#include <stdio.h>
#include "Eigen/Dense"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>

const float RADIUS_OF_EARTH = 6378388.0;

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd calculate_joacobian(const VectorXd& v);
VectorXd convert_NED_to_ECEF()


#endif /* tools_hpp */
