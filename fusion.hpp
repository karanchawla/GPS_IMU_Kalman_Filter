//
//  fusion.hpp
//  Fusion
//
//  Created by Karan on 4/9/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#ifndef fusion_hpp
#define fusion_hpp

#include <stdio.h>
#include "Eigen/Dense"
#include "ekf.hpp"
#include "geo_ned.hpp"


using namespace geodectic_converter;

class Fusion
{
public:
    const int _n = 6;
    bool _initialized;
    long long timestamp;
    Eigen::MatrixXd _P;
    Eigen::MatrixXd _Q;
    EKF _KF;
    
    
    
private:
    
    
    
    
};

#endif /* fusion_hpp */
