//
//  main.cpp
//  ExtendedKalmanFilter
//
//  Created by Karan on 4/6/18.
//  Copyright © 2018 Karan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "fusion.hpp"

int main(int argc, const char * argv[])
{
    std::ifstream ip("/Users/karan/Desktop/EKF/EKF/Data.csv");
    if(!ip.is_open())
    {
        std::cerr << "Failed to open the data file";
        std::exit(EXIT_FAILURE);
    }
    
    std::string timestamp;
    std::string ax;
    std::string yaw_rate;
    std::string yaw;
    std::string course;
    std::string lat;
    std::string lon;
    
    std::string str;
    std::getline(ip, str); // Skip the first line
    
    while(std::getline(ip, str))
    {
        std::istringstream iss(str);
        std::string token;
        
        while (std::getline(iss, token, ','))
        {
            // process each token
            std::cout << token.size() << " ";
        }
        std::cout << std::endl;
    }
    
    ip.close();
    
    return 0;
}
