#include "fusion.hpp"
#include "lib/Eigen/Dense"
#include <mutex>
#include <time.h>
#include "paramReader.hpp"

// Read these parameters from the file
struct EKFParams
{
    double varGPS, varSpeed, varYaw, varAcc, maxAccel, maxTurnRate, maxYawAccel, xOff, yOff;
};

inline static EKFParams getDefaultParams()
{
    ParameterReader pd;
    EKFParams params;
    params.varGPS = atof(pd.getData( "vargps" ).c_str());
    params.varSpeed = atof(pd.getData( "varspeed" ).c_str());
    params.varYaw = atof(pd.getData( "varyaw" ).c_str());
    params.varAcc = atof(pd.getData( "varaccel" ).c_str());
    params.maxAccel = atof(pd.getData( "maxaccel" ).c_str());
    params.maxTurnRate = atof(pd.getData( "maxturnrate" ).c_str());
    params.maxYawAccel = atof(pd.getData( "maxyawaccel" ).c_str());
    params.xOff = atof(pd.getData("xOff").c_str());
    params.yOff = atof(pd.getData("yOff").c_str());
    return params;
}

/**
 * @brief Class with the highest level of abstraction for the GPS INS estimation system. Uses D      ataPoint and Fusion
 * classes to perform filtering and state estimation on the robot.
 */
class GpsIns
{
    public:
        GpsIns(bool verbose = false);
        ~GpsIns();
         void read_gps_data(double lat, double lon, double alt);
         void read_imu_data(double yaw_rate, double long_accel);
         void read_encoders(double vel);
         void loop();
         void set_data();
         Eigen::VectorXd get_estimated_state() const;
         DataPoint get_sensor_data();
     private:
         Fusion* filter;
         DataPoint* _sensor_data;
	    EKFParams params;
         clock_t _dt;
	    clock_t _cur_time;
	    clock_t _prev_time;

         int _imucounter, _gpscounter, _enccounter;
         int _prev_imu_counter, _prev_gps_counter, _prev_enc_counter;
         std::mutex m;
         DataPointType _data_type;
         Eigen::VectorXd _raw_data;
	    bool verbose;
};
