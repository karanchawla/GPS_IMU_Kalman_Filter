# ExtendedKalmanFilter
EKF to fuse GPS, IMU and encoder readings to estimate the pose of a ground robot in the navigation frame. 


Wikipedia writes: In the extended Kalman filter, the state transition and observation models need not be linear functions of the state but may instead be differentiable functions.

```
x_k = g(x_k), u_k-1 + w_k-1
z_k = h(x_k) + v_k
```

![EKF step](https://github.com/karanchawla/GPS_INS_Fusion/blob/master/figure/Extended-Kalman-Filter-Step.png)

Where `w_k` and `v_k` are the process and observation noises which are both assumed to be zero mean Multivariate Gaussian noises with covariance matrix `Q` and `R` respectively.

The function `g` can be used to compute the predicted state from the previous estimate and similarly the function `h` can be used to compute the predicted measurement from the predicted state. However, `g` and `h` cannot be applied to the covariance directly. Instead a matrix of partial derivatives (the Jacobian matrix) is computed.

At each time step, the Jacobian is evaluated with current predicted states. These matrices can be used in the Kalman filter equations. This process essentially linearizes the non-linear function around the current estimate.

Here we have a velocity sensor (encoders/GPS velocity), which measures the vehicle speed (`v`) in heading direction (`psi`), a yaw rate sensor (`psi_dot`) and an accelerometer which measures longitudinal velocity which both have to fused with the position (`x` & `y`) from the GPS sensor.

![Ground robot model](https://github.com/karanchawla/GPS_INS_Fusion/blob/master/figure/CTRV-Model.png)

## References
1. @balzer82 for his tutorials on Kalman Filters.
2. Probabilistic Robotics by Thrun, Burgard, and Fox.
