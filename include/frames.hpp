// include/uav/frames.hpp
#pragma once
#include <Eigen/Dense>

namespace uav {
using Eigen::Matrix3d;
// returns a rotation matrix that rotates a coordinate system
// by the three angles (euler angles). basically any vector multiplied
// by this matrix is re-expressed in that rotated frame
//
// based on the content of chapter 2 (2.4 and 2.5 in particular)
Matrix3d rotation(double phi, double theta, double psi);

// rotation matrix bringing body angular velocity to the inertial frame of reference
// from eqn 3.16
Matrix3d body_rates_to_euler_rates(double phi, double theta);

} // namespace uav
