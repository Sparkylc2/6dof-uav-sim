// include/uav/frames.hpp
#pragma once
#include <Eigen/Dense>

namespace uav {
// returns a rotation matrix that rotates a coordinate system
// by the three angles (euler angles). basically any vector multiplied
// by this matrix is re-expressed in that rotated frame
//
// based on the content of chapter 2 (2.4 and 2.5 in particular)
Eigen::Matrix3d rotation(double phi, double theta, double psi);

} // namespace uav
