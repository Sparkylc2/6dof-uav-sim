#include "uav/frames.hpp"
#include <cmath>

namespace uav {

Eigen::Matrix3d rotation(double phi, double theta, double psi) {
    const double c_phi = std::cos(phi), s_phi = std::sin(phi);
    const double c_theta = std::cos(theta), s_theta = std::sin(theta);
    const double c_psi = std::cos(psi), s_psi = std::sin(psi);

    Eigen::Matrix3d R;
    // this comes from multiplying the three rotation matrices
    // for the three different rotations
    R << c_theta * c_psi, c_theta * s_psi, -s_theta, s_phi * s_theta * c_psi - c_phi * s_psi,
        s_phi * s_theta * s_psi + c_phi * c_psi, s_phi * c_theta,
        c_phi * s_theta * c_psi + s_phi * s_psi, c_phi * s_theta * s_psi - s_phi * c_psi,
        c_phi * c_theta;

    return R;
}

} // namespace uav
