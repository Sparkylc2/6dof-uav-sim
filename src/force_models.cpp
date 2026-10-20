#include "force_models.hpp"
#include "frames.hpp"

namespace uav {

Wrench gravity(double m, const Vector4d &e) {
    Wrench w;

    Eigen::Vector3d g = {0, 0, m * 9.81};
    // gravity is along +d in NED, rotation(e) is body -> vehicle so transpose to go vehicle ->
    // body
    Eigen::Matrix3d R_grav = rotation(e).transpose();
    w.F = R_grav * g;

    return w;
}

} // namespace uav
