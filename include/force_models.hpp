#pragma once
#include <Eigen/Dense>
namespace uav {
using Eigen::Vector3d;
using Eigen::Vector4d;

// a force and moment acting on a body, both in the body frame
struct Wrench {
    Vector3d F = Vector3d::Zero(); // force [fx, fy, fz]
    Vector3d M = Vector3d::Zero(); // moment [l, m, n]

    Wrench &operator+=(const Wrench &o) {
        F += o.F;
        M += o.M;
        return *this;
    }
};

// gravity on a body of mass m with attitude e (F3 / Q11)
// acts at the center of mass so there's no moment
Wrench gravity(double m, const Vector4d &e);

} // namespace uav
