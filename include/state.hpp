#pragma once
#include <Eigen/Dense>

namespace uav {
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

using StateVector = Eigen::Matrix<double, 13, 1>;

enum IDX {
    POS = 0,
    VEL = 3,
    QUAT = 6,
    RATES = 10

};

struct State {
    // state layout is the same as B.6
    // [pn pe pd | u v w | e0 e1 e2 e3 | p q r]
    StateVector state = StateVector::Zero();

    State() {
        state.setZero();
        quat() << 1.0, 0.0, 0.0, 0.0;
    }

    // returns a mutable view into the state vector
    // position (inertial FoR) [x, y, z]
    Ref<Vector3d> pos() { return state.segment<3>(IDX::POS); }
    // velocity (ground velocity in body FoR) [u, v, w]
    Ref<Vector3d> vel() { return state.segment<3>(IDX::VEL); }
    // attitude of body relative to vehicle frame [e0, e1, e2, e3]
    Ref<Vector4d> quat() { return state.segment<4>(IDX::QUAT); }
    // angular velocity (body FoR) [p, q, r]
    Ref<Vector3d> rates() { return state.segment<3>(IDX::RATES); }

    // const versions
    Ref<const Vector3d> pos() const { return state.segment<3>(0); }
    Ref<const Vector3d> vel() const { return state.segment<3>(3); }
    Ref<const Vector4d> quat() const { return state.segment<4>(6); }
    Ref<const Vector3d> rates() const { return state.segment<3>(10); }
};
} // namespace uav
