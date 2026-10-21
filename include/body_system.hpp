#pragma once
#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "body.hpp"
#include "body_state.hpp"
#include "force_models.hpp"

namespace uav {
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

// stores all the simulated bodies and does the physics updates
struct BodySystem {
    //  uses pointers so references to a body don't invalidate when more are added
    std::vector<std::unique_ptr<Body>> bodies;

    // steps all bodies in the system forward
    void step(double dt);

    // computes the x_dot state vector
    StateVector derivatives(const Body &b, const BodyState &x, const Vector3d &F,
                            const Vector3d &T);
};

} // namespace uav
