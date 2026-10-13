#pragma once
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace uav {
using Eigen::Vector3d;
using Eigen::Vector4d;
// stores the simulated bodies state
// everything is in SI
struct Body {
    double m = 0.0;                // mass
    Vector4d I = Vector4d::Zero(); // moment of inertia (the four components left after assumptions,
                                   // Ix, Iy, Iz, Ixz)

    Vector3d p = Vector3d::Zero(); // position (inertial FoR) [x, y, z]
    Vector3d r = Vector3d::Zero(); // rotation (inertial FoR) [phi, theta, psi]
    Vector3d w = Vector3d::Zero(); // rotational velocity (body FoR) [p, q, r]
    Vector3d v = Vector3d::Zero(); // velocity (ground velocity in body FoR) [u, v, w]

    Vector3d F = Vector3d::Zero(); // forces on the body (gets cleared after each step)
    Vector3d M = Vector3d::Zero(); // moments on the body (also gets cleared after each step)

    // easier than some complicated system so i can just set it in main
    std::function<Vector3d(Body *, double)> get_forces;  // body, dt
    std::function<Vector3d(Body *, double)> get_moments; // body, dt

    void reset() {
        p = Vector3d::Zero();
        r = Vector3d::Zero();
        v = Vector3d::Zero();
        w = Vector3d::Zero();
    }
};

// stores all the simulated bodies and does the physics updates
struct BodySystem {
    std::vector<Body *> bodies;

    // steps the simulation using the simplified eqn. discussed in chapter 3
    void step(double dt);

    // steps the simulation just using the raw vector operations rather than an algebraic result
    void step_full(double dt);
};

// generates a body with the same characteristics as the zagi flying wing
// everything is from Table E.1 in appendix E
Body *generate_zagi_flying_wing();

} // namespace uav
