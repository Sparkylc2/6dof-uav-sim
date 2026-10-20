#pragma once
#include <Eigen/Dense>

#include "force_models.hpp"
#include "state.hpp"

namespace uav {
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

// a simulated body
// everything is in SI
struct Body {
    double m = 0.0;                // mass
    Vector4d I = Vector4d::Zero(); // moment of inertia (the four components left after the
                                   // symmetry assumptions, Ix, Iy, Iz, Ixz)

    State state; // carries the full 13 dim state

    // this will be overwritten by subclasses to implement different craft and their specific forces
    virtual Wrench wrench(const State &x) const = 0;
    virtual ~Body() = default;
};

// a basic drone, for now just gravity plus the forces and moments from the keyboard
// uses the zagi's mass properties for now (Table E.1 in appendix E)
struct Drone : Body {
    // the forces and moments the keyboard is asking for (body frame)
    Wrench input;

    Drone();

    // reads the keyboard into input, call once per frame before stepping
    void read_keyboard();

    // the keyboard interaction, body frame force and moment
    Vector3d keyboard_force() const;
    Vector3d keyboard_moment() const;

    Wrench wrench(const State &x) const override;
};
} // namespace uav
