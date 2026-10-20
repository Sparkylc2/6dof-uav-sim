#include "body.hpp"
#include <raylib.h>

namespace uav {

Drone::Drone() {
    m = 1.56;

    const double Ix = 0.1147, Iy = 0.0576, Iz = 0.1712, Ixz = 0.0015;

    // following the symmetry assumption that airplanes use
    I << Ix, Iy, Iz, Ixz;
}

void Drone::read_keyboard() {
    input.F = keyboard_force();
    input.M = keyboard_moment();
}

Vector3d Drone::keyboard_force() const {
    Vector3d f = Vector3d::Zero();

    double mult = 2.0 * m * 9.81;
    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }

    if (IsKeyDown(KEY_SPACE)) {
        f += Vector3d(0.0, 0.0, -mult);
    }
    if (IsKeyDown(KEY_LEFT_SHIFT)) {
        f += Vector3d(0.0, 0.0, mult);
    }

    return f;
}

Vector3d Drone::keyboard_moment() const {
    Vector3d T = Vector3d::Zero();

    double mult = 0.5;

    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }
    if (IsKeyDown(KEY_W)) {
        T.y() -= mult;
    }
    if (IsKeyDown(KEY_S)) {
        T.y() += mult;
    }
    if (IsKeyDown(KEY_A)) {
        T.x() -= mult;
    }
    if (IsKeyDown(KEY_D)) {
        T.x() += mult;
    }
    if (IsKeyDown(KEY_Q)) {
        T.z() -= mult;
    }
    if (IsKeyDown(KEY_E)) {
        T.z() += mult;
    }

    return T;
}

Wrench Drone::wrench(const State &x) const {
    Wrench w = gravity(m, x.quat());
    w += input;
    return w;
}

} // namespace uav
