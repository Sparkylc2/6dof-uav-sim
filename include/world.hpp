#pragma once
#include "air_system.hpp"
#include "body_system.hpp"

namespace uav {

// describes the world the bodies live in
struct World {
    AirSystem atmosphere; // atmospheric model for wind and turbulence
    BodySystem bodies;    // dynamics model for the bodies
    double t = 0.0;

    // steps the atmosphere and the bodies
    void step(double dt);
};
} // namespace uav
