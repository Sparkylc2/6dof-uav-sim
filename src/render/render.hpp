#pragma once
#include <Eigen/Dense>
#include <raylib.h>

namespace uav::rendering {

extern Camera3D CAMERA; // global camera

// raylib uses a different coord convention (y is up)
inline Vector3 ned_to_rl(const Eigen::Vector3d &p) {
    return {static_cast<float>(p.y()), static_cast<float>(-p.z()), static_cast<float>(-p.x())};
}

// wraps the intial window sizing and framerate
void init_rendering();

// initializes the global camera
void init_camera(const Eigen::Vector3d &pos, const Eigen::Vector3d &target,
                 const Eigen::Vector3d &up, double fovy);

// begin and end for rendering, wraps BeginMode, clearing background etc.
void begin();
void end();

} // namespace uav::rendering
