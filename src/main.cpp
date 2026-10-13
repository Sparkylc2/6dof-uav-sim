// src/main.cpp
#include "render/render.hpp"
#include "uav/frames.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <raylib.h>

int main() {
    uav::rendering::init_rendering();

    uav::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);

    double t = 0.0;
    while (!WindowShouldClose()) {
        t += GetFrameTime();

        const Eigen::Matrix3d R = uav::rotation(0.5 * std::sin(t), 0.0, t).transpose();
        const Eigen::Vector3d o = Eigen::Vector3d::Zero();

        uav::rendering::begin();
        DrawLine3D(uav::rendering::ned_to_rl(o), uav::rendering::ned_to_rl(R.col(0) * 3.0),
                   RED); // x_b (nose)
        DrawLine3D(uav::rendering::ned_to_rl(o), uav::rendering::ned_to_rl(R.col(1) * 3.0),
                   GREEN); // y_b (right wing)
        DrawLine3D(uav::rendering::ned_to_rl(o), uav::rendering::ned_to_rl(R.col(2) * 3.0),
                   BLUE); // z_b (down)
        uav::rendering::end();
    }
    CloseWindow();
}
