#include "render/render.hpp"

namespace uav::rendering {

Camera3D CAMERA = {};

void init_rendering() {
    InitWindow(1280, 720, "uav_sim");
    SetTargetFPS(60);
}

void init_camera(const Eigen::Vector3d &pos, const Eigen::Vector3d &target,
                 const Eigen::Vector3d &up, double fovy) {
    CAMERA.position = ned_to_rl(pos);
    CAMERA.target = ned_to_rl(target);
    CAMERA.up = ned_to_rl(up);
    CAMERA.fovy = static_cast<float>(fovy);
    CAMERA.projection = CAMERA_PERSPECTIVE;
}

void begin() {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    BeginMode3D(CAMERA);
    DrawGrid(200, 5.0f); // 1 km across so there's still ground to see at flying speed
}

void end() {
    EndMode3D();
    DrawFPS(10, 10);
    EndDrawing();
}

} // namespace uav::rendering
