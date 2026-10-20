// src/main.cpp
#include "body.hpp"
#include "body_system.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "state.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <raylib.h>

using namespace uav;
int main() {
    bool first_person = false;

    uav::rendering::init_rendering();
    uav::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);
    Mesh mesh = GenMeshCube(0.25, 0.25, 0.25f);

    Model model = LoadModelFromMesh(mesh);

    uav::BodySystem sys;
    auto drone = std::make_unique<uav::Drone>();
    uav::Drone &body = *drone;
    sys.bodies.push_back(std::move(drone));

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            body.state = uav::State();
        }
        if (IsKeyPressed(KEY_F)) {
            first_person = !first_person;
        }
        body.read_keyboard();
        const double dt = GetFrameTime();
        sys.step(dt);
        const Eigen::Matrix3d R_b2i = rotation(body.state.quat());

        Matrix transform = rendering::to_rl_transform(R_b2i, body.state.pos());
        model.transform = transform;

        if (first_person) {
            rendering::CAMERA.position = rendering::ned_to_rl(body.state.pos());
            rendering::CAMERA.target =
                rendering::ned_to_rl(body.state.pos() + R_b2i * Vector3d(1.0, 0.0, 0.0));
            rendering::CAMERA.up = rendering::ned_to_rl(R_b2i * Vector3d(0.0, 0.0, -1.0));
            rendering::CAMERA.fovy = 90.0f;
        } else {
            uav::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);
        }

        uav::rendering::begin();
        DrawModel(model, {0, 0, 0}, 1.0f, GRAY);
        uav::rendering::end();
    }
    CloseWindow();
}
