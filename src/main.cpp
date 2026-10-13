// src/main.cpp
#include "dynamics.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <raylib.h>

using namespace uav;
int main() {
    bool first_person = false;

    uav::rendering::init_rendering();
    uav::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);
    Mesh mesh = GenMeshCube(0.5f, 0.5f, 0.5f);
    Model model = LoadModelFromMesh(mesh);

    uav::Body *body = generate_zagi_flying_wing();

    auto forces = [](uav::Body *b, double) -> Eigen::Vector3d {
        Eigen::Vector3d g = {0, 0, b->m * 9.81};
        Eigen::Matrix3d R_grav = rotation(b->r.x(), b->r.y(), b->r.z());
        Vector3d f = R_grav * g;

        double mult = 2.0 * b->m * 9.81;
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
    };

    auto moments = [](uav::Body *, double) -> Eigen::Vector3d {
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
    };

    body->get_forces = forces;
    body->get_moments = moments;

    uav::BodySystem sys;
    sys.bodies.push_back(body);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            body->reset();
        }
        if (IsKeyPressed(KEY_F)) {
            first_person = !first_person;
        }
        const double dt = GetFrameTime();
        sys.step(dt);
        const Eigen::Matrix3d R_b2i =
            uav::rotation(body->r.x(), body->r.y(), body->r.z()).transpose();

        Matrix transform = rendering::to_rl_transform(R_b2i, body->p);
        model.transform = transform;

        if (first_person) {
            rendering::CAMERA.position = rendering::ned_to_rl(body->p);
            rendering::CAMERA.target =
                rendering::ned_to_rl(body->p + R_b2i * Vector3d(1.0, 0.0, 0.0));
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
