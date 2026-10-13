#include "dynamics.hpp"
#include "frames.hpp"
#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>
#include <random>

static int g_checks = 0;
static int g_failures = 0;

// these i had claude generate for me
#define CHECK(cond)                                                                                \
    do {                                                                                           \
        ++g_checks;                                                                                \
        if (!(cond)) {                                                                             \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                            \
        }                                                                                          \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                      \
    do {                                                                                           \
        ++g_checks;                                                                                \
        const double a_ = (a), b_ = (b);                                                           \
        if (!(std::abs(a_ - b_) <= (tol))) {                                                       \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s = %.17g, expected %.17g\n", __FILE__, __LINE__, #a, a_,    \
                        b_);                                                                       \
        }                                                                                          \
    } while (0)

#define CHECK_VEC_NEAR(a, b, tol)                                                                  \
    do {                                                                                           \
        ++g_checks;                                                                                \
        const Eigen::Vector3d a_ = (a), b_ = (b);                                                  \
        if (!((a_ - b_).norm() <= (tol))) {                                                        \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s = [%.12g %.12g %.12g], expected [%.12g %.12g %.12g]\n",    \
                        __FILE__, __LINE__, #a, a_.x(), a_.y(), a_.z(), b_.x(), b_.y(), b_.z());   \
        }                                                                                          \
    } while (0)

constexpr double pi = std::numbers::pi;

// the full inertia matrix from the packed [Ix, Iy, Iz, Ixz]
// with Ixy = Iyz = 0 from symmetry about the x-z plane
static Eigen::Matrix3d inertia_matrix(const Eigen::Vector4d &I) {
    Eigen::Matrix3d J;
    J << I(0), 0, -I(3), 0, I(1), 0, -I(3), 0, I(2);
    return J;
}

static std::unique_ptr<uav::Body> make_zagi() {
    std::unique_ptr<uav::Body> b(uav::generate_zagi_flying_wing());
    b->p.setZero();
    b->r.setZero();
    b->v.setZero();
    b->w.setZero();
    return b;
}

// ---- frames ----
static void test_rotation_is_orthonormal() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> ang(-pi, pi);
    double worst_orth = 0.0, worst_det = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Eigen::Matrix3d R = uav::rotation(ang(rng), ang(rng) / 2.0, ang(rng));
        worst_orth = std::max(worst_orth, (R * R.transpose() - Eigen::Matrix3d::Identity()).norm());
        worst_det = std::max(worst_det, std::abs(R.determinant() - 1.0));
    }
    CHECK_NEAR(worst_orth, 0.0, 1e-12);
    CHECK_NEAR(worst_det, 0.0, 1e-12);
}

// rotation() is R_v^b (NED -> body) so the body axes in NED are the columns of the transpose
static void test_rotation_sign_conventions() {
    const Eigen::Vector3d nose = uav::rotation(0, 0, pi / 2).transpose().col(0);
    CHECK_VEC_NEAR(nose, Eigen::Vector3d(0, 1, 0), 1e-12); // psi = 90 deg means nose points east

    const Eigen::Vector3d nose_pitched = uav::rotation(0, pi / 6, 0).transpose().col(0);
    CHECK(nose_pitched.z() < 0.0); // theta > 0 means nose up (NED d is down)

    const Eigen::Vector3d right_wing = uav::rotation(pi / 6, 0, 0).transpose().col(1);
    CHECK(right_wing.z() > 0.0); // phi > 0 meaning the right wing down
}

// at phi = theta = 0 the Euler rates equal the body rates
static void test_euler_rates_identity_when_level() {
    const Eigen::Matrix3d E = uav::body_rates_to_euler_rates(0.0, 0.0);
    CHECK_NEAR((E - Eigen::Matrix3d::Identity()).norm(), 0.0, 1e-15);
}

// ---- dynamics ----
// if u have zero body rates, v_dot = f / m
static void test_constant_force_gives_linear_velocity() {
    auto b = make_zagi();
    const double a = 2.0;
    b->get_forces = [a](uav::Body *body, double) { return Eigen::Vector3d(body->m * a, 0, 0); };
    uav::BodySystem sys;
    sys.bodies.push_back(b.get());

    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i)
        sys.step(dt);

    CHECK_VEC_NEAR(b->v, Eigen::Vector3d(a * 1.0, 0, 0), 1e-12); // u = a t after 1 s
    CHECK(b->p.x() > 0.0);                                       // meaning it moved north
    CHECK_NEAR(b->p.y(), 0.0, 1e-15);
    CHECK_NEAR(b->p.z(), 0.0, 1e-15);
}

// position rates are R_b^v [u v w]^T, and its yawed 90 deg and flying forward
// along the body x axis, the craft must move east
static void test_yawed_forward_flight_moves_east() {
    auto b = make_zagi();
    b->r = Eigen::Vector3d(0, 0, pi / 2);
    b->v = Eigen::Vector3d(10.0, 0, 0);
    uav::BodySystem sys;
    sys.bodies.push_back(b.get());

    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i)
        sys.step(dt);

    CHECK_VEC_NEAR(b->p, Eigen::Vector3d(0, 10.0, 0), 1e-9); // 10 m east after 1 s
}

// if you start from rest the gyro terms should vanish, so w_dot = I^-1*tau.
// this checks gamma, gamma_3, gamma_4, gamma_8 (eqn. 3.13) and the 1/Iy term against
// I.inverse()
static void test_angular_accel_from_rest_matches_Jinv() {
    auto b = make_zagi();
    const Eigen::Vector3d moment(0.3, -0.2, 0.1);
    b->get_moments = [moment](uav::Body *, double) { return moment; };
    uav::BodySystem sys;
    sys.bodies.push_back(b.get());

    const double dt = 1e-3;
    sys.step(dt);

    const Eigen::Vector3d expected = dt * inertia_matrix(b->I).inverse() * moment;
    CHECK_VEC_NEAR(b->w, expected, 1e-12);
}

// eqn 3.10 is I * w_dot = m - w_dot x (I*w).
// with no applied moment an euler step should equal w0 + dt * I^-1 (-w0 x I*w0).
// this checks gamma_1, gamma_2, gamma_5, gamma_6, gamma_7 (eqn. 3.13)
static void test_gyroscopic_terms_match_vector_form() {
    auto b = make_zagi();
    const Eigen::Vector3d w0(0.5, -0.3, 0.8);
    b->w = w0;
    uav::BodySystem sys;
    sys.bodies.push_back(b.get());

    const double dt = 1e-3;
    sys.step(dt);

    const Eigen::Matrix3d J = inertia_matrix(b->I);
    const Eigen::Vector3d expected = w0 + dt * J.inverse() * (-w0.cross(J * w0));
    CHECK_VEC_NEAR(b->w, expected, 1e-12);
}

// with Ixz 0 a pure roll moment shouldn't give any yaw rates
static void test_roll_yaw_coupling() {
    for (const double Ixz : {0.0, 0.0015}) {
        auto b = make_zagi();
        b->I.w() = Ixz;
        b->get_moments = [](uav::Body *, double) { return Eigen::Vector3d(0.1, 0, 0); };
        uav::BodySystem sys;
        sys.bodies.push_back(b.get());

        for (int i = 0; i < 100; ++i)
            sys.step(1e-3);

        if (Ixz == 0.0)
            CHECK(b->w.z() == 0.0);
        else
            CHECK(std::abs(b->w.z()) > 1e-6);
    }
}

int main() {
    test_rotation_is_orthonormal();
    test_rotation_sign_conventions();
    test_euler_rates_identity_when_level();
    test_constant_force_gives_linear_velocity();
    test_yawed_forward_flight_moves_east();
    test_angular_accel_from_rest_matches_Jinv();
    test_gyroscopic_terms_match_vector_form();
    test_roll_yaw_coupling();

    std::printf("%d/%d passed\n (pray its all right xD)", g_checks - g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
