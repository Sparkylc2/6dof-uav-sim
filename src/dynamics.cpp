#include "dynamics.hpp"
#include "frames.hpp"

namespace uav {

Body *generate_zagi_flying_wing() {
    Body *b = new Body();
    b->m = 1.56;

    const double Ix = 0.1147, Iy = 0.0576, Iz = 0.1712, Ixz = 0.0015;

    // following the symmetry assumption that airplanes use
    b->I << Ix, Iy, Iz, Ixz;

    // i'll add the other coefficients for the flight model later (if i make it that far xD)
    return b;
}

void BodySystem::step(double dt) {
    for (Body *b : bodies) {

        // the forces on the body
        Vector3d F = Vector3d::Zero(), T = Vector3d::Zero();
        if (b->get_forces) {
            F = b->get_forces(b, dt);
        }
        if (b->get_moments) {
            T = b->get_moments(b, dt);
        }

        // velocity
        Vector3d vel_body = b->v;
        // angular velocity
        Vector3d w_vel_body = b->w;

        // the velocity components
        const double u = vel_body.x(), v = vel_body.y(), w = vel_body.z();
        const double p = w_vel_body.x(), q = w_vel_body.y(), r = w_vel_body.z();

        // the moments of inertia we need
        const double Ix = b->I.x(), Iy = b->I.y(), Iz = b->I.z();
        const double Ixz = b->I.w();

        // the moments
        const double L = T.x();
        const double M = T.y();
        const double N = T.z();

        // current rotation relative to inertial
        const double phi = b->r.x();
        const double theta = b->r.y();
        const double psi = b->r.z();

        // factors in the final expression (3.13) for angular accel
        const double g = 1.0 / (Ix * Iz - Ixz * Ixz);
        const double g1 = g * Ixz * (Ix - Iy + Iz);
        const double g2 = g * (Iz * Iz - Iz * Iy + Ixz * Ixz);
        const double g3 = g * Iz;
        const double g4 = g * Ixz;
        const double g5 = (Iz - Ix) / Iy;
        const double g6 = Ixz / Iy;
        const double g7 = g * (Ix * Ix - Iy * Ix + Ixz * Ixz);
        const double g8 = g * Ix;

        // --- accel ---
        Vector3d v_dot_body;

        // from eqn 3.15 chapter 3
        v_dot_body.x() = r * v - q * w;
        v_dot_body.y() = p * w - r * u;
        v_dot_body.z() = q * u - p * v;

        // adding the force term
        v_dot_body += 1 / b->m * F;

        // --- rotational accel ---
        Vector3d w_dot_body;

        // from eqn 3.17 chapter 3
        w_dot_body.x() = g1 * p * q - g2 * q * r;
        w_dot_body.y() = g5 * p * r - g6 * (p * p - r * r);
        w_dot_body.z() = g7 * p * q - g1 * q * r;

        w_dot_body.x() += g3 * L + g4 * N;
        w_dot_body.y() += M / Iy;
        w_dot_body.z() += g4 * L + g8 * N;

        // the following can most definitely be done with a formulation with one matrix but for the
        // sake of simplicity i'm splitting it up

        // once we have v_dot and w_dot we need to(euler for now since no forces are
        // dependent on current state [for now])
        //
        // both in body frame
        b->v += v_dot_body * dt;
        b->w += w_dot_body * dt;

        // using eqn 3.14 you can take the body velocities and reframe them in the inertial
        // transpose since we are going the other way
        const Vector3d vel_inertial = rotation(phi, theta, psi).transpose() * b->v;

        // using eqn 3.16 you can also take the body angular velocities and reframe them in the
        // inertial too
        const Vector3d w_vel_inertial = body_rates_to_euler_rates(phi, theta) * b->w;

        // updating the body with the new integrated values for the inertial frame
        b->p += vel_inertial * dt;
        b->r += w_vel_inertial * dt;
    };
}
} // namespace uav
