#include <V_spring_particle_particle.h>

//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  V - potential energy for the spring
// 
//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0


void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    // see the lagrangian formula in lecture 3 page 33
    // 
    // Set up the Matrix B formed by combine -I and I
    // No need to set up B since it's just q1-q0 3D distance

    /*Eigen::MatrixXd B(3, 6);
    B << -1., 0., 0., 1., 0., 0.,
        0., -1., 0., 0., 1., 0.,
        0., 0., -1., 0., 0., 1.;
    Eigen::VectorXd q(6, 1);
    q << q0, q1;
    V = 0.5 * stiffness * pow((sqrt(q.transpose() * B.transpose() * B * q)), 2);*/

    double deformed_length;
    deformed_length = (q1 - q0).norm();

    V = 0.5 * stiffness * std::pow(deformed_length - l0, 2);
}