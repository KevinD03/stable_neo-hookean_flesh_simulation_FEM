#include <assemble_forces.h>
#include <iostream>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system


void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
    f.resize(q.rows());
    f.setZero();

    for (int i = 0; i < T.rows(); i++) {

        Eigen::Vector12d dV;
        Eigen::RowVector4i element = T.row(i);
        dV_linear_tetrahedron_dq(dV, q, V, element, v0(i), C, D);
        
        // set up for each vertices, x0,x1,x2,x3
        // 
        // x0  0   ......
        // y0  0   ......
        // z0  0   ......
        //  0  x1  ......
        //  0  y1  ......
        //  0  z1  ......


        // setup x0
        f(3 * T(i, 0)) -= dV(0);
        f(3 * T(i, 0) + 1) -= dV(1);
        f(3 * T(i, 0) + 2) -= dV(2);

        // setup x1
        f(3 * T(i, 1)) -= dV(3);
        f(3 * T(i, 1) + 1) -= dV(4);
        f(3 * T(i, 1) + 2) -= dV(5);

        // set up x2
        f(3 * T(i, 2)) -= dV(6);
        f(3 * T(i, 2) + 1) -= dV(7);
        f(3 * T(i, 2) + 2) -= dV(8);

        // set up x3
        f(3 * T(i, 3)) -= dV(9);
        f(3 * T(i, 3) + 1) -= dV(10);
        f(3 * T(i, 3) + 2) -= dV(11);
    
    }
        
};