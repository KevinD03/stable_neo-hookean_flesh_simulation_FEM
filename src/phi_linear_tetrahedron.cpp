#include <phi_linear_tetrahedron.h>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the energy density
//Output:
//  phi - the 4x1 values the basis functions

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
    // this is for a single tetrahedron

    phi.setZero();

    // each vector is a x,y,z coordinate
    Eigen::Vector3d underformed_X_0 = V.row(element(0));
    Eigen::Vector3d underformed_X_1 = V.row(element(1));
    Eigen::Vector3d underformed_X_2 = V.row(element(2));
    Eigen::Vector3d underformed_X_3 = V.row(element(3));

    Eigen::Matrix3d T;
    T.block<3, 1>(0, 0) = underformed_X_1 - underformed_X_0;
    T.block<3, 1>(0, 1) = underformed_X_2 - underformed_X_0;
    T.block<3, 1>(0, 2) = underformed_X_3 - underformed_X_0;

    // transfer to barycentric corrdinate
    Eigen::Vector3d vector_phi123 = T.inverse() * (x - underformed_X_0);
    double phi_0 = 1 - vector_phi123(0) - vector_phi123(1) - vector_phi123(2);

    phi(0) = phi_0;
    phi.tail<3>() = vector_phi123;
}