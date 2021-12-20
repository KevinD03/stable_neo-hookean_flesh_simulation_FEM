#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the energy density
//Output:
//  dphi - the 4x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX

void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    dphi.setZero();

    // each vector is a x,y,z coordinate
    Eigen::Vector3d underformed_X_0 = V.row(element(0));
    Eigen::Vector3d underformed_X_1 = V.row(element(1));
    Eigen::Vector3d underformed_X_2 = V.row(element(2));
    Eigen::Vector3d underformed_X_3 = V.row(element(3));

    Eigen::Matrix3d T;
    T.block<3, 1>(0, 0) = underformed_X_1 - underformed_X_0;
    T.block<3, 1>(0, 1) = underformed_X_2 - underformed_X_0;
    T.block<3, 1>(0, 2) = underformed_X_3 - underformed_X_0;
 
    dphi.block<1, 3>(0, 0) = -1.0 * T.inverse().colwise().sum();
    dphi.block<3, 3>(1, 0) = T.inverse();

}