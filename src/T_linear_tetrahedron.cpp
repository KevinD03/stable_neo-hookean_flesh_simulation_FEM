#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

//Input:
// qdot - generalied velocity of FEM system
// element - vertex indices of the element
// density - material density
// volume - volume of tetrahedron
//Output:
//  T - kinetic energy of tetrahedron

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::Matrix1212d M;
    mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

    Eigen::Vector12d element_qdot;
    for (int i = 0; i < element.size(); i++) {
        element_qdot.segment(i * 3, 3) = qdot.segment(3 * element(i), 3);
    }

    // lecture 4 slide 57
    T = 0.5 * element_qdot.transpose() * M * qdot;
}