#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron


void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
    // column is x0,x1,x2,x3
    // row is x,y,z coordinate for each x 
    Eigen::Matrix34d x;
    x.setZero();
    for (int i = 0; i < 4; i++) {
        x.block<3, 1>(0, i) = q.segment<3>(3 * element(i));
    }

    Eigen::Matrix43d dphi;
    Eigen::Matrix3d F;

    // lecture 4 slide 78 - 79
    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        dphi_linear_tetrahedron_dX(dphi, V, element, X);

        // per-tetrahedron deformation gradient
        F = x * dphi;
        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}