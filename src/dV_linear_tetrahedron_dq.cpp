#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  v0 - the undeformed tetrahedron volume
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  dV - the 12x1 gradient of the potential energy for a single tetrahedron


void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    Eigen::Matrix34d x;
    x.setZero();
    for (int i = 0; i < 4; i++) {
        x.block<3, 1>(0, i) = q.segment<3>(3 * element(i));
    }

    Eigen::Matrix43d dphi;
    Eigen::Matrix3d F;
    Eigen::Vector9d dpsi;
    Eigen::MatrixXd B(9, 12);
    B.setZero();

    auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        F = x * dphi;

        // dpsi is gradient of F 
        dpsi_neo_hookean_dF(dpsi, F, C, D);
        
        // set up matrix B  9*12
        // lecture slide page 85

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 3; k++) {
                    B(i + k * 3, j * 3 + k) = dphi(j, i);
                }
            }
        }

        dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}