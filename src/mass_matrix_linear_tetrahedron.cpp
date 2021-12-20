//#include <unsupported/Eigen/SparseExtra>
#include <mass_matrix_linear_tetrahedron.h>

//Input:
//  qdot - generalized velocity for the FEM system
//  element - the 1x4 vertex indices for this tetrahedron
//  density - density of material
//  volume - the undeformed tetrahedron volume
//Output:
//  M - dense 12x12 per-tetrahedron mass matrix

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
     M.setZero();

     // understand M as a 4*4 but aech component is a 3*3
     // 4*4 because x0, x1, x2, x3; 4 vertices on tetrahedron
     // each component is a 3*3 because x,y,z; each vertice contain a x,y,z coordinate
     // x phi 00,     0   ,     0   , ...... same for phi01, phi02, phi03
     //     0     y phi 00,     0   , ......
     //     0   ,     0   , z phi 00, ......
     //    same for phi10,
     //         phi20,
     //         phi30

     Eigen::Matrix4d M0;
     M0(0, 0) = 1.0 / 10.0 * density * volume;
     M0(0, 1) = 1.0 / 20.0 * density * volume;
     M0(0, 2) = 1.0 / 20.0 * density * volume;
     M0(0, 3) = 1.0 / 20.0 * density * volume;
     M0(1, 0) = 1.0 / 20.0 * density * volume;
     M0(1, 1) = 1.0 / 10.0 * density * volume;
     M0(1, 2) = 1.0 / 20.0 * density * volume;
     M0(1, 3) = 1.0 / 20.0 * density * volume;
     M0(2, 0) = 1.0 / 20.0 * density * volume;
     M0(2, 1) = 1.0 / 20.0 * density * volume;
     M0(2, 2) = 1.0 / 10.0 * density * volume;
     M0(2, 3) = 1.0 / 20.0 * density * volume;
     M0(3, 0) = 1.0 / 20.0 * density * volume;
     M0(3, 1) = 1.0 / 20.0 * density * volume;
     M0(3, 2) = 1.0 / 20.0 * density * volume;
     M0(3, 3) = 1.0 / 10.0 * density * volume;

     Eigen::Matrix3d I;
     I.setIdentity();


     for (int i = 0; i < 4; i++) {
         for (int j = 0; j < 4; j++) {
             M.block<3, 3>(i * 3, j * 3) = I * M0(i, j);
         }
     }

     //Eigen::saveMarket(M, "single_mass_mesh.txt");
 }