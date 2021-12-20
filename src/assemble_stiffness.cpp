#include <assemble_stiffness.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix


void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
    K.resize(q.rows(), q.rows());
    K.setZero();
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> coordinate_triplets;
    coordinate_triplets.reserve(q.rows() * q.rows());

    for (int tetrahedron_i = 0; tetrahedron_i < T.rows(); tetrahedron_i++) {
        Eigen::Matrix1212d d2V;
        Eigen::RowVector4i element = T.row(tetrahedron_i);
        d2V_linear_tetrahedron_dq2(d2V, q, V, element, v0(tetrahedron_i), C, D);
        
        // describe the relation between a selected vertice and other 3 vertices
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i), 3 * T(tetrahedron_i, j), -d2V(3 * i, 3 * j)) );
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i), 3 * T(tetrahedron_i, j) + 1, -d2V(3 * i, 3 * j + 1)) );
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i), 3 * T(tetrahedron_i, j) + 2, -d2V(3 * i, 3 * j + 2)) );

                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 1, 3 * T(tetrahedron_i, j), -d2V(3 * i + 1, 3 * j)));
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 1, 3 * T(tetrahedron_i, j) + 1, -d2V(3 * i + 1, 3 * j + 1)) );
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 1, 3 * T(tetrahedron_i, j) + 2, -d2V(3 * i + 1, 3 * j + 2)) );

                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 2, 3 * T(tetrahedron_i, j), -d2V(3 * i + 2, 3 * j)));
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 2, 3 * T(tetrahedron_i, j) + 1, -d2V(3 * i + 2, 3 * j + 1)) );
                coordinate_triplets.push_back( Triplet(3 * T(tetrahedron_i, i) + 2, 3 * T(tetrahedron_i, j) + 2, -d2V(3 * i + 2, 3 * j + 2)) );

            }
        }
    }
    K.setFromTriplets(coordinate_triplets.begin(), coordinate_triplets.end());
        
};
