#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

// this is to visualize the matrix mesh
//#include <unsupported/Eigen/SparseExtra>

//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
    // for visual reference, see the stiffness matrix in lecture 3 page 70
    // it describe relationships between x0's phi and x1's phi and so on

    M.resize(qdot.rows(), qdot.rows());
    M.setZero();

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> tripletList;

    for (int tetrahedra_i = 0; tetrahedra_i < T.rows(); tetrahedra_i++) {
    // now select each tetrahedra and calculate each mass matrix
        Eigen::Matrix1212d M_i;
        Eigen::RowVector4i element_i;
        element_i = T.row(tetrahedra_i);

        // where calculation for each tetrahedra happens
        mass_matrix_linear_tetrahedron(M_i, qdot, element_i, density, v0(tetrahedra_i));

        // now distrubute the sparse matrix
        // see explain in mass_matrix_linear_tetrahedron for clarity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 3; k++) {
                    tripletList.push_back(Triplet(3 * T(tetrahedra_i, i) + k, 3 * T(tetrahedra_i, j) + k, M_i(3 * i + k, 3 * j + k) ));
                }
            }
        }
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());

    // this is to visualize the matrix mesh
    //Eigen::saveMarket(M, "mass_mesh.txt");
}