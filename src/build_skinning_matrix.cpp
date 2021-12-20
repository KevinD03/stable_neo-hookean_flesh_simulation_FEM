#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix 

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> tripletList;

    N.resize(V_skin.rows(), V.rows());

    // looping through each row of V_skin, each row contain the coordinate of one skin vertice, (x,y,z)
    for (int mesh_vertice = 0; mesh_vertice < V_skin.rows(); mesh_vertice++) {
        Eigen::Vector3d current_skin = V_skin.row(mesh_vertice);
        int wanted_vertice_on_tetrahedron;
        Eigen::Vector4d wanted_phi;

        double phi_norm = 9999999999999;

        // loop through all tetrahedral and find which tetrahedral our skinning vertice is in
        for (int i = 0; i < T.rows(); i++) {
            Eigen::Vector4d cur_phi;
            Eigen::RowVectorXi element = T.row(i);
            phi_linear_tetrahedron(cur_phi, V, element, current_skin);

            // find the tetrahedral our skinning vertice is in and the tetrahedral's phi
            if (cur_phi.norm() < phi_norm) {
                phi_norm = cur_phi.norm();
                wanted_phi = cur_phi;
                wanted_vertice_on_tetrahedron = i;
            }
        }
        for (int k = 0; k < 4; k++){
            tripletList.push_back(Triplet(mesh_vertice, T(wanted_vertice_on_tetrahedron, k), wanted_phi(k)));
        }
    }

    N.setFromTriplets(tripletList.begin(), tripletList.end());

}