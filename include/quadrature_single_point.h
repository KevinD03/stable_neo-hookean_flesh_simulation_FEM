#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    Eigen::Vector3d X0, X1, X2, X3;
    X0 = q.segment<3>(3 * element(0));
    X1 = q.segment<3>(3 * element(1));
    X2 = q.segment<3>(3 * element(2));
    X3 = q.segment<3>(3 * element(3));

    // centroid
    Eigen::Vector3d X;
    X[0] = (X0(0) + X1(0) + X2(0) + X3(0)) * 0.25;
    X[1] = (X0(1) + X1(1) + X2(1) + X3(1)) * 0.25;
    X[2] = (X0(2) + X1(2) + X2(0) + X3(2)) * 0.25;

    integrand(integrated, q, element, X);

    integrated *= volume;
}

