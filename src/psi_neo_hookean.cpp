#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>

//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the neohookean energy


void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    // reference wiki https://en.wikipedia.org/wiki/Neo-Hookean_solid neo-hookean solid
    // code from prof David's matlab code demo

    double J = F.determinant();
    double I = (F.transpose() * F).trace();
    
    /*psi = C * (I / pow(J, 2.0 / 3.0) - 3) + D * pow(J - 1, 2);*/


    // to use the stable neo-hookean we will use the formula in the research paper page 4
    // equation 11,12,13

    double mu = C * 2;
    //lambda = D * 2;

    psi = mu / 2 * (I - 3) - mu * (J - 1) + D * pow(J - 1, 2);


}