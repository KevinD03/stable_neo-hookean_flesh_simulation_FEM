#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    
   dw.setZero();

    double mu = C * 2;

    double F1_1, F1_2, F1_3, F2_1, F2_2, F2_3, F3_1, F3_2, F3_3;
    F1_1 = F(0, 0);
    F1_2 = F(0, 1);
    F1_3 = F(0, 2);
    F2_1 = F(1, 0);
    F2_2 = F(1, 1);
    F2_3 = F(1, 2);
    F3_1 = F(2, 0);
    F3_2 = F(2, 1);
    F3_3 = F(2, 2);


    /*dw(0) = C * (F1_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 - (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) - D * (F2_2 * F3_3 - F2_3 * F3_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(1) = C * (F2_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) + D * (F1_2 * F3_3 - F1_3 * F3_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(2) = C * (F3_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 - (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(3) = C * (F1_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) + D * (F2_1 * F3_3 - F2_3 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(4) = C * (F2_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 - (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(5) = C * (F3_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(6) = C * (F1_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 - (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) - D * (F2_1 * F3_2 - F2_2 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(7) = C * (F2_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(8) = C * (F3_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 - (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    */

    dw(0) = -mu * (F2_2 * F3_3 - F2_3 * F3_2) + F1_1 * mu - D * (F2_2 * F3_3 - F2_3 * F3_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(1) = mu * (F1_2 * F3_3 - F1_3 * F3_2) + F2_1 * mu + D * (F1_2 * F3_3 - F1_3 * F3_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(2) = -mu * (F1_2 * F2_3 - F1_3 * F2_2) + F3_1 * mu - D * (F1_2 * F2_3 - F1_3 * F2_2) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(3) = mu * (F2_1 * F3_3 - F2_3 * F3_1) + F1_2 * mu + D * (F2_1 * F3_3 - F2_3 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(4) = -mu * (F1_1 * F3_3 - F1_3 * F3_1) + F2_2 * mu - D * (F1_1 * F3_3 - F1_3 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(5) = mu * (F1_1 * F2_3 - F1_3 * F2_1) + F3_2 * mu + D * (F1_1 * F2_3 - F1_3 * F2_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(6) = -mu * (F2_1 * F3_2 - F2_2 * F3_1) + F1_3 * mu - D * (F2_1 * F3_2 - F2_2 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(7) = mu * (F1_1 * F3_2 - F1_2 * F3_1) + F2_3 * mu + D * (F1_1 * F3_2 - F1_2 * F3_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    dw(8) = -mu * (F1_1 * F2_2 - F1_2 * F2_1) + F3_3 * mu - D * (F1_1 * F2_2 - F1_2 * F2_1) * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0; 
    
}