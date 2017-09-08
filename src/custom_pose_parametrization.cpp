#include "ceres_analytical/custom_pose_parametrization.h"

namespace ceres{

  bool CustomPoseParameterization::Plus(const double* x,
                                        const double* delta,
                                        double* x_plus_delta) const {

    //for the first part of x (corresponding to the quaternion)
    const double norm_delta =
        sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
    if (norm_delta > 0.0) {
      const double sin_delta_by_delta = (sin(norm_delta) / norm_delta);
      double q_delta[4];
      q_delta[0] = cos(norm_delta);
      q_delta[1] = sin_delta_by_delta * delta[0];
      q_delta[2] = sin_delta_by_delta * delta[1];
      q_delta[3] = sin_delta_by_delta * delta[2];
      QuaternionProduct(q_delta, x, x_plus_delta);
    } else {
      for (int i = 0; i < 4; ++i) {
        x_plus_delta[i] = x[i];
      }
    }


    //second part of x corresponding to the translation
    const double* x_translation= x+4;
    double* x_plus_delta_translation= x_plus_delta+4;
    const double* delta_translation = delta + 3;
    VectorRef(x_plus_delta_translation, 3) = ConstVectorRef(x_translation, 3) + ConstVectorRef(delta_translation, 3);


    return true;
  }

  bool CustomPoseParameterization::ComputeJacobian(const double* x,
                                                   double* jacobian) const {

    // MatrixRef jacobian(jacobian, 7, 6);
    // jacobian.setZero();

    //quat
    jacobian[0] =  x[3]; jacobian[1]  =  x[2]; jacobian[2]  = -x[1];  // NOLINT
    jacobian[6] = -x[2]; jacobian[7]  =  x[3]; jacobian[8]  =  x[0];  // NOLINT
    jacobian[12] =  x[1]; jacobian[13]  = -x[0]; jacobian[14]  =  x[3];  // NOLINT
    jacobian[18] = -x[0]; jacobian[19] = -x[1]; jacobian[20] = -x[2]; // NOLINT


    //translation
    jacobian[27] =  1; jacobian[28]  =  0; jacobian[29]  = 0;  // NOLINT
    jacobian[33] = 0; jacobian[34]  =  1; jacobian[35]  =  0;  // NOLINT
    jacobian[39] =  0; jacobian[40]  = 0; jacobian[41]  =  1;  // NOLINT



    //rest is zeros
    jacobian[3] =  0; jacobian[4]  =  0; jacobian[5]  = 0;  // NOLINT
    jacobian[9] = 0; jacobian[10]  =  0; jacobian[11]  =  0;  // NOLINT
    jacobian[15] =  0; jacobian[16]  = 0; jacobian[17]  =  0;  // NOLINT
    jacobian[21] = 0; jacobian[22] = 0; jacobian[23] = 0; // NOLINT
    jacobian[24] =  0; jacobian[25]  =  0; jacobian[26]  = 0;  // NOLINT
    jacobian[30] = 0; jacobian[31]  =  0; jacobian[32]  =  0;  // NOLINT
    jacobian[36] =  0; jacobian[37]  = 0; jacobian[38]  =  0;  // NOLINT



    // Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor> > M(jacobian);
    // std::cout << "M is " << std::endl << M << '\n';




    //the 7x6 jacobian that maps from the global 7 to the local 6

    // jacobian[0] = 1; jacobian[1]  = 0; jacobian[2]  = 0;  // NOLINT
    // jacobian[3] = 0; jacobian[4]  =  1; jacobian[5]  = 0;  // NOLINT
    // jacobian[6] = 0; jacobian[7]  =  0; jacobian[8]  =  1;  // NOLINT
    // jacobian[9] =  0; jacobian[10] = 0; jacobian[11] =  0;  // NOLINT
    return true;
  }

}
