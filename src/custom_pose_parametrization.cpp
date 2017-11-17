#include "ceres_analytical/custom_pose_parametrization.h"

namespace ceres{

  bool CustomPoseParameterization::Plus(const double* x,
                                        const double* delta,
                                        double* x_plus_delta) const {

    //for the first part of x (corresponding to the quaternion), 
    // taken from QuaternionLocalParametrization
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

    // //JUST IDENTITY for use with ErrorAnalytical
    MatrixRef j_eigen(jacobian, 7, 6);
    j_eigen.setZero();
    j_eigen.block(0,0,6,6).setIdentity(); //=Eigen::MatrixXd::Identity(6,6);

    return true;
  }


}
