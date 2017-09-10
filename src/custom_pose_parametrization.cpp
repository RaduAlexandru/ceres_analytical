#include "ceres_analytical/custom_pose_parametrization.h"

namespace ceres{

  bool CustomPoseParameterization::Plus(const double* x,
                                        const double* delta,
                                        double* x_plus_delta) const {

    //for the first part of x (corresponding to the quaternion)
    // std::cout << "custom: delta rot is " << delta[0] << " " << delta[1] << " "<< delta[2] << '\n';
    // std::cout << "custom: delta trans is " << delta[3] << " " << delta[4] << " "<< delta[5] << '\n';
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


    // //new implementation of the plus of a quaternion http://citeseerx.ist.psu.edu/viewdoc/download?rep=rep1&type=pdf&doi=10.1.1.132.20
    // //also the same used in line 330 of Sophus https://github.com/strasdat/Sophus/blob/master/sophus/so3.hpp
    // const double theta = sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
    // // std::cout << "delta is " << delta[0] << " " << delta[1] << " "<< delta[2] << '\n';
    // if (theta >std::numeric_limits<double>::epsilon()){
    //   double q_delta[4];
    //   q_delta[0]=cos( (1.0/2.0) * theta);
    //   q_delta[1]=sin((1.0/2.0)*theta)/theta*delta[1];
    //   q_delta[2]=sin((1.0/2.0)*theta)/theta*delta[2];
    //   q_delta[3]=sin((1.0/2.0)*theta)/theta*delta[3];
    //   QuaternionProduct(q_delta, x, x_plus_delta);
    // } else {
    //   for (int i = 0; i < 4; ++i) {
    //     x_plus_delta[i] = x[i];
    //   }
    // }


    //my own imlementation
    // double q_delta[4];
    // AngleAxisToQuaternion(delta,q_delta);
    // QuaternionProduct(q_delta, x, x_plus_delta);




    //second part of x corresponding to the translation
    const double* x_translation= x+4;
    double* x_plus_delta_translation= x_plus_delta+4;
    const double* delta_translation = delta + 3;
    VectorRef(x_plus_delta_translation, 3) = ConstVectorRef(x_translation, 3) + ConstVectorRef(delta_translation, 3);




    return true;
  }

  bool CustomPoseParameterization::ComputeJacobian(const double* x,
                                                   double* jacobian) const {


    //  std::cout << "x is " << '\n';
    //  for (size_t i = 0; i < 7; i++) {
    //    std::cout << "x: " << x[i] << '\n';
    //  }
    //  std::cout  << '\n';


    //This case is used with the ReprojectionErrorWithoutIntrinsics cost function.
    // MatrixRef j_eigen(jacobian, 7, 6);
    // j_eigen.setZero();
    // j_eigen.block(4,3,3,3)=Eigen::MatrixXd::Identity(3,3);
    //
    // jacobian[0] =  x[3]; jacobian[1]  =  x[2]; jacobian[2]  = -x[1];  // NOLINT
    // jacobian[6] = -x[2]; jacobian[7]  =  x[3]; jacobian[8]  =  x[0];  // NOLINT
    // jacobian[12] =  x[1]; jacobian[13]  = -x[0]; jacobian[14]  =  x[3];  // NOLINT
    // jacobian[18] = -x[0]; jacobian[19] = -x[1]; jacobian[20] = -x[2]; // NOLINT
    //
    // Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor> > M(jacobian);
    // std::cout << "local parametrization jacobian is " << std::endl << M << '\n';




    // //JUST IDENTITY for use with ErrorAnalytical
    MatrixRef j_eigen(jacobian, 7, 6);
    j_eigen.setZero();
    j_eigen.block(0,0,6,6)=Eigen::MatrixXd::Identity(6,6);
    //
    // std::cout << "j eigen is " << j_eigen << '\n';
    //





    //the 7x6 jacobian that maps from the global 7 to the local 6

    // jacobian[0] = 1; jacobian[1]  = 0; jacobian[2]  = 0;  // NOLINT
    // jacobian[3] = 0; jacobian[4]  =  1; jacobian[5]  = 0;  // NOLINT
    // jacobian[6] = 0; jacobian[7]  =  0; jacobian[8]  =  1;  // NOLINT
    // jacobian[9] =  0; jacobian[10] = 0; jacobian[11] =  0;  // NOLINT
    return true;
  }


}
