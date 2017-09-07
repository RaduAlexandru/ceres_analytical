#pragma once

#include "ceres/rotation.h"

#define CAMERA_NR_INTRINSICS 3  //BAL has 3 intrinsics, focal and 2 distorsion params

namespace ceres {
namespace examples {


// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 10 parameters. 4 for rotation, 3 for
// translation, 1 for focal length and 2 for radial distortion. The
// principal point is not modeled (i.e. it is assumed be located at
// the image center).
struct ReprojectionErrorWithoutIntrinsics {
  // (u, v): the position of the observation with respect to the image
  // center point.
  ReprojectionErrorWithoutIntrinsics(const Eigen::Vector2d & point2D)
      : observed_x(point2D (0)), observed_y(point2D (1)) {}

  template <typename T>
  bool operator()(const T * const qtvec,
                  const T * const point3D,
                  const T * const camera_params,  //will be kept fixed
                  T * residuals) const {



    //colmap

    const T * const tvec = qtvec + 4;  //pointer to the translation part

    //swap the qvec to w,x,y,z as ceres wants it   //TODO we are not using eigen for the moment, uncomment when we do
    // const T qvec[4] = {T ( qtvec[3] ), T ( qtvec[0] ), T ( qtvec[1] ), T ( qtvec[2] ) };
    const T qvec[4] = {T ( qtvec[0] ), T ( qtvec[1] ), T ( qtvec[2] ), T ( qtvec[3] ) };

    // Rotate
    T point3D_local[3];
    const T scale = T(1) / sqrt(qvec[0] * qvec[0] +
                                qvec[1] * qvec[1] +
                                qvec[2] * qvec[2] +
                                qvec[3] * qvec[3]);
    const T unit_q[4] = {
      scale * qvec[0],
      scale * qvec[1],
      scale * qvec[2],
      scale * qvec[3],
    };
    const T t2 =  unit_q[0] * unit_q[1];
    const T t3 =  unit_q[0] * unit_q[2];
    const T t4 =  unit_q[0] * unit_q[3];
    const T t5 = -unit_q[1] * unit_q[1];
    const T t6 =  unit_q[1] * unit_q[2];
    const T t7 =  unit_q[1] * unit_q[3];
    const T t8 = -unit_q[2] * unit_q[2];
    const T t9 =  unit_q[2] * unit_q[3];
    const T t1 = -unit_q[3] * unit_q[3];
    point3D_local[0] = T(2) * ((t8 + t1) * point3D[0] + (t6 - t4) * point3D[1] + (t3 + t7) * point3D[2]) + point3D[0];  // NOLINT
    point3D_local[1] = T(2) * ((t4 + t6) * point3D[0] + (t5 + t1) * point3D[1] + (t9 - t2) * point3D[2]) + point3D[1];  // NOLINT
    point3D_local[2] = T(2) * ((t7 - t3) * point3D[0] + (t2 + t9) * point3D[1] + (t5 + t8) * point3D[2]) + point3D[2];  // NOLINT
    // ceres::UnitQuaternionRotatePoint ( qvec, point3D, point3D_local );


    //translate
    point3D_local[0] += tvec[0];
    point3D_local[1] += tvec[1];
    point3D_local[2] += tvec[2];

    // std::cout << "point local_x is " << point3D_local[2]  << '\n';
    // std::cout << "point clocal is " << point3D_local[0] << " " << point3D_local[1] << " " << point3D_local[2]  << '\n';

    //TODO why is this needed
    // if ( point3D_local[2] < T ( 0 ) )
    // {
    //   std::cout << "nothing" << '\n';
    //   return false;
    // }


    // Normalize to image plane.

    // The sign change comes from
    // the camera model that Noah Snavely's Bundler assumes, whereby
    // the camera coordinate system has a negative z axis.
    point3D_local[0] /= -point3D_local[2];
    point3D_local[1] /= -point3D_local[2];

    // Distort.
    //TODO
    const T distortion=T(1.0);


    //project
    const T focal = camera_params[0];
    const T predicted_x = focal * distortion * point3D_local[0];
    const T predicted_y = focal * distortion * point3D_local[1];


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;


    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const Eigen::Vector2d & point2D) {
    return (new ceres::AutoDiffCostFunction<
            ReprojectionErrorWithoutIntrinsics,
              2, //2 residuals
              7,  //7 for rotation and translation
              3,  //3 for the 3D point
              CAMERA_NR_INTRINSICS  //intrinsics of the camera (will be kept fixed)
              >(
                new ReprojectionErrorWithoutIntrinsics( point2D)));
  }

  double observed_x;
  double observed_y;
};

}  // namespace examples
}  // namespace ceres
