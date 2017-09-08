#pragma once

#include "ceres/rotation.h"

#define CAMERA_NR_INTRINSICS 3  //BAL has 3 intrinsics, focal and 2 distorsion params





namespace ceres{

class ErrorAnalytical : public SizedCostFunction<2, /* number of residuals */
                             7, /* size of first parameter, corresponding to the camera pose */
                             3, /* size of the second parameter corresponding to the point*/
                             CAMERA_NR_INTRINSICS /* size of the second parameter corresponding to the intrinsics*/> {
 public:
   ErrorAnalytical(const Eigen::Vector2d & point2D) : observed_x(point2D(0)), observed_y(point2D(1)) {}
  virtual ~ErrorAnalytical() {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    // double w_x = parameters[0][0];
    // double w_x = parameters[0][0];
    // double w_x = parameters[0][0];

    double p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points that needs to be rotated
    double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion

    //rotate
    double scale = 1.0 / sqrt(qvec[0] * qvec[0] +
                                qvec[1] * qvec[1] +
                                qvec[2] * qvec[2] +
                                qvec[3] * qvec[3]);
    double unit_q[4] = {
      scale * qvec[0],
      scale * qvec[1],
      scale * qvec[2],
      scale * qvec[3],
    };
    double t2 =  unit_q[0] * unit_q[1];
    double t3 =  unit_q[0] * unit_q[2];
    double t4 =  unit_q[0] * unit_q[3];
    double t5 = -unit_q[1] * unit_q[1];
    double t6 =  unit_q[1] * unit_q[2];
    double t7 =  unit_q[1] * unit_q[3];
    double t8 = -unit_q[2] * unit_q[2];
    double t9 =  unit_q[2] * unit_q[3];
    double t1 = -unit_q[3] * unit_q[3];
    p3d_local[0] = 2 * ((t8 + t1) * p3d_global[0] + (t6 - t4) * p3d_global[1] + (t3 + t7) * p3d_global[2]) + p3d_global[0];  // NOLINT
    p3d_local[1] = 2 * ((t4 + t6) * p3d_global[0] + (t5 + t1) * p3d_global[1] + (t9 - t2) * p3d_global[2]) + p3d_global[1];  // NOLINT
    p3d_local[2] = 2 * ((t7 - t3) * p3d_global[0] + (t2 + t9) * p3d_global[1] + (t5 + t8) * p3d_global[2]) + p3d_global[2];  // NOLINT




    // camera[3,4,5] are the translation.
    p3d_local[0] += parameters[0][4];
    p3d_local[1] += parameters[0][5];
    p3d_local[2] += parameters[0][6];


    //dehomogeneous
    double xp = - p3d_local[0] / p3d_local[2];
    double yp = - p3d_local[1] / p3d_local[2];


    //NO DISTORSION
    //TODO
    double distortion=1.0;

    // Compute final projected point position.
    double focal =parameters[2][0];
    double predicted_x = focal * distortion * xp;
    double predicted_y = focal * distortion * yp;


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;
    std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';


    //Blanco's notation
    double f=focal;
    double gx=p3d_global[0];
    double gy=p3d_global[1];
    double gz=p3d_global[2];
    double gx2=gx*gx;
    double gy2=gy*gy;
    double gz2=gz*gz;
    Eigen::Matrix3d R;
    double R_array[9];
    QuaternionToRotation(qvec,R_array);
    //Copy into eigen matrix
    R(0,0)=R_array[0];
    R(0,1)=R_array[1];
    R(0,2)=R_array[2];
    R(1,0)=R_array[3];
    R(1,1)=R_array[4];
    R(1,2)=R_array[5];
    R(2,0)=R_array[6];
    R(2,1)=R_array[7];
    R(2,2)=R_array[8];


    Eigen::MatrixXd jacobian_point(2,3);
    jacobian_point(0,0)=f/gz;
    jacobian_point(0,0)=0.0;
    jacobian_point(0,0)=-f*(gx/gz2);

    jacobian_point(0,0)=0.0;
    jacobian_point(0,0)=f/gz;
    jacobian_point(0,0)=-f*(gx/gz2);
    jacobian_point=jacobian_point*R;

    std::cout << "jacobian point is " << std::endl << jacobian_point << '\n';


    Eigen::MatrixXd local_j(7,6);
    local_j.block<6,6>(0,0)=Eigen::MatrixXd::Identity(6,6);
    std::cout << "local j"  << std::endl << local_j << '\n';


    if (jacobians != NULL && jacobians[0] != NULL) {
      std::cout << "computing jacobian" << '\n';

      //jacobian of camera_pose is of size 2x7 , fill in a row majow order jacobians[0][0...13]
      jacobians[0][0] = f/gz;
      jacobians[0][1] = 0.0;
      jacobians[0][2] = -f*(gx/gz2);
      jacobians[0][3] = -f*( (gx*gy) /gz2 );
      jacobians[0][4] = f*( 1.0 +  gx2/gz2 );
      jacobians[0][5] = -f*( gy/gz );

      jacobians[0][6] = 0.0;
      jacobians[0][7] = f/gz;
      jacobians[0][8] = -f*(gy/gz2);
      jacobians[0][9] = -f*(1.0+ (gy2/gz2) );
      jacobians[0][10] = f*( (gx*gy)/gz2) ;
      jacobians[0][11] = f* (gx/gz);
    }

    if (jacobians != NULL && jacobians[1] != NULL) {

      //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[0][0...5]
      jacobians[1][0] = jacobian_point(0,0);
      jacobians[1][1] = jacobian_point(0,1);
      jacobians[1][2] = jacobian_point(0,2);

      jacobians[1][3] = jacobian_point(1,0);
      jacobians[1][4] = jacobian_point(1,1);
      jacobians[1][5] = jacobian_point(1,2);


    }
    return true;
  }

private:
  double observed_x;
  double observed_y;
};

} //namespace ceres
