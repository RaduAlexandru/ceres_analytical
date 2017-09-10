#pragma once

#include "ceres/rotation.h"
#include <limits>

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
    // std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';


    //Blanco's notation
    double f=focal;
    double gx=p3d_local[0];
    double gy=p3d_local[1];
    double gz=p3d_local[2];
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


    Eigen::MatrixXd jacobian_projection(2,3);
    Eigen::MatrixXd jacobian_point3d(2,3);
    //Blanco's one is wrong because the whole thing has to be negated
    // jacobian_point(0,0)=f/gz;
    // jacobian_point(0,1)=0.0;
    // jacobian_point(0,2)=-f*gx/gz2;
    //
    // jacobian_point(1,0)=0.0;
    // jacobian_point(1,1)=f/gz;
    // jacobian_point(1,2)=-f*gx/gz2;
    // jacobian_point=jacobian_point*R;


    // //orbslam one  (works if the local coordinated are left as they are)
    // jacobian_point(0,0)=-f/gz;
    // jacobian_point(0,1)=0.0;
    // jacobian_point(0,2)=gx/gz2*f;
    //
    // jacobian_point(1,0)=0.0;
    // jacobian_point(1,1)=-f/gz;
    // jacobian_point(1,2)=gy/gz2*f;
    // jacobian_point= jacobian_point*R;



    // //orbslam one (copied directly)
    jacobian_projection(0,0)=-f/gz;
    jacobian_projection(0,1)=0.0;
    jacobian_projection(0,2)=gx/gz2*f;

    jacobian_projection(1,0)=0.0;
    jacobian_projection(1,1)=-f/gz;
    jacobian_projection(1,2)=gy/gz2*f;
    jacobian_point3d= jacobian_projection*R;
    // std::cout << "jacobian_point3d is " << std::endl << jacobian_point3d << '\n';




    //sympy one (fucking works!!)
//     double p3d_x=p3d_global[0];
//     double p3d_y=p3d_global[1];
//     double p3d_z=p3d_global[2];
//     double t_x=parameters[0][4];
//     double t_y=parameters[0][5];
//     double t_z=parameters[0][6];
//     //x
//     jacobian_point(0,0)=-R(2,0)*f*(-p3d_x - R(0,1)*p3d_y - R(0,2)*p3d_z - t_x)/std::pow((R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z),2) - f/(R(2,0)*p3d_x + R(2,1)*p3d_y + 9*p3d_z + t_z);
//     jacobian_point(1,0)=-R(2,0)*f*(-R(1,0)*p3d_x - R(1,1)*p3d_y - R(1,2)*p3d_z - t_y)/std::pow((R(1,0)*p3d_x + R(2,1)*p3d_y +
// R(2,2)*p3d_z + t_z),2) - R(1,0)*f/(R(2.0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z);
//     //y
//     jacobian_point(0,1)=-R(2,1)*f*(-p3d_x - R(0,1)*p3d_y - R(0,2)*p3d_z - t_x)/std::pow((R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z),2) - R(0,1)*f/(R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z);
//     jacobian_point(1,1)=-R(2,1)*f*(-R(1,0)*p3d_x - R(1,1)*p3d_y - R(1,2)*p3d_z - t_y)/std::pow((R(2,0)*p3d_x + R(2,1)*p3d_y
// + R(2,2)*p3d_z + t_z),2) - R(1,1)*f/(R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z);
//     //z
//     jacobian_point(0,2)=-R(2,2)*f*(-p3d_x - R(0,1)*p3d_y - R(0,2)*p3d_z - t_x)/std::pow((R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z),2) - R(0,2)*f/(R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z);
//     jacobian_point(1,2)=-R(2,2)*f*(-R(1,0)*p3d_x - R(1,1)*p3d_y - R(1,2)*p3d_z - t_y)/std::pow((R(2,0)*p3d_x + R(2,1)*p3d_y
// + R(2,2)*p3d_z + t_z),2) - R(1,2)*f/(R(2,0)*p3d_x + R(2,1)*p3d_y + R(2,2)*p3d_z + t_z);


      //blanco one but with the local point  (also works but the last cost is not as good as the sympy one)
      // jacobian_point(0,0)=f/p3d_local[2];
      // jacobian_point(0,1)=0.0;
      // jacobian_point(0,2)=-f*p3d_local[0]/(p3d_local[2]*p3d_local[2]);
      //
      // jacobian_point(1,0)=0.0;
      // jacobian_point(1,1)=f/p3d_local[2];
      // jacobian_point(1,2)=-f*p3d_local[0]/(p3d_local[2]*p3d_local[2]);
      // jacobian_point=jacobian_point*R;





    // std::cout << "jacobian point is " << std::endl << jacobian_point << '\n';


    if (jacobians != NULL && jacobians[0] != NULL) {

      // std::cout << "computing jacobian" << '\n';

      //jacobian of camera_pose is of size 2x7 , fill in a row majow order jacobians[0][0...13]
      // jacobians[0][0] = f/gz;
      // jacobians[0][1] = 0.0;
      // jacobians[0][2] = -f*(gx/gz2);
      // jacobians[0][3] = -f*( (gx*gy) /gz2 );
      // jacobians[0][4] = f*( 1.0 +  gx2/gz2 );
      // jacobians[0][5] = -f*( gy/gz );
      // jacobians[0][6] = 0.0;
      //
      // jacobians[0][7] = 0.0;
      // jacobians[0][8] = f/gz;
      // jacobians[0][9] = -f*(gy/gz2);
      // jacobians[0][10] = -f*(1.0+ (gy2/gz2) );
      // jacobians[0][11] = f*( (gx*gy)/gz2) ;
      // jacobians[0][12] = f* (gx/gz);
      // jacobians[0][13] = 0.0;


      // //the orblam one
      // jacobians[0][0] = gx*gy/gz2 *f;
      // jacobians[0][1] = -(1+(gx*gx/gz2)) *f;
      // jacobians[0][2] = gy/gz *f;
      // jacobians[0][3] = -1./gz *f;
      // jacobians[0][4] = 0;
      // jacobians[0][5] = gx/gz2 *f;
      // jacobians[0][6] = 0.0;
      //
      // jacobians[0][7] = (1+gy*gy/gz2) *f;
      // jacobians[0][8] = -gx*gy/gz2 *f;
      // jacobians[0][9] = -gx/gz *f;
      // jacobians[0][10] = 0;
      // jacobians[0][11] = -1./gz *f ;
      // jacobians[0][12] = gy/gz *f;
      // jacobians[0][13] = 0.0;

      //my own test one (good one)
      // jacobians[0][0] = -f/gz;
      // jacobians[0][1] = 0.0;
      // jacobians[0][2] = gx/gz2*f;
      // jacobians[0][3] = -f*( (gx*gy) /gz2 );
      // jacobians[0][4] = f*( 1.0 +  gx2/gz2 );
      // jacobians[0][5] = f*( gy/gz );
      // jacobians[0][6] = 0.0;
      //
      // jacobians[0][7] = 0.0;
      // jacobians[0][8] = -f/gz;
      // jacobians[0][9] = gy/gz2*f;
      // jacobians[0][10] = -f*(1.0+ (gy2/gz2) );
      // jacobians[0][11] = f*( (gx*gy)/gz2) ;
      // jacobians[0][12] = -f* (gx/gz);
      // jacobians[0][13] = 0.0;

      //my own test one ()2
      // jacobians[0][0] = gx*gy/gz2 *f;
      // jacobians[0][1] = -(1+(gx*gx/gz2)) *f;
      // jacobians[0][2] = -gy/gz *f;
      // jacobians[0][3] = 1./gz *f;
      // jacobians[0][4] = 0;
      // jacobians[0][5] = gx/gz2 *f;
      // jacobians[0][6] = 0.0;
      //
      // jacobians[0][7] = (1+gy*gy/gz2) *f;
      // jacobians[0][8] = -gx*gy/gz2 *f;
      // jacobians[0][9] = gx/gz *f;
      // jacobians[0][10] = 0;
      // jacobians[0][11] = 1./gz *f ;
      // jacobians[0][12] = -gy/gz *f;
      // jacobians[0][13] = 0.0;

      //Eigen
      MatrixRef j_eigen(jacobians[0], 2, 7);
      j_eigen.setZero();
      // j_eigen.block(0,0,2,3)=jacobian_point;
      // p3d_local[0] = 2 * ((t8 + t1) * p3d_global[0] + (t6 - t4) * p3d_global[1] + (t3 + t7) * p3d_global[2]) + p3d_global[0];  // NOLINT
      // p3d_local[1] = 2 * ((t4 + t6) * p3d_global[0] + (t5 + t1) * p3d_global[1] + (t9 - t2) * p3d_global[2]) + p3d_global[1];  // NOLINT
      // p3d_local[2] = 2 * ((t7 - t3) * p3d_global[0] + (t2 + t9) * p3d_global[1] + (t5 + t8) * p3d_global[2]) + p3d_global[2];  // NOLINT
      // Eigen::Vector3d g={p3d_global[0],p3d_global[1],p3d_global[2]};
      Eigen::Vector3d camera_trans={parameters[0][4],parameters[0][5],parameters[0][6]};
      p3d_local[0] -= camera_trans[0];
      p3d_local[1] -= camera_trans[1];
      p3d_local[2] -= camera_trans[2];
      Eigen::Vector3d g={p3d_local[0],p3d_local[1],p3d_local[2]};

      // std::cout << "g local is " << g << '\n';
      // std::cout << "camera transi is" << camera_trans << '\n';
      Eigen::Matrix3d g_hat;
      g_hat << 0.0, -g(2), g(1),
              g(2), 0.0, -g(0),
              -g(1), g(0), 0.0;
      // g_hat << -g(1), g(2), 0,
      //         g(0), 0.0, g(2),
      //          0, -g(0), g(1);
      // std::cout << "g_hat is" << std::endl  << g_hat << '\n';
      Eigen::MatrixXd right_side(3,7);
      right_side.setZero();
      // right_side.block(0,0,3,3)=Eigen::Matrix3d::Identity();
      // right_side.block(0,3,3,3)=g_hat;

      right_side.block(0,0,3,3)=-g_hat*2;
      right_side.block(0,3,3,3)=Eigen::Matrix3d::Identity();
      // std::cout << "jacobian point is " << std::endl << jacobian_point << '\n';
      // std::cout << "righ side" << std::endl << right_side << '\n';

      j_eigen=jacobian_projection*right_side;
      //TODO remove this
      // j_eigen.block(0,0,3,3)=Eigen::Matrix3d::Zero();

      j_eigen(0,6)=0.0;
      j_eigen(1,6)=0.0;


      //testint th multiplication of the jacobian point with jhat
      // Eigen::MatrixXd result= jacobian_point*g_hat;
      // std::cout << "jacobian point is " << std::endl << jacobian_point << '\n';
      // std::cout << "g_hat is " << std::endl << g_hat << '\n';
      // std::cout << "result is " << std::endl << result << '\n';

      //testing guillermos thingy https://arxiv.org/pdf/1312.0788.pdf
      // double angle_axis_arr[3];
      // RotationMatrixToAngleAxis(R_array,angle_axis_arr);
      // Eigen::Vector3d v={angle_axis_arr[0],angle_axis_arr[1],angle_axis_arr[2]};
      // Eigen::Matrix3d v_hat;
      // v_hat << 0.0, -v(2), v(1),
      //         v(2), 0.0, -v(0),
      //         -v(1), v(0), 0.0;
      // Eigen::MatrixXd result= -R*g_hat * (v*v.transpose() + (R.transpose() - Eigen::Matrix3d::Identity())*v_hat )/(v.norm());
      // result=jacobian_projection*result;
      // // std::cout << "result is " << std::endl << result << '\n';
      // j_eigen.block(0,0,3,3)=result;







      //derived by hand (the last block of 3x3 is correct)
      // MatrixRef j_eigen(jacobians[0], 2, 7);
      // j_eigen.setZero();
      // j_eigen.block(0,3,2,3)=jacobian_point;
      // j_eigen(0,0)=f*gy*gx/gz2;
      // j_eigen(0,1)=f+f*gx*gx/gz2;
      // j_eigen(0,2)=-f*gy/gz;
      // j_eigen(1,0)=-f-f*gy*gy/gz2;
      // j_eigen(1,1)=f*gx*gy/gz2;
      // j_eigen(1,2)=f*gx/gz;


      // if (j_eigen(0,0)<135 && j_eigen(0,0)>133){
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      //   std::cout << "FUCK YES" << '\n';
      // }

      // Eigen::Map<Eigen::Matrix<double,2,7,Eigen::RowMajor> > M(jacobians[0]);
      // std::cout << "jacobian of quaternion is " << std::endl << M << '\n';
      //
      //
      // Eigen::MatrixXd parametrization(7,6);
      // parametrization.setZero();
      // parametrization.block(0,0,6,6)=Eigen::MatrixXd::Identity(6,6);
      // Eigen::MatrixXd result=j_eigen*parametrization;
      //
      //
      // std::cout << "parametrized is " << std::endl << result << '\n';



      //test of multiplication
      // Eigen::MatrixXd J(2,7);
      // J << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0;
      // Eigen::MatrixXd parametrization(7,6);
      // parametrization.setZero();
      // parametrization.block(0,0,6,6)=Eigen::MatrixXd::Identity(6,6);
      // Eigen::MatrixXd result=J*parametrization;
      //
      //
      // std::cout << "J is " << std::endl << J << '\n';
      // std::cout << "parametrization is " << std::endl << parametrization << '\n';
      // std::cout << "result is " << std::endl << result << '\n';




    }

    if (jacobians != NULL && jacobians[1] != NULL) {
      // std::cout << "computing jacobian point" << '\n';

      //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[0][0...5]
      jacobians[1][0] = jacobian_point3d(0,0);
      jacobians[1][1] = jacobian_point3d(0,1);
      jacobians[1][2] = jacobian_point3d(0,2);

      jacobians[1][3] = jacobian_point3d(1,0);
      jacobians[1][4] = jacobian_point3d(1,1);
      jacobians[1][5] = jacobian_point3d(1,2);


    }

    if (jacobians != NULL && jacobians[2] != NULL) {

      std::cout << "requiring jacobians of camera instrinsics" << '\n';

    }

    // exit(1);
    return true;
  }

private:
  double observed_x;
  double observed_y;
};

} //namespace ceres
