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

    double p3d_only_rot[3]; //resulting point after only rotation but NOT translating
    double p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points in world frame
    double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion

    //rotate
    QuaternionRotatePoint(qvec,p3d_global,p3d_only_rot);

    // camera[3,4,5] are the translation.
    p3d_local[0] = p3d_only_rot[0]+ parameters[0][4];
    p3d_local[1] = p3d_only_rot[1]+ parameters[0][5];
    p3d_local[2] = p3d_only_rot[2]+ parameters[0][6];


    //Normalize to image plane.
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
    R=Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(R_array);


    Eigen::MatrixXd jacobian_projection(2,3);
    Eigen::MatrixXd jacobian_point3d(2,3);

    jacobian_projection(0,0)=-f/gz;
    jacobian_projection(0,1)=0.0;
    jacobian_projection(0,2)=gx/gz2*f;

    jacobian_projection(1,0)=0.0;
    jacobian_projection(1,1)=-f/gz;
    jacobian_projection(1,2)=gy/gz2*f;



    if (jacobians != NULL && jacobians[0] != NULL) {


      MatrixRef j_eigen(jacobians[0], 2, 7);
      j_eigen.setZero();
      Eigen::Vector3d g={p3d_only_rot[0],p3d_only_rot[1],p3d_only_rot[2]};
      Eigen::Matrix3d g_hat;
      g_hat << 0.0, -g(2), g(1),
              g(2), 0.0, -g(0),
              -g(1), g(0), 0.0;

      Eigen::MatrixXd right_side(3,7);
      right_side.setZero();
      right_side.block(0,0,3,3)=-g_hat*2;  //TODO why is it scaled down by a half?
      right_side.block(0,3,3,3)=Eigen::Matrix3d::Identity();

      j_eigen=jacobian_projection*right_side;

    }

    if (jacobians != NULL && jacobians[1] != NULL) {

      //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[0][0...5]
      jacobian_point3d= jacobian_projection*R;
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




class ErrorAnalyticalOpt : public SizedCostFunction<2, /* number of residuals */
                             7, /* size of first parameter, corresponding to the camera pose */
                             3, /* size of the second parameter corresponding to the point*/
                             CAMERA_NR_INTRINSICS /* size of the second parameter corresponding to the intrinsics*/> {
 public:
   ErrorAnalyticalOpt(const Eigen::Vector2d & point2D) : observed_x(point2D(0)), observed_y(point2D(1)) {}
  virtual ~ErrorAnalyticalOpt() {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {

    double p3d_only_rot[3]; //resulting point after only rotation but NOT translating
    double p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points in world frame
    double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion

    //rotate
    QuaternionRotatePoint(qvec,p3d_global,p3d_only_rot);

    // camera[3,4,5] are the translation.
    p3d_local[0] = p3d_only_rot[0]+ parameters[0][4];
    p3d_local[1] = p3d_only_rot[1]+ parameters[0][5];
    p3d_local[2] = p3d_only_rot[2]+ parameters[0][6];


    //Normalize to image plane.
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
    double R_array[9];
    QuaternionToRotation(qvec,R_array);


    //precompute stuff
    double _z=-f/ p3d_local[2];
    double z_2=p3d_local[2]*p3d_local[2];
    double x_z=f*p3d_local[0]/ z_2;
    double y_z=f*p3d_local[1]/ z_2;



    if (jacobians != NULL && jacobians[0] != NULL) {


      //precompute some more stuff
      double x_rot_2=p3d_only_rot[0]*2;
      double y_rot_2=p3d_only_rot[1]*2;
      double z_rot_2=p3d_only_rot[2]*2;

      //block of rotation
      jacobians[0][0]=x_z*y_rot_2;
      jacobians[0][1]=_z* z_rot_2 + x_z*-x_rot_2;
      jacobians[0][2]=_z*-y_rot_2;

      jacobians[0][7]= _z*-z_rot_2 + y_z*y_rot_2;
      jacobians[0][8]= y_z*-x_rot_2;
      jacobians[0][9]= _z* x_rot_2;


      //block of translation
      jacobians[0][3]=_z;
      jacobians[0][4]=0.0;
      jacobians[0][5]=x_z;

      jacobians[0][10]=0.0;
      jacobians[0][11]=_z;
      jacobians[0][12]=y_z;


      //last two zeros
      jacobians[0][6]=0.0;
      jacobians[0][13]=0.0;

    }

    if (jacobians != NULL && jacobians[1] != NULL) {

      //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[1][0...5]
      jacobians[1][0]= _z*R_array[0] + x_z*R_array[6];
      jacobians[1][1]= _z*R_array[1] + x_z*R_array[7];
      jacobians[1][2]= _z*R_array[2] + x_z*R_array[8];

      jacobians[1][3]= _z*R_array[3] + y_z*R_array[6];
      jacobians[1][4]= _z*R_array[4] + y_z*R_array[7];
      jacobians[1][5]= _z*R_array[5] + y_z*R_array[8];


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




class ToImagePlane : public SizedCostFunction<2, /* number of residuals */
                             7, /* size of first parameter, corresponding to the camera pose */
                             3 /* size of the second parameter corresponding to the point*/> {
 public:
   ToImagePlane() {}
  virtual ~ToImagePlane() {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {

    double p3d_only_rot[3]; //resulting point after only rotation but NOT translating
    double p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points in world frame
    double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion

    //rotate
    QuaternionRotatePoint(qvec,p3d_global,p3d_only_rot);

    // camera[3,4,5] are the translation.
    p3d_local[0] = p3d_only_rot[0]+ parameters[0][4];
    p3d_local[1] = p3d_only_rot[1]+ parameters[0][5];
    p3d_local[2] = p3d_only_rot[2]+ parameters[0][6];


    //Normalize to image plane.
    double xp = - p3d_local[0] / p3d_local[2];
    double yp = - p3d_local[1] / p3d_local[2];

    //The output will be just the points in the image plane
    residuals[0] = xp;
    residuals[1] = yp;
    // std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';



    //Compute jacobians
    Eigen::Matrix3d R;
    double R_array[9];
    QuaternionToRotation(qvec,R_array);


    //precompute stuff
    double _z=-1.0/ p3d_local[2];
    double z_2=p3d_local[2]*p3d_local[2];
    double x_z=p3d_local[0]/ z_2;
    double y_z=p3d_local[1]/ z_2;


    if (jacobians != NULL && jacobians[0] != NULL) {

      // MatrixRef j_eigen(jacobians[0], 2, 7);
      // j_eigen.setZero();
      // j_eigen.block(0,0,2,3)= jacobian_camera_rotation;
      // j_eigen.block(0,3,2,3)= jacobian_camera_translation;

      //precompute some more stuff
      double x_rot_2=p3d_only_rot[0]*2;
      double y_rot_2=p3d_only_rot[1]*2;
      double z_rot_2=p3d_only_rot[2]*2;

      //block of rotation
      jacobians[0][0]=x_z*y_rot_2;
      jacobians[0][1]=_z* z_rot_2 + x_z*-x_rot_2;
      jacobians[0][2]=_z*-y_rot_2;

      jacobians[0][7]= _z*-z_rot_2 + y_z*y_rot_2;
      jacobians[0][8]= y_z*-x_rot_2;
      jacobians[0][9]= _z* x_rot_2;


      //block of translation
      jacobians[0][3]=_z;
      jacobians[0][4]=0.0;
      jacobians[0][5]=x_z;

      jacobians[0][10]=0.0;
      jacobians[0][11]=_z;
      jacobians[0][12]=y_z;


      //last two zeros
      jacobians[0][6]=0.0;
      jacobians[0][13]=0.0;


    }

    if (jacobians != NULL && jacobians[1] != NULL) {

      //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[1][0...5]
      jacobians[1][0]= _z*R_array[0] + x_z*R_array[6];
      jacobians[1][1]= _z*R_array[1] + x_z*R_array[7];
      jacobians[1][2]= _z*R_array[2] + x_z*R_array[8];

      jacobians[1][3]= _z*R_array[3] + y_z*R_array[6];
      jacobians[1][4]= _z*R_array[4] + y_z*R_array[7];
      jacobians[1][5]= _z*R_array[5] + y_z*R_array[8];



    }


    // exit(1);
    return true;
  }

};





struct ErrorAnalyticalWithDistortion {
  ErrorAnalyticalWithDistortion(const Eigen::Vector2d & point2D) : observed_x(point2D(0)), observed_y(point2D(1)) {
    compute_normalization_to_image_place.reset(new ceres::CostFunctionToFunctor<2, 7, 3>(new ToImagePlane));
  }

  template <typename T>
  bool operator()(const T* qtvec,
                  const T* p3d_global,
                  const T* cam_intrinsics,
                  T* residuals) const {

    T p_img_place[2];
    (*compute_normalization_to_image_place)(qtvec, p3d_global, p_img_place);

    //Distorsion
    const T& l1 = cam_intrinsics[1];
    const T& l2 = cam_intrinsics[2];
    const T r2 = p_img_place[0]*p_img_place[0] + p_img_place[1]*p_img_place[1];
    const T distortion = 1.0 + r2  * (l1 + l2  * r2);

    //project
    const T focal =cam_intrinsics[0];
    const T predicted_x = focal * distortion * p_img_place[0];
    const T predicted_y = focal * distortion * p_img_place[1];


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;

    return true;
  }

  double observed_x;
  double observed_y;
  std::unique_ptr<ceres::CostFunctionToFunctor<2, 7, 3> > compute_normalization_to_image_place;


  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const Eigen::Vector2d & point2D) {
    return (new ceres::AutoDiffCostFunction<
            ErrorAnalyticalWithDistortion,
              2, //2 residuals
              7,  //7 for rotation and translation
              3,  //3 for the 3D point
              CAMERA_NR_INTRINSICS  //intrinsics of the camera 
              >(
                new ErrorAnalyticalWithDistortion( point2D)));
  }
};





} //namespace ceres
