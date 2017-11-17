#pragma once

#include "ceres/rotation.h"
#include <limits>

namespace ceres{

// TODO: ErrorAnalyticalOptFullRelPose
// projection function with fx,fy,cx,cy
class ErrorAnalyticalOptK : public SizedCostFunction<2, /* number of residuals */
                             7, /* size of first parameter, corresponding to the camera pose */
                             3 /* point size*/ >{
 public:
   ErrorAnalyticalOptK(const double* point2D, const double * intrinsics)
   : observed_x(point2D[0]), observed_y(point2D[1]), m_intrinsics ( intrinsics )
   {}
  virtual ~ErrorAnalyticalOptK() {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {

    double p3d_only_rot[3]; //resulting point after only rotation but NOT translating
    double p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    const double * p3d_global = parameters[1];//double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points in world frame
    const double * qvec = parameters[0];// double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion

    //rotate
    QuaternionRotatePoint(qvec,p3d_global,p3d_only_rot);

    // camera[4,5,6] are the translation, first part is quaternion
    p3d_local[0] = p3d_only_rot[0] + parameters[0][4];
    p3d_local[1] = p3d_only_rot[1] + parameters[0][5];
    p3d_local[2] = p3d_only_rot[2] + parameters[0][6];

    //Normalize to image plane.
    const double xp = p3d_local[0] / p3d_local[2];
    const double yp = p3d_local[1] / p3d_local[2];

    //NO DISTORTION

    // Compute final projected point position.
    const double & focal_x = m_intrinsics[0];
    const double & focal_y = m_intrinsics[1];
    const double & c_x = m_intrinsics[2];
    const double & c_y = m_intrinsics[3];
    const double predicted_x = focal_x * xp + c_x;
    const double predicted_y = focal_y * yp + c_y;


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;
    //std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';



    //Blanco's notation
    double fx=focal_x;
    double fy=focal_y;
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

    jacobian_projection(0,0)=fx/gz;
    jacobian_projection(0,1)=0.0;
    jacobian_projection(0,2)=-gx/gz2*fx;

    jacobian_projection(1,0)=0.0;
    jacobian_projection(1,1)=fy/gz;
    jacobian_projection(1,2)=-gy/gz2*fy;


    if ( jacobians != NULL )
    {


      const double z_2 = p3d_local[2]*p3d_local[2];
      if ( jacobians[0] != NULL) {

        // //jacobian of the pose (2x7)
        // MatrixRef j_eigen(jacobians[0], 2, 7);
        // j_eigen.setZero();
        // Eigen::Vector3d g={p3d_only_rot[0],p3d_only_rot[1],p3d_only_rot[2]};
        // Eigen::Matrix3d g_hat;
        // g_hat << 0.0, -g(2), g(1),
        //         g(2), 0.0, -g(0),
        //         -g(1), g(0), 0.0;
        //
        // Eigen::MatrixXd right_side(3,7);
        // right_side.setZero();
        // right_side.block(0,0,3,3)=-g_hat*2;  //TODO why is it scaled down by a half?
        // right_side.block(0,3,3,3)=Eigen::Matrix3d::Identity();
        //
        // j_eigen=jacobian_projection*right_side;



        // //Without eigen matrices
        // block of rotation x
        jacobians[0][0]= -focal_x *   p3d_local[0]*p3d_only_rot[1]/z_2*2; // -fx*gx*gy/gz2
        jacobians[0][1]=  focal_x/p3d_local[2]*p3d_only_rot[2]*2 + focal_x*p3d_local[0]/z_2*p3d_only_rot[0]*2;
        jacobians[0][2]= -focal_x/p3d_local[2]*p3d_only_rot[1]*2;

        // translation x
        jacobians[0][3]=  focal_x /   p3d_local[2];
        jacobians[0][4]= 0.0;
        jacobians[0][5]= -focal_x *   p3d_local[0] / z_2;

        // block or rotation y
        jacobians[0][7]= -focal_y/p3d_local[2]*p3d_only_rot[2]*2 - focal_y*p3d_local[1]/z_2*p3d_only_rot[1]*2;
        jacobians[0][8]= focal_y *   p3d_local[1]*p3d_only_rot[0]/z_2*2;//fy*gx*gy/gz2
        jacobians[0][9]=  focal_y *   p3d_only_rot[0]/p3d_local[2]*2; // fy*gx/gz

        //block of translation y
        jacobians[0][10]= 0.0;
        jacobians[0][11]= focal_y /   p3d_local[2];
        jacobians[0][12]=-focal_y *   p3d_local[1] / z_2;

        //last two zeros
        jacobians[0][6] = 0.0;
        jacobians[0][13]= 0.0;




      }

      if ( jacobians[1] != NULL) {
         //precompute stuff
         //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[0][0...5]
       jacobian_point3d= jacobian_projection*R;
       jacobians[1][0] = jacobian_point3d(0,0);
       jacobians[1][1] = jacobian_point3d(0,1);
       jacobians[1][2] = jacobian_point3d(0,2);

       jacobians[1][3] = jacobian_point3d(1,0);
       jacobians[1][4] = jacobian_point3d(1,1);
       jacobians[1][5] = jacobian_point3d(1,2);





      }
    }

    // exit(1);
    return true;
  }

private:
  double observed_x;
  double observed_y;
  const double * m_intrinsics;
};

// Automatic Diff:
class ErrorAutoDiffOptK {//: public SizedCostFunction<2, /* number of residuals */
                         //    7, /* size of first parameter, corresponding to the camera pose */
                         //    3 /* point size */ >{

 public:
   ErrorAutoDiffOptK(const double* point2D, const double * intrinsics)
   : observed_x(point2D[0]), observed_y(point2D[1]), m_intrinsics ( intrinsics )
   {}
  virtual ~ErrorAutoDiffOptK() {}

  template <typename T>
  bool operator()(const T* qtvec, const T* p3d_global, T* residuals) const {
    T p3d_only_rot[3]; //resulting point after only rotation but NOT translating
    T p3d_local[3]; //resulting point after rotation and translation in the local camera frame
    //const T * p3d_global = parameters[1];//double p3d_global[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points in world frame
    //const T * qvec = parameters[0];// double qvec[4] = { parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3] };  //Swap the w if it's an eigen quaternion
    const T* qvec = qtvec;

    //rotate
    QuaternionRotatePoint(qvec,p3d_global,p3d_only_rot);

    // camera[4,5,6] are the translation, first part is quaternion
    p3d_local[0] = p3d_only_rot[0]+ qtvec[4];
    p3d_local[1] = p3d_only_rot[1]+ qtvec[5];
    p3d_local[2] = p3d_only_rot[2]+ qtvec[6];

    //Normalize to image plane.
    const T xp = p3d_local[0] / p3d_local[2];
    const T yp = p3d_local[1] / p3d_local[2];

    //NO DISTORTION
    // const double distortion=1.0;

    // Compute final projected point position.
    const T focal_x = T(m_intrinsics[0]);
    const T focal_y = T(m_intrinsics[1]);
    const T c_x = T(m_intrinsics[2]);
    const T c_y = T(m_intrinsics[3]);
    const T predicted_x = focal_x * xp + c_x;
    const T predicted_y = focal_y * yp + c_y;


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - T(observed_x);
    residuals[1] = predicted_y - T(observed_y);
    //std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';

    // exit(1);
    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const double* point2D, const double* intrinsics) {
    return (new ceres::AutoDiffCostFunction<ErrorAutoDiffOptK,
              2, //2 residuals
              7,  //7 for rotation and translation
              3  //3 for the 3D point
              >( new ErrorAutoDiffOptK( point2D, intrinsics)));
  }
private:
  double observed_x;
  double observed_y;
  const double * m_intrinsics;
};



} //namespace ceres
