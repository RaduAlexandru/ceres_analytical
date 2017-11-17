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
    const double & fx = m_intrinsics[0];
    const double & fy = m_intrinsics[1];
    const double & c_x = m_intrinsics[2];
    const double & c_y = m_intrinsics[3];
    const double predicted_x = fx * xp + c_x;
    const double predicted_y = fy * yp + c_y;


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;
    //std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';


    //precompute stuff (this are the entries of the jacobian of the projection)
    double z_2=p3d_local[2]*p3d_local[2];
    double j_proj_00=fx/ p3d_local[2];
    double j_proj_11=fy/ p3d_local[2];
    double j_proj_02=-fx*p3d_local[0]/ z_2;
    double j_proj_12=-fy*p3d_local[1]/ z_2;


    if ( jacobians != NULL )
    {


      if ( jacobians[0] != NULL) {

        //precompute some more stuff
        double x_rot_2=p3d_only_rot[0]*2;
        double y_rot_2=p3d_only_rot[1]*2;
        double z_rot_2=p3d_only_rot[2]*2;

        //block of rotation x
        jacobians[0][0]=j_proj_02*y_rot_2;
        jacobians[0][1]=j_proj_00*z_rot_2 - j_proj_02*x_rot_2;
        jacobians[0][2]=-j_proj_00*y_rot_2;

        // translation x
        jacobians[0][3]= j_proj_00;
        jacobians[0][4]= 0.0;
        jacobians[0][5]= j_proj_02;


        // block of rotation y
        jacobians[0][7]= -j_proj_11*z_rot_2 + j_proj_12*y_rot_2;
        jacobians[0][8]= -j_proj_12*x_rot_2;
        jacobians[0][9]=  j_proj_11*x_rot_2;

        //block of translation y
        jacobians[0][10]= 0.0;
        jacobians[0][11]= j_proj_11;
        jacobians[0][12]= j_proj_12;

        //last two zeros
        jacobians[0][6] = 0.0;
        jacobians[0][13]= 0.0;


      }

      if ( jacobians[1] != NULL) {

         //precompute stuff
         double R_array[9];
         QuaternionToRotation(qvec,R_array);

         //jacobian of p3d is of size 2x3 , fill in a row majow order jacobians[0][0...5]
        jacobians[1][0]= j_proj_00*R_array[0] + j_proj_02*R_array[6];
        jacobians[1][1]= j_proj_00*R_array[1] + j_proj_02*R_array[7];
        jacobians[1][2]= j_proj_00*R_array[2] + j_proj_02*R_array[8];

        jacobians[1][3]= j_proj_11*R_array[3] + j_proj_12*R_array[6];
        jacobians[1][4]= j_proj_11*R_array[4] + j_proj_12*R_array[7];
        jacobians[1][5]= j_proj_11*R_array[5] + j_proj_12*R_array[8];





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
