#include "ceres/sized_cost_function.h"


namespace ceres{

class ProjectionAnalyticalError : public SizedCostFunction<2 /* number of residuals */,
                             9 /* size of first parameter, corresponding to the camera */,
                             3 /* size of the second parameter corresponding to the point*/  > {
 public:
   ProjectionAnalyticalError(double observed_x, double observed_y) : observed_x(observed_x), observed_y(observed_y) {}
  virtual ~ProjectionAnalyticalError() {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    // double w_x = parameters[0][0];
    // double w_x = parameters[0][0];
    // double w_x = parameters[0][0];

    double p[3]; //resulting point after rotation
    double pt[3]={parameters[1][0],parameters[1][1],parameters[1][2]};  //3D points that needs to be rotated
    // AngleAxisRotatePoint(parameters[0], parameters[1], p);

    //rotate
    double angle_axis[3];
    angle_axis[0]=parameters[0][0];
    angle_axis[1]=parameters[0][1];
    angle_axis[2]=parameters[0][2];
    double theta2 = angle_axis[0]*angle_axis[0] + angle_axis[1]*angle_axis[1] + angle_axis[2]*angle_axis[2];

    if (theta2 > (std::numeric_limits<double>::epsilon())) {
      double theta = sqrt(theta2);
      double costheta = cos(theta);
      double sintheta = sin(theta);
      double theta_inverse = 1.0 / theta;
      double w[3] = { angle_axis[0] * theta_inverse,
                      angle_axis[1] * theta_inverse,
                      angle_axis[2] * theta_inverse };
      double w_cross_pt[3] = { w[1] * pt[2] - w[2] * pt[1],
                               w[2] * pt[0] - w[0] * pt[2],
                               w[0] * pt[1] - w[1] * pt[0] };
      double tmp = (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (1.0 - costheta);
      p[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
      p[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
      p[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
    }else{
      std::cout << "angle is very small " << '\n';
      std::cout << "Implement the Taylor expansion near zero" << '\n';
    }



    // camera[3,4,5] are the translation.
    p[0] += parameters[0][3];
    p[1] += parameters[0][4];
    p[2] += parameters[0][5];


    //dehomogeneous
    double xp = - p[0] / p[2];
    double yp = - p[1] / p[2];


    //NO DISTORSION
    //TODO
    double distortion=1.0;

    // Compute final projected point position.
    double focal =parameters[0][6];
    double predicted_x = focal * distortion * xp;
    double predicted_y = focal * distortion * yp;


    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;
    std::cout << " residual is " << residuals[0] << ", " << residuals[1] << '\n';


    if (jacobians != NULL && jacobians[0] != NULL) {






      //jacobian of camera is of size 2x9 , fill in a row majow order jacobians[0][0...17]
      jacobians[0][0] = 0;

      //jacobian of point is of size 2x3 , fill in a row majow order jacobians[0][0...7]
    }
    return true;
  }

private:
  double observed_x;
  double observed_y;
};

} //namespace ceres
