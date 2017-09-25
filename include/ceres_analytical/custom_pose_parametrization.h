#pragma once

#include <vector>
#include <iostream>
#include "ceres/internal/port.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/internal/disable_warnings.h"
#include "ceres/local_parameterization.h"


#include <algorithm>
#include "Eigen/Geometry"
// #include "ceres/householder_vector.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/fixed_array.h"
#include "ceres/rotation.h"
#include "glog/logging.h"

namespace ceres{

  class  CustomPoseParameterization : public LocalParameterization {
   public:
    virtual ~CustomPoseParameterization() {}
    virtual bool Plus(const double* x,
                      const double* delta,
                      double* x_plus_delta) const;
    virtual bool ComputeJacobian(const double* x,
                                 double* jacobian) const;
    virtual int GlobalSize() const { return 7; }
    virtual int LocalSize() const { return 6; }
  };

}
