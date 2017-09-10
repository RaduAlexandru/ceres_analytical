#pragma once

#include <vector>
#include <iostream>

#include "ceres/local_parameterization.h"
#include "ceres/internal/port.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/internal/disable_warnings.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/fixed_array.h"
#include "ceres/rotation.h"
#include "glog/logging.h"


namespace ceres{

class CERES_EXPORT CustomLocalParameterization  {
 public:
  virtual ~CustomLocalParameterization();

  // Generalization of the addition operation,
  //
  //   x_plus_delta = Plus(x, delta)
  //
  // with the condition that Plus(x, 0) = x.
  virtual bool Plus(const double* x,
                    const double* delta,
                    double* x_plus_delta) const = 0;

  // The jacobian of Plus(x, delta) w.r.t delta at delta = 0.
  //
  // jacobian is a row-major GlobalSize() x LocalSize() matrix.
  virtual bool ComputeJacobian(const double* x, double* jacobian) const = 0;

  // local_matrix = global_matrix * jacobian
  //
  // global_matrix is a num_rows x GlobalSize  row major matrix.
  // local_matrix is a num_rows x LocalSize row major matrix.
  // jacobian(x) is the matrix returned by ComputeJacobian at x.
  //
  // This is only used by GradientProblem. For most normal uses, it is
  // okay to use the default implementation.
  virtual bool MultiplyByJacobian(const double* x,
                                  const int num_rows,
                                  const double* global_matrix,
                                  double* local_matrix) const;

  // Size of x.
  virtual int GlobalSize() const = 0;

  // Size of delta.
  virtual int LocalSize() const = 0;
};

}
