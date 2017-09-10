#include "ceres_analytical/custom_local_parametrization.h"

namespace ceres {

  using std::vector;

  CustomLocalParameterization::~CustomLocalParameterization() {
  }

  bool CustomLocalParameterization::MultiplyByJacobian(const double* x,
                                                 const int num_rows,
                                                 const double* global_matrix,
                                                 double* local_matrix) const {

    std::cout << "lol" << '\n';
    exit(1);
    Matrix jacobian(GlobalSize(), LocalSize());
    if (!ComputeJacobian(x, jacobian.data())) {
      return false;
    }

    MatrixRef(local_matrix, num_rows, LocalSize()) = ConstMatrixRef(global_matrix, num_rows, GlobalSize()) * jacobian;

    std::cout << "lol" << '\n';

    return true;
  }



}
