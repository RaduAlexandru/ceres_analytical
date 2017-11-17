// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres_analytical/ba_problem.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Core"
#include "ceres/rotation.h"
#include "glog/logging.h"
#include "ceres_analytical/random.h"

namespace ceres {
namespace examples {
namespace {
typedef Eigen::Map<Eigen::VectorXd> VectorRef;
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;

void PerturbPoint3(const double sigma, double* point) {
  for (int i = 0; i < 3; ++i) {
    point[i] += RandNormal() * sigma;
  }
}

double Median(std::vector<double>* data) {
  int n = data->size();
  std::vector<double>::iterator mid_point = data->begin() + n / 2;
  std::nth_element(data->begin(), mid_point, data->end());
  return *mid_point;
}

}  // namespace

BAProblem::BAProblem( bool use_quaternions ) {
  use_quaternions_ = use_quaternions;
  
  num_cameras_ = 50;
  num_points_ = 100;
  num_observations_ = num_cameras_ * num_points_;
  
  VLOG(1) << "Header: " << num_cameras_
          << " " << num_points_
          << " " << num_observations_;

  point_index_ = new int[num_observations_];
  camera_index_ = new int[num_observations_];
  observations_ = new double[2 * num_observations_];

  num_parameters_ = camera_block_size() * num_cameras_ + point_block_size() * num_points_;
  parameters_ = new double[num_parameters_];
  
  double * points = parameters_ + camera_block_size() * num_cameras_;
  for ( int i = 0; i < num_points_; ++i )
  {
     double * point = points + point_block_size() * i;
     for ( int j = 0; j < point_block_size(); ++j)
       point[j] = 5*RandNormal();
  }
  double * cams = parameters_;
  for ( int i = 0; i < num_cameras_; ++i )
  {
     double * camera = cams + camera_block_size() * i;
     double * rot = camera;
     double n = 0; // normalizer
     for ( int j = 0; j < rotation_block_size(); ++j )
     {
        rot[j] = RandNormal();
        n += rot[j]*rot[j];
     }
     n = sqrt( n );
     for ( int j = 0; j < rotation_block_size(); ++j )
        rot[j] /= n;
     
     double * trans = camera + rotation_block_size();
     for ( int j = 0; j < NUM_PARAMS_FOR_TRANSLATION; ++j )
     {
        trans[j] = 2*RandNormal();
     }
     double * intrinsics = trans + NUM_PARAMS_FOR_TRANSLATION;
     intrinsics[0] = 640 + RandNormal();
     intrinsics[1] = 640 + RandNormal();
     intrinsics[2] = 320 + RandNormal();
     intrinsics[3] = 240 + RandNormal();
  }

  for (int i = 0; i < num_observations_; ++i) {
    camera_index_[i] = i % num_cameras_;
    point_index_[i] = i % num_points_;
    
    double * point = points + point_block_size() * point_index_[i];
    double * camera = cams + camera_block_size() * camera_index_[i];
    double * trans = camera + rotation_block_size();
    double * intrinsics = trans + NUM_PARAMS_FOR_TRANSLATION;
    
    double p3d_only_rot[3];
    double p3d_local[3];
    QuaternionRotatePoint(camera,point,p3d_only_rot);
    p3d_local[0] = p3d_only_rot[0] + camera[4];
    p3d_local[1] = p3d_only_rot[1] + camera[5];
    p3d_local[2] = p3d_only_rot[2] + camera[6];

    //Normalize to image plane.
    const double xp = p3d_local[0] / p3d_local[2];
    const double yp = p3d_local[1] / p3d_local[2];
    const double predicted_x = intrinsics[0] * xp + intrinsics[2];
    const double predicted_y = intrinsics[1] * yp + intrinsics[3];
    
    double obs [2] = { predicted_x, predicted_y};
    for (int j = 0; j < 2; ++j) {
      observations_ [2*i + j] = obs[j]+RandNormal();
    }
  }
}

void BAProblem::CameraToAngleAxisAndCenter(const double* camera,
                                            double* angle_axis,
                                            double* center) const {
  VectorRef angle_axis_ref(angle_axis, 3);
  if (use_quaternions_) {
    QuaternionToAngleAxis(camera, angle_axis);
  } else {
    angle_axis_ref = ConstVectorRef(camera, 3);
  }

  // c = -R't
  Eigen::VectorXd inverse_rotation = -angle_axis_ref;
  AngleAxisRotatePoint(inverse_rotation.data(),
                       camera + rotation_block_size(),
                       center);
  VectorRef(center, 3) *= -1.0;
}

void BAProblem::AngleAxisAndCenterToCamera(const double* angle_axis,
                                            const double* center,
                                            double* camera) const {
  ConstVectorRef angle_axis_ref(angle_axis, 3);
  if (use_quaternions_) {
    AngleAxisToQuaternion(angle_axis, camera);
  } else {
    VectorRef(camera, 3) = angle_axis_ref;
  }

  // t = -R * c
  AngleAxisRotatePoint(angle_axis,
                       center,
                       camera + rotation_block_size() );
  VectorRef(camera + rotation_block_size(), 3) *= -1.0;
}


void BAProblem::Normalize() {
  // Compute the marginal median of the geometry.
  std::vector<double> tmp(num_points_);
  Eigen::Vector3d median;
  double* points = mutable_points();
  for (int i = 0; i < point_block_size(); ++i) {
    for (int j = 0; j < num_points_; ++j) {
      tmp[j] = points[point_block_size() * j + i];
    }
    median(i) = Median(&tmp);
  }

  for (int i = 0; i < num_points_; ++i) {
    VectorRef point(points + point_block_size() * i, 3);
    tmp[i] = (point - median).lpNorm<1>();
  }

  const double median_absolute_deviation = Median(&tmp);

  // Scale so that the median absolute deviation of the resulting
  // reconstruction is 100.
  const double scale = 100.0 / median_absolute_deviation;

  VLOG(2) << "median: " << median.transpose();
  VLOG(2) << "median absolute deviation: " << median_absolute_deviation;
  VLOG(2) << "scale: " << scale;

  // X = scale * (X - median)
  for (int i = 0; i < num_points_; ++i) {
    VectorRef point(points + point_block_size() * i, 3);
    point = scale * (point - median);
  }

  double* cameras = mutable_cameras();
  double angle_axis[3];
  double center[3];
  for (int i = 0; i < num_cameras_; ++i) {
    double* camera = cameras + camera_block_size() * i;
    CameraToAngleAxisAndCenter(camera, angle_axis, center);
    // center = scale * (center - median)
    VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
    AngleAxisAndCenterToCamera(angle_axis, center, camera);
  }
}

void BAProblem::Perturb(const double rotation_sigma,
                         const double translation_sigma,
                         const double point_sigma) {
  CHECK_GE(point_sigma, 0.0);
  CHECK_GE(rotation_sigma, 0.0);
  CHECK_GE(translation_sigma, 0.0);

  double* points = mutable_points();
  if (point_sigma > 0) {
    for (int i = 0; i < num_points_; ++i) {
      PerturbPoint3(point_sigma, points + point_block_size() * i);
    }
  }

  for (int i = 0; i < num_cameras_; ++i) {
    double* camera = mutable_cameras() + camera_block_size() * i;

    double angle_axis[3];
    double center[3];
    // Perturb in the rotation of the camera in the angle-axis
    // representation.
    CameraToAngleAxisAndCenter(camera, angle_axis, center);
    if (rotation_sigma > 0.0) {
      PerturbPoint3(rotation_sigma, angle_axis);
    }
    AngleAxisAndCenterToCamera(angle_axis, center, camera);

    if (translation_sigma > 0.0) {
      PerturbPoint3(translation_sigma, camera + rotation_block_size() );
    }
  }
}

BAProblem::~BAProblem() {
  delete []point_index_;
  delete []camera_index_;
  delete []observations_;
  delete []parameters_;
}

}  // namespace examples
}  // namespace ceres
