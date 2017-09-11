#include <cmath>
#include <cstdio>
#include <iostream>
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres_analytical/bal_problem.h"
#include "ceres_analytical/snavely_reprojection_error.h"

#include "ceres_analytical/projection_analytical_error.h"
#include "ceres_analytical/reprojection_error_without_intrinsics.h"
#include "ceres_analytical/error_analytical.h"

#include "ceres_analytical/custom_pose_parametrization.h"
#include "ceres_analytical/custom_local_parametrization.h"

#include <chrono>

#include <Eigen/Dense>

using namespace ceres::examples;
using namespace ceres;


void BuildProblem(BALProblem* bal_problem, Problem* problem);
void BuildProblemAnalytical(BALProblem* bal_problem, Problem* problem);

void BuildProblemWithoutIntrinsics(BALProblem* bal_problem, Problem* problem);
void BuildProblemWithoutIntrinsicsAnalytical(BALProblem* bal_problem, Problem* problem);

void report_camera_parameters(BALProblem& bal_problem);




int main(int argc, char** argv) {
  ceres::Solver::Summary summary;
  google::InitGoogleLogging(argv[0]);
  if (argc != 2) { std::cerr << "usage: simple_bundle_adjuster <bal_problem>\n"; return 1; }

  //Solver options
  ceres::Solver::Options solver_options;
  solver_options.linear_solver_type = ceres::DENSE_SCHUR;
  solver_options.minimizer_progress_to_stdout = true;

  //problem options
  Problem::Options problem_options;
  // problem_options.disable_all_safety_checks=true; //Gain 5% in speed


  //READ file
  BALProblem bal_problem(argv[1],true); //Second argument is for using quaternions or not


  //create problem
  Problem problem(problem_options);
  bal_problem.Normalize(); //TODO check if necesary
  // BuildProblem(&bal_problem, &problem);
  // BuildProblemWithoutIntrinsics(&bal_problem, &problem);  //optimizes a problem only over poses and points, without the instrinsics of the camera

  BuildProblemWithoutIntrinsicsAnalytical(&bal_problem, &problem);

  // BuildProblemWithoutIntrinsicsSophus(&bal_problem, &problem);   //optimizes only poses and points, using sophus and eigen
  // BuildProblemWithoutIntrinsicsSophusAnalytical(&bal_problem, &problem);  //optimizes poses, points, using sophus with analytical jacobian


  // BuildProblemSophus(&bal_problem, &problem);  //builds problem using sophus
  // BuildProblemAnalytical(&bal_problem, &problem);


  //solve
  std::chrono::high_resolution_clock::time_point t1_opt = std::chrono::high_resolution_clock::now();
  ceres::Solve(solver_options, &problem, &summary);
  std::chrono::high_resolution_clock::time_point t2_opt = std::chrono::high_resolution_clock::now();


  //report summaries and time
  auto duration_opt = std::chrono::duration_cast<std::chrono::microseconds>( t2_opt - t1_opt ).count();
  std::cout << "total took " << duration_opt * 1e-6<< " s " << std::endl;
  std::cout << summary.FullReport() << "\n";


  //report camera parameters
  // report_camera_parameters(bal_problem);


  //Write to ply file
  bal_problem.WriteToPLYFile("result.ply");



  return 0;
}



void BuildProblem(BALProblem* bal_problem, Problem* problem) {
  const int point_block_size = bal_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = bal_problem->camera_block_size();  //camera is 9 parameters
  double* points = bal_problem->mutable_points();
  double* cameras = bal_problem->mutable_cameras();

  // Observations is 2*num_observations long array observations =
  // [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
  // and y positions of the observation.
  const double* observations = bal_problem->observations();
  for (int i = 0; i < bal_problem->num_observations(); ++i) {
    CostFunction* cost_function;
    // Each Residual block takes a point and a camera as input and
    // outputs a 2 dimensional residual.
    cost_function = SnavelyReprojectionErrorWithQuaternions::Create(observations[2 * i + 0],
                                                                    observations[2 * i + 1]);


    // Each observation correponds to a pair of a camera and a point
    // which are identified by camera_index()[i] and point_index()[i]
    // respectively.
    double* camera = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    problem->AddResidualBlock(cost_function, NULL, camera, point);


    // std::cout << "adding residual block" << '\n';


  }


  //Set parametrization for quaternion
  LocalParameterization* camera_parameterization =
  new ProductParameterization( new QuaternionParameterization(),
                               new IdentityParameterization(6));
  for (int i = 0; i < bal_problem->num_cameras(); ++i) {
    problem->SetParameterization(cameras + camera_block_size * i,
                                 camera_parameterization);

  }


}

void BuildProblemAnalytical(BALProblem* bal_problem, Problem* problem){

  const int point_block_size = bal_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = bal_problem->camera_block_size();  //camera is 9 parameters
  double* points = bal_problem->mutable_points();
  double* cameras = bal_problem->mutable_cameras();
  const double* observations = bal_problem->observations();

  for (int i = 0; i < bal_problem->num_observations(); ++i) {
    CostFunction* cost_function;
    cost_function = new ProjectionAnalyticalError(observations[2 * i + 0],
                                                  observations[2 * i + 1]);


    // // Each observation correponds to a pair of a camera and a point
    // // which are identified by camera_index()[i] and point_index()[i]
    // // respectively.
    double* camera = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    problem->AddResidualBlock(cost_function, NULL, camera, point);
  }

}

void BuildProblemWithoutIntrinsics(BALProblem* bal_problem, Problem* problem){

  const int point_block_size = bal_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = bal_problem->camera_block_size();  //camera is 9 parameters
  double* points = bal_problem->mutable_points();
  double* cameras = bal_problem->mutable_cameras();
  const double* observations = bal_problem->observations();

  for (int i = 0; i < bal_problem->num_observations(); ++i) {
    Eigen::Vector2d obs={observations[2 * i + 0], observations[2 * i + 1]};

    CostFunction* cost_function;
    cost_function = ReprojectionErrorWithoutIntrinsics::Create(obs);


    // // Each observation correponds to a pair of a camera and a point
    // // which are identified by camera_index()[i] and point_index()[i]
    // // respectively.
    double* cam_pose = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    double* cam_intrinsics = cameras + camera_block_size * bal_problem->camera_index()[i] + 7;  //move +7 because of the quat paramerization
    problem->AddResidualBlock(cost_function, NULL, cam_pose, point, cam_intrinsics);

    problem->SetParameterBlockConstant(cam_intrinsics);


    // std::cout << "num residuals " << cost_function->num_residuals() << '\n';
    // std::cout << "number of parameter blocks " << cost_function->parameter_block_sizes().size() << '\n';  //Is 3 because we have 3 parameter blocks (pose,point3d,intrisnics)
    // std::cout << "size of first parameter block " << cost_function->parameter_block_sizes()[0] << "\n";
    // std::cout << "size of second parameter block " << cost_function->parameter_block_sizes()[1] << "\n";
    // std::cout << "size of third parameter block " << cost_function->parameter_block_sizes()[2] << "\n";
    // std::cout  << '\n';

  }


  //TODO recheck if its correct
  //Set parametrization for quaternion
  LocalParameterization* camera_parameterization =
  new ProductParameterization( new QuaternionParameterization(),
                               new IdentityParameterization(3));
  for (int i = 0; i < bal_problem->num_cameras(); ++i) {
    problem->SetParameterization(cameras + camera_block_size * i,
                                 camera_parameterization);

  }

}

void BuildProblemWithoutIntrinsicsAnalytical(BALProblem* bal_problem, Problem* problem){

  const int point_block_size = bal_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = bal_problem->camera_block_size();  //camera is 9 parameters
  double* points = bal_problem->mutable_points();
  double* cameras = bal_problem->mutable_cameras();
  const double* observations = bal_problem->observations();

  for (int i = 0; i < bal_problem->num_observations(); ++i) {
    Eigen::Vector2d obs={observations[2 * i + 0], observations[2 * i + 1]};

    CostFunction* cost_function;
    // cost_function = ReprojectionErrorWithoutIntrinsics::Create(obs);
    // cost_function = new ErrorAnalytical(obs);
    cost_function = ErrorAnalyticalWithDistortion::Create(obs);


    // // Each observation correponds to a pair of a camera and a point
    // // which are identified by camera_index()[i] and point_index()[i]
    // // respectively.
    double* cam_pose = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    double* cam_intrinsics = cameras + camera_block_size * bal_problem->camera_index()[i] + 7;  //move +7 because of the quat paramerization
    problem->AddResidualBlock(cost_function, NULL, cam_pose, point, cam_intrinsics);

    // problem->SetParameterBlockConstant(point);
    // problem->SetParameterBlockConstant(cam_intrinsics);
    // problem->SetParameterBlockConstant(cam_pose);


    // std::cout << "num residuals " << cost_function->num_residuals() << '\n';
    // std::cout << "number of parameter blocks " << cost_function->parameter_block_sizes().size() << '\n';  //Is 3 because we have 3 parameter blocks (pose,point3d,intrisnics)
    // std::cout << "size of first parameter block " << cost_function->parameter_block_sizes()[0] << "\n";
    // std::cout << "size of second parameter block " << cost_function->parameter_block_sizes()[1] << "\n";
    // std::cout << "size of third parameter block " << cost_function->parameter_block_sizes()[2] << "\n";
    // std::cout  << '\n';

  }


  //Set parametrization for quaternion
  // LocalParameterization* camera_parameterization= new ProductParameterization( new QuaternionParameterization(),
  //                              new IdentityParameterization(3));
  LocalParameterization* camera_parameterization =  new CustomPoseParameterization( );
  for (int i = 0; i < bal_problem->num_cameras(); ++i) {
    problem->SetParameterization(cameras + camera_block_size * i,
                                 camera_parameterization);

  }

}



void report_camera_parameters(BALProblem& bal_problem){

  //report final parameters
  double* cameras = bal_problem.mutable_cameras();
  const int camera_block_size = bal_problem.camera_block_size();  //camera is 9 parameters
  for (size_t i = 0; i < bal_problem.num_cameras(); i++) {
    double* camera = cameras + camera_block_size * i;


    //get rotation part
    double rotation_mat_flattened[9];
    QuaternionToRotation(camera,rotation_mat_flattened);
    Eigen::Matrix3d R= Eigen::Map<Eigen::Matrix3d >(rotation_mat_flattened);

    //get translation part
    Eigen::Vector3d T ( camera[4], camera[5], camera[6]);

    // //combine them
    Eigen::Matrix4d Trans; // Your Transformation Matrix
    Trans.setIdentity();   // Set to Identity to make bottom row of Matrix 0,0,0,1
    Trans.block<3,3>(0,0) = R;
    Trans.block<3,1>(0,3)= T;

    std::cout << "Transformation matrix is" << std::endl << Trans << '\n';
  }



}
