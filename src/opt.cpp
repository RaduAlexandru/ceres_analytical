#include <cmath>
#include <cstdio>
#include <iostream>
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "bal_problem.h"
#include "ceres_analytical/snavely_reprojection_error.h"

#include <chrono>

using namespace ceres::examples;
using namespace ceres;


void BuildProblem(BALProblem* bal_problem, Problem* problem);




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
  BALProblem bal_problem(argv[1],false); //Second argument is for using quaternions or not


  //create problem
  Problem problem(problem_options);
  bal_problem.Normalize(); //TODO check if necesary
  BuildProblem(&bal_problem, &problem);





  //solve
  std::chrono::high_resolution_clock::time_point t1_opt = std::chrono::high_resolution_clock::now();
  ceres::Solve(solver_options, &problem, &summary);
  std::chrono::high_resolution_clock::time_point t2_opt = std::chrono::high_resolution_clock::now();


  //report summaries and time
  auto duration_opt = std::chrono::duration_cast<std::chrono::microseconds>( t2_opt - t1_opt ).count();
  std::cout << "total took " << duration_opt * 1e-6<< " s " << std::endl;
  // std::cout << summary.FullReport() << "\n";

  //report final parameters


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
    cost_function = SnavelyReprojectionError::Create(observations[2 * i + 0],
                                                     observations[2 * i + 1]);


    // Each observation correponds to a pair of a camera and a point
    // which are identified by camera_index()[i] and point_index()[i]
    // respectively.
    double* camera = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    problem->AddResidualBlock(cost_function, NULL, camera, point);
  }
}
