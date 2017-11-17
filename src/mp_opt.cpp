#include <cmath>
#include <cstdio>
#include <iostream>
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres_analytical/ba_problem.h"

#include "ceres_analytical/mp_error_analytical.h"
#include "ceres_analytical/custom_pose_parametrization.h"

#include <chrono>

#include <Eigen/Dense>

using namespace ceres::examples;
using namespace ceres;


void BuildProblem(BAProblem* bal_problem, Problem* problem);
void BuildProblemAnalytical(BAProblem* bal_problem, Problem* problem);
void CompareErrorFunctions ( BAProblem * ba_AD, BAProblem * ba_AN );
void CompareResidualBlocks ( Problem * p_AD, Problem * p_AN );

int main(int argc, char** argv) {
  ceres::Solver::Summary summary;
  google::InitGoogleLogging(argv[0]);

  //Solver options
  ceres::Solver::Options solver_options;
  solver_options.linear_solver_type = ceres::DENSE_SCHUR;//DENSE_QR;
  solver_options.minimizer_progress_to_stdout = true;
  solver_options.max_num_iterations = 1000;
  // solver_options.use_explicit_schur_complement = true;
  // solver_options.num_threads=4;

  //problem options
  Problem::Options problem_options;
  problem_options.disable_all_safety_checks=true; //Gain 5% in speed
  
  std::cout << "creating new BA problem." << std::endl;
  //READ file
  BAProblem ba_problem(true); //Second argument is for using quaternions or not
  std::cout << "created new BA problem." << std::endl;
  ba_problem.Normalize(); //TODO check
  std::cout << "normalized new BA problem." << std::endl;
  BAProblem ba_problem_analytic ( ba_problem ); //Second argument is for using quaternions or not
  std::cout << "deepCopied new BA problem." << std::endl;
  //ba_problem_analytic.Normalize();
  std::cout << "normalized copied BA problem." << std::endl;
  
  //create ceres problem
  Problem problem(problem_options);
  Problem problem_analytic(problem_options);
  BuildProblem(&ba_problem, &problem);
  BuildProblemAnalytical(&ba_problem_analytic, &problem_analytic);
  
  //CompareErrorFunctions ( &ba_problem, &ba_problem_analytic);
  
  CompareResidualBlocks( &problem, & problem_analytic);
  
  //solve
  std::chrono::high_resolution_clock::time_point t1_opt = std::chrono::high_resolution_clock::now();
  ceres::Solve(solver_options, &problem, &summary);
  std::chrono::high_resolution_clock::time_point t2_opt = std::chrono::high_resolution_clock::now();

  //report summaries and time
  auto duration_opt = std::chrono::duration<double>( t2_opt - t1_opt ).count();
  std::cout << summary.FullReport() << "\n";
  std::cout << "total took " << duration_opt << " s " << std::endl;
  
  //solve
  t1_opt = std::chrono::high_resolution_clock::now();
  ceres::Solve(solver_options, &problem_analytic, &summary);
  t2_opt = std::chrono::high_resolution_clock::now();

  //report summaries and time
  duration_opt = std::chrono::duration<double>( t2_opt - t1_opt ).count();
  std::cout << summary.FullReport() << "\n";
  std::cout << "total took " << duration_opt << " s " << std::endl;
  
  return 0;
}

/**
 * Computes the function `Y = f(X)` with `X = [x, y, z]` and `Y = [u, v]`.
 *
 * This function also computes the Jacobian `J` of `[u, v] = f(x, y, z)` as
 *  J = | du/dx  du/dy  du/dz|
 *      | dv/dx  dv/dy  dv/dz|.
 */
inline bool autoDiff(const CostFunction * func,
                     const double * cams,
                     const double * pt,
                     Eigen::Vector2d & Y,
                     Eigen::Matrix<double, 2, 7>& J, Eigen::Matrix<double,2,3> & Jp)
{
    // Wrap input point double** as required by AutoDiff<> class
    const double* parametersX = cams;
    const double ** parameters = new const double*[2];
    parameters[0] = cams;
    parameters[1] = pt;
    

    // Init memory for Jacobian
    double** jacobians = new double*[2];
    jacobians[0] = new double[14];
    jacobians[1] = new double[6];
    
    func->Evaluate( parameters, Y.data(), jacobians);
    
    // Copy to output variable J
    {
       double * jd = J.data();
       for(size_t j=0; j<14; ++j)
          jd[j] = jacobians[0][j];
       double * jp = Jp.data();
       for(size_t j=0; j<6; ++j)
          jp[j] = jacobians[1][j];
    }

    for(size_t i=0; i<2; i++)
        delete[] jacobians[i];
    delete[] jacobians;
    delete[] parameters;
    
    return true;
}

std::vector<CostFunction * > autodiffCF;
std::vector<ResidualBlockId> autodiffRB;
void BuildProblem(BAProblem* bal_problem, Problem* problem) {
  autodiffCF.clear();
  autodiffRB.clear();
  const int point_block_size = bal_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = bal_problem->camera_block_size();  //camera is 7+4 parameters
  double* points = bal_problem->mutable_points();
  double* cameras = bal_problem->mutable_cameras();
  const double* observations = bal_problem->observations();
  LossFunction * loss_function = NULL; // new HuberLoss(1);
  for (int i = 0; i < bal_problem->num_observations(); ++i) {
    
    const double * cam_intrinsics = cameras + camera_block_size * bal_problem->camera_index()[i] + 7;  //move +7 because of the quat paramerization     
    const double * cam_observation = observations + 2 * i;
    
    CostFunction* cost_function;
    // Each Residual block takes a point and a camera as input and
    // outputs a 2 dimensional residual.
    cost_function = ErrorAutoDiffOptK::Create(cam_observation, cam_intrinsics);

    autodiffCF.push_back( cost_function );
    
    // Each observation correponds to a pair of a camera and a point
    // which are identified by camera_index()[i] and point_index()[i]
    // respectively.
    double* camera = cameras + camera_block_size * bal_problem->camera_index()[i];
    double* point = points + point_block_size * bal_problem->point_index()[i];
    ResidualBlockId rb = problem->AddResidualBlock(cost_function, loss_function, camera, point);
    autodiffRB.push_back(rb);
    
    problem->SetParameterBlockConstant( point );
  }

  //Set parametrization for quaternion
  LocalParameterization* camera_parameterization =
    new ProductParameterization( new QuaternionParameterization(), new IdentityParameterization(3) );
  for (int i = 0; i < bal_problem->num_cameras(); ++i) {
    problem->SetParameterization(cameras + camera_block_size * i, camera_parameterization);
  }
}
std::vector<CostFunction * > analyticalCF;
std::vector<ResidualBlockId> analyticalRB;
void BuildProblemAnalytical(BAProblem* ba_problem, Problem* problem){
  analyticalCF.clear();
  analyticalRB.clear();
  const int point_block_size = ba_problem->point_block_size();    //point is 3 parameters
  const int camera_block_size = ba_problem->camera_block_size();  //camera is 7+4 parameters
  double* points = ba_problem->mutable_points();
  double* cameras = ba_problem->mutable_cameras();
  const double* observations = ba_problem->observations();
  LossFunction * loss_function = NULL; //new HuberLoss(1);
  for (int i = 0; i < ba_problem->num_observations(); ++i) {
    const double* cam_intrinsics = cameras + camera_block_size * ba_problem->camera_index()[i] + 7;  //move +7 because of the quat paramerization
    const double * cam_observation = observations + 2 * i;
    
    CostFunction* cost_function;
    cost_function = new ErrorAnalyticalOptK( cam_observation, cam_intrinsics ); // Optimized code for analytical error without intrinsics
    //cost_function = ErrorAnalyticalWithDistortion::Create(obs);  //distortion is autodiff, rest is analytical

    analyticalCF.push_back( cost_function );
    // // Each observation correponds to a pair of a camera and a point
    // // which are identified by camera_index()[i] and point_index()[i]
    // // respectively.
    double* camera = cameras + camera_block_size * ba_problem->camera_index()[i];
    double* point = points + point_block_size * ba_problem->point_index()[i];
    ResidualBlockId rb = problem->AddResidualBlock(cost_function, loss_function, camera, point);
    analyticalRB.push_back(rb);

    problem->SetParameterBlockConstant(point);
    // problem->SetParameterBlockConstant(cam_pose);
  }

  //Set parametrization for the pose
  LocalParameterization* camera_parameterization =  new CustomPoseParameterization( );
  for (int i = 0; i < ba_problem->num_cameras(); ++i) {
    problem->SetParameterization(cameras + camera_block_size * i, camera_parameterization);
  }
}

void CompareErrorFunctions ( BAProblem * ba_AD, BAProblem * ba_AN )
{
   if ( analyticalCF.size() != autodiffCF.size() )
   {
      std::cerr << "sizes are not equal!"<< std::endl;
   }
   for ( int i = 0; i < analyticalCF.size() && i < autodiffCF.size(); ++i)
   {
      const double* cams1 = ba_AD->mutable_cameras() + ba_AD->camera_block_size() * ba_AD->camera_index()[i];
      const double* pt1 = ba_AD->mutable_points() + ba_AD->point_block_size() * ba_AD->point_index()[i];
      const double* cams2 = ba_AN->mutable_cameras() + ba_AN->camera_block_size() * ba_AN->camera_index()[i];
      const double* pt2 = ba_AN->mutable_points() + ba_AN->point_block_size() * ba_AN->point_index()[i];
    
      Eigen::Vector2d resAD, resAN;
      Eigen::Matrix<double,2,7> JAD, JAN;
      Eigen::Matrix<double,2,3> JpAD, JpAN;
      autoDiff(autodiffCF[i], cams1, pt1, resAD,JAD, JpAD );
      autoDiff(analyticalCF[i], cams2, pt2, resAN, JAN, JpAN);
      //std::cout << "resAD=["<<resAD.transpose()<<"] resAN=[" << resAN.transpose() << "]" << std::endl;
      std::cout << "JAD=["<< JAD.row(0) <<";"<<JAD.row(1)<<"]"<<std::endl<<"JAN=[" << JAN.row(0) <<";"<<JAN.row(1) << "]" << std::endl;
      std::cout << "JADp=["<< JAD.data()[0]<<";"<< JAD.data()[1]<<";"<< JAD.data()[2]<<";"<< JAD.data()[3]<<";"<< JAD.data()[4]<<";"<< JAD.data()[5]<<";"<< JAD.data()[6]<<";"<<JAD.data()[7]<<";"<< JAD.data()[8]<<";"<< JAD.data()[9]<<";"<< JAD.data()[10]<<";"<< JAD.data()[11]<<";"<< JAD.data()[12]<<";"<< JAD.data()[13]<<"];"<<std::endl
       << "JANp=["<< JAN.data()[0]<<";"<< JAN.data()[1]<<";"<< JAN.data()[2]<<";"<< JAN.data()[3]<<";"<< JAN.data()[4]<<";"<< JAN.data()[5]<<";"<< JAN.data()[6]<<";"<<JAN.data()[7]<<";"<< JAN.data()[8]<<";"<< JAN.data()[9]<<";"<< JAN.data()[10]<<";"<< JAN.data()[11]<<";"<< JAN.data()[12]<<";"<< JAN.data()[13]<<"];"<<std::endl;
      //std::cout << "JpAD=["<< JpAD.row(0) <<";"<<JpAD.row(1)<<"]"<<std::endl<<"JpAN=[" << JpAN.row(0) <<";"<<JpAN.row(1) << "]" << std::endl;
      
      std::cout << "dR="<<(resAD-resAN).norm() << " dJp=" <<(JpAD-JpAN).norm()<< std::endl;
   }
}

void CompareResidualBlocks ( Problem * p_AD, Problem * p_AN )
{
   if ( analyticalRB.size() != autodiffRB.size() )
   {
      std::cerr << "sizes are not equal!"<< std::endl;
   }
   std::vector<double> residualAD, residualAN;
   std::vector<double> gradientAD, gradientAN;
   CRSMatrix jacobianAD, jacobianAN;
   Problem::EvaluateOptions eo_AD;
   eo_AD.residual_blocks = autodiffRB;
   p_AD->Evaluate(eo_AD, NULL, &residualAD, &gradientAD, &jacobianAD);
   Problem::EvaluateOptions eo_AN;
   eo_AN.residual_blocks = analyticalRB;
   p_AN->Evaluate(eo_AN, NULL, &residualAN, &gradientAN, &jacobianAN);
   jacobianAD.values;
   for ( int i = 0; i < analyticalRB.size() && i < autodiffRB.size() && i < 10; ++i)
   {
      std::cout << "i=" << i<<std::endl;
      
      std::stringstream ss;
      ss << "AD=" << residualAD.size() << " [" << residualAD[2*i] <<";" << residualAD[2*i+1] << "] ";
      ss << "gAD=" << gradientAD.size() << " [" << gradientAD[i] <<"];";
      std::cout << ss.str() << std::endl;
      std::stringstream ss1;
      ss1<<"AN=" << residualAN.size() << " [" << residualAN[2*i] <<";" << residualAN[2*i+1] << "] ";
      ss1<<"gAN=" << gradientAN.size() << " [" << gradientAN[i] << "];";
      std::cout << ss1.str() << std::endl;
      
      std::stringstream ss2,ss3;
      ss2<<"c=["; ss3 <<"v=[";
      for ( int idx = jacobianAD.rows[2*i]; idx < jacobianAD.rows[2*(i+1)]; ++idx)
      {
         ss2 <<jacobianAD.cols[idx];
         ss3 << jacobianAD.values[idx];
         if ( idx < jacobianAD.rows[2*(i+1)] -1 )
         {
            ss2<<";"; ss3<<";";
         }
      }
      ss2 << "];"; ss3<<"];";
      std::cout <<"jAD "<< jacobianAD.rows[2*i]<<" "<<jacobianAD.rows[2*(i+1)]<<" " << ss2.str() <<" "<< ss3.str()   << std::endl;

      std::stringstream ss4,ss5;
      ss4<<"c=["; ss5 <<"v=[";
      for ( int idx = jacobianAN.rows[2*i]; idx < jacobianAN.rows[2*(i+1)]; ++idx)
      {
         ss4 <<jacobianAN.cols[idx];
         ss5 << jacobianAN.values[idx];
         if ( idx < jacobianAN.rows[2*(i+1)] -1 )
         {
            ss4<<";"; ss5<<";";
         }
      }
      ss4 << "];"; ss5<<"];";
      std::cout <<"jAN= "<< jacobianAN.rows[2*i]<<" "<<jacobianAN.rows[2*(i+1)]<<" " << ss4.str() <<" "<< ss5.str()   << std::endl;

      
//       Eigen::Vector2d resAD, resAN;
//       Eigen::Matrix<double,2,7> JAN, JAD;
//       Eigen::Matrix<double,2,3> JpAD, JpAN;
//       
//       //std::cout << "resAD=["<<resAD.transpose()<<"] resAN=[" << resAN.transpose() << "]" << std::endl;
//       std::cout << "JAD=["<< JAD.row(0) <<";"<<JAD.row(1)<<"]"<<std::endl<<"JAN=[" << JAN.row(0) <<";"<<JAN.row(1) << "]" << std::endl;
//       //std::cout << "JpAD=["<< JpAD.row(0) <<";"<<JpAD.row(1)<<"]"<<std::endl<<"JpAN=[" << JpAN.row(0) <<";"<<JpAN.row(1) << "]" << std::endl;
//       
//       std::cout << "dR="<<(resAD-resAN).norm() << " dJp=" <<(JpAD-JpAN).norm() << " dJ="<<(JAD-JAN).norm()<< std::endl;
   }
   std::cout <<"J_AD=["<< jacobianAD.num_rows << "x" << jacobianAD.num_cols<<"] vals="<<jacobianAD.values.size()
      << " J_AN=["<< jacobianAN.num_rows << "x" << jacobianAN.num_cols<<"] vals="<<jacobianAN.values.size() << std::endl;      
}


// JAD=[ 651.767 -7214.61 -38.6311  83.6083 -18111.4   -26294 -38.7575;-2810.28 -13253.3        0  11765.3 -14397.2        0  215.778]
// JAN=[ -7690.5 -3553.38        0        0  7715.66        0  215.778; 3627.86 -38.6311  83.6083 -20488.2  1385.87 -38.7575        0]
