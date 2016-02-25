/**
 * @file Predictor.h
 * @brief Dynamics Predictor
 * @author Paul Drews
 */

#pragma once

#include <gtsam_unstable/base/dllexport.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/base/ProductLieGroup.h>

#define X_DIM 4
#define U_DIM 2
#define THETA_DIM 100

namespace gtsam {

typedef Matrix<double,4,100> MatrixBasis;
typedef Matrix<double,100,1> MatrixTheta;
typedef Matrix<double,4,1> MatrixX;
typedef Matrix<double,2,1> MatrixU;

struct Predictor {
  Predictor(Matrix31 x,double dt);
  Matrix31 operator()(Matrix51 theta, OptionalJacobian<3,5> H);
private:
  //std:vector<RBF> rbf;
  Matrix35 dt_rbf_;
  Matrix31 x_;

};

} // \namespace gtsam
