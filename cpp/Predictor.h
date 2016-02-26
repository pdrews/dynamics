/**
 * @file Predictor.h
 * @brief Dynamics Predictor
 * @author Paul Drews
 */

#pragma once

#include <gtsam_unstable/base/dllexport.h>
#include <gtsam/geometry/Pose3.h>

#define X_DIM 4
#define U_DIM 2
#define THETA_DIM 100

namespace gtsam {

typedef Matrix<double,X_DIM,THETA_DIM> MatrixBasis;
typedef Matrix<double,THETA_DIM,1> MatrixTheta;
typedef Matrix<double,X_DIM,1> MatrixX;
typedef Matrix<double,U_DIM,1> MatrixU;

struct Predictor {
  Predictor(MatrixX x, MatrixU u, double heading, double dt);
  MatrixX operator()(MatrixTheta theta, OptionalJacobian<X_DIM,THETA_DIM> H);
private:
  MatrixBasis dt_basis_;
  MatrixX x_;

};

} // \namespace gtsam
