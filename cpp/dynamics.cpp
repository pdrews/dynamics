#include "dynamics.h"
#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam_unstable/dynamics/Predictor.h>

#include <iostream>

using namespace gtsam;

int main()
{
  NonlinearFactorGraph graph;
  Vector3 sigmas;
  sigmas << 0.3,0.3,0.3;
  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigmas);

  Values vals;
  Matrix51 Theta;
  Theta << 0.1,0.1,0.1,0.1,0.1;

  vals.insert(50, Theta);
  Expression<Matrix51> theta_expr(50);
  Matrix31 X, X2;
  X << 0,0,0;
  X2 << 0.3,0.3,0.3;

  Predictor predict(X, 0.1);
  Expression<Matrix31> predict_expr(predict, theta_expr);
  graph.addExpressionFactor(R, X2, predict_expr);
  graph.print("Factor Graph:\n");
  LevenbergMarquardtOptimizer optimizer(graph, vals);
  Values result = optimizer.optimize();
  result.print("Final Result:\n");

  Matrix31 Xpredict;
  Matrix51 ThetaOptimized = result.at<Matrix51>(50);
  Matrix35 J;
  Xpredict = predict(ThetaOptimized, J);
  std::cout << "Predicted X: \n" << Xpredict << std::endl;
}
