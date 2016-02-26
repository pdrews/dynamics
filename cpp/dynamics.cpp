#include "dynamics.h"
#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Values.h>

#include <iostream>

using namespace gtsam;

int main()
{

  //Load up all of our data
  int numSamples;
  Matrix<double,X_DIM,Dynamic> X;
  Matrix<double,U_DIM,Dynamic> U;
  Matrix<double,1,Dynamic> heading;
  double dt = 0.02;

  //Iterate through each pair in our data, putting it in our FG
  Expression<MatrixTheta> theta_expr(50);
  Vector3 sigmas;
  sigmas << 0.3,0.3,0.3;
  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigmas);

  for(int i=0; i<numSamples-1; i++) {
    Predictor predict(X.col(i),U.col(i),heading(i), dt);
    Expression<MatrixX> predict_expr(predict, theta_expr);
    graph.addExpressionFactor(R, X.col(i+1), predict_expr);
  }
  //Solve the FG and spit out the thetas
  NonlinearFactorGraph graph;

  Values vals;
  MatrixTheta Theta = MatrixTheta::Zeros();
  vals.insert(50, Theta);


//  graph.print("Factor Graph:\n");
  LevenbergMarquardtOptimizer optimizer(graph, vals);
  Values result = optimizer.optimize();
  result.print("Final Result:\n");

  // Get the average and mean squared error
  MatrixX squared_error = MatrixX::Zeros();
  for(int i=0; i<numSamples-1; i++) {
    Predictor predict(X.col(i),U.col(i),heading(i), dt);
    MatrixX predicted = predict(theta_expr);
    squared_error += (predicted - X.col(i+1)).cwiseProduct(predicted - X.col(i+1));
  }
  squared_error = squared_error.array() / (numSamples-1);
  std::cout << "Squared error is \n" << squared_error << std::endl;

}
