#include "dynamics.h"
#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/nonlinear/Values.h>

#include <boost/tokenizer.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace gtsam;
using namespace std;

int main()
{

  //Load up all of our data
  ifstream data;
  data.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/data.txt");
  string line;
  int numSamples = 0;
  while(getline(data, line)) {numSamples++;}
  cout << "Found " << numSamples << " Samples" << endl;
  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> X(X_DIM,numSamples);
  Eigen::Matrix<double,U_DIM,Eigen::Dynamic> U(U_DIM,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> heading(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> beta(1,numSamples);
  data.close();
  data.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/data.txt");
  int sampleNum = 0;
  while(getline(data, line)) {
    double left, right;
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char>> toks(line,sep);
    boost::tokenizer<boost::char_separator<char>>::iterator tok = toks.begin();
    //Roll
    X(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Ux
    X(1,sampleNum) = atof(tok->c_str());
    tok++;
    //Uy
    X(2,sampleNum) = atof(tok->c_str());
    tok++;
    //Yaw rate
    X(3,sampleNum) = atof(tok->c_str());
    tok++;
    //Yaw
    heading(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Throttle
    U(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Steering
    U(1,sampleNum) = atof(tok->c_str());
    left = atof(tok->c_str());
    tok++;
    right = atof(tok->c_str());
    X(4,sampleNum) = (left+right) / 2.0;

    beta(0,sampleNum) = atan2(X(2,sampleNum),X(1,sampleNum));

    sampleNum++;
  }
  cout << "Finished reading " << sampleNum << " samples" << endl;
  double dt = 1.0/40.0;

  //Iterate through each pair in our data, putting it in our FG
  Expression<MatrixTheta> theta_expr(50);
  MatrixX sigmas = MatrixX::Constant(0.5);

  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigmas);

  NonlinearFactorGraph graph;
  for(int i=0; i<numSamples-1; i++) {
    if(std::abs(X(1,i) - X(1,i+1))/dt > 15) {
      cout << "Skipping a sample between bag files" << endl;
      continue;
    }
    // First add states
    Predictor predict(X.col(i),U.col(i),heading(0,i), dt);
    Expression<MatrixX> predict_expr(predict, theta_expr);
    graph.addExpressionFactor<MatrixX>(R, X.col(i+1), predict_expr);
  }

  //Solve the FG and spit out the thetas
  Values vals;
  MatrixTheta Theta = MatrixTheta::Zero();
  vals.insert(50, Theta);
  //  graph.print("Factor Graph:\n");
  LevenbergMarquardtOptimizer optimizer(graph, vals);
  Values result = optimizer.optimize();
  result.print("Final Result:\n");

  // Get the average and mean squared error
  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_pred(X_DIM,numSamples);
  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_actual(X_DIM,numSamples);

  MatrixX squared_error = MatrixX::Zero();
  MatrixX squared_error_der = MatrixX::Zero();
  MatrixBasis J;
  MatrixTheta thetaEst = result.at<MatrixTheta>(50);
  ofstream sed, sedp, cov;
  sed.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/sed.txt");
  sedp.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/sedp.txt");
  int skips = 0;
  for(int i=0; i<numSamples-1; i++) {
    if(std::abs(X(1,i) - X(1,i+1))/dt > 15) {
      cout << "Skipping an error between bag files" << endl;
      continue;

    }
    Predictor predict(X.col(i),U.col(i),heading(i), dt);
    MatrixX predicted = predict(thetaEst,J);
    squared_error += (predicted - X.col(i+1)).cwiseProduct(predicted - X.col(i+1));

    //Calc der stuff
    MatrixX der = (X.col(i) - X.col(i+1)) / dt;
    MatrixX derPred = (X.col(i) - predicted) / dt;
    sed << der(0) << "," << der(1) << "," << der(2) << "," << der(3) << endl;
    sedp << derPred(0) << "," << derPred(1) << "," << derPred(2) << "," << derPred(3) << endl;
//    der_pred.col(i) = derPred.transpose();
//    der.col(i) = der.transpose();
    squared_error_der += (der - derPred).cwiseProduct(der - derPred);
  }
  squared_error = squared_error.array() / ((numSamples-skips)-1);
  squared_error_der = squared_error_der.array() / ((numSamples-skips)-1);
  cout << "Squared error is \n" << squared_error << endl;
  cout << "Squared error of derivative is \n" << squared_error_der << endl;

  //Print all the marginals
  Marginals marginals(graph, result);
  Matrix covMatrix = marginals.marginalCovariance(50);
  cov.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/cov.txt");
  for(int i=0; i<covMatrix.rows(); i++) {
    for(int j=0; j<covMatrix.cols(); j++){
      cov << covMatrix(i,j) << ",";
    }
    cov << endl;
  }
  sed.close();
  sedp.close();
  cov.close();
}
