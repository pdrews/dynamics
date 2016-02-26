#include "dynamics.h"
#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Values.h>

#include <boost/tokenizer.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

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
  data.close();
  data.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/data.txt");
  int sampleNum = 0;
  while(getline(data, line)) {
//    cout << line << endl;
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char>> toks(line,sep);
    boost::tokenizer<boost::char_separator<char>>::iterator tok = toks.begin();
//    cout << *tok << " ";
    X(0,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    X(1,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    X(2,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    X(3,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    heading(0,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    U(0,sampleNum) = atof(tok->c_str());
    tok++;
//    cout << *tok << " ";
    U(1,sampleNum) = atof(tok->c_str());
//    cout << "X " << X.col(sampleNum) << endl;
//    cout << "U " << U.col(sampleNum) << endl;
//    cout << "heading " << heading(sampleNum) << endl;
    sampleNum++;
  }
  cout << "Finished reading " << sampleNum << " samples" << endl;
  double dt = 0.02;

  //Iterate through each pair in our data, putting it in our FG
  Expression<MatrixTheta> theta_expr(50);
  MatrixX sigmas = MatrixX::Constant(0.5);

  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigmas);

  NonlinearFactorGraph graph;
  for(int i=0; i<numSamples-1; i++) {
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
  MatrixX squared_error = MatrixX::Zero();
  MatrixBasis J;
  MatrixTheta thetaEst = result.at<MatrixTheta>(50);
  for(int i=0; i<numSamples-1; i++) {
    Predictor predict(X.col(i),U.col(i),heading(i), dt);
    MatrixX predicted = predict(thetaEst,J);
    squared_error += (predicted - X.col(i+1)).cwiseProduct(predicted - X.col(i+1));
  }
  squared_error = squared_error.array() / (numSamples-1);
  cout << "Squared error is \n" << squared_error << endl;

}
