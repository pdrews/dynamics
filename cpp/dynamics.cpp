#include "dynamics.h"
#include "PredictorExpression.h"
//#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>


#include <boost/tokenizer.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace gtsam;
using namespace std;
using symbol_shorthand::T;
using symbol_shorthand::V;
using symbol_shorthand::B;
using symbol_shorthand::P;
int main()
{

  //Load up all of our data
  ifstream data;
  data.open("../data.txt");
  string line;
  int numSamples = 0;
  while(getline(data, line)) {numSamples++;}
  cout << "Found " << numSamples << " Samples" << endl;
  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> X(X_DIM,numSamples);
  Eigen::Matrix<double,U_DIM,Eigen::Dynamic> U(U_DIM,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> heading(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> betam(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> Vxm(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> ThetaDm(1,numSamples);

  data.close();
  data.open("../data.txt");
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

    betam(0,sampleNum) = atan2(X(2,sampleNum),X(1,sampleNum));
    ThetaDm(0,sampleNum) = X(3,sampleNum);
    Vxm(0,sampleNum) = X(1,sampleNum);

    sampleNum++;
  }
  cout << "Finished reading " << sampleNum << " samples" << endl;
  double dt = 1.0/40.0;

  //Define expressions for all of our parameters
  double_ Calpha('p',0);    //About 10000???
  double_ mu('p',1);        //About 0.5??
  double_ mass('p',2);      //About 20kg
  double_ a('p',3);         //About .45m
  double_ b('p',4);         //About .35m
  double_ Iz('p', 5);       //About 2??????
  //Parameters for fitting steering and throttle
  double_ c1('p',6); //throttle
  double_ c2('p',7); //steering
//  double_ c3('p',8);
  //Derived parameters
  double_ Fzf = (9.8 * 0.5) * MulExpression(DivExpression(b,(a+b)), mass);
  double_ Fzr = (9.8 * 0.5) * MulExpression(DivExpression(a,(a+b)), mass);

  NonlinearFactorGraph graph;
  Values vals;

  //Put in initial guesses
  vals.insert(P(0), 10000.0);
  vals.insert(P(1), 0.5);
  vals.insert(P(2), 20.0);
  vals.insert(P(3), 0.45);
  vals.insert(P(4), 0.35);
  vals.insert(P(5), 2.0);

  vals.insert(P(6), 18.0);
  vals.insert(P(7), 70.0);
//  vals.insert(P(8), 15.0);


  //Iterate through each pair in our data, putting it in our FG
  Matrix1 sigma = Matrix1::Constant(0.1);

  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigma);


  for(int i=0; i<numSamples-2; i++) {
//  for(int i=0; i<2; i++) {
    //put vx, beta, thetadot in the graph
//    vals.insert(B(i), betam(0,i));
//    vals.insert(T(i), ThetaDm(0,i));
//    vals.insert(V(i), Vxm(0,i));
    //Prior factors
//    PriorFactor<double> priorB(B(i), betam(0,i), R);
//    PriorFactor<double> priorT(T(i), ThetaDm(0,i), R);
//    PriorFactor<double> priorV(V(i), Vxm(0,i), R);
//    graph.add(priorB);
//    graph.add(priorT);
//    graph.add(priorV);
    double_ beta_(betam(0,i)), thetad_(ThetaDm(0,i)), vx_(Vxm(0,i));
//    double_ beta_('b',i), thetad_('t',i), vx_('v',i);
    double_ steering_ = U(1,i)*c2;
    double_ throttle_ = U(0,i)*c1;

    if(std::abs(X(1,i) - X(1,i+1))/dt > 15) {
      cout << "Skipping a sample between bag files" << endl;
      continue;
    }
    //Add prediction factor requirements
    double_ gamma = gammaExpression(mu, Fzr, throttle_);
    double_ alphaf = alphafExpression(beta_, a, vx_, thetad_, steering_);
    double_ alphar = alpharExpression(beta_, a, vx_, thetad_);
    double_ Fyr = FyrExpression(mu, Fzr, gamma, Calpha, alphar);
    double_ Fyf = FyfExpression(mu, Fzr, Calpha, alphaf);

    //Add prediction factors
    double_ betadot = betadotExpression(mu, mass, vx_, thetad_, Fyf, Fyr);
    double_ vxdot = VxdExpression(steering_, Fyf, throttle_, mass, vx_, beta_, thetad_);
    double_ thetaddot = thetaddExpression(a, Fyf, b, Fyr, Iz);

    graph.addExpressionFactor<double>(R, betam(0,i+1), dt*betadot);
    graph.addExpressionFactor<double>(R, Vxm(0,i+1), dt*vxdot);
    graph.addExpressionFactor<double>(R, ThetaDm(0,i+1), dt*thetaddot);

//    cout << "Expression for Beta " << i << endl;
//    betadot.print("Beta: ");

  }
//  graph.print("Our factor graph:\n");
//  vals.print("InitialValues\n");


  //Solve the FG and spit out the thetas
//  MatrixTheta Theta = MatrixTheta::Zero();
//  vals.insert(50, Theta);
//  //  graph.print("Factor Graph:\n");
  LevenbergMarquardtParams params;
//  params.verbosity = LevenbergMarquardtParams::Verbosity::DELTA;
//  params.verbosityLM = LevenbergMarquardtParams::VerbosityLM::TRYDELTA;
  LevenbergMarquardtOptimizer optimizer(graph, vals, params);
//  optimizer.iterate();
  Values result = optimizer.optimize();
  cout << "We did " << optimizer.iterations() << " Iterations" << endl;
//  optimizer.print("Optimizer:");
  cout << "mu: " << result.at<double>(P(1)) << endl;
  cout << "mass: " << result.at<double>(P(2)) << endl;
  cout << "a: " << result.at<double>(P(3)) << endl;
  cout << "b: " << result.at<double>(P(4)) << endl;
      result.print("Final Result:\n");
//
//  // Get the average and mean squared error
//  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_pred(X_DIM,numSamples);
//  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_actual(X_DIM,numSamples);
//
//  MatrixX squared_error = MatrixX::Zero();
//  MatrixX squared_error_der = MatrixX::Zero();
//  MatrixBasis J;
//  MatrixTheta thetaEst = result.at<MatrixTheta>(50);
//  ofstream sed, sedp, cov;
//  sed.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/sed.txt");
//  sedp.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/sedp.txt");
//  int skips = 0;
//  for(int i=0; i<numSamples-1; i++) {
//    if(std::abs(X(1,i) - X(1,i+1))/dt > 15) {
//      cout << "Skipping an error between bag files" << endl;
//      continue;
//
//    }
//    Predictor predict(X.col(i),U.col(i),heading(i), dt);
//    MatrixX predicted = predict(thetaEst,J);
//    squared_error += (predicted - X.col(i+1)).cwiseProduct(predicted - X.col(i+1));
//
//    //Calc der stuff
//    MatrixX der = (X.col(i) - X.col(i+1)) / dt;
//    MatrixX derPred = (X.col(i) - predicted) / dt;
//    sed << der(0) << "," << der(1) << "," << der(2) << "," << der(3) << endl;
//    sedp << derPred(0) << "," << derPred(1) << "," << derPred(2) << "," << derPred(3) << endl;
////    der_pred.col(i) = derPred.transpose();
////    der.col(i) = der.transpose();
//    squared_error_der += (der - derPred).cwiseProduct(der - derPred);
//  }
//  squared_error = squared_error.array() / ((numSamples-skips)-1);
//  squared_error_der = squared_error_der.array() / ((numSamples-skips)-1);
//  cout << "Squared error is \n" << squared_error << endl;
//  cout << "Squared error of derivative is \n" << squared_error_der << endl;
//
//  //Print all the marginals
//  Marginals marginals(graph, result);
//  Matrix covMatrix = marginals.marginalCovariance(50);
//  cov.open("/media/data/logs/MPPI_SystemID/Current/Model_Parameters/cov.txt");
//  for(int i=0; i<covMatrix.rows(); i++) {
//    for(int j=0; j<covMatrix.cols(); j++){
//      cov << covMatrix(i,j) << ",";
//    }
//    cov << endl;
//  }
//  sed.close();
//  sedp.close();
//  cov.close();
}
